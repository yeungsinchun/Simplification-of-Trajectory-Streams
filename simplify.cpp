#include <atomic>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <iomanip>
#include <memory>
#include <mutex>
#include <thread>

#include "simplify_geometry.h"
#include "simplify_io.h"
#include "timer.h"

// ===========================================================================
//  Core algorithm (opt-in TIMER sites; active only when --time is set)
// ===========================================================================

// Per-thread state for one anchor update: the wedge F, clip buffers, and the
// per-step cache for anchors whose F is the whole bbox (they all clip the
// same polygon against the same Gi, so the result is identical).
struct AnchorWorker {
    std::vector<Point> F, bbox_result;
    sh_double::FastClipBuffers clip_buffers;
    AxisBounds bbox_bounds;
    bool bbox_cached = false, bbox_hit = false;
};

// Persistent helpers for the anchor loop. Within one stream step every live
// anchor's find_F + clip reads only shared step data (P, Gi) and writes only
// its own S/bounds/latch, so anchors can be updated in any order on any
// thread. Survivors are compacted afterwards in index order, so the output
// is identical to the sequential loop.
//
// Work is claimed in small chunks through one ticket holding (generation,
// next index); a claim is valid only while the generation still matches, and
// run() waits for finished items rather than for every helper, so a helper
// the OS has descheduled never stalls a step. Idle helpers spin for up to
// kIdleSpin and then sleep on a condition variable.
class AnchorPool {
  public:
    explicit AnchorPool(int participants) : workers_(participants) {
        for (int t = 1; t < participants; ++t)
            threads_.emplace_back([this, t] { helper_loop(t); });
    }
    ~AnchorPool() {
        stop_.store(true);
        publish(0, nullptr, nullptr);
        for (auto &thread : threads_)
            thread.join();
    }
    AnchorPool(const AnchorPool &) = delete;
    AnchorPool &operator=(const AnchorPool &) = delete;

    int participants() const { return int(workers_.size()); }
    AnchorWorker &worker(int t) { return workers_[t]; }

    // Calls fn(begin, end, worker) on chunks covering [0, count) across all
    // participants (the caller included); returns once every call finished.
    template <class Fn> void run(int count, Fn &fn) {
        const uint64_t generation = publish(
            count, &fn, [](void *job, int begin, int end, AnchorWorker &w) {
                (*static_cast<Fn *>(job))(begin, end, w);
            });
        drain(0, generation);
        while (done_.load(std::memory_order_acquire) != count)
            cpu_relax();
    }

  private:
    using Call = void (*)(void *, int, int, AnchorWorker &);
    struct Job {
        std::atomic<void *> ctx{nullptr};
        std::atomic<Call> call{nullptr};
        std::atomic<int> count{0};
    };
    static constexpr int kChunk = 4;

    // Steps arrive every few microseconds while a run is busy; helpers idle
    // this long (e.g. after the run) go to sleep.
    static constexpr auto kIdleSpin = std::chrono::milliseconds(1);

    static void cpu_relax() {
#if defined(__aarch64__)
        asm volatile("yield");
#elif defined(__x86_64__) || defined(__i386__)
        __builtin_ia32_pause();
#endif
    }

    // Jobs alternate between two slots; a slot is rewritten only after the
    // generation using it has finished and the next one has been published.
    uint64_t publish(int count, void *ctx, Call call) {
        const uint64_t generation = (ticket_.load() >> 32) + 1;
        Job &job = jobs_[generation & 1];
        job.ctx.store(ctx, std::memory_order_relaxed);
        job.call.store(call, std::memory_order_relaxed);
        job.count.store(count, std::memory_order_relaxed);
        done_.store(0, std::memory_order_relaxed);
        ticket_.store(generation << 32);
        if (sleepers_.load() > 0) {
            std::lock_guard<std::mutex> lock(mutex_);
            wake_.notify_all();
        }
        return generation;
    }

    void drain(int t, uint64_t generation) {
        AnchorWorker &w = workers_[t];
        const Job &job = jobs_[generation & 1];
        uint64_t ticket = ticket_.load(std::memory_order_acquire);
        for (;;) {
            if ((ticket >> 32) != generation)
                return;
            const int begin = int(uint32_t(ticket));
            const int count = job.count.load(std::memory_order_relaxed);
            if (begin >= count)
                return;
            if (!ticket_.compare_exchange_weak(ticket, ticket + kChunk,
                                               std::memory_order_acq_rel,
                                               std::memory_order_acquire))
                continue;
            const int end = std::min(begin + kChunk, count);
            void *ctx = job.ctx.load(std::memory_order_relaxed);
            const Call call = job.call.load(std::memory_order_relaxed);
            call(ctx, begin, end, w);
            done_.fetch_add(end - begin, std::memory_order_release);
            ticket = ticket_.load(std::memory_order_acquire);
        }
    }

    void helper_loop(int t) {
        uint64_t seen = 0;
        for (;;) {
            uint64_t generation;
            auto idle_since = std::chrono::steady_clock::now();
            for (int spins = 0; (generation = ticket_.load() >> 32) == seen;
                 ++spins) {
                cpu_relax();
                if (spins % 64 != 63 ||
                    std::chrono::steady_clock::now() - idle_since < kIdleSpin)
                    continue;
                std::unique_lock<std::mutex> lock(mutex_);
                sleepers_.fetch_add(1);
                wake_.wait(lock, [&] { return (ticket_.load() >> 32) != seen; });
                sleepers_.fetch_sub(1);
            }
            seen = generation;
            if (stop_.load())
                return;
            drain(t, generation);
        }
    }

    std::vector<AnchorWorker> workers_;
    std::vector<std::thread> threads_;
    Job jobs_[2];
    std::atomic<uint64_t> ticket_{0};
    std::atomic<int> done_{0}, sleepers_{0};
    std::atomic<bool> stop_{false};
    std::mutex mutex_;
    std::condition_variable wake_;
};

// Parallel updates pay a hand-off and a join per step, so only steps with at
// least this much work (live anchors x Gi edges) use the pool; small steps
// (coarse epsilon, few live anchors) stay on the calling thread.
constexpr size_t kParallelMinWork = 256;

struct StabScratch {
    std::vector<Point> P, Gi;
    std::vector<std::vector<Point>> S;
    std::vector<int> active;
    std::vector<AxisBounds> stab_bounds;
    std::vector<uint8_t> anchor_outside;
    PreparedClipPolygon prepared_Gi;
    AnchorWorker serial;
    int threads = 1;
    std::unique_ptr<AnchorPool> pool;
    std::vector<uint8_t> keep;
};

// Steps anchors active[begin..end): S[i] = F(S[i], P[i]) ∩ Gi, and keep[k] is
// 0 when anchor active[k] dies. Out of line so the sequential loop and the
// pool share a single copy of the inlined geometry (a second copy pushes GCC
// past its inlining budget and slows the sequential loop by ~20%); one call
// covers a whole step or chunk, so the call itself costs nothing measurable.
__attribute__((noinline)) void update_anchors(StabScratch &scratch, int begin,
                                              int end, AnchorWorker &w) {
    auto &S = scratch.S;
    auto &stab_bounds = scratch.stab_bounds;
    auto &anchor_outside = scratch.anchor_outside;
    for (int k = begin; k < end; ++k) {
        const int i = scratch.active[k];
        bool full_bbox, disjoint;
        {
            TIMER("find_F");
            full_bbox = find_F(scratch.P[i], S[i], w.F, &scratch.prepared_Gi,
                               &disjoint, &stab_bounds[i], &anchor_outside[i]);
        }
        if (disjoint) {
            scratch.keep[k] = false;
            continue;
        }
        TIMER("intersect");
        if (full_bbox && w.bbox_cached) {
            S[i] = w.bbox_result;
            stab_bounds[i] = w.bbox_bounds;
            scratch.keep[k] = w.bbox_hit;
            continue;
        }
        const bool hit = intersect_prepared(
            w.F, scratch.prepared_Gi, S[i],
            anchor_outside[i] && !full_bbox ? nullptr : &stab_bounds[i],
            w.clip_buffers);
        if (full_bbox) {
            w.bbox_result = S[i];
            w.bbox_bounds = stab_bounds[i];
            w.bbox_hit = hit;
            w.bbox_cached = true;
        }
        scratch.keep[k] = hit;
    }
}

int get_longest_stab(const std::vector<Point> &stream, int cur,
                     std::vector<Point> &simplified, double EPSILON,
                     double DELTA, StabScratch &scratch) {
    TIMER("get_longest_stab");
    const Point& p0 = stream[cur];
    auto &P = scratch.P;
    auto &Gi = scratch.Gi;
    auto &S = scratch.S;
    auto &active = scratch.active;
    auto &stab_bounds = scratch.stab_bounds;
    auto &anchor_outside = scratch.anchor_outside;
    {
        TIMER("boundary_P");
        P = get_boundary_points_from_grid(p0, EPSILON, DELTA);
    }
    std::array<Point, 2> buffer = {p0, p0};
    const int Pn = (int)P.size();
    S.resize(Pn);
    stab_bounds.resize(Pn);
    anchor_outside.assign(Pn, 0);
    active.clear();
    active.reserve(Pn);
    for (int i = 0; i < Pn; ++i) {
        S[i].clear();
        S[i].push_back(P[i]);
        active.push_back(i);
    }

    cur++;
    while (cur < int(stream.size())) {
        {
            TIMER("hull_Gi");
            Gi = get_conv_from_grid(stream[cur], EPSILON, DELTA);
        }
        prepare_clip_polygon(Gi, scratch.prepared_Gi);
        const int count = int(active.size());
        scratch.keep.resize(count);
        // Timers keep global state, so --time runs stay sequential.
        if (scratch.threads > 1 && !timer_detail::enabled() &&
            size_t(count) * scratch.prepared_Gi.edges.size() >=
                kParallelMinWork) {
            if (!scratch.pool)
                scratch.pool = std::make_unique<AnchorPool>(scratch.threads);
            for (int t = 0; t < scratch.pool->participants(); ++t)
                scratch.pool->worker(t).bbox_cached = false;
            auto job = [&](int begin, int end, AnchorWorker &w) {
                update_anchors(scratch, begin, end, w);
            };
            scratch.pool->run(count, job);
        } else {
            scratch.serial.bbox_cached = false;
            update_anchors(scratch, 0, count, scratch.serial);
        }
        size_t surviving = 0;
        for (int k = 0; k < count; ++k)
            if (scratch.keep[k])
                active[surviving++] = active[k];
        active.resize(surviving);
        if (active.empty())
            break;
        const int chosen = active.back();
        buffer[0] = P[chosen];
        buffer[1] = S[chosen].front();
        cur++;
    }
    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    return cur;
}

// ===========================================================================
//  Web-server trace mode (--web-server)
// ===========================================================================
//
// Mirrors get_longest_stab/simplify exactly, but instead of only emitting the
// final two-point segment per prefix, it records every intermediate value the
// paper's construction produces (the boundary anchors P, the delta-disk hull
// Gi, the free-space wedge F(S,p), and the resulting stab region S) at every
// step of every prefix.  With --json-stream (used by the Flask server), stdout
// is NDJSON: header line, one prefix per line, then done.  Without it, the
// whole trace is one JSON object.  No human-readable text is ever written in
// this mode so stdout stays machine-readable.
namespace webtrace {

struct Candidate {
    int grid_pt_idx = 0;
    bool alive = true;
    std::vector<Point> F;      // free-space wedge F(S_{i-1}[p0], pi) fed to intersect() this step
    std::vector<Point> F_Si;   // free-space wedge F(S_i[p0], pi) computed after intersect()
    std::vector<Point> S;      // resulting stab region (new_S[i]) after this step, if alive
};

struct StepTrace {
    int stream_idx = 0;        // index into the full input stream of the point consumed this step
    Point pi{0, 0};
    std::vector<Point> Gi;     // conv(G_i): convex hull of the delta-disk grid samples around pi
    std::vector<Candidate> candidates;
    std::array<Point, 2> buffer{Point(0, 0), Point(0, 0)};
};

struct PrefixTrace {
    Point p0{0, 0};
    int p0_idx = 0;            // index of p0 in the full stream
    int end_idx = 0;           // stream index of the last point consumed (output[1] comes from here)
    std::vector<Point> P;      // boundary (hull) anchors of the delta-disk grid samples for this prefix
    std::vector<StepTrace> steps;
    std::array<Point, 2> output{Point(0, 0), Point(0, 0)};
};

// ---------------------------------------------------------------------------
//  Minimal dependency-free JSON writer
// ---------------------------------------------------------------------------

inline void write_num(std::ostream& os, double v) {
    if (!std::isfinite(v)) { os << "null"; return; }
    os << v;
}

inline void write_point(std::ostream& os, const Point& p) {
    os << '[';
    write_num(os, CGAL::to_double(p.x()));
    os << ',';
    write_num(os, CGAL::to_double(p.y()));
    os << ']';
}

inline void write_points(std::ostream& os, const std::vector<Point>& pts) {
    os << '[';
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (i) os << ',';
        write_point(os, pts[i]);
    }
    os << ']';
}

inline void write_prefix(std::ostream& os, const PrefixTrace& p) {
    os << "{\"p0\":";       write_point(os, p.p0); os << ',';
    os << "\"p0_idx\":"  << p.p0_idx << ',';
    os << "\"end_idx\":" << p.end_idx << ',';
    os << "\"P\":";    write_points(os, p.P);   os << ',';
    os << "\"output\":["; write_point(os, p.output[0]); os << ','; write_point(os, p.output[1]); os << "],";
    os << "\"steps\":[";
    for (std::size_t si = 0; si < p.steps.size(); ++si) {
        if (si) os << ',';
        const StepTrace& s = p.steps[si];
        os << "{\"stream_idx\":" << s.stream_idx << ',';
        os << "\"pi\":"; write_point(os, s.pi); os << ',';
        os << "\"Gi\":"; write_points(os, s.Gi); os << ',';
        os << "\"buffer\":["; write_point(os, s.buffer[0]); os << ','; write_point(os, s.buffer[1]); os << "],";
        os << "\"candidates\":[";
        for (std::size_t ci = 0; ci < s.candidates.size(); ++ci) {
            if (ci) os << ',';
            const Candidate& c = s.candidates[ci];
            os << "{\"idx\":" << c.grid_pt_idx
               << ",\"alive\":" << (c.alive ? "true" : "false") << ',';
            os << "\"F\":"; write_points(os, c.F); os << ',';
            os << "\"F_Si\":"; write_points(os, c.F_Si); os << ',';
            os << "\"S\":"; write_points(os, c.S);
            os << "}";
        }
        os << "]}";
    }
    os << "]}";
}

inline void write_stream_header(std::ostream& os, double EPSILON, double DELTA,
                                const std::vector<Point>& stream) {
    os << std::setprecision(17);
    os << "{\"type\":\"header\",";
    os << "\"eps\":";       write_num(os, EPSILON);                       os << ',';
    os << "\"delta\":";     write_num(os, DELTA);                         os << ',';
    os << "\"grid_val\":";  write_num(os, GRID_val(EPSILON, DELTA));      os << ',';
    os << "\"r_val\":";     write_num(os, R_val(EPSILON, DELTA));         os << ',';
    os << "\"expected_frechet\":"; write_num(os, std::sqrt(expected_frechet_squared)); os << ',';
    os << "\"bbox\":[";
    write_num(os, BMIN); os << ','; write_num(os, BMIN); os << ',';
    write_num(os, BMAX); os << ','; write_num(os, BMAX);
    os << "],";
    os << "\"stream\":"; write_points(os, stream);
    os << "}\n";
}

inline void write_stream_prefix(std::ostream& os, const PrefixTrace& p) {
    os << std::setprecision(17);
    os << "{\"type\":\"prefix\",\"data\":";
    write_prefix(os, p);
    os << "}\n";
}

inline void write_stream_done(std::ostream& os, double time_ms,
                              const std::vector<Point>& simplified) {
    os << std::setprecision(17);
    os << "{\"type\":\"done\",";
    os << "\"time_ms\":"; write_num(os, time_ms); os << ',';
    os << "\"simplified\":"; write_points(os, simplified); os << ',';
    os << "\"frechet_distance\":null";
    os << "}\n";
}

inline void write_json(std::ostream& os, double EPSILON, double DELTA, double time_ms,
                       const std::vector<Point>& stream,
                       const std::vector<Point>& simplified,
                       const std::vector<PrefixTrace>& prefixes) {
    os << std::setprecision(17);
    os << "{";
    os << "\"eps\":";       write_num(os, EPSILON);                       os << ',';
    os << "\"delta\":";     write_num(os, DELTA);                         os << ',';
    os << "\"time_ms\":";   write_num(os, time_ms);                       os << ',';
    os << "\"grid_val\":";  write_num(os, GRID_val(EPSILON, DELTA));      os << ',';
    os << "\"r_val\":";     write_num(os, R_val(EPSILON, DELTA));         os << ',';
    os << "\"expected_frechet\":"; write_num(os, std::sqrt(expected_frechet_squared)); os << ',';
    os << "\"bbox\":[";
    write_num(os, BMIN); os << ','; write_num(os, BMIN); os << ',';
    write_num(os, BMAX); os << ','; write_num(os, BMAX);
    os << "],";
    os << "\"stream\":";      write_points(os, stream);      os << ',';
    os << "\"simplified\":"; write_points(os, simplified);   os << ',';
    os << "\"prefixes\":[";
    for (std::size_t pi = 0; pi < prefixes.size(); ++pi) {
        if (pi) os << ',';
        write_prefix(os, prefixes[pi]);
    }
    os << "]}\n";
}

}  // namespace webtrace

// Web-trace twin of get_longest_stab: identical control flow, additionally
// records P, Gi, F[i], new_S[i], alive/dead, and buffer at every step.
int get_longest_stab_web(const std::vector<Point>& stream, int cur,
                         std::vector<Point>& simplified,
                         double EPSILON, double DELTA,
                         std::vector<webtrace::PrefixTrace>& prefixes) {
    const Point& p0 = stream[cur];
    std::vector<Point> P = get_boundary_points_from_grid(p0, EPSILON, DELTA);
    std::array<Point, 2> buffer = {p0, p0};
    const int Pn = (int)P.size();
    std::vector<std::vector<Point>> S(Pn);
    for (int i = 0; i < Pn; ++i) S[i] = {P[i]};
    int dead_cnt = 0;
    std::vector<int> dead(Pn);
    std::vector<std::vector<Point>> new_S(Pn);
    std::vector<std::vector<Point>> F(Pn);
    std::vector<Point> Gi;

    webtrace::PrefixTrace trace;
    trace.p0 = p0;
    trace.p0_idx = cur;   // p0 = stream[cur] (cur not yet incremented)
    trace.P = P;

    cur++;
    while (cur < int(stream.size())) {
        Gi = get_conv_from_grid(stream[cur], EPSILON, DELTA);

        webtrace::StepTrace step;
        step.stream_idx = cur;
        step.pi = stream[cur];
        step.Gi = Gi;
        step.candidates.reserve(Pn);

        for (int i = 0; i < Pn; ++i) {
            webtrace::Candidate cand;
            cand.grid_pt_idx = i;
            if (dead[i]) {
                cand.alive = false;
                step.candidates.push_back(std::move(cand));
                continue;
            }
            
            find_F(P[i], S[i], F[i]);
            cand.F = F[i];
            
            if (!intersect(F[i], Gi, new_S[i])) {
                dead[i] = true;
                dead_cnt++;
                cand.alive = false;
            } else {
                cand.alive = true;
                cand.S = new_S[i];
                std::vector<Point> F_Si_temp;
                find_F(P[i], new_S[i], F_Si_temp);
                cand.F_Si = F_Si_temp;
            }
            step.candidates.push_back(std::move(cand));
        }

        bool has_candidate = false;
        for (int i = Pn - 1; i >= 0 && !has_candidate; --i) {
            if (dead[i] || new_S[i].empty()) continue;
            buffer[0] = P[i];
            buffer[1] = new_S[i].front();
            has_candidate = true;
        }
        step.buffer = buffer;
        trace.steps.push_back(std::move(step));

        if (!has_candidate || dead_cnt == Pn) break;
        for (int i = 0; i < Pn; ++i)
            if (!dead[i]) S[i].swap(new_S[i]);
        cur++;
    }
    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    trace.end_idx = cur;
    trace.output = buffer;
    prefixes.push_back(std::move(trace));
    return cur;
}

std::vector<Point> simplify_web(const std::vector<Point>& stream,
                                double EPSILON, double DELTA) {
    std::vector<Point> simplified;
    std::vector<webtrace::PrefixTrace> prefixes;
    configure_bbox(stream, EPSILON, DELTA);

    std::ostream* stream_out = nullptr;
    if (json_stream_flag && json_output_path.empty()) {
        stream_out = &std::cout;
        if (!stream.empty())
            get_boundary_points_from_grid(stream[0], EPSILON, DELTA);
        webtrace::write_stream_header(*stream_out, EPSILON, DELTA, stream);
        stream_out->flush();
    }

    double core_ms = 0;
    int cur = 0;
    while (cur != int(stream.size())) {
        auto t0 = std::chrono::high_resolution_clock::now();
        cur = get_longest_stab_web(stream, cur, simplified, EPSILON, DELTA, prefixes);
        core_ms += std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - t0).count();
        if (stream_out) {
            webtrace::write_stream_prefix(*stream_out, prefixes.back());
            stream_out->flush();
        }
    }
    
    std::cerr << "SIMPLIFY_CORE_MS: " << std::fixed << std::setprecision(4) << core_ms << '\n';
    
    if (stream_out) {
        webtrace::write_stream_done(*stream_out, core_ms, simplified);
        stream_out->flush();
    } else if (!json_output_path.empty()) {
        std::ofstream ofs(json_output_path);
        if (!ofs) {
            std::cerr << "Failed to open output file: " << json_output_path << '\n';
            return simplified;
        }
        webtrace::write_json(ofs, EPSILON, DELTA, core_ms, stream, simplified, prefixes);
        ofs.close();
    } else {
        webtrace::write_json(std::cout, EPSILON, DELTA, core_ms, stream, simplified, prefixes);
    }
    return simplified;
}

std::vector<Point> simplify(const std::vector<Point>& stream,
                            double EPSILON, double DELTA) {
    std::vector<Point> simplified;
    configure_bbox(stream, EPSILON, DELTA);
    if (time_flag) {
        reset_timing();
        timer_detail::enabled() = true;
    }
    auto t0 = std::chrono::high_resolution_clock::now();
    {
        TIMER("total");
        StabScratch scratch;
        scratch.threads = simplify_threads;
        int cur = 0;
        while (cur != int(stream.size()))
            cur = get_longest_stab(stream, cur, simplified, EPSILON, DELTA,
                                   scratch);
    }
    double ms = std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t0).count();

    std::cerr << "SIMPLIFY_CORE_MS: " << std::fixed << std::setprecision(4) << ms << '\n';
    if (time_flag) {
        print_timing_summary();
        print_timing_machine();
        timer_detail::enabled() = false;
    }

    return simplified;
}

// ===========================================================================
//  main
// ===========================================================================

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    
    int test_case_no = -1;
    int code = get_repo_root(argv, repo_root);
    if (code != 0) {
        return code;
    }

    code = parse_arguments(argc, argv, test_case_no);
    if (code != 0) {
        return code;
    }

    std::vector<Point> stream;
    if (help_flag) {
        return 0;
    }
    code = read_stream(test_case_no, argv, stream);
    if (code != 0) {
        return code;
    }

    std::vector<Point> simplified = web_server_flag
        ? simplify_web(stream, EPSILON, DELTA)
        : simplify(stream, EPSILON, DELTA);
    stream = std::move(simplified);

    if (out_flag) {
        out_stream(test_case_no, argv, stream);
    }

    if (!web_server_flag) maybe_run_frechet(test_case_no);

    return 0;
}
