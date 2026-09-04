#include <chrono>
#include <fstream>
#include <iomanip>

#include "simplify_geometry.h"
#include "simplify_io.h"
#include "timer.h"

// ===========================================================================
//  Core algorithm
// ===========================================================================


int get_longest_stab(const std::vector<Point>& stream, int cur,
                     std::vector<Point>& simplified,
                     double EPSILON, double DELTA) {
    TIMER("get_longest_stab");
    const Point& p0 = stream[cur];
    std::vector<Point> P;
    {
        TIMER("get_boundary_points");
        P = get_boundary_points_from_grid(p0, EPSILON, DELTA);
    }
    std::array<Point, 2> buffer = {p0, p0};
    const int Pn = (int)P.size();
    if (timer_detail::enabled())
        timer_detail::counters()["boundary_candidates"] += Pn;
    std::vector<std::vector<Point>> S(Pn);
    for (int i = 0; i < Pn; ++i) S[i] = {P[i]};
    int dead_cnt = 0;
    std::vector<int> dead(Pn);
    std::vector<std::vector<Point>> new_S(Pn);
    std::vector<std::vector<Point>> F(Pn);
    std::vector<Point> Gi;
    std::vector<std::array<double, 2>> Gi_xy;
    std::vector<int> tangents;
    tangents.reserve(2);

    cur++;
    while (cur < int(stream.size())) {
        if (timer_detail::enabled())
            timer_detail::counters()["stab_steps"]++;
        {
            TIMER("get_conv_from_grid");
            Gi = get_conv_from_grid(stream[cur], EPSILON, DELTA);
            sh_double::prepare_convex_xy(Gi, Gi_xy);
        }
        for (int i = 0; i < Pn; ++i) {
            if (dead[i]) continue;
            if (timer_detail::enabled())
                timer_detail::counters()["alive_candidate_iters"]++;

            // Conservative Frechet-safe prune: if Gi is robustly outside the
            // tangent cone of S at P[i], F n Gi is empty and the candidate dies
            // without building F or running the clip.  On a miss, reuse the
            // tangents already computed for find_F.
            tangents.clear();
            if (wedge_gi_disjoint(P[i], S[i], Gi, &tangents)) {
                if (timer_detail::enabled())
                    timer_detail::counters()["wedge_prune_hits"]++;
                dead[i] = true;
                dead_cnt++;
                continue;
            }

            {
                TIMER("find_F");
                const std::vector<int>* tan_ptr =
                    (tangents.size() == 2) ? &tangents : nullptr;
                find_F(P[i], S[i], F[i], tan_ptr);
            }

            bool hit;
            {
                TIMER("intersect");
                hit = intersect_prepared(F[i], Gi_xy, new_S[i]);
            }
            if (!hit) {
                dead[i] = true;
                dead_cnt++;
            }
        }
        bool has_candidate = false;
        for (int i = Pn - 1; i >= 0 && !has_candidate; --i) {
            if (dead[i] || new_S[i].empty()) continue;
            buffer[0] = P[i];
            buffer[1] = new_S[i].front();
            has_candidate = true;
        }
        if (!has_candidate || dead_cnt == Pn) break;
        {
            TIMER("update_S");
            for (int i = 0; i < Pn; ++i)
                if (!dead[i]) S[i].swap(new_S[i]);
        }
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
    std::vector<std::array<double, 2>> Gi_xy;
    std::vector<int> tangents;
    tangents.reserve(2);

    webtrace::PrefixTrace trace;
    trace.p0 = p0;
    trace.p0_idx = cur;   // p0 = stream[cur] (cur not yet incremented)
    trace.P = P;

    cur++;
    while (cur < int(stream.size())) {
        Gi = get_conv_from_grid(stream[cur], EPSILON, DELTA);
        sh_double::prepare_convex_xy(Gi, Gi_xy);

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

            tangents.clear();
            if (wedge_gi_disjoint(P[i], S[i], Gi, &tangents)) {
                dead[i] = true;
                dead_cnt++;
                cand.alive = false;
                step.candidates.push_back(std::move(cand));
                continue;
            }

            find_F(P[i], S[i], F[i],
                   (tangents.size() == 2) ? &tangents : nullptr);
            cand.F = F[i];

            if (!intersect_prepared(F[i], Gi_xy, new_S[i])) {
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
    double ms = 0.0;
    {
        TIMER("total");
        configure_bbox(stream, EPSILON, DELTA);
        auto t0 = std::chrono::high_resolution_clock::now();
        {
            TIMER("simplify");
            int cur = 0;
            while (cur != int(stream.size()))
                cur = get_longest_stab(stream, cur, simplified, EPSILON, DELTA);
        }
        ms = std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - t0).count();
    }

    std::cerr << "SIMPLIFY_CORE_MS: " << std::fixed << std::setprecision(4) << ms << '\n';
    if (time_flag) {
        // Extra counters that are not TIMER scopes: emit as KEY: value lines
        // so the benchmark harness can scrape them without parsing the tree.
        std::cerr << "BOUNDARY_CANDIDATES_SUM: "
                  << timer_detail::counters()["boundary_candidates"] << '\n';
        std::cerr << "STAB_STEPS: "
                  << timer_detail::counters()["stab_steps"] << '\n';
        std::cerr << "ALIVE_CANDIDATE_ITERS: "
                  << timer_detail::counters()["alive_candidate_iters"] << '\n';
        std::cerr << "WEDGE_PRUNE_HITS: "
                  << timer_detail::counters()["wedge_prune_hits"] << '\n';
        std::cerr << "SIMPLIFIED_POINTS: " << simplified.size() << '\n';
        print_timing_summary();
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
    if (time_flag) timer_detail::enabled() = true;

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
