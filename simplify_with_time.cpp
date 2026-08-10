#include <chrono>

#include "simplify_geometry.h"
#include "simplify_io.h"
#include "timer.h"

// ===========================================================================
//  Core algorithm — timed (TIMER instrumentation + --time breakdown)
// ===========================================================================
//
// This is the instrumented twin of simplify.cpp.  It is built as the separate
// `simplify_time` executable so the production `simplify` binary keeps a hot
// path with zero timing hooks.  Timing is always enabled here; the breakdown
// prints to stderr at exit.

int get_longest_stab_time(const std::vector<Point>& stream, int cur,
                          std::vector<Point>& simplified,
                          double EPSILON, double DELTA) {
    const Point& p0 = stream[cur];
    std::vector<Point> P;
    { TIMER("get_points_from_grid"); P = get_points_from_grid(p0, EPSILON, DELTA); }
    std::array<Point, 2> buffer = {p0, p0};
    const int Pn = (int)P.size();
    std::vector<std::vector<Point>> S(Pn);
    for (int i = 0; i < Pn; ++i) S[i] = {P[i]};
    int dead_cnt = 0;
    std::vector<int> dead(Pn);
    std::vector<std::vector<Point>> new_S(Pn);
    std::vector<std::vector<Point>> F(Pn);
    std::vector<Point> Gi;

    cur++;
    while (cur < int(stream.size())) {
        TIMER("stab.step");
        { TIMER("get_conv_from_grid"); Gi = get_conv_from_grid(stream[cur], EPSILON, DELTA); }
        {
            TIMER("stab.clip_loop");
            for (int i = 0; i < Pn; ++i) {
                if (dead[i]) continue;
                { TIMER("find_F"); find_F(P[i], S[i], F[i]); }
                bool hit;
                { TIMER("intersect"); hit = intersect(F[i], Gi, new_S[i]); }
                if (!hit) { dead[i] = true; dead_cnt++; }
            }
        }
        bool has_candidate = false;
        {
            TIMER("stab.scan");
            for (int i = Pn - 1; i >= 0 && !has_candidate; --i) {
                if (dead[i] || new_S[i].empty()) continue;
                buffer[0] = P[i];
                buffer[1] = new_S[i].front();
                has_candidate = true;
            }
        }
        if (!has_candidate || dead_cnt == Pn) break;
        {
            TIMER("stab.commit");
            for (int i = 0; i < Pn; ++i)
                if (!dead[i]) S[i].swap(new_S[i]);
        }
        cur++;
    }
    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    return cur;
}

std::vector<Point> simplify_time(const std::vector<Point>& stream,
                                 double EPSILON, double DELTA) {
    std::vector<Point> simplified;
    std::cout << "Running Simplification (EPSILON: " << EPSILON
              << ", DELTA: " << DELTA << ")...\n";
    {
        TIMER("total");
        configure_bbox(stream, EPSILON, DELTA);
        auto t0 = std::chrono::high_resolution_clock::now();
        int cur = 0;
        while (cur != int(stream.size())) {
            TIMER("prefix");
            cur = get_longest_stab_time(stream, cur, simplified, EPSILON, DELTA);
        }
        double ms = std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - t0).count();
        std::cout << "SIMPLIFY_CORE_MS: " << std::fixed << std::setprecision(4) << ms << '\n';
    }
    std::cout << "Simplified to " << simplified.size() << " points ("
              << (100.0 * simplified.size() / stream.size()) << "%).\n";
    return simplified;
}

// ===========================================================================
//  main
// ===========================================================================

int main(int argc, char** argv) {
    int test_case_no = -1;
    int code = get_repo_root(argv, repo_root);
    if (code != 0) {
        return code;
    }

    code = parse_arguments(argc, argv, test_case_no);
    if (code != 0) {
        return code;
    }

    // simplify_time always produces the breakdown; --time is accepted for
    // compatibility but the instrumentation is on regardless.
    timer_detail::enabled() = true;

    std::vector<Point> stream;
    code = read_stream(test_case_no, argv, stream);
    if (code != 0) {
        return code;
    }

    std::vector<Point> simplified = simplify_time(stream, EPSILON, DELTA);
    print_timing_summary();
    stream = std::move(simplified);

    if (out_flag) {
        out_stream(test_case_no, argv, stream);
    }

    maybe_run_frechet(test_case_no);

    return 0;
}
