#include <chrono>

#include "simplify_geometry.h"
#include "simplify_io.h"

// ===========================================================================
//  Core algorithm — bare (no instrumentation)
// ===========================================================================
//
// This is the production front-end: no timing hooks in the hot path.  The
// instrumented twin lives in simplify_with_time.cpp and is built as the separate
// `simplify_time` executable.

int get_longest_stab(const std::vector<Point>& stream, int cur,
                     std::vector<Point>& simplified,
                     double EPSILON, double DELTA) {
    const Point& p0 = stream[cur];
    std::vector<Point> P = get_points_from_grid(p0, EPSILON, DELTA);
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
        Gi = get_conv_from_grid(stream[cur], EPSILON, DELTA);
        for (int i = 0; i < Pn; ++i) {
            if (dead[i]) continue;
            find_F(P[i], S[i], F[i]);
            if (!intersect(F[i], Gi, new_S[i])) {
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
        for (int i = 0; i < Pn; ++i)
            if (!dead[i]) S[i].swap(new_S[i]);
        cur++;
    }
    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    return cur;
}

std::vector<Point> simplify(const std::vector<Point>& stream,
                            double EPSILON, double DELTA) {
    std::vector<Point> simplified;
    std::cout << "Running Simplification (EPSILON: " << EPSILON
              << ", DELTA: " << DELTA << ")...\n";
    configure_bbox(stream, EPSILON, DELTA);
    auto t0 = std::chrono::high_resolution_clock::now();
    int cur = 0;
    while (cur != int(stream.size()))
        cur = get_longest_stab(stream, cur, simplified, EPSILON, DELTA);
    double ms = std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t0).count();
    std::cout << "SIMPLIFY_CORE_MS: " << std::fixed << std::setprecision(4) << ms << '\n';
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

    std::vector<Point> stream;
    code = read_stream(test_case_no, argv, stream);
    if (code != 0) {
        return code;
    }

    std::vector<Point> simplified = simplify(stream, EPSILON, DELTA);
    stream = std::move(simplified);

    if (out_flag) {
        out_stream(test_case_no, argv, stream);
    }

    maybe_run_frechet(test_case_no);

    return 0;
}
