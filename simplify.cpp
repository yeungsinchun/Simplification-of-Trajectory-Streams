#include <chrono>
#include <format>
#include <iostream>
#if __has_include(<print>)
#include <print>
#endif

#include "simplify_core.h"
#include "simplify_io.h"
#include "timer.h"
#include "web_trace.h"

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
        simplified = Simplifier(EPSILON, DELTA).simplify(stream);
    }
    double ms = std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t0).count();

#if __has_include(<print>)
    std::println(stderr, "SIMPLIFY_CORE_MS: {:.4f}", ms);
#else
    std::cerr << std::format("SIMPLIFY_CORE_MS: {:.4f}\n", ms);
#endif
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
