#ifndef SIMPLIFY_IO_H
#define SIMPLIFY_IO_H

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#if __has_include(<print>)
#include <print>
#endif

#include "simplify_geometry.h"

// ===========================================================================
//  Global parameters
// ===========================================================================

inline double DELTA = 200;
inline double EPSILON = 0.5;

inline std::filesystem::path repo_root;
inline bool out_flag = false;
inline bool dist_flag = false;
inline bool web_server_flag = false;
inline bool json_stream_flag = false;
inline bool help_flag = false;
inline bool time_flag = false;
inline std::string json_output_path = "";

// ===========================================================================
//  Help
// ===========================================================================

inline void print_help(const char* prog) {
#if __has_include(<print>)
    std::println("Usage: {} [options]", prog);
    std::println("  --in <id>        Read input from data/taxi/<id>.txt "
                 "(resolved absolutely)");
    std::println("  --out            Write output to data/<id>/original.txt & "
                 "simplify.txt (resolved absolutely; requires --in <id>)");
    std::println(
        "  --dist           After output, compute Frechet distance by invoking "
        "'julia scripts/frechet.jl' with --in <id> --path <simplify.txt>");
    std::println("  -d <delta>       Override DELTA (default {})", DELTA);
    std::println("  -e <epsilon>     Override EPSILON (default {})", EPSILON);
    std::println("  --dump-intersect Dump every (F_poly, Gi_poly) pair fed to "
                 "intersect() to data/<id>/intersect_pairs.txt");
    std::println("  --web-server     Emit a machine-readable JSON trace of the "
                 "algorithm to stdout for the web visualizer (suppresses all "
                 "other stdout text)");
    std::println("  --json-stream    With --web-server, emit NDJSON (header, "
                 "one prefix per line, done)");
    std::println("  --json-output <path>  Write JSON trace to file instead of "
                 "stdout (use with --web-server)");
    std::println("  --time           Opt-in phase timers (stderr TIMING "
                 "SUMMARY + TIMER_MS lines)");
    std::println("  -h               Show this help and exit");
    std::println("");
    std::println(
        "Shorthand: {} <id> [flags] is equivalent to '--in <id> --out [flags]'",
        prog);
#else
    std::cout << std::format("Usage: {} [options]\n", prog)
              << "  --in <id>        Read input from data/taxi/<id>.txt "
                 "(resolved absolutely)\n"
              << "  --out            Write output to data/<id>/original.txt & "
                 "simplify.txt (resolved absolutely; requires --in <id>)\n"
              << "  --dist           After output, compute Frechet distance by "
                 "invoking 'julia scripts/frechet.jl' with --in <id> --path "
                 "<simplify.txt>\n"
              << std::format("  -d <delta>       Override DELTA (default {})\n",
                             DELTA)
              << std::format(
                     "  -e <epsilon>     Override EPSILON (default {})\n",
                     EPSILON)
              << "  --dump-intersect Dump every (F_poly, Gi_poly) pair fed to "
                 "intersect() to data/<id>/intersect_pairs.txt\n"
              << "  --web-server     Emit a machine-readable JSON trace of the "
                 "algorithm to stdout for the web visualizer (suppresses all "
                 "other stdout text)\n"
              << "  --json-stream    With --web-server, emit NDJSON (header, "
                 "one prefix per line, done)\n"
              << "  --json-output <path>  Write JSON trace to file instead of "
                 "stdout (use with --web-server)\n"
              << "  --time           Opt-in phase timers (stderr TIMING "
                 "SUMMARY + TIMER_MS lines)\n"
              << "  -h               Show this help and exit\n"
              << "\n"
              << std::format("Shorthand: {} <id> [flags] is equivalent to "
                             "'--in <id> --out [flags]'\n",
                             prog);
#endif
}

inline int parse_arguments(int argc, char** argv, int& test_case_no) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i],"--out") == 0) out_flag = true;
        else if (strcmp(argv[i],"--dist") == 0) dist_flag = true;
        else if (strcmp(argv[i],"--web-server") == 0) web_server_flag = true;
        else if (strcmp(argv[i],"--json-stream") == 0) json_stream_flag = true;
        else if (strcmp(argv[i],"--json-output") == 0 && i+1 < argc) {
            json_output_path = argv[++i];
        }
        else if (strcmp(argv[i],"--time") == 0) time_flag = true;
        else if (strcmp(argv[i],"--gui") == 0 || strcmp(argv[i],"-F") == 0 ||
                 strcmp(argv[i],"-G") == 0 || strcmp(argv[i],"-S") == 0) {
#if __has_include(<print>)
            std::println(stderr, "GUI options require simplify_with_gui");
#else
            std::cerr << std::format("GUI options require simplify_with_gui\n");
#endif
            return 1;
        }
        else if (strcmp(argv[i],"-d") == 0 && i+1 < argc) {
            try {
                DELTA = std::stod(argv[++i]);
            } catch (...) {
#if __has_include(<print>)
                std::println(stderr, "Invalid -d value");
#else
                std::cerr << std::format("Invalid -d value\n");
#endif
                return 1;
            }
        }
        else if (strcmp(argv[i],"-e") == 0 && i+1 < argc) {
            try {
                EPSILON = std::stod(argv[++i]);
            } catch (...) {
#if __has_include(<print>)
                std::println(stderr, "Invalid -e value");
#else
                std::cerr << std::format("Invalid -e value\n");
#endif
                return 1;
            }
        }
        else if (strcmp(argv[i],"-h") == 0) { print_help(argv[0]); help_flag = true; return 0; }
        else if (strcmp(argv[i],"--in") == 0 && i+1 < argc) {
            try {
                test_case_no = std::stoi(argv[++i]);
            } catch (...) {
#if __has_include(<print>)
                std::println(stderr, "Invalid --in argument");
#else
                std::cerr << std::format("Invalid --in argument\n");
#endif
                return 1;
            }
        }
    }

    if (dist_flag) out_flag = true;

    if (test_case_no == -1 && argc >= 2 && argv[1][0] != '-') {
        try {
            test_case_no = std::stoi(argv[1]);
            out_flag = true;
        } catch (...) {
#if __has_include(<print>)
            std::println("Command parse error");
#else
            std::cout << std::format("Command parse error\n");
#endif
            return 1;
        }
    }

    if (argc == 1) {
        print_help(argv[0]);
        help_flag = true;
        return 0;
    }
    return 0;
}

inline int get_repo_root(char** argv, std::filesystem::path& repo_root) {
    auto find_repo_root = [](const std::filesystem::path& start, int max_levels) -> std::filesystem::path {
        auto dir = std::filesystem::weakly_canonical(start);
        for (int i = 0; i < max_levels && !dir.empty(); ++i) {
            if (std::filesystem::is_directory(dir / "data")) return dir;
            const auto parent = dir.parent_path();
            if (parent == dir) break;
            dir = parent;
        }
        return {};
    };
    try {
        // search with reference to the path passed by argv[0]
        repo_root = find_repo_root(argv[0], 5);
    } catch (const std::filesystem::filesystem_error& e) {
        // search with reference to the shell location
        try {
            repo_root = find_repo_root(std::filesystem::current_path(), 5);
        } catch (const std::filesystem::filesystem_error& fallback_error) {
#if __has_include(<print>)
            std::println(stderr,
                         "Error: could not resolve the data directory: {}",
                         fallback_error.what());
#else
            std::cerr << std::format(
                "Error: could not resolve the data directory: {}\n",
                fallback_error.what());
#endif
            return 1;
        }
    }
    return 0;
}

inline int read_stream(int test_case_no, char** argv, std::vector<Point>& stream) {
    if (test_case_no != -1) {
        auto simp_orig = repo_root / "data" / std::to_string(test_case_no) / "original.txt";
        [[maybe_unused]] auto simp_output =
            repo_root / "data" / std::to_string(test_case_no) / "simplify.txt";
        if (!web_server_flag) {
#if __has_include(<print>)
            std::println("Input file: {}", simp_orig.string());
            if (out_flag)
                std::println("Output file: {}", simp_output.string());
#else
            std::cout << std::format("Input file: {}\n", simp_orig.string());
            if (out_flag)
                std::cout << std::format("Output file: {}\n",
                                         simp_output.string());
#endif
        }
        std::ifstream fin(simp_orig.string());
        if (!fin) {
#if __has_include(<print>)
            std::println(stderr, "Cannot open {}", simp_orig.string());
#else
            std::cerr << std::format("Cannot open {}\n", simp_orig.string());
#endif
            return 1;
        }

        // Optimize IO performance
        fin.sync_with_stdio(false);
        fin.tie(nullptr);

        int N = 0;
        if (!(fin >> N)) {
#if __has_include(<print>)
            std::println(stderr, "Empty or invalid input in {}",
                         simp_orig.string());
#else
            std::cerr << std::format("Empty or invalid input in {}\n",
                                     simp_orig.string());
#endif
            return 1;
        }
        stream.clear(); stream.reserve(N);

        // Batch read points
        double x, y;
        for (int i = 0; i < N; ++i) {
            if (!(fin >> x >> y)) {
#if __has_include(<print>)
                std::println(stderr, "Malformed pair at index {} in {}", i,
                             simp_orig.string());
#else
                std::cerr << std::format("Malformed pair at index {} in {}\n",
                                         i, simp_orig.string());
#endif
                return 1;
            }
            stream.emplace_back(x, y);
        }
        if (!web_server_flag) {
#if __has_include(<print>)
            std::println("Loaded {} points.", stream.size());
#else
            std::cout << std::format("Loaded {} points.\n", stream.size());
#endif
        }
    }
    return 0;
}

inline int out_stream(int test_case_no, char** argv, const std::vector<Point>& stream) {
    std::filesystem::path dir = repo_root / "data" / std::to_string(test_case_no);
    std::filesystem::create_directories(dir);
    std::ofstream simp(dir / "simplify.txt");

    // Optimize IO performance
    simp.sync_with_stdio(false);
    simp.tie(nullptr);

    simp << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::size_t N = stream.size();
    simp << N << '\n';
    for (const auto& p : stream) {
        simp << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()) << '\n';
    }
    simp.close();
    return 0;
}

// The Frechet-distance post-step shared by both front-ends.
inline void maybe_run_frechet(int test_case_no) {
    if (dist_flag && test_case_no != -1) {
        if (!web_server_flag) {
#if __has_include(<print>)
            std::println("Calculating Frechet distance...");
#else
            std::cout << std::format("Calculating Frechet distance...\n")
                      << std::flush;
#endif
            std::fflush(stdout);
        }
        std::filesystem::path frechet_path = repo_root / "scripts" / "frechet.jl";
        std::string cmd1 = std::string("julia \"") + frechet_path.string() + "\" --in " + std::to_string(test_case_no) + " --path \"" + (repo_root / "data" / std::to_string(test_case_no) / "simplify.txt").string() + "\" --raw";
        int rc = std::system(cmd1.c_str());
        (void)rc;
    }
}

#endif // SIMPLIFY_IO_H
