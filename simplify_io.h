#ifndef SIMPLIFY_IO_H
#define SIMPLIFY_IO_H

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

#if defined(__linux__)
#include <sched.h>
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
// Default threads for the per-step anchor updates: the usable CPUs, capped
// at 8. On Linux, containers (Cloud Run, docker --cpus) and taskset limit CPUs
// through the affinity mask or a cgroup quota that hardware_concurrency()
// does not see, and oversubscribed helpers would only slow the run down.
inline int default_simplify_threads() {
    long n = long(std::thread::hardware_concurrency());
#if defined(__linux__)
    cpu_set_t mask;
    if (sched_getaffinity(0, sizeof(mask), &mask) == 0)
        n = std::min<long>(n, CPU_COUNT(&mask));
    auto limit_by_quota = [&n](long quota, long period) {
        if (quota > 0 && period > 0)
            n = std::min(n, std::max(1L, (quota + period - 1) / period));
    };
    std::string quota;
    long period = 0;
    if (std::ifstream v2("/sys/fs/cgroup/cpu.max"); v2 >> quota >> period) {
        if (quota != "max")
            limit_by_quota(std::atol(quota.c_str()), period);
    } else {
        long v1_quota = 0;
        std::ifstream q("/sys/fs/cgroup/cpu/cpu.cfs_quota_us");
        std::ifstream p("/sys/fs/cgroup/cpu/cpu.cfs_period_us");
        if (q >> v1_quota && p >> period)
            limit_by_quota(v1_quota, period);
    }
#endif
    return int(std::clamp(n, 1L, 8L));
}

// Threads for the per-step anchor updates (--threads). Output is identical
// for every value; 0 on the command line restores the default.
inline int simplify_threads = default_simplify_threads();
inline std::string json_output_path = "";

// ===========================================================================
//  Help
// ===========================================================================

inline void print_help(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "  --in <id>        Read input from data/taxi/<id>.txt (resolved absolutely)\n"
              << "  --out            Write output to data/<id>/original.txt & simplify.txt (resolved absolutely; requires --in <id>)\n"
              << "  --dist           After output, compute Frechet distance by invoking 'julia scripts/frechet.jl' with --in <id> --path <simplify.txt>" << '\n'
              << "  -d <delta>       Override DELTA (default " << DELTA << ")\n"
              << "  -e <epsilon>     Override EPSILON (default " << EPSILON << ")\n"
              << "  --dump-intersect Dump every (F_poly, Gi_poly) pair fed to "
                 "intersect() to data/<id>/intersect_pairs.txt\n"
              << "  --web-server     Emit a machine-readable JSON trace of the algorithm to "
                 "stdout for the web visualizer (suppresses all other stdout text)\n"
              << "  --json-stream    With --web-server, emit NDJSON (header, one prefix per line, done)\n"
              << "  --json-output <path>  Write JSON trace to file instead of stdout (use with --web-server)\n"
              << "  --time           Opt-in phase timers (stderr TIMING SUMMARY + TIMER_MS lines; runs single-threaded)\n"
              << "  --threads <n>    Threads for per-step anchor updates (default " << simplify_threads << "; output is identical for any n)\n"
              << "  -h               Show this help and exit\n"
              << "\n"
              << "Shorthand: " << prog << " <id> [flags] is equivalent to '--in <id> --out [flags]'\n";
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
        else if (strcmp(argv[i],"--threads") == 0 && i+1 < argc) {
            try { simplify_threads = std::stoi(argv[++i]); } catch(...) { std::cerr << "Invalid --threads value\n"; return 1; }
            if (simplify_threads < 0) { std::cerr << "Invalid --threads value\n"; return 1; }
            if (simplify_threads == 0)
                simplify_threads = default_simplify_threads();
        }
        else if (strcmp(argv[i],"--gui") == 0 || strcmp(argv[i],"-F") == 0 ||
                 strcmp(argv[i],"-G") == 0 || strcmp(argv[i],"-S") == 0) {
            std::cerr << "GUI options require simplify_with_gui\n";
            return 1;
        }
        else if (strcmp(argv[i],"-d") == 0 && i+1 < argc) {
            try { DELTA = std::stod(argv[++i]); } catch(...) { std::cerr << "Invalid -d value\n"; return 1; }
        }
        else if (strcmp(argv[i],"-e") == 0 && i+1 < argc) {
            try { EPSILON = std::stod(argv[++i]); } catch(...) { std::cerr << "Invalid -e value\n"; return 1; }
        }
        else if (strcmp(argv[i],"-h") == 0) { print_help(argv[0]); help_flag = true; return 0; }
        else if (strcmp(argv[i],"--in") == 0 && i+1 < argc) {
            try { test_case_no = std::stoi(argv[++i]); }
            catch(...) { std::cerr << "Invalid --in argument\n"; return 1; }
        }
    }

    if (dist_flag) out_flag = true;

    if (test_case_no == -1 && argc >= 2 && argv[1][0] != '-') {
        try {
            test_case_no = std::stoi(argv[1]);
            out_flag = true;
        } catch (...) {
            std::cout << "Command parse error\n";
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
            std::cerr << "Error: could not resolve the data directory: " << fallback_error.what() << "\n";
            return 1;
        }
    }
    return 0;
}

inline int read_stream(int test_case_no, char** argv, std::vector<Point>& stream) {
    if (test_case_no != -1) {
        auto simp_orig = repo_root / "data" / std::to_string(test_case_no) / "original.txt";
        auto simp_output = repo_root / "data" / std::to_string(test_case_no) / "simplify.txt";
        if (!web_server_flag) {
            std::cout << "Input file: " << simp_orig.string() << '\n';
            if (out_flag) std::cout << "Output file: " << simp_output.string() << '\n';
        }
        std::ifstream fin(simp_orig.string());
        if (!fin) { std::cerr << "Cannot open " << simp_orig.string() << "\n"; return 1; }
        
        // Optimize IO performance
        fin.sync_with_stdio(false);
        fin.tie(nullptr);
        
        int N = 0;
        if (!(fin >> N)) { std::cerr << "Empty or invalid input in " << simp_orig.string() << "\n"; return 1; }
        stream.clear(); stream.reserve(N);
        
        // Batch read points
        double x, y;
        for (int i = 0; i < N; ++i) {
            if (!(fin >> x >> y)) { std::cerr << "Malformed pair at index " << i << " in " << simp_orig.string() << "\n"; return 1; }
            stream.emplace_back(x, y);
        }
        if (!web_server_flag) std::cout << "Loaded " << stream.size() << " points." << "\n";
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
        if (!web_server_flag) std::cout << "Calculating Frechet distance...\n" << std::flush;
        std::filesystem::path frechet_path = repo_root / "scripts" / "frechet.jl";
        std::string cmd1 = std::string("julia \"") + frechet_path.string() + "\" --in " + std::to_string(test_case_no) + " --path \"" + (repo_root / "data" / std::to_string(test_case_no) / "simplify.txt").string() + "\" --raw";
        int rc = std::system(cmd1.c_str());
        (void)rc;
    }
}

#endif // SIMPLIFY_IO_H
