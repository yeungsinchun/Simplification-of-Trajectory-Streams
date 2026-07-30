#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "simplify_geometry.h"
#include "timer.h"

// ===========================================================================
//  Global parameters
// ===========================================================================

double DELTA = 200;
double EPSILON = 0.5;

std::filesystem::path repo_root;
bool out_flag = false;
bool dist_flag = false;

// ===========================================================================
//  Help
// ===========================================================================

static void print_help() {
    std::cout << "Usage: simplify [options]\n"
              << "  --in <id>        Read input from data/taxi/<id>.txt (resolved absolutely)\n"
              << "  --out            Write output to data/<id>/original.txt & simplify.txt (resolved absolutely; requires --in <id>)\n"
              << "  --dist           After output, compute Frechet distance by invoking ./frechet (Julia wrapper) with --in <id> --path <simplify.txt>\n"
              << "  -d <delta>       Override DELTA (default " << DELTA << ")\n"
              << "  -e <epsilon>     Override EPSILON (default " << EPSILON << ")\n"
              << "  --dump-intersect Dump every (F_poly, Gi_poly) pair fed to "
                 "intersect() to data/<id>/intersect_pairs.txt\n"
              << "  -h               Show this help and exit\n"
              << "\n"
              << "Shorthand: simplify <id> [flags] is equivalent to '--in <id> --out [flags]'\n";
}

int parse_arguments(int argc, char** argv, int& test_case_no) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i],"--out") == 0) out_flag = true;
        else if (strcmp(argv[i],"--dist") == 0) dist_flag = true;
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
        else if (strcmp(argv[i],"-h") == 0) { print_help(); return 0; }
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
        print_help();
        return 0;
    }
    return 0;
}

int get_repo_root(char** argv, std::filesystem::path& repo_root) {
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

int read_stream(int test_case_no, char** argv, std::vector<Point>& stream) {
    if (test_case_no != -1) {
        auto simp_orig = repo_root / "data" / std::to_string(test_case_no) / "original.txt";
        std::cout << simp_orig.string() << '\n';
        std::ifstream fin(simp_orig.string());
        if (!fin) { std::cerr << "Cannot open " << simp_orig.string() << "\n"; return 1; }
        int N = 0;
        if (!(fin >> N)) { std::cerr << "Empty or invalid input in " << simp_orig.string() << "\n"; return 1; }
        stream.clear(); stream.reserve(N);
        for (int i = 0; i < N; ++i) {
            double x,y; if (!(fin >> x >> y)) { std::cerr << "Malformed pair at index " << i << " in " << simp_orig.string() << "\n"; return 1; }
            stream.emplace_back(x, y);
        }
    }
    return 0;
}

int out_stream(int test_case_no, char** argv, const std::vector<Point>& stream) {
    std::filesystem::path dir = repo_root / "data" / std::to_string(test_case_no);
    std::filesystem::create_directories(dir);

    std::ofstream simp(dir / "simplify.txt");
    simp << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::size_t N = stream.size();
    simp << N << '\n';
    for (const auto& p : stream) {
        simp << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()) << '\n';
    }
    simp.close();
    std::cout << "Output Written\n";
    return 0;
}

// ===========================================================================
//  Core algorithm
// ===========================================================================

int get_longest_stab(const std::vector<Point>& stream, int cur,
                     std::vector<Point>& simplified,
                     double EPSILON, double DELTA) {
    SIGNPOST_EVENT(stab.start, "cur=%d", cur);
    const Point& p0 = stream[cur];
    std::vector<Point> P;
    {
        SIGNPOST_BEGIN(get_points_from_grid);
        P = get_points_from_grid(p0, EPSILON, DELTA);
        SIGNPOST_END(get_points_from_grid);
        SIGNPOST_EVENT(stab.grid_points, "Pn=%lu", (unsigned long)P.size());
    }
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
        const Point& pi = stream[cur];
        SIGNPOST_BEGIN(stab.step);
        SIGNPOST_EVENT(stab.consume, "cur=%d", cur);
        Gi = get_conv_from_grid(pi, EPSILON, DELTA);

        for (int i = 0; i < Pn; ++i) {
            if (dead[i]) continue;
            {
                SIGNPOST_BEGIN(find_F);
                find_F(P[i], S[i], F[i]);
                SIGNPOST_END(find_F);
            }
            bool hit;
            {
                SIGNPOST_BEGIN(intersect);
                hit = intersect(F[i], Gi, new_S[i]);
                SIGNPOST_END(intersect);
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
        SIGNPOST_END(stab.step);
        if (!has_candidate || dead_cnt == Pn) break;

        for (int i = 0; i < Pn; ++i) {
            if (!dead[i]) S[i].swap(new_S[i]);
        }
        cur++;
    }
    SIGNPOST_EVENT(stab.emit, "buffer=(%{public}lf,%{public}lf)->(%{public}lf,%{public}lf)",
                   (double)CGAL::to_double(buffer[0].x()), (double)CGAL::to_double(buffer[0].y()),
                   (double)CGAL::to_double(buffer[1].x()), (double)CGAL::to_double(buffer[1].y()));

    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    return cur;
}

std::vector<Point> simplify(const std::vector<Point>& stream,
                            double EPSILON, double DELTA) {
    SIGNPOST_BEGIN(simplify.run);
    configure_bbox(stream, EPSILON, DELTA);
    std::vector<Point> simplified;
    std::cout << "Simplifying...\n";
    
    // Start core algorithm timing (excluding configure_bbox setup)
    auto _core_start = std::chrono::high_resolution_clock::now();
    
    int cur = 0;
    int prefix = 0;
    while (cur != int(stream.size())) {
        SIGNPOST_EVENT(simplify.prefix, "from=%d size=%lu", cur, (unsigned long)simplified.size());
        cur = get_longest_stab(stream, cur, simplified, EPSILON, DELTA);
        prefix++;
    }
    
    auto _core_end = std::chrono::high_resolution_clock::now();
    double _core_ms = std::chrono::duration<double, std::milli>(_core_end - _core_start).count();
    fprintf(stderr, "SIMPLIFY_CORE_MS %.4f\n", _core_ms);
    
    SIGNPOST_END(simplify.run);
    std::cout << "Expected Frechet distance: " << std::sqrt(expected_frechet_squared) << '\n';
    std::cout << "The original stream of size " << stream.size() << " is simplified to " << simplified.size() << " points.\n";
    SIGNPOST_EVENT(simplify.done, "n=%lu m=%lu prefixes=%d",
                   (unsigned long)stream.size(), (unsigned long)simplified.size(), prefix);
    return simplified;
}

// ===========================================================================
//  main
// ===========================================================================

int main(int argc, char** argv) {
    int test_case_no = -1;
    int code = get_repo_root(argv, repo_root);
    std::cout << repo_root.string() << std::endl;
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

    SIGNPOST_BEGIN(app.run);
    std::vector<Point> simplified = simplify(stream, EPSILON, DELTA);
    SIGNPOST_END(app.run);
    stream = std::move(simplified);

    if (out_flag) {
        out_stream(test_case_no, argv, stream);
    }

    if (dist_flag && test_case_no != -1) {
        std::filesystem::path frechet_path = repo_root / "scripts" / "frechet";
        std::string cmd1 = std::string("\"") + frechet_path.string() + "\" --in " + std::to_string(test_case_no) + " --path \"" + (repo_root / "data" / std::to_string(test_case_no) / "simplify.txt").string() + "\"";
        int rc = std::system(cmd1.c_str());
    }

    return 0;
}
