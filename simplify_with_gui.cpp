#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <QApplication>
#include "drawing.h"
#include "simplify_geometry.h"

// ===========================================================================
//  Global parameters
// ===========================================================================

double DELTA = 200;
double EPSILON = 0.5;

std::filesystem::path repo_root;
bool out_flag = false;
bool dist_flag = false;
bool showF = false;
bool showG = false;
bool showS = false;
bool keep_polygons = false;
bool show_simplified = true;
bool show_labels = false;
bool show_dots = true;
bool help_flag = false;

// ===========================================================================
//  Help
// ===========================================================================

static void print_help() {
    std::cout << "Usage: simplify_with_gui [options]\n"
              << "  --in <id>        Read input from data/<id>/original.txt (resolved absolutely)\n"
              << "  --out            Write output to data/<id>/simplify.txt (requires --in <id>)\n"
              << "  --dist           After output, compute Frechet distance with --in <id>\n"
              << "  --keep           Do not clear F/G/S polygons between steps\n"
              << "  --no-simp        Do not show the simplified curve in GUI\n"
              << "  --labels         Show vertex labels (indices) in GUI\n"
              << "  --no-dots        Do not overlay data/<id>/dots_simplified.txt when present\n"
              << "  -F/-G/-S         Show F, G, or S debug polygons in GUI\n"
              << "  -d <delta>       Override DELTA (default " << DELTA << ")\n"
              << "  -e <epsilon>     Override EPSILON (default " << EPSILON << ")\n"
              << "  -h               Show this help and exit\n"
              << "\n"
              << "Shorthand: simplify_with_gui <id> [flags] is equivalent to '--in <id> --out [flags]'\n"
              << "If data/<id>/dots_simplified.txt (or DOTS.txt) exists, it is overlaid as a DOTS curve.\n";
}

int parse_arguments(int argc, char** argv, int& test_case_no) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--out") == 0) out_flag = true;
        else if (std::strcmp(argv[i], "--dist") == 0) dist_flag = true;
        else if (std::strcmp(argv[i], "--keep") == 0) keep_polygons = true;
        else if (std::strcmp(argv[i], "--no-simp") == 0) show_simplified = false;
        else if (std::strcmp(argv[i], "--labels") == 0) show_labels = true;
        else if (std::strcmp(argv[i], "--no-dots") == 0) show_dots = false;
        else if (std::strcmp(argv[i], "-F") == 0) showF = true;
        else if (std::strcmp(argv[i], "-G") == 0) showG = true;
        else if (std::strcmp(argv[i], "-S") == 0) showS = true;
        else if (std::strcmp(argv[i], "-d") == 0 && i + 1 < argc) {
            try { DELTA = std::stod(argv[++i]); }
            catch (...) { std::cerr << "Invalid -d value\n"; return 1; }
        }
        else if (std::strcmp(argv[i], "-e") == 0 && i + 1 < argc) {
            try { EPSILON = std::stod(argv[++i]); }
            catch (...) { std::cerr << "Invalid -e value\n"; return 1; }
        }
        else if (std::strcmp(argv[i], "-h") == 0) {
            print_help();
            help_flag = true;
            return 0;
        }
        else if (std::strcmp(argv[i], "--in") == 0 && i + 1 < argc) {
            try { test_case_no = std::stoi(argv[++i]); }
            catch (...) { std::cerr << "Invalid --in argument\n"; return 1; }
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
        help_flag = true;
        return 0;
    }
    return 0;
}

int get_repo_root(char** argv, std::filesystem::path& root) {
    auto find_repo_root = [](const std::filesystem::path& start, int max_levels) {
        auto dir = std::filesystem::weakly_canonical(start);
        for (int i = 0; i < max_levels && !dir.empty(); ++i) {
            if (std::filesystem::is_directory(dir / "data")) return dir;
            const auto parent = dir.parent_path();
            if (parent == dir) break;
            dir = parent;
        }
        return std::filesystem::path{};
    };

    try {
        root = find_repo_root(argv[0], 5);
    } catch (const std::filesystem::filesystem_error&) {
        try {
            root = find_repo_root(std::filesystem::current_path(), 5);
        } catch (const std::filesystem::filesystem_error& error) {
            std::cerr << "Error: could not resolve the data directory: " << error.what() << "\n";
            return 1;
        }
    }
    return 0;
}

int read_stream(int test_case_no, std::vector<Point>& stream) {
    if (test_case_no != -1) {
        auto input_path = repo_root / "data" / std::to_string(test_case_no) / "original.txt";
        std::cout << input_path.string() << '\n';
        std::ifstream fin(input_path);
        if (!fin) {
            std::cerr << "Cannot open " << input_path.string() << "\n";
            return 1;
        }
        int count = 0;
        if (!(fin >> count)) {
            std::cerr << "Empty or invalid input in " << input_path.string() << "\n";
            return 1;
        }
        stream.clear();
        stream.reserve(count);
        for (int i = 0; i < count; ++i) {
            double x, y;
            if (!(fin >> x >> y)) {
                std::cerr << "Malformed pair at index " << i << " in " << input_path.string() << "\n";
                return 1;
            }
            stream.emplace_back(x, y);
        }
    }
    return 0;
}

int out_stream(int test_case_no, const std::vector<Point>& stream) {
    std::filesystem::path dir = repo_root / "data" / std::to_string(test_case_no);
    std::filesystem::create_directories(dir);

    std::ofstream output(dir / "simplify.txt");
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << stream.size() << '\n';
    for (const auto& point : stream) {
        output << CGAL::to_double(point.x()) << ' ' << CGAL::to_double(point.y()) << '\n';
    }
    std::cout << "Output Written\n";
    return 0;
}

// Load N\\n x y polyline files written by dots / simplify / benchmarks.
bool load_polyline_file(const std::filesystem::path& path, std::vector<Point>& out) {
    std::ifstream fin(path);
    if (!fin) return false;
    int count = 0;
    if (!(fin >> count) || count < 0) return false;
    out.clear();
    out.reserve(count);
    for (int i = 0; i < count; ++i) {
        double x, y;
        if (!(fin >> x >> y)) return false;
        out.emplace_back(x, y);
    }
    return !out.empty();
}

// Overlay DOTS (or other named baseline) when its output file is present.
void maybe_overlay_baseline(MultiViewer& viewer, int test_case_no) {
    if (!show_dots || test_case_no == -1) return;
    const std::filesystem::path dir = repo_root / "data" / std::to_string(test_case_no);
    const std::filesystem::path candidates[] = {
        dir / "dots_simplified.txt",
        dir / "DOTS.txt",
    };
    std::vector<Point> dots;
    std::filesystem::path used;
    for (const auto& path : candidates) {
        if (load_polyline_file(path, dots)) {
            used = path;
            break;
        }
    }
    if (used.empty()) {
        std::cout << "No DOTS overlay found under " << dir.string() << '\n';
        return;
    }
    // Hot pink matches the web compare layer (#f472b6).
    viewer.addCurve(dots, QColor(244, 114, 182), QStringLiteral("DOTS"));
    std::cout << "Overlaid DOTS (" << dots.size() << " points) from " << used.string() << '\n';
}

// ===========================================================================
//  Core algorithm
// ===========================================================================

int get_longest_stab(const std::vector<Point>& stream, int cur,
                     std::vector<Point>& simplified,
                     double epsilon, double delta,
                     MultiViewer* viewer = nullptr) {
    const Point& p0 = stream[cur];
    std::vector<Point> P;
    P = get_boundary_points_from_grid(p0, epsilon, delta);

    if (viewer) {
        viewer->markP0(p0);
        viewer->addOriginalPoint(p0);
    }

    std::array<Point, 2> buffer = {p0, p0};
    const int Pn = static_cast<int>(P.size());
    std::vector<std::vector<Point>> S(Pn);
    for (int i = 0; i < Pn; ++i) S[i] = {P[i]};
    int dead_cnt = 0;
    std::vector<int> dead(Pn);
    std::vector<std::vector<Point>> new_S(Pn);
    std::vector<std::vector<Point>> F(Pn);
    std::vector<Point> Gi;

    cur++;
    while (cur < static_cast<int>(stream.size())) {
        const Point& pi = stream[cur];
        Gi = get_conv_from_grid(pi, epsilon, delta);

        for (int i = 0; i < Pn; ++i) {
            if (dead[i]) continue;
            find_F(P[i], S[i], F[i]);
            bool hit;
            hit = intersect(F[i], Gi, new_S[i]);
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

        if (viewer) {
            const QColor step_colors[] = {Qt::red, Qt::blue, Qt::green, Qt::magenta, Qt::cyan};
            const QColor color = step_colors[cur % 5];
            for (int i = 0; i < Pn; ++i) {
                if (dead[i]) continue;
                if (showF) viewer->addPolygon(Polygon(F[i].begin(), F[i].end()), color);
                if (showG) viewer->addPolygon(Polygon(Gi.begin(), Gi.end()), color);
                if (showS) viewer->addPolygon(Polygon(S[i].begin(), S[i].end()), color);
                break;
            }
        }

        for (int i = 0; i < Pn; ++i) {
            if (!dead[i]) S[i].swap(new_S[i]);
        }
        if (viewer) {
            viewer->addOriginalPoint(pi);
            viewer->markPi(pi);
            viewer_process_events();
            if (!keep_polygons) viewer->clearPolygons();
        }
        cur++;
    }

    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    if (viewer) {
        viewer->addSimplifiedPoint(buffer[0]);
        viewer->addSimplifiedPoint(buffer[1]);
        viewer->clearMarkedP0();
        viewer->clearMarkedPi();
        viewer_process_events();
    }
    return cur;
}

std::vector<Point> simplify(const std::vector<Point>& stream,
                            double epsilon, double delta,
                            MultiViewer* viewer = nullptr) {
    configure_bbox(stream, epsilon, delta);
    std::vector<Point> simplified;
    std::cout << "Simplifying...\n";

    auto core_start = std::chrono::high_resolution_clock::now();
    int cur = 0;
    while (cur != static_cast<int>(stream.size())) {
        cur = get_longest_stab(stream, cur, simplified, epsilon, delta, viewer);
    }
    auto core_end = std::chrono::high_resolution_clock::now();
    double core_ms = std::chrono::duration<double, std::milli>(core_end - core_start).count();
    std::fprintf(stderr, "SIMPLIFY_CORE_MS %.4f\n", core_ms);

    std::cout << "Expected Frechet distance: " << std::sqrt(expected_frechet_squared) << '\n';
    std::cout << "The original stream of size " << stream.size()
              << " is simplified to " << simplified.size() << " points.\n";
    return simplified;
}

// ===========================================================================
//  main
// ===========================================================================

int main(int argc, char** argv) {
    int test_case_no = -1;
    int code = get_repo_root(argv, repo_root);
    std::cout << repo_root.string() << std::endl;
    if (code != 0) return code;

    code = parse_arguments(argc, argv, test_case_no);
    if (code != 0) return code;
    if (help_flag) return 0;

    std::vector<Point> stream;
    code = read_stream(test_case_no, stream);
    if (code != 0) return code;

    QApplication app(argc, argv);
    MultiViewer viewer;
    viewer.setParameters(DELTA, EPSILON);
    viewer.setShowLabels(show_labels);
    viewer.show();
    viewer.addOriginalPoints(stream);
    viewer_process_events();

    std::vector<Point> simplified = simplify(stream, EPSILON, DELTA, &viewer);
    stream = std::move(simplified);

    if (show_simplified) {
        viewer.clearSimplified();
        viewer.addSimplifiedPoints(stream);
        viewer_process_events();
    }

    maybe_overlay_baseline(viewer, test_case_no);
    viewer_process_events();

    if (out_flag) {
        code = out_stream(test_case_no, stream);
        if (code != 0) return code;
    }

    if (dist_flag && test_case_no != -1) {
        std::filesystem::path frechet_path = repo_root / "scripts" / "frechet";
        std::string command = std::string("\"") + frechet_path.string()
            + "\" --in " + std::to_string(test_case_no)
            + " --path \"" + (repo_root / "data" / std::to_string(test_case_no) / "simplify.txt").string() + "\"";
        std::system(command.c_str());
    }

    return app.exec();
}
