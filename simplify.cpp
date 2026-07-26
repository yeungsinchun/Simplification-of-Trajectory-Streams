#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <algorithm>
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

bool out_flag = false;
bool dist_flag = false;

// ===========================================================================
//  Help & debug printing
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

// Copied from examples/Boolean_set_operations/print_utils.cpp
template<class Kernel, class Container>
void print_polygon (const CGAL::Polygon_2<Kernel, Container>& P)
{
  typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator  vit;

  std::cout << "[ " << P.size() << " vertices:";
  for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    std::cout << " (" << *vit << ')';
  std::cout << " ]" << std::endl;

  return;
}

static void print_poly_info(const Polygon& P, const char* name = "poly") {
    std::cerr << name << ": size=" << P.size();
    bool simple = false;
    try {
        simple = P.is_simple();
        std::cerr << ", simple=" << (simple ? "yes" : "no");
    } catch (...) {
        std::cerr << ", simple=err\n";
        return;
    }
    if (!simple) {
        std::cerr << "\n";
        return;
    }
    try {
        std::cerr << ", orient="
                  << (P.is_counterclockwise_oriented() ? "CCW" : "CW");
    } catch (...) {
        std::cerr << ", orient=err";
    }
    try {
        std::cerr << ", area=" << CGAL::to_double(P.area()) << "\n";
    } catch (...) {
        std::cerr << ", area=err\n";
    }
}

template<class Kernel, class Container>
void print_polygon_with_holes
    (const CGAL::Polygon_with_holes_2<Kernel, Container>& pwh)
{
  if (! pwh.is_unbounded())
  {
    std::cout << "{ Outer boundary = ";
    print_polygon (pwh.outer_boundary());
  }
  else
    std::cout << "{ Unbounded polygon." << std::endl;

  typename CGAL::Polygon_with_holes_2<Kernel,Container>::
                                             Hole_const_iterator  hit;
  unsigned int                                                     k = 1;

  std::cout << "  " << pwh.number_of_holes() << " holes:" << std::endl;
  for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k)
  {
    std::cout << "    Hole #" << k << " = ";
    print_polygon (*hit);
  }
  std::cout << " }" << std::endl;

  return;
}

// ===========================================================================
//  Exact intersection (ground truth) — EPECK CGAL::intersection
// ===========================================================================

// Computes the intersection of two convex polygons using EPECK exact
// constructions. Mirrors the semantics of the fast intersect() in
// simplify_geometry.h: returns true and fills `result` (CCW, no duplicate
// closing vertex, >= 3 vertices) when the intersection has positive area;
// otherwise clears `result` and returns false.
inline bool intersect_exact(const std::vector<Point>& P_verts,
                            const std::vector<Point>& Q_verts,
                            std::vector<Point>& result) {
    using EPolygon = CGAL::Polygon_2<Epeck>;
    using EPwh     = CGAL::Polygon_with_holes_2<Epeck>;

    result.clear();
    if (P_verts.size() < 3 || Q_verts.size() < 3) return false;

    EPolygon P, Q;
    for (const auto& p : P_verts) P.push_back(conv_to_exact(p));
    for (const auto& q : Q_verts) Q.push_back(conv_to_exact(q));
    if (P.orientation() == CGAL::CLOCKWISE) P.reverse_orientation();
    if (Q.orientation() == CGAL::CLOCKWISE) Q.reverse_orientation();

    std::vector<EPwh> out;
    CGAL::intersection(P, Q, std::back_inserter(out));
    if (out.empty()) return false;

    const EPolygon& ob = out.front().outer_boundary();
    if (ob.size() < 3) return false;

    result.reserve(ob.size());
    for (auto it = ob.vertices_begin(); it != ob.vertices_end(); ++it) {
        result.emplace_back(CGAL::to_double(it->x()), CGAL::to_double(it->y()));
    }
    return true;
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
#ifdef DUMP_DIVERGENCE
            {
                // Does RAW O'Rourke output (pre-guard) actually differ from the
                // exact result, or only the guard verdict?  Compare the two
                // polygons cyclically within a coordinate tolerance.
                std::vector<Point> Fd = orourke_cgal::dedup_consecutive(F[i]);
                std::vector<Point> Gd = orourke_cgal::dedup_consecutive(Gi);
                std::vector<Point> raw = orourke_cgal::convex_intersect_robust(Fd, Gd);
                if (!raw.empty() && raw.front() == raw.back()) raw.pop_back();
                std::vector<Point> ex;
                bool ex_hit = intersect_exact_convex(F[i], Gi, ex);

                // Length-independent shape comparison: symmetric Hausdorff
                // distance between the RAW and EXACT vertex sets.  Immune to
                // degenerate edges, vertex count, and cyclic ordering, so it
                // measures GENUINE geometric divergence only.
                auto directed_hausdorff = [](const std::vector<Point>& A,
                                             const std::vector<Point>& B) -> double {
                    double worst = 0.0;
                    for (const auto& a : A) {
                        double best = 1e300;
                        for (const auto& b : B) {
                            double dx = CGAL::to_double(a.x()) - CGAL::to_double(b.x());
                            double dy = CGAL::to_double(a.y()) - CGAL::to_double(b.y());
                            best = std::min(best, dx*dx + dy*dy);
                        }
                        worst = std::max(worst, best);
                    }
                    return std::sqrt(worst);
                };

                static long n_match = 0;      // RAW ~= EXACT geometrically
                static long n_genuine = 0;     // RAW genuinely differs
                static long n_both_empty = 0;
                static double max_hd = 0.0;
                struct Reporter {
                    long *m, *g, *be; double* mh;
                    ~Reporter() {
                        std::ofstream f("hausdorff_tally.txt");
                        f << "match=" << *m << " genuine_diff=" << *g
                          << " both_empty=" << *be << " max_hausdorff=" << *mh << "\n";
                    }
                };
                static Reporter rep{&n_match, &n_genuine, &n_both_empty, &max_hd};

                bool raw_ok = raw.size() >= 3;
                if (!raw_ok && !ex_hit) { n_both_empty++; }
                else if (raw_ok && ex_hit) {
                    double hd = std::max(directed_hausdorff(raw, ex),
                                         directed_hausdorff(ex, raw));
                    max_hd = std::max(max_hd, hd);
                    if (hd < 1e-3) n_match++;
                    else {
                        n_genuine++;
                        static int dumped = 0;
                        if (dumped < 5) {
                            dumped++;
                            std::ofstream d("divergence_case_" + std::to_string(dumped) + ".txt");
                            d << std::setprecision(std::numeric_limits<double>::max_digits10);
                            d << "# hausdorff=" << hd << " raw=" << raw.size()
                              << " exact=" << ex.size() << "\n";
                            d << "F " << F[i].size() << "\n";
                            for (auto& p : F[i]) d << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
                            d << "Gi " << Gi.size() << "\n";
                            for (auto& p : Gi) d << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
                            d << "RAW " << raw.size() << "\n";
                            for (auto& p : raw) d << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
                            d << "EXACT " << ex.size() << "\n";
                            for (auto& p : ex) d << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
                            d.close();
                        }
                    }
                } else {
                    // one empty, one not — a genuine hit/miss disagreement
                    n_genuine++;
                }
            }
#endif
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
    int cur = 0;
    int prefix = 0;
    while (cur != int(stream.size())) {
        SIGNPOST_EVENT(simplify.prefix, "from=%d size=%lu", cur, (unsigned long)simplified.size());
        cur = get_longest_stab(stream, cur, simplified, EPSILON, DELTA);
        prefix++;
    }
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

// TODO: clean up this mess of flags
int main(int argc, char** argv) {
    int test_case_no = -1;

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

    std::vector<Point> stream;
    std::filesystem::path repo_root;
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
        repo_root = find_repo_root(argv[0], 5);
    } catch (const std::filesystem::filesystem_error& e) {
        try {
            repo_root = find_repo_root(std::filesystem::current_path(), 5);
        } catch (const std::filesystem::filesystem_error& fallback_error) {
            std::cerr << "Error: could not resolve the data directory: " << fallback_error.what() << "\n";
            return 1;
        }
    }
    if (test_case_no != -1) {
        auto simp_orig = repo_root / "data" / std::to_string(test_case_no) / "original.txt";
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

    SIGNPOST_BEGIN(app.run);
    std::vector<Point> simplified = simplify(stream, EPSILON, DELTA);
    stream = std::move(simplified);
    SIGNPOST_END(app.run);

    if (out_flag) {
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
    }

    if (dist_flag && test_case_no != -1) {
        std::filesystem::path frechet_path = repo_root / "scripts" / "frechet";
        std::string cmd1 = std::string("\"") + frechet_path.string() + "\" --in " + std::to_string(test_case_no) + " --path \"" + (repo_root / "data" / std::to_string(test_case_no) / "simplify.txt").string() + "\"";
        int rc = std::system(cmd1.c_str());
    }

    return 0;
}
