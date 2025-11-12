#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <cgal/simple_cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <fstream>
#include <cstdlib>
#include <filesystem>
#include <cstring>
#include <QApplication>
#include "drawing.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = Kernel::Point_2;
using Segment = Kernel::Segment_2;
using Ray = Kernel::Ray_2;
using Bbox = CGAL::Bbox_2;
using Rect = CGAL::Iso_rectangle_2<Kernel>;
using Polygon = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes =  CGAL::Polygon_with_holes_2<Kernel>;

bool showF = false; // true if -F is passed
bool showG = false; // true if -G is passed
bool showS = false; // true if -S is passed

constexpr double TOL = 1e-6;
constexpr double SQRT2 =
    1.41421356237; // sqrt in STL does not have constexpr version !?

// Bounding square [BMIN, BMAX]^2
extern const double BMIN = -10000;
extern const double BMAX = 10000;
// Bounding data points (convex hull should not exceed bounding box)
const double DATAMIN = -8000;
const double DATAMAX= 8000;
// TODO: DELTA should not be too high that convex hull goes out of bounding square, which may cause the program to crash
constexpr double DELTA = 200;
constexpr double EPSILON = 0.5;
// Runtime-tunable parameters (overridable via -d and -e flags)
static double G_EPSILON = EPSILON;
static double G_DELTA = DELTA;
static inline double GRID_val() { return G_EPSILON * G_DELTA / (2 * SQRT2); }
constexpr double GRID = EPSILON * DELTA / (2 * SQRT2);

static void normalize_stream(std::vector<Point> &stream) {
    if (stream.empty()) return;
    double minx = std::numeric_limits<double>::infinity();
    double miny = std::numeric_limits<double>::infinity();
    double maxx = -std::numeric_limits<double>::infinity();
    double maxy = -std::numeric_limits<double>::infinity();
    for (const auto &p : stream) {
        double x = CGAL::to_double(p.x());
        double y = CGAL::to_double(p.y());
        minx = std::min(minx, x);
        miny = std::min(miny, y);
        maxx = std::max(maxx, x);
        maxy = std::max(maxy, y);
    }
    double data_min = std::min(minx, miny);
    double data_max = std::max(maxx, maxy);
    if (data_max <= data_min) return; // nothing to do
    double scale = (DATAMAX - DATAMIN) / (data_max - data_min);
    for (auto &p : stream) {
        double x = CGAL::to_double(p.x());
        double y = CGAL::to_double(p.y());
        double nx = DATAMIN + (x - data_min) * scale;
        double ny = DATAMIN + (y - data_min) * scale;
        p = Point(nx, ny);
    }
}

// Copied from examples/Boolean_set_operations/print_utils.cpp
// Pretty-print a CGAL polygon.
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

// TO BE REMOVED: print simple/orientation/area
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

// Pretty-print a polygon with holes.
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


bool point_in_convex(const Point &p, const std::vector<Point> &poly, bool ccw = true) {
    for (size_t i = 0, n = poly.size(); i < n; ++i) {
        const auto &a = poly[i];
        const auto &b = poly[(i + 1) % n];
        if (CGAL::orientation(a, b, p) == (ccw ? CGAL::RIGHT_TURN : CGAL::LEFT_TURN))
            return false;
    }
    return true;
}

std::optional<Point> ray_hit_bbox(const Point &p, const Point &dir) {
    double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y());
    double dx = CGAL::to_double(dir.x()) - px,
           dy = CGAL::to_double(dir.y()) - py;
    double best = std::numeric_limits<double>::infinity(), hx = 0, hy = 0;
    auto consider = [&best, &hx, &hy, dx, dy, px, py](double t) {
        if (t <= 0)
            return;
        double x = px + t * dx, y = py + t * dy;
        if (x < BMIN - 1e-8 || x > BMAX + 1e-8 || y < BMIN - 1e-8 ||
            y > BMAX + 1e-8)
            return;
        if (t < best) {
            best = t;
            hx = x;
            hy = y;
        }
    };
    if (std::abs(dx) > 1e-18) {
        consider((BMIN - px) / dx);
        consider((BMAX - px) / dx);
    }
    if (std::abs(dy) > 1e-18) {
        consider((BMIN - py) / dy);
        consider((BMAX - py) / dy);
    }
    if (!std::isfinite(best))
        return std::nullopt;
    return Point(hx, hy);
}

enum class Bbox_edge {
    BL = 0,
    BOTTOM = 1,
    BR = 2,
    RIGHT = 3,
    TR = 4,
    TOP = 5,
    TL = 6,
    LEFT = 7,
};

// Return current bbox corners (clockwise starting from bottom-left)
static std::array<Point, 4> current_bbox_corner() {
    return {Point(BMIN, BMIN), Point(BMAX, BMIN), Point(BMAX, BMAX),
            Point(BMIN, BMAX)};
}

// Return current bbox as vector of points
static std::vector<Point> current_bbox() {
    auto c = current_bbox_corner();
    return std::vector<Point>(c.begin(), c.end());
}

std::optional<Bbox_edge> which_edge(const Point &s) {
    double x = CGAL::to_double(s.x()), y = CGAL::to_double(s.y());
    bool on_left   = std::abs(x - BMIN) < TOL;
    bool on_right  = std::abs(x - BMAX) < TOL;
    bool on_bottom = std::abs(y - BMIN) < TOL;
    bool on_top    = std::abs(y - BMAX) < TOL;

    // Corners first
    if (on_left && on_bottom) return Bbox_edge::BL;
    if (on_right && on_bottom) return Bbox_edge::BR;
    if (on_right && on_top) return Bbox_edge::TR;
    if (on_left && on_top) return Bbox_edge::TL;

    // Edges
    if (on_bottom) return Bbox_edge::BOTTOM;
    if (on_right)  return Bbox_edge::RIGHT;
    if (on_top)    return Bbox_edge::TOP;
    if (on_left)   return Bbox_edge::LEFT;
    return std::nullopt;
}

void append_rect_pts(std::vector<Point> &out, Bbox_edge from, Bbox_edge to,
                     bool ccw) {
    auto corners = current_bbox_corner(); // order: BL(0), BR(1), TR(2), TL(3)

    auto next = [&](int idx) {
        return (idx + (ccw ? 1 : 7)) % 8; // modulo 8
    };

    int i = static_cast<int>(from);
    int j = static_cast<int>(to);
    if (i == j) return;
    for (int k = next(i); k != j; k = next(k)) {
        if ((k & 1) == 0) { // corner indices: 0,2,4,6
            int ci = k / 2; // map 0->0(BL),2->1(BR),4->2(TR),6->3(TL)
            out.push_back(corners[ci]);
        }
    }
}

std::vector<Point> get_conv_from_grid(const Point &p) {
    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    const double r = (1.0 + G_EPSILON / 2.0) * G_DELTA;
    const double GRID = GRID_val();
    const double r2 = r * r;
    // treat p as (0, 0), find the topmost index of grid point that is contained
    // in the G_p
    int y_min = -(r / GRID);
    int y_max = -y_min;
    std::vector<Point> boundaries;
    boundaries.reserve(2 * y_max + 1);
    for (int y = y_min; y <= y_max; y++) {
        const double y_actual = y * GRID;
        const int x_min = sqrt(r2 - y_actual * y_actual) / GRID;
        boundaries.emplace_back(px + x_min * GRID, py + y_actual);
    }
    for (int y = y_min; y <= y_max; y++) {
        const double y_actual = y * GRID;
        const int x_max = -(sqrt(r2 - y_actual * y_actual) / GRID);
        if (x_max != 0) // to avoid duplicates
            boundaries.emplace_back(px + x_max * GRID, py + y_actual);
    }
    std::vector<Point> conv;
    // This uses the O(n log n) implementation
    // cuz h = O(n) in this case
    CGAL::convex_hull_2(boundaries.begin(), boundaries.end(),
                        std::back_inserter(conv));
    return conv;
}

std::vector<Point> get_points_from_grid(const Point &p) {
    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    const double r = (1.0 + G_EPSILON / 2.0) * G_DELTA;
    const double GRID = GRID_val();
    const double r2 = r * r;
    // treat p as (0, 0), find the topmost index of grid point that is contained
    // in the G_p
    int y_min = -(r / GRID);
    int y_max = -y_min;
    std::vector<Point> points;
    points.reserve(2 * y_max + 1);
    for (int y = y_min; y <= y_max; y++) {
        const double y_actual = y * GRID;
        const int x_min = sqrt(r2 - y_actual * y_actual) / GRID;
        const int x_max = -x_min;
        for (int x = x_min; x <= x_max; x++) {
            points.emplace_back(px + x * GRID, py + y_actual);
        }
    }
    return points;
}

std::vector<int> find_tangent_idx(const Point &p,
                                 const std::vector<Point> &S) {
    int n = S.size();
    std::vector<int> tangent;
    for (int i = 0; i < n; i++) {
        Point pred = S[(i - 1 + n) % n];
        Point now = S[i];
        Point succ = S[(i + 1 + n) % n];
        if (CGAL::orientation(p, pred, now) // be careful of the collinear case
                != CGAL::orientation(p, now, succ) &&
            CGAL::orientation(p, now, succ) != CGAL::COLLINEAR) {
            // std::cerr << "[find tangent idx]: " << to_string(CGAL::orientation(p, pred, now)) << ' ' << to_string(CGAL::orientation(p, now, succ)) << std::endl;
            tangent.push_back(i);
        }
    }
    return tangent;
}

std::optional<Point> intersect_ray_with_rect(const Point& p, const Point& direction) { // This is ok?
    // degenerate ray direction
    if (CGAL::to_double(CGAL::squared_distance(p, direction)) < TOL*TOL) return std::nullopt;

    Ray ray(p, direction);
    Bbox box(BMIN, BMIN, BMAX, BMAX);

    if (auto obj = CGAL::intersection(Rect(box), ray)) {
        if (const Point* ip = std::get_if<Point>(&*obj)) {
            // ignore self-hit at the origin of the ray
            if (CGAL::to_double(CGAL::squared_distance(*ip, p)) < TOL*TOL) return std::nullopt;
            return *ip;
        }
        if (const Segment* seg = std::get_if<Segment>(&*obj)) {
            // Ray overlaps an edge of the rectangle: pick the correct endpoint
            const Point& a = seg->source();
            const Point& b = seg->target();
            double da = CGAL::to_double(CGAL::squared_distance(a, p));
            double db = CGAL::to_double(CGAL::squared_distance(b, p));
            if (da < TOL*TOL) return b;           // p is at a -> take the far endpoint
            if (db < TOL*TOL) return a;           // p is at b -> take the far endpoint
            return (da < db) ? a : b;             // otherwise take the nearer endpoint
        }
    }
    return std::nullopt;
}

std::vector<Point> find_F(const Point& p, const std::vector<Point>& S) {
    // assert(S.size() != 2); // wait why this check??
    if (S.size() == 1) {
        auto F = current_bbox();
        return F;
    }
    if (point_in_convex(p, S)) {
        auto F = current_bbox();
        return F;
    }
    std::vector<int> tangent = find_tangent_idx(p, S);
    assert(tangent.size() == 2);

    std::optional<Point> hit1 = intersect_ray_with_rect(p, S[tangent[0]]);
    std::optional<Point> hit2 = intersect_ray_with_rect(p, S[tangent[1]]);
    if (!hit1 || !hit2) {
        std::cerr << "Ray doesn't intersect with bounding box!\n";
        auto F = current_bbox();
        return F;
    }
    std::optional<Bbox_edge> e1 = which_edge(hit1.value());
    std::optional<Bbox_edge> e2 = which_edge(hit2.value());

    int n = int(S.size());
    assert(n >= 3);
    CGAL_precondition(Polygon(S.begin(), S.end()).is_counterclockwise_oriented());
    assert(tangent[1] - tangent[0] - 1 >= 1 || tangent[0] + n - tangent[1] - 1 >= 1);
    std::vector<Point> F;
    if (CGAL::right_turn(p, S[tangent[0]], S[tangent[1]]))  {
        // [i..j] (inclusive)
        std::copy(S.begin() + tangent[0], S.begin() + tangent[1] + 1,
                  std::back_inserter(F));
        // walk from hit2 to hit1 ccw to close
        F.push_back(hit2.value());
        append_rect_pts(F, e2.value(), e1.value(), true);
        F.push_back(hit1.value());
    } else {
        // [j..n-1] + [0..i] (inclusive)
        std::copy(S.begin() + tangent[1], S.end(),
                  std::back_inserter(F));
        std::copy(S.begin(), S.begin() + tangent[0] + 1,
                  std::back_inserter(F));
        F.push_back(hit1.value());
        append_rect_pts(F, e1.value(), e2.value(), true);
        F.push_back(hit2.value());
    }
    if (!e1 || !e2) {
        std::cerr << "Cannot determine which Bbox edge the ray intersect with.\n";
        auto F = current_bbox();
        return F;
    }
    return F;
}

int get_longest_stab(const std::vector<Point> &stream, int cur,
                     std::vector<Point> &simplified, MultiViewer* viewer = nullptr) {
    const Point& p0 = stream[cur];
    int p0cur = cur;
    std::cerr << "[p0]: " << cur << std::endl;
    if (viewer) viewer->addOriginalPoint(p0);
    if (viewer) viewer->markP0(p0);
    std::vector<Point> P = get_points_from_grid(p0);
    std::array<Point, 2> buffer = {p0};
    std::vector<std::vector<Point>> S(P.size(), std::vector<Point>{p0});
    int dead_cnt = 0;
    std::vector<int> dead(P.size());
    while (cur < int(stream.size())) {
        const Point& pi  = stream[cur];
        if (viewer) viewer->addOriginalPoint(pi);
        std::vector<std::vector<Polygon_with_holes>> new_S(P.size());
        for (int i = 0; i < int(P.size()); i++) {
            if (dead[i]) {
                continue;
            }
            std::vector<Point> F = find_F(p0, S[i]);
            std::vector<Point> Gi = get_conv_from_grid(pi);
            Polygon F_poly(F.begin(), F.end());
            Polygon Gi_poly(Gi.begin(), Gi.end());
            Polygon S_poly(S[i].begin(), S[i].end());

            if (showF && viewer) 
                viewer->addPolygon(F_poly);
            else if (showG && viewer)
                viewer->addPolygon(Gi_poly);
            else if (showS && viewer)
                viewer->addPolygon(S_poly);

            CGAL::intersection(F_poly, Gi_poly, back_inserter(new_S[i]));

            if (new_S[i].size() == 0) {
                dead[i] = true;
                dead_cnt++;
                continue;
            }
            assert(new_S[i].size() == 1);
            buffer[1] = *new_S[i].begin()->outer_boundary().vertices_begin();
        }
        if (dead_cnt == int(P.size())) {
            break;
        }
        for (int i = 0; i < int(P.size()); i++) {
            if (dead[i]) continue;
            S[i].clear();
            std::copy(new_S.begin()->begin()->outer_boundary().vertices_begin(),
                      new_S.begin()->begin()->outer_boundary().vertices_end(),
                      std::back_inserter(S[i]));
        }
        cur++;
        if (viewer) viewer_process_events();
    }
    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    if (viewer) viewer->addSimplifiedPoint(buffer[0]);
    if (viewer) viewer->addSimplifiedPoint(buffer[1]);
    if (viewer) viewer_process_events();
    return cur;
}

std::vector<Point> simplify(const std::vector<Point> &stream, MultiViewer* viewer = nullptr) {
    std::vector<Point> simplified;
    int cur = 0;
    while (cur != int(stream.size())) {
        cur = get_longest_stab(stream, cur, simplified, viewer);
    }
    return simplified;
}

static void print_help() {
    std::cout << "Usage: simplify [options]\n"
              << "  --in <id>        Read input from ../data/taxi/<id>.txt\n"
              << "  --out            Write output to ../data/taxi_simplified/<id>/original.txt & simplified.txt (requires --in)\n"
              << "  --dist           After output, compute Frechet distance via ./frechet -n <id>\n"
              << "  --gui            Show GUI viewer (omit for headless run)\n"
              << "  -d <delta>       Override DELTA (default " << DELTA << ")\n"
              << "  -e <epsilon>     Override EPSILON (default " << EPSILON << ")\n"
              << "  -F/-G/-S         Debug polygon display modes\n"
              << "  -h               Show this help and exit\n";
}

int main(int argc, char** argv) {
    bool out_flag = false;      // write outputs if true
    bool gui_flag = false;      // show viewer if true
    bool dist_flag = false;     // compute frechet if true
    int test_case_no = -1;      // provided by --in <id>

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i],"-F") == 0) showF = true;
        else if (strcmp(argv[i],"-S") == 0) showS = true;
        else if (strcmp(argv[i],"-G") == 0) showG = true;
        else if (strcmp(argv[i],"--gui") == 0) gui_flag = true;
        else if (strcmp(argv[i],"--out") == 0) out_flag = true;
        else if (strcmp(argv[i],"--dist") == 0) dist_flag = true;
        else if (strcmp(argv[i],"-d") == 0 && i+1 < argc) {
            try { G_DELTA = std::stod(argv[++i]); } catch(...) { std::cerr << "Invalid -d value\n"; return 1; }
        }
        else if (strcmp(argv[i],"-e") == 0 && i+1 < argc) {
            try { G_EPSILON = std::stod(argv[++i]); } catch(...) { std::cerr << "Invalid -e value\n"; return 1; }
        }
        else if (strcmp(argv[i],"-h") == 0) { print_help(); return 0; }
        else if (strcmp(argv[i],"--in") == 0 && i+1 < argc) {
            try { test_case_no = std::stoi(argv[++i]); }
            catch(...) { std::cerr << "Invalid --in argument\n"; return 1; }
        }
    }

    if (dist_flag) out_flag = true;

    // Load input stream
    std::vector<Point> stream;
    if (test_case_no != -1) {
        std::string input_file = std::string("../data/taxi/") + std::to_string(test_case_no) + ".txt";
        std::ifstream fin(input_file);
        if (!fin) { std::cerr << "Cannot open " << input_file << "\n"; return 1; }
        int N = 0;
        if (!(fin >> N)) { std::cerr << "Empty or invalid input in " << input_file << "\n"; return 1; }
        stream.resize(N);
        for (int i = 0; i < N; ++i) {
            if (!(fin >> stream[i])) { std::cerr << "Malformed point at index " << i << " in " << input_file << "\n"; return 1; }
        }
    } else {
        int N = 0;
        if (!(std::cin >> N)) { std::cerr << "Expected N from stdin\n"; return 1; }
        stream.resize(N);
        for (int i = 0; i < N; ++i) { if (!(std::cin >> stream[i])) { std::cerr << "Malformed stdin point at " << i << "\n"; return 1; } }
    }

    normalize_stream(stream);

    // Optional GUI
    QApplication app(argc, argv);
    MultiViewer viewer;
    MultiViewer* vptr = nullptr;
    if (gui_flag) {
        vptr = &viewer;
        viewer.setParameters(G_DELTA, G_EPSILON);
        viewer.show();
    }

    // Simplify
    std::vector<Point> simplified = simplify(stream, vptr);

    // Optional output
    if (out_flag) {
        if (test_case_no == -1) {
            std::cerr << "--out requires --in <id> to determine output location\n";
        } else {
            std::filesystem::path dir = std::filesystem::path("../data/taxi_simplified") / std::to_string(test_case_no);
            std::filesystem::create_directories(dir);
            // original
            std::ofstream orig(dir / "original.txt");
            for (const auto& p : stream) orig << CGAL::to_double(p.x()) << ' ';
            orig << '\n';
            for (const auto& p : stream) orig << CGAL::to_double(p.y()) << ' ';
            orig << '\n';
            // simplified
            std::ofstream simp(dir / "simplified.txt");
            for (const auto& p : simplified) simp << CGAL::to_double(p.x()) << ' ';
            simp << '\n';
            for (const auto& p : simplified) simp << CGAL::to_double(p.y()) << ' ';
            simp << '\n';
        }
        std::cout << "Output Written\n";
    }

    // Optional distance computation
    if (dist_flag && test_case_no != -1) {
        // Assume running from build/; frechet script in project root
        // Try ./frechet first, then python3 frechet, then python3 frechet.py
        std::string cmd1 = std::string("../frechet -id ") + std::to_string(test_case_no);
        int rc = std::system(cmd1.c_str());
    }

    if (gui_flag) return app.exec();
    return 0;
}