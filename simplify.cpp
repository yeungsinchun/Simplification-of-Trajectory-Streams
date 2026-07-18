#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>
#include <QApplication>
#include "drawing.h"
#include "simplify_geometry.h"
#include "timer.h"

// ===========================================================================
//  Global parameters
// ===========================================================================

double DELTA = 200;
double EPSILON = 0.5;

bool showF = false; // true if -F is passed
bool showG = false; // true if -G is passed
bool showS = false; // true if -S is passed
bool out_flag = false;
bool gui_flag = false;
bool dist_flag = false;

// ===========================================================================
//  Trajectory-simplification geometry
// ===========================================================================

static constexpr double BMIN = -10000;
static constexpr double BMAX =  10000;
static constexpr double TOL  = 1e-6;

inline bool point_in_convex(const Point& p, const std::vector<Point>& poly, bool ccw = true) {
    const CGAL::Orientation outside = ccw ? CGAL::RIGHT_TURN : CGAL::LEFT_TURN;
    for (size_t i = 0, n = poly.size(); i < n; ++i) {
        if (CGAL::orientation(poly[i], poly[(i + 1) % n], p) == outside) {
            return false;
        }
    }
    return true;
}

inline std::optional<Point> ray_hit_bbox(const Point& p, const Point& dir) {
    double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y());
    double dx = CGAL::to_double(dir.x()) - px,
           dy = CGAL::to_double(dir.y()) - py;
    double best = std::numeric_limits<double>::infinity(), hx = 0, hy = 0;
    auto consider = [&best, &hx, &hy, dx, dy, px, py](double t) {
        if (t <= 0) return;
        double x = px + t * dx, y = py + t * dy;
        if (x < BMIN - 1e-8 || x > BMAX + 1e-8 ||
            y < BMIN - 1e-8 || y > BMAX + 1e-8) return;
        if (t < best) { best = t; hx = x; hy = y; }
    };
    if (std::abs(dx) > 1e-18) {
        consider((BMIN - px) / dx);
        consider((BMAX - px) / dx);
    }
    if (std::abs(dy) > 1e-18) {
        consider((BMIN - py) / dy);
        consider((BMAX - py) / dy);
    }
    if (!std::isfinite(best)) return std::nullopt;
    return Point(hx, hy);
}

enum class Bbox_edge {
    BL = 0, BOTTOM = 1, BR = 2, RIGHT = 3,
    TR = 4, TOP = 5,    TL = 6, LEFT = 7
};

inline const std::array<Point, 4>& current_bbox_corner_static() {
    static const std::array<Point, 4> corners = {
        Point(BMIN, BMIN), Point(BMAX, BMIN), Point(BMAX, BMAX), Point(BMIN, BMAX)
    };
    return corners;
}

inline std::array<Point, 4> current_bbox_corner() {
    return current_bbox_corner_static();
}

inline const std::vector<Point>& current_bbox_static() {
    static const std::vector<Point> bbox(
        current_bbox_corner_static().begin(), current_bbox_corner_static().end());
    return bbox;
}

inline std::vector<Point> current_bbox() {
    return std::vector<Point>(current_bbox_corner_static().begin(),
                              current_bbox_corner_static().end());
}

inline std::optional<Bbox_edge> which_edge(const Point& s) {
    double x = CGAL::to_double(s.x()), y = CGAL::to_double(s.y());
    bool on_left   = std::abs(x - BMIN) < TOL;
    bool on_right  = std::abs(x - BMAX) < TOL;
    bool on_bottom = std::abs(y - BMIN) < TOL;
    bool on_top    = std::abs(y - BMAX) < TOL;

    if (on_left && on_bottom) return Bbox_edge::BL;
    if (on_right && on_bottom) return Bbox_edge::BR;
    if (on_right && on_top) return Bbox_edge::TR;
    if (on_left && on_top) return Bbox_edge::TL;
    if (on_bottom) return Bbox_edge::BOTTOM;
    if (on_right)  return Bbox_edge::RIGHT;
    if (on_top)    return Bbox_edge::TOP;
    if (on_left)   return Bbox_edge::LEFT;
    return std::nullopt;
}

inline void append_rect_pts(std::vector<Point>& out, Bbox_edge from, Bbox_edge to, bool ccw) {
    const auto& corners = current_bbox_corner_static();
    auto next = [&](int idx) { return (idx + (ccw ? 1 : 7)) % 8; };

    int i = static_cast<int>(from);
    int j = static_cast<int>(to);
    if (i == j) return;
    for (int k = next(i); k != j; k = next(k)) {
        if ((k & 1) == 0) {
            int ci = k / 2;
            out.push_back(corners[ci]);
        }
    }
}

template <class F>
inline void for_each_grid_row(double px, double py, double r, double GRID, F&& f) {
    const double r2 = r * r;
    const int y_min = -(r / GRID) - 1;
    const int y_max = -y_min;
    for (int y = y_min; y <= y_max; y++) {
        const double y_actual = y * GRID;
        const double dy2 = y_actual * y_actual;
        if (dy2 > r2) continue;
        const double dx = std::sqrt(r2 - dy2);
        const int x_min = static_cast<int>(-dx / GRID - 1);
        const int x_max = static_cast<int>(+dx / GRID + 1);
        f(px, py, y_actual, x_min, x_max);
    }
}

inline double GRID_val(double EPSILON, double DELTA) {
    return EPSILON * DELTA / (2.0 * std::sqrt(2.0));
}

inline double R_val(double EPSILON, double DELTA) {
    return (1.0 + EPSILON / 2.0) * DELTA;
}

inline double& expected_frechet_ref() {
    static double value = 0;
    return value;
}

inline std::vector<Point> get_points_from_grid(const Point& p, double EPSILON, double DELTA) {
    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    const double GRID = GRID_val(EPSILON, DELTA);
    if (DELTA == 0) return std::vector<Point>{p};

    const double r = R_val(EPSILON, DELTA);
    const double r2 = r * r;

    const int j_min = static_cast<int>(std::floor(-r / GRID));
    const int j_max = static_cast<int>(std::ceil(r / GRID));
    const int cell_count = j_max - j_min + 1;
    const int corner_count = cell_count + 1;

    // Adjacent cells share corners, so deduplicate them in a corner-sized
    // buffer. Cell indices end at j_max, but their upper corners reach
    // j_max + 1 on each axis.
    std::vector<Point> points;
    points.reserve(size_t(corner_count) * corner_count);
    std::vector<uint8_t> seen(size_t(corner_count) * corner_count, 0);

    auto add_corner = [&](int ji, int ki) {
        const int ix = ji - j_min, iy = ki - j_min;
        const size_t index = size_t(ix) * corner_count + iy;
        if (seen[index]) return;
        seen[index] = 1;
        points.emplace_back(px + ji * GRID, py + ki * GRID);
    };

    for (int j = j_min; j <= j_max; ++j) {
        const double x0 = j * GRID, x1 = (j + 1) * GRID;
        for (int k = j_min; k <= j_max; ++k) {
            const double y0 = k * GRID, y1 = (k + 1) * GRID;
            if (x0 * x0 + y0 * y0 > r2 &&
                x1 * x1 + y0 * y0 > r2 &&
                x0 * x0 + y1 * y1 > r2 &&
                x1 * x1 + y1 * y1 > r2) continue;
            add_corner(j,     k);
            add_corner(j + 1, k);
            add_corner(j + 1, k + 1);
            add_corner(j,     k + 1);
        }
    }

    // Diagonal of the boundary cell: r + GRID * sqrt(2) = (1 + EPSILON) * DELTA.
    // Tight Fréchet upper bound for this construction.
    const double diagonal_sq = (r + GRID * std::sqrt(2.0)) * (r + GRID * std::sqrt(2.0));
    if (diagonal_sq > expected_frechet_ref()) {
        expected_frechet_ref() = diagonal_sq;
    }
    return points;
}

inline std::vector<Point> get_conv_from_grid(const Point& p, double EPSILON, double DELTA) {
    std::vector<Point> points = get_points_from_grid(p, EPSILON, DELTA);
    std::vector<Point> conv;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(conv));
    return conv;
}

inline std::vector<int> find_tangent_idx(const Point& p, const std::vector<Point>& S) {
    const int n = static_cast<int>(S.size());
    std::vector<int> tangent;
    tangent.reserve(2);

    auto orient_at = [&](int i) {
        return CGAL::orientation(S[i], S[(i + 1) % n], p);
    };

    CGAL::Orientation prev = orient_at(n - 1);
    for (int i = 0; i < n; ++i) {
        const CGAL::Orientation cur = orient_at(i);
        if (cur != prev && cur != CGAL::COLLINEAR) {
            tangent.push_back(i);
            if (tangent.size() == 2) break;
        }
        prev = cur;
    }
    return tangent;
}

inline std::optional<Point> intersect_ray_with_rect(const Point& p, const Point& direction) {
    if (CGAL::to_double(CGAL::squared_distance(p, direction)) < TOL * TOL)
        return std::nullopt;

    Ray ray(p, direction);
    Bbox box(BMIN, BMIN, BMAX, BMAX);

    if (auto obj = CGAL::intersection(Rect(box), ray)) {
        if (const Point* ip = std::get_if<Point>(&*obj)) {
            if (CGAL::to_double(CGAL::squared_distance(*ip, p)) < TOL * TOL)
                return std::nullopt;
            return *ip;
        }
        if (const Segment* seg = std::get_if<Segment>(&*obj)) {
            const Point& a = seg->source();
            const Point& b = seg->target();
            double da = CGAL::to_double(CGAL::squared_distance(a, p));
            double db = CGAL::to_double(CGAL::squared_distance(b, p));
            if (da < TOL * TOL) return b;
            if (db < TOL * TOL) return a;
            return (da < db) ? a : b;
        }
    }
    return std::nullopt;
}

inline void find_F(const Point& p, const std::vector<Point>& S,
                   std::vector<Point>& F) {
    F.clear();
    assert(S.size() != 2);
    if (S.size() == 1 || point_in_convex(p, S)) {
        const auto& bbox = current_bbox_static();
        F.assign(bbox.begin(), bbox.end());
        return;
    }

    std::vector<int> tangent = find_tangent_idx(p, S);
    if (tangent.size() != 2) {
        return;
    }

    auto hit1 = ray_hit_bbox(p, S[tangent[0]]);
    auto hit2 = ray_hit_bbox(p, S[tangent[1]]);
    if (!hit1 || !hit2) {
        std::cerr << "Ray doesn't intersect with bounding box!\n";
        const auto& bbox = current_bbox_static();
        F.assign(bbox.begin(), bbox.end());
        return;
    }

    auto e1 = which_edge(hit1.value());
    auto e2 = which_edge(hit2.value());
    if (!e1 || !e2) {
        std::cerr << "Cannot determine which Bbox edge the ray intersect with.\n";
        const auto& bbox = current_bbox_static();
        F.assign(bbox.begin(), bbox.end());
        return;
    }

    int n = int(S.size());
    assert(n >= 3);
    assert(tangent[1] - tangent[0] - 1 >= 1 || tangent[0] + n - tangent[1] - 1 >= 1);

    F.reserve(n + 4);
    // Raw-double right_turn: sign((S[t1] - p) x (S[t2] - p))
    const double ax = CGAL::to_double(S[tangent[0]].x() - p.x()),
                  ay = CGAL::to_double(S[tangent[0]].y() - p.y());
    const double bx = CGAL::to_double(S[tangent[1]].x() - p.x()),
                  by = CGAL::to_double(S[tangent[1]].y() - p.y());
    const bool is_right_turn = (ax * by - ay * bx) < 0;
    if (is_right_turn) {
        std::copy(S.begin() + tangent[0], S.begin() + tangent[1] + 1, std::back_inserter(F));
        F.push_back(hit2.value());
        append_rect_pts(F, e2.value(), e1.value(), true);
        F.push_back(hit1.value());
    } else {
        std::copy(S.begin() + tangent[1], S.end(), std::back_inserter(F));
        std::copy(S.begin(), S.begin() + tangent[0] + 1, std::back_inserter(F));
        F.push_back(hit1.value());
        append_rect_pts(F, e1.value(), e2.value(), true);
        F.push_back(hit2.value());
    }

}

// ===========================================================================
//  Help & debug printing
// ===========================================================================

static void print_help() {
    std::cout << "Usage: simplify [options]\n"
              << "  --in <id>        Read input from data/taxi/<id>.txt (resolved absolutely)\n"
              << "  --out            Write output to data/<id>/original.txt & simplify.txt (resolved absolutely; requires --in <id>)\n"
              << "  --dist           After output, compute Frechet distance by invoking ./frechet (Julia wrapper) with --in <id> --path <simplify.txt>\n"
              << "  --gui            Show GUI viewer \n"
              << "  --keep           Do not clear F/G/S polygons in GUI between steps\n"
              << "  --no-simp        Do not show the simplified curve in GUI\n"
              << "  --labels         Show vertex labels (indices) in GUI\n"
              << "  -d <delta>       Override DELTA (default " << DELTA << ")\n"
              << "  -e <epsilon>     Override EPSILON (default " << EPSILON << ")\n"
              << "  -F/-G/-S         Debug polygon display modes\n"
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
//  Core algorithm
// ===========================================================================

inline bool segment_within_frechet(const std::vector<Point>& stream,
                                   int first, int last,
                                   const Point& a, const Point& b,
                                   double bound) {
    const double ax = CGAL::to_double(a.x());
    const double ay = CGAL::to_double(a.y());
    const double dx = CGAL::to_double(b.x()) - ax;
    const double dy = CGAL::to_double(b.y()) - ay;
    const double length_sq = dx * dx + dy * dy;
    const double safe_bound = bound - std::max(1.0, bound) * 1e-9;
    const double bound_sq = safe_bound * safe_bound;
    const double tolerance = std::numeric_limits<double>::epsilon() *
                             std::max(1.0, bound_sq) * 64.0;

    auto squared_distance = [](double x1, double y1, double x2, double y2) {
        const double x = x1 - x2;
        const double y = y1 - y2;
        return x * x + y * y;
    };

    const double first_x = CGAL::to_double(stream[first].x());
    const double first_y = CGAL::to_double(stream[first].y());
    const double last_x = CGAL::to_double(stream[last].x());
    const double last_y = CGAL::to_double(stream[last].y());
    if (squared_distance(first_x, first_y, ax, ay) > bound_sq + tolerance ||
        squared_distance(last_x, last_y, ax + dx, ay + dy) > bound_sq + tolerance) {
        return false;
    }

    if (length_sq <= std::numeric_limits<double>::epsilon()) {
        for (int i = first; i <= last; ++i) {
            const double px = CGAL::to_double(stream[i].x());
            const double py = CGAL::to_double(stream[i].y());
            if (squared_distance(px, py, ax, ay) > bound_sq + tolerance) return false;
        }
        return true;
    }

    double reachable = 0.0;
    for (int i = first; i <= last; ++i) {
        const double px = CGAL::to_double(stream[i].x());
        const double py = CGAL::to_double(stream[i].y());
        const double rx = ax - px;
        const double ry = ay - py;
        const double linear = 2.0 * (rx * dx + ry * dy);
        const double constant = rx * rx + ry * ry - bound_sq;
        const double discriminant = linear * linear - 4.0 * length_sq * constant;
        if (discriminant < -tolerance) return false;

        const double root = std::sqrt(std::max(0.0, discriminant));
        const double lo = std::max(0.0, (-linear - root) / (2.0 * length_sq));
        const double hi = std::min(1.0, (-linear + root) / (2.0 * length_sq));
        if (lo > hi + 1e-12 || reachable > hi + 1e-12) return false;
        reachable = std::max(reachable, lo);
    }
    return true;
}

int get_longest_stab(const std::vector<Point>& stream, int cur,
                     std::vector<Point>& simplified,
                     double EPSILON, double DELTA,
                     MultiViewer* viewer = nullptr) {
    TIMER("get_longest_stab");
    const Point& p0 = stream[cur];
    int p0cur = cur;
    std::vector<Point> P;
    {
        TIMER("get_points_from_grid");
        P = get_points_from_grid(p0, EPSILON, DELTA);
    }
    if (viewer) {
        TIMER("gui_update");
        viewer->markP0(p0);
        viewer->addOriginalPoint(p0);
    }
    std::array<Point, 2> buffer = {p0, p0};
    const int Pn = (int)P.size();
    std::vector<std::vector<Point>> S(Pn, std::vector<Point>{p0});
    int dead_cnt = 0;
    std::vector<int> dead(Pn);
    std::vector<std::vector<Point>> new_S(Pn);
    std::vector<char> new_dead(Pn, 0);
    std::vector<std::vector<Point>> F_work(Pn);
    std::vector<char> F_cached(Pn, 0);
    std::vector<Point> Gi;

    cur++;
    while (cur < int(stream.size())) {
        const Point& pi  = stream[cur];
        {
            TIMER("get_conv_from_grid");
            Gi = get_conv_from_grid(pi, EPSILON, DELTA);
        }
        std::fill(new_dead.begin(), new_dead.end(), 0);
        if (viewer) std::fill(F_cached.begin(), F_cached.end(), 0);
        for (int i = 0; i < Pn; ++i) {
            if (dead[i]) continue;
            find_F(P[i], S[i], F_work[i]);
            if (!intersect(F_work[i], Gi, new_S[i])) {
                new_dead[i] = 1;
            } else if (viewer) {
                F_cached[i] = 1;
            }
        }

        bool verified_candidate = false;
        for (int i = Pn - 1; i >= 0 && !verified_candidate; --i) {
            if (dead[i] || new_dead[i]) continue;
            for (const Point& endpoint : new_S[i]) {
                if (segment_within_frechet(stream, p0cur, cur, P[i], endpoint,
                                           (1.0 + EPSILON) * DELTA)) {
                    buffer[0] = P[i];
                    buffer[1] = endpoint;
                    verified_candidate = true;
                    break;
                }
            }
        }
        if (!verified_candidate) break;

        // GUI updates (sequential, main thread only).
        if (viewer) {
            TIMER("gui_update");
            const QColor stepColors[] = {
                Qt::red, Qt::blue, Qt::green, Qt::magenta, Qt::cyan
            };
            QColor c = stepColors[cur % 6];
            for (int i = 0; i < Pn; ++i) {
                if (dead[i] || new_dead[i]) continue;
                if (showF && F_cached[i]) {
                    Polygon F_poly(F_work[i].begin(), F_work[i].end());
                    viewer->addPolygon(F_poly, c);
                }
                if (showG) {
                    Polygon Gi_poly(Gi.begin(), Gi.end());
                    viewer->addPolygon(Gi_poly, c);
                }
                if (showS) {
                    Polygon S_poly(S[i].begin(), S[i].end());
                    viewer->addPolygon(S_poly, c);
                }
                break;
            }
        }
        for (int i = 0; i < Pn; ++i) {
            if (new_dead[i]) {
                dead[i] = true;
                dead_cnt++;
            }
        }
        if (dead_cnt == int(P.size())) {
            break;
        }
        if (viewer) {
            TIMER("gui_update");
            viewer->addOriginalPoint(pi);
            viewer->markPi(pi);
        }
        {
            TIMER("update_S");
            for (int i = 0; i < int(P.size()); i++) {
                if (dead[i]) continue;
                S[i].swap(new_S[i]);
            }
        }
        cur++;
        if (viewer) {
            TIMER("gui_update");
            viewer_process_events();
        }
    }

    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    if (viewer) {
        TIMER("gui_update");
        viewer->addSimplifiedPoint(buffer[0]);
        viewer->addSimplifiedPoint(buffer[1]);
        viewer->clearPolygons();
        viewer->clearMarkedP0();
        viewer->clearMarkedPi();
        viewer_process_events();
    }
    return cur;
}

std::vector<Point> simplify(const std::vector<Point>& stream,
                            double EPSILON, double DELTA,
                            MultiViewer* viewer = nullptr) {
    TIMER("simplify");
    std::vector<Point> simplified;
    std::cout << "Simplifying...\n";
    int cur = 0;
    while (cur != int(stream.size())) {
        cur = get_longest_stab(stream, cur, simplified, EPSILON, DELTA, viewer);
    }
    std::cout << "Expected Frechet distance: " << std::sqrt(expected_frechet_ref()) << '\n';
    std::cout << "The original stream of size " << stream.size() << " is simplified to " << simplified.size() << " points.\n";
    return simplified;
}

// ===========================================================================
//  main
// ===========================================================================

// TODO: clean up this mess of flags
int main(int argc, char** argv) {
    int test_case_no = -1;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i],"-F") == 0) showF = true;
        else if (strcmp(argv[i],"-S") == 0) showS = true;
        else if (strcmp(argv[i],"-G") == 0) showG = true;
        else if (strcmp(argv[i],"--gui") == 0) gui_flag = true;
        else if (strcmp(argv[i],"--out") == 0) out_flag = true;
        else if (strcmp(argv[i],"--dist") == 0) dist_flag = true;
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

    QApplication app(argc, argv);
    MultiViewer viewer;
    MultiViewer* vptr = nullptr;
    if (gui_flag) {
        vptr = &viewer;
        viewer.setParameters(DELTA, EPSILON);
        viewer.setShowLabels(false);
        viewer.show();
    }

    std::vector<Point> simplified;
    {
        TIMER("total");
        simplified = simplify(stream, EPSILON, DELTA, vptr);
        stream = std::move(simplified);
    }
    print_timing_summary();

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

    int gui_result = 0;
    if (gui_flag) {
        gui_result = app.exec();
    }

    if (dist_flag && test_case_no != -1) {
        std::filesystem::path frechet_path = repo_root / "scripts" / "frechet";
        std::string cmd1 = std::string("\"") + frechet_path.string() + "\" --in " + std::to_string(test_case_no) + " --path \"" + (repo_root / "data" / std::to_string(test_case_no) / "simplify.txt").string() + "\"";
        int rc = std::system(cmd1.c_str());
    }

    return gui_result;
}
