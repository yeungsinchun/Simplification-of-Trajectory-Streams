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
#include <exception>
#include <fstream>
#include <filesystem>
#include <QApplication>
#include <thread>
#include <chrono>
#include "drawing.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = Kernel::Point_2;
using Segment = Kernel::Segment_2;
using Ray = Kernel::Ray_2;
using Bbox = CGAL::Bbox_2;
using Rect = CGAL::Iso_rectangle_2<Kernel>;
using Polygon = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes =  CGAL::Polygon_with_holes_2<Kernel>;

constexpr double TOL = 1e-6;
constexpr double EPSILON = 0.5;
constexpr double SQRT2 =
    1.41421356237; // sqrt in STL does not have constexpr version !?

// Bounding square [BMIN, BMAX]^2
// Compute from input data in main() so the display fits the points.
double BMIN;
double BMAX;
constexpr double DELTA = 0.2;
constexpr double GRID = EPSILON * DELTA / (2 * SQRT2);

/*
void draw_points_svg(const std::vector<Point> &pts,
                     const std::string &filename = "../data/points.svg",
                     int W = 800, int H = 800, double margin = 20.0,
                     double point_radius = 3.5) {
    double minx = 1e300, miny = 1e300, maxx = -1e300, maxy = -1e300;
    for (const auto &p : pts) {
        double x = CGAL::to_double(p.x()), y = CGAL::to_double(p.y());
        minx = std::min(minx, x);
        miny = std::min(miny, y);
        maxx = std::max(maxx, x);
        maxy = std::max(maxy, y);
    }
    if (pts.empty() || !(minx <= maxx)) {
        minx = -1.0;
        miny = -1.0;
        maxx = 1.0;
        maxy = 1.0;
    }

    double dx = maxx - minx;
    double dy = maxy - miny;
    if (dx == 0)
        dx = 1.0;
    if (dy == 0)
        dy = 1.0;
    double scale = std::min((W - 2 * margin) / dx, (H - 2 * margin) / dy);

    auto mapx = [&](double x) { return margin + (x - minx) * scale; };
    auto mapy = [&](double y) { return H - (margin + (y - miny) * scale); };

    std::filesystem::create_directories(
        std::filesystem::path(filename).parent_path());
    std::ofstream os(filename);
    if (!os)
        return;

    os << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
       << "width=\"" << W << "\" height=\"" << H << "\" viewBox=\"0 0 " << W
       << " " << H << "\">\n";

    // background
    os << "<rect x=\"0\" y=\"0\" width=\"" << W << "\" height=\"" << H
       << "\" fill=\"white\" stroke=\"none\" />\n";

    // points as filled circles
    os.setf(std::ios::fixed);
    os << std::setprecision(6);
    for (const auto &p : pts) {
        double x = mapx(CGAL::to_double(p.x()));
        double y = mapy(CGAL::to_double(p.y()));
        os << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\""
           << point_radius << "\" fill=\"black\" stroke=\"none\" />\n";
    }

    os << "</svg>\n";
    os.close();
}
*/

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

// Write an SVG snapshot visualizing p (green dot), polygon S (blue outline),
// and polygon F (red outline, translucent fill). Overwrites the same file
// each call so the last write reflects the most recent state before a crash.
static void write_F_svg(const Point& p,
                        const std::vector<Point>& S,
                        const std::vector<Point>& F,
                        const std::string& filename = "../data/F.svg",
                        int W = 800, int H = 800, double margin = 20.0)
{
    // Gather bounds from p, S, and F
    double minx = std::numeric_limits<double>::infinity();
    double miny = std::numeric_limits<double>::infinity();
    double maxx = -std::numeric_limits<double>::infinity();
    double maxy = -std::numeric_limits<double>::infinity();
    auto consider = [&](const Point& q){
        double x = CGAL::to_double(q.x());
        double y = CGAL::to_double(q.y());
        minx = std::min(minx, x); miny = std::min(miny, y);
        maxx = std::max(maxx, x); maxy = std::max(maxy, y);
    };
    consider(p);
    for (const auto& q : S) consider(q);
    for (const auto& q : F) consider(q);
    if (!std::isfinite(minx) || !std::isfinite(miny) || !std::isfinite(maxx) || !std::isfinite(maxy)) {
        minx = -1.0; miny = -1.0; maxx = 1.0; maxy = 1.0;
    }
    double dx = maxx - minx, dy = maxy - miny;
    if (dx <= 0) dx = 1.0;
    if (dy <= 0) dy = 1.0;
    // add 5% padding
    double pad_x = 0.05 * dx, pad_y = 0.05 * dy;
    minx -= pad_x; maxx += pad_x; miny -= pad_y; maxy += pad_y;
    dx = maxx - minx; dy = maxy - miny;
    double scale = std::min((W - 2 * margin) / dx, (H - 2 * margin) / dy);
    auto mapx = [&](double x) { return margin + (x - minx) * scale; };
    auto mapy = [&](double y) { return H - (margin + (y - miny) * scale); };

    std::filesystem::create_directories(std::filesystem::path(filename).parent_path());
    std::ofstream os(filename);
    if (!os) return;
    os.setf(std::ios::fixed);
    os << std::setprecision(6);
    os << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W
       << "\" height=\"" << H << "\" viewBox=\"0 0 " << W << " " << H
       << "\">\n";
    os << "<rect x=\"0\" y=\"0\" width=\"" << W << "\" height=\"" << H
       << "\" fill=\"white\" stroke=\"none\"/>\n";

    auto draw_poly = [&](const std::vector<Point>& poly,
                         const char* stroke, const char* fill,
                         double stroke_width, double fill_opacity) {
        if (poly.empty()) return;
        os << "<polygon points=\"";
        for (const auto& q : poly) {
            os << mapx(CGAL::to_double(q.x())) << "," << mapy(CGAL::to_double(q.y())) << " ";
        }
        os << "\" stroke=\"" << stroke << "\" stroke-width=\"" << stroke_width
           << "\" fill=\"" << fill << "\" fill-opacity=\"" << fill_opacity << "\"/>\n";
    };

    // S: blue outline, no fill
    draw_poly(S, "#1f77b4", "none", 2.0, 0.0);
    // F: green outline, no fill
    draw_poly(F, "#3f1398ff", "none", 2.0, 0.0);
    // p: green point
    os << "<circle cx=\"" << mapx(CGAL::to_double(p.x()))
       << "\" cy=\"" << mapy(CGAL::to_double(p.y()))
       << "\" r=\"3\" fill=\"#2ca02c\" stroke=\"none\"/>\n";

    os << "</svg>\n";
}

// REMOVE LATER: print simple/orientation/area
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

// REMOVE LATER: Sanitize polygon: drop duplicate and collinear vertices, enforce CCW
static void sanitize_polygon(Polygon &P, const char* name = "poly") {
    auto safe_is_simple = [](const Polygon& poly) -> bool {
        try { return poly.is_simple(); } catch (...) { return false; }
    };
    auto safe_is_clockwise = [](const Polygon& poly) -> bool {
        try { return poly.is_clockwise_oriented(); } catch (...) { return false; }
    };
    // gather vertices, remove consecutive duplicates
    std::vector<Point> pts;
    pts.reserve(P.size());
    for (auto it = P.vertices_begin(); it != P.vertices_end(); ++it) {
        if (pts.empty() || *it != pts.back()) pts.push_back(*it);
    }
    if (pts.size() > 1 && pts.front() == pts.back()) pts.pop_back();

    // remove collinear triples
    auto is_collinear = [](const Point &a, const Point &b, const Point &c) {
        return CGAL::orientation(a, b, c) == CGAL::COLLINEAR;
    };
    std::vector<Point> cleaned;
    cleaned.reserve(pts.size());
    for (size_t i = 0; i < pts.size(); ++i) {
        if (pts.size() <= 2) { cleaned.push_back(pts[i]); continue; }
        const Point &a = pts[(i + pts.size() - 1) % pts.size()];
        const Point &b = pts[i];
        const Point &c = pts[(i + 1) % pts.size()];
        if (!is_collinear(a, b, c)) cleaned.push_back(b);
    }

    Polygon Q(cleaned.begin(), cleaned.end());
    if (Q.size() >= 3 && safe_is_clockwise(Q)) Q.reverse_orientation();
    if (!safe_is_simple(Q)) {
        std::cerr << "sanitize_polygon: still not simple for '" << name << "'\n";
        print_polygon(Q);
    }
    P = std::move(Q);
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
    const double r = (1.0 + EPSILON / 2.0) * DELTA;
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
    const double r = (1.0 + EPSILON / 2.0) * DELTA;
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
    // CGAL_precondition(Polygon(S.begin(), S.end()).is_counterclockwise_oriented());
    assert(S.size() != 2); // this could happen?
    if (S.size() == 1 || point_in_convex(p, S)) {
        auto F = current_bbox();
        write_F_svg(p, S, F);
        std::cerr << "Case " << 1 << std::endl;
        return F;
    }
    std::vector<int> tangent = find_tangent_idx(p, S);
    // std::cerr << "tangent size: " << tangent.size() << std::endl;
    assert(tangent.size() == 2);

    std::optional<Point> hit1 = intersect_ray_with_rect(p, S[tangent[0]]);
    std::optional<Point> hit2 = intersect_ray_with_rect(p, S[tangent[1]]);
    if (!hit1 || !hit2) {
        std::cerr << "Ray doesn't intersect with bounding box!\n";
        auto F = current_bbox();
        write_F_svg(p, S, F);
        std::cerr << "Case " << 2 << std::endl;
        return F;
    }
    std::optional<Bbox_edge> e1 = which_edge(hit1.value());
    std::optional<Bbox_edge> e2 = which_edge(hit2.value());

    int n = int(S.size());
    assert(n >= 3);
    // std::cerr << n << ' ' << tangent[1] << ' ' << tangent[0] << std::endl;
    // This is false
    assert(tangent[1] - tangent[0] - 1 >= 1 || tangent[0] + n - tangent[1] - 1 >= 1);
    std::vector<Point> F;
    if (CGAL::squared_distance(S[(tangent[0] + 1) % n], p) <
        CGAL::squared_distance(S[(tangent[1] + 1) % n], p)) {
        // [i..j] (inclusive)
        std::cerr << "First\n";
        std::copy(S.begin() + tangent[0], S.begin() + tangent[1] + 1,
                  std::back_inserter(F));
        // walk from hit2 to hit1 ccw to close
        F.push_back(hit2.value());
        append_rect_pts(F, e2.value(), e1.value(), true);
        F.push_back(hit1.value());
    } else {
        // reverse([j..n-1] + [0..i]) (inclusive)
        std::cerr << "Second\n";
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
        write_F_svg(p, S, F);
        std::cerr << "Case " << 3 << std::endl;
        return F;
    }
    write_F_svg(p, S, F);
    for (size_t i = 0; i < F.size(); i++) {
        std::cerr << "(" << F[i].x() << ',' << F[i].y() << ") ";
    }
    std::cerr << '\n';
    std::cerr << hit1.value().x() << ' ' << hit1.value().y() << '\n';
    std::cerr << hit2.value().x() << ' ' << hit2.value().y() << '\n';
    std::cerr << "Case " << 4 << std::endl;
    return F;
}

int get_longest_stab(const std::vector<Point> &stream, int cur,
                     std::vector<Point> &simplified, MultiViewer& viewer) {
    const Point& p0 = stream[cur];
    viewer.addOriginalPoint(p0);
    cur++; // What if stream ends here?
    std::vector<Point> P = get_points_from_grid(p0);
    std::array<Point, 2> buffer = {p0};
    std::vector<std::vector<Point>> S(P.size());
    for (int i = 0; i < int(P.size()); i++) {
        S[i].push_back(P[i]);
    }
    int dead_cnt = 0;
    std::vector<int> dead(P.size());
    while (cur < int(stream.size())) {
        if (dead_cnt == int(P.size())) {
            break;
        }
        const Point& pi  = stream[cur];
        viewer.addOriginalPoint(pi);
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
            if (i == 0) {
                // viewer.addPolygon(F_poly);
                // viewer.addPolygon(Gi_poly);
            }
            /*
            std::cerr << "F\n";
            viewer.addPolygon(F_poly);
            print_polygon(F_poly);
            poly_test(F_poly);
            std::cerr << "Gi\n";
            viewer.addPolygon(Gi_poly);
            print_polygon(Gi_poly);
            poly_test(Gi_poly);
            std::cerr << "done\n";
            std::cerr.flush();
            */
            // Validate and sanitize polygons before boolean ops
            // Only print minimal info before sanitize to avoid triggering
            // CGAL preconditions on invalid polygons.
            print_poly_info(S_poly, "S_poly(before)");
            print_poly_info(F_poly, "F_poly(before)");
            print_poly_info(Gi_poly, "Gi_poly(before)");
            // sanitize_polygon(F_poly, "F_poly");
            // sanitize_polygon(Gi_poly, "Gi_poly");
            print_poly_info(F_poly, "F_poly(after)");
            print_poly_info(Gi_poly, "Gi_poly(after)");
            if (F_poly.is_clockwise_oriented()) F_poly.reverse_orientation();
            if (Gi_poly.is_clockwise_oriented()) Gi_poly.reverse_orientation();
            // Early reject invalid/degenerate polygons to avoid CGAL preconditions
            auto safe_is_simple = [](const Polygon& poly) -> bool {
                try { return poly.is_simple(); } catch (...) { return false; }
            };
            if (F_poly.size() < 3 || !safe_is_simple(F_poly)) {
                std::cerr << "Rejecting F_poly (degenerate/invalid) at i=" << i << "\n";
                dead[i] = true; dead_cnt++; continue;
            }
            if (Gi_poly.size() < 3 || !safe_is_simple(Gi_poly)) {
                std::cerr << "Rejecting Gi_poly (degenerate/invalid) at i=" << i << "\n";
                dead[i] = true; dead_cnt++; continue;
            }
            // Log the exact source location of the next intersection call.
            std::cerr << "[DBG] Next line will call CGAL::intersection at "
                      << __FILE__ << ":" << (__LINE__ + 1) << "\n";
            try {
                CGAL::intersection(F_poly, Gi_poly, back_inserter(new_S[i]));
            } catch (const std::exception &e) {
                std::cerr << "CGAL::intersection threw: " << e.what() << "\n";
                std::cerr << "At grid idx i=" << i << "\n";
                print_polygon(F_poly);
                print_polygon(Gi_poly);
                throw;
            }
            if (new_S[i].size() == 0) {
                dead[i] = true;
                dead_cnt++;
                continue;
            }
            assert(new_S[i].size() == 1);
            buffer[1] = *new_S[i].begin()->outer_boundary().vertices_begin();
            /*
            std::cout << "throw\n";
            std::cout << new_S[i].size() << '\n';
            */
        }
        for (int i = 0; i < int(P.size()); i++) {
            if (dead[i]) continue;
            S[i].clear();
            std::copy(new_S.begin()->begin()->outer_boundary().vertices_begin(),
                      new_S.begin()->begin()->outer_boundary().vertices_end(),
                      std::back_inserter(S[i]));
        }
        /*
        Polygon S_poly(S[0].begin(), S[0].end());
        std::cout << "S\n";
        viewer.addPolygon(S_poly);
        poly_test(S_poly);
        */
        viewer_process_events();
        // std::this_thread::sleep_for(std::chrono::seconds(1));
        cur++;
    }
    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    viewer.addSimplifiedPoint(buffer[0]);
    viewer.addSimplifiedPoint(buffer[1]);
    viewer_process_events();
    return cur;
}

std::vector<Point> simplify(const std::vector<Point> &stream, MultiViewer& viewer) {
    std::vector<Point> simplified;
    int cur = 0;
    while (cur != int(stream.size())) {
        cur = get_longest_stab(stream, cur, simplified, viewer);
    }
    return simplified;
}

int main(int argc, char** argv) {
    int N;
    std::cin >> N;
    std::vector<Point> stream(N);
    for (int i = 0; i < N; i++) {
        std::cin >> stream[i];
    }
    // Compute bounding square
    if (!stream.empty()) {
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
        double dx = maxx - minx;
        double dy = maxy - miny;
        double span = std::max(dx, dy);
        double pad = span * 0.05;
        double overall_min = std::min(minx, miny) - pad;
        double overall_max = std::max(maxx, maxy) + pad;
        BMIN = overall_min;
        BMAX = overall_max;
    }
    QApplication app(argc, argv);
    MultiViewer viewer;
    viewer.show();
     simplify(stream, viewer);
    return app.exec();
}