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
constexpr double DELTA = 1e4;
constexpr double SQRT2 =
    1.41421356237; // sqrt in STL does not have constexpr version !?
constexpr double GRID = EPSILON * DELTA / (2 * SQRT2);

constexpr double BMIN = -2e6;
constexpr double BMAX = 2e6;

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
    BOTTOM, // 0
    RIGHT,  // 1
    TOP,    // 2
    LEFT,   // 3
};

const std::array<Point, 4> bbox_corner = {
    Point(BMIN, BMIN), Point(BMAX, BMIN), Point(BMAX, BMAX),
    Point(BMIN, BMAX)};

const std::vector<Point> bbox =
    std::vector<Point>(bbox_corner.begin(), bbox_corner.end());

std::optional<Bbox_edge> which_edge(const Point &s) {
    double x = CGAL::to_double(s.x()), y = CGAL::to_double(s.y());
    if (std::abs(y - BMIN) < TOL)
        return Bbox_edge::BOTTOM;
    if (std::abs(x - BMAX) < TOL)
        return Bbox_edge::RIGHT;
    if (std::abs(y - BMAX) < TOL)
        return Bbox_edge::TOP;
    if (std::abs(x - BMIN) < TOL)
        return Bbox_edge::LEFT;
    return std::nullopt;
}

void append_rect_pts(std::vector<Point> &out, Bbox_edge e1, Bbox_edge e2,
                     bool ccw) {
    if (ccw) { // e.g., if a_edge = bottom (0), b_edge = top (2): need to append
               // bottom right (1) and top right (2)
        int begin = (static_cast<int>(e1) + 1) % 4;
        int end = (static_cast<int>(e2) + 1) % 4;
        while (begin != end) {
            out.push_back(bbox_corner[begin]);
            begin = (begin + 1) % 4;
        }
    } else {
        int begin = (static_cast<int>(e1) + 3) % 4;
        int end = (static_cast<int>(e2) + 3) % 4;
        while (begin != end) {
            out.push_back(bbox_corner[begin]);
            begin = (begin + 3) % 4;
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
    auto to_string = [=](CGAL::Epeck::Orientation e) {
        if (e == CGAL::COLLINEAR) {
            return "collinear";
        } else if (e == CGAL::LEFT_TURN) {
            return "left";
        } else if (e == CGAL::RIGHT_TURN) {
            return "right";
        } else {
            return "what??";
        }
    };
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
    // Precondition: S should be counterclockwise
    assert(S.size() != 2); // this could happen?
    if (S.size() == 1 || point_in_convex(p, S)) return bbox;
    std::vector<int> tangent = find_tangent_idx(p, S);
    // std::cerr << "tangent size: " << tangent.size() << std::endl;
    assert(tangent.size() == 2);

    std::optional<Point> hit1 = intersect_ray_with_rect(p, S[tangent[0]]);
    std::optional<Point> hit2 = intersect_ray_with_rect(p, S[tangent[1]]);
    if (!hit1 || !hit2) {
        std::cerr << "Ray doesn't intersect with bounding box!\n";
        return bbox;
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
        std::copy(S.begin() + tangent[0], S.begin() + tangent[1] + 1,
                  std::back_inserter(F));
        // walk from hit2 to hit1 ccw to close
        F.push_back(hit2.value());
        append_rect_pts(F, e2.value(), e1.value(), true);
        F.push_back(hit1.value());
    } else {
        // reverse([j..n-1] + [0..i]) (inclusive)
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
        return bbox;
    }
    return F;
}

void poly_test(Polygon p) {
    // std::cerr << p.is_simple() << ' ' << p.is_counterclockwise_oriented() << ' ' << p.is_convex() << '\n';
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
            CGAL::intersection(F_poly, Gi_poly, back_inserter(new_S[i]));
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
        std::cout << "[cur]: " << cur << '\n';
    }
    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    viewer.addSimplifiedPoint(buffer[0]);
    viewer.addSimplifiedPoint(buffer[1]);
    viewer_process_events();
    std::cout << "simplified!\n";
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
    QApplication app(argc, argv);
    MultiViewer viewer;
    viewer.show();
    simplify(stream, viewer);
    return app.exec();
}