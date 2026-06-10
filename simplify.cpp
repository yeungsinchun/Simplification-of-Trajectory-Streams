#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/intersections.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <fstream>
#include <cstdlib>
#include <filesystem>
#include <cstring>
#include <chrono>
#include <array>
#include <iomanip>
#include <algorithm>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <QApplication>
#include "drawing.h"
#include "geometry_kernel.h"

using Segment = Kernel::Segment_2;
using Ray = Kernel::Ray_2;
using Bbox = CGAL::Bbox_2;
using Rect = CGAL::Iso_rectangle_2<Kernel>;
using Polygon = CGAL::Polygon_2<Kernel>;
using BgPoint = boost::geometry::model::d2::point_xy<double>;
using BgPolygon = boost::geometry::model::polygon<BgPoint, true, true>;

bool showF = false; // true if -F is passed
bool showG = false; // true if -G is passed
bool showS = false; // true if -S is passed
bool keepPolygons = false; // true if --keep is passed
bool showSimp = true; // false if --no-simp is passed
bool showLabels = false; // true if --labels is passed

constexpr double TOL = 1e-6;
constexpr double SQRT2 =
    1.41421356237; // sqrt in STL does not have constexpr version !?

// Bounding square [BMIN, BMAX]^2
extern const double BMIN = -10000;
extern const double BMAX = 10000;
// Bounding data points (convex hull should not exceed bounding box)
const double DATAMIN = -8000;
const double DATAMAX = 8000;
// TODO: DELTA should not be too high that convex hull goes out of bounding square, which may cause the program to crash
double DELTA = 200;
double EPSILON = 0.6;

double GRID_val() { return EPSILON * DELTA / (2 * SQRT2); }

double R_val() { return (1.0 + EPSILON / 2.0) * DELTA; }

struct TimingStat {
    const char* name;
    long long total_ns = 0;
    long long calls = 0;
};

enum class TimedOp : std::size_t {
    GetPointsFromGrid = 0,
    FindF,
    PointInConvex,
    FindTangentIdx,
    CgalIntersection,
    WhichEdge,
    AppendRectPts,
    GetConvFromGrid,
    BoostPolygonIntersection,
    Count
};

static std::array<TimingStat, static_cast<std::size_t>(TimedOp::Count)> g_timing_stats{{
    {"get_points_from_grid"},
    {"find_F"},
    {"point_in_convex"},
    {"find_tangent_idx"},
    {"CGAL::intersection"},
    {"which_edge"},
    {"append_rect_pts"},
    {"get_conv_from_grid"},
    {"boost::geometry::intersection"},
}};

struct ScopedTimer {
    TimedOp op;
    std::chrono::steady_clock::time_point start;

    explicit ScopedTimer(TimedOp timed_op)
        : op(timed_op), start(std::chrono::steady_clock::now()) {}

    ~ScopedTimer() {
        const auto end = std::chrono::steady_clock::now();
        auto& stat = g_timing_stats[static_cast<std::size_t>(op)];
        stat.total_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        stat.calls += 1;
    }
};

static void print_timing_report() {
    long long grand_total_ns = 0;
    for (const auto& stat : g_timing_stats) {
        grand_total_ns += stat.total_ns;
    }

    std::cout << "Timing report:\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::left << std::setw(28) << "Operation"
              << std::right << std::setw(14) << "Total (ms)"
              << std::setw(10) << "% Time"
              << std::setw(12) << "Calls"
              << std::setw(14) << "Avg (ms)"
              << '\n';
    std::cout << std::string(78, '-') << '\n';

    for (const auto& stat : g_timing_stats) {
        const double total_ms = static_cast<double>(stat.total_ns) / 1'000'000.0;
        const double avg_ms = stat.calls == 0 ? 0.0 : total_ms / static_cast<double>(stat.calls);
        const double pct = grand_total_ns == 0 ? 0.0 : (100.0 * static_cast<double>(stat.total_ns) / static_cast<double>(grand_total_ns));

        std::cout << std::left << std::setw(28) << stat.name
                  << std::right << std::setw(14) << total_ms
                  << std::setw(10) << pct
                  << std::setw(12) << stat.calls
                  << std::setw(14) << avg_ms
                  << '\n';
    }

    std::cout << std::string(78, '-') << '\n';
    std::cout << std::left << std::setw(28) << "Total"
              << std::right << std::setw(14) << (static_cast<double>(grand_total_ns) / 1'000'000.0)
              << std::setw(10) << 100.0
              << std::setw(12) << "-"
              << std::setw(14) << "-"
              << '\n';
}

static void reset_timing_stats() {
    for (auto& stat : g_timing_stats) {
        stat.total_ns = 0;
        stat.calls = 0;
    }
}

static BgPolygon to_bg_polygon(const std::vector<Point>& points) {
    BgPolygon polygon;
    auto& outer = polygon.outer();
    outer.reserve(points.size() + 1);
    for (const auto& point : points) {
        outer.emplace_back(CGAL::to_double(point.x()), CGAL::to_double(point.y()));
    }
    if (!points.empty()) {
        outer.emplace_back(CGAL::to_double(points.front().x()), CGAL::to_double(points.front().y()));
    }
    boost::geometry::correct(polygon);
    return polygon;
}

static std::vector<Point> from_bg_polygon(const BgPolygon& polygon) {
    std::vector<Point> points;
    const auto& outer = polygon.outer();
    if (outer.size() <= 1) {
        return points;
    }

    points.reserve(outer.size() - 1);
    for (std::size_t i = 0; i + 1 < outer.size(); ++i) {
        points.emplace_back(outer[i].x(), outer[i].y());
    }

    return points;
}

static bool points_equal_eps(const Point& lhs, const Point& rhs) {
    return CGAL::to_double(CGAL::squared_distance(lhs, rhs)) < TOL * TOL;
}

static double cross2d(const Point& a, const Point& b, const Point& c) {
    const double abx = CGAL::to_double(b.x()) - CGAL::to_double(a.x());
    const double aby = CGAL::to_double(b.y()) - CGAL::to_double(a.y());
    const double bcx = CGAL::to_double(c.x()) - CGAL::to_double(b.x());
    const double bcy = CGAL::to_double(c.y()) - CGAL::to_double(b.y());
    return abx * bcy - aby * bcx;
}

static std::vector<Point> normalize_convex_polygon(std::vector<Point> points) {
    if (points.empty()) {
        return {};
    }

    std::vector<Point> normalized;
    normalized.reserve(points.size());
    for (const auto& point : points) {
        if (normalized.empty() || !points_equal_eps(normalized.back(), point)) {
            normalized.push_back(point);
        }
    }

    if (normalized.size() >= 2 && points_equal_eps(normalized.front(), normalized.back())) {
        normalized.pop_back();
    }

    if (normalized.size() < 3) {
        return {};
    }

    BgPolygon polygon = to_bg_polygon(normalized);
    boost::geometry::unique(polygon);
    boost::geometry::correct(polygon);

    std::vector<Point> cleaned = from_bg_polygon(polygon);
    if (cleaned.size() < 3) {
        return {};
    }

    std::vector<Point> result;
    result.reserve(cleaned.size());
    for (const auto& point : cleaned) {
        while (result.size() >= 2) {
            const Point& a = result[result.size() - 2];
            const Point& b = result[result.size() - 1];
            if (std::abs(cross2d(a, b, point)) < TOL) {
                result.pop_back();
                continue;
            }
            break;
        }
        if (result.empty() || !points_equal_eps(result.back(), point)) {
            result.push_back(point);
        }
    }

    while (result.size() >= 2 && points_equal_eps(result.front(), result.back())) {
        result.pop_back();
    }
    while (result.size() >= 3) {
        const Point& a = result[result.size() - 2];
        const Point& b = result[result.size() - 1];
        const Point& c = result.front();
        if (std::abs(cross2d(a, b, c)) < TOL) {
            result.pop_back();
            continue;
        }
        break;
    }
    while (result.size() >= 3) {
        const Point& a = result.back();
        const Point& b = result.front();
        const Point& c = result[1];
        if (std::abs(cross2d(a, b, c)) < TOL) {
            result.erase(result.begin());
            continue;
        }
        break;
    }

    if (result.size() < 3) {
        return {};
    }

    Polygon polygon_check(result.begin(), result.end());
    if (!polygon_check.is_simple()) {
        return {};
    }
    if (polygon_check.orientation() == CGAL::CLOCKWISE) {
        std::reverse(result.begin(), result.end());
        polygon_check = Polygon(result.begin(), result.end());
    }
    if (polygon_check.orientation() != CGAL::COUNTERCLOCKWISE) {
        return {};
    }

    return result;
}

static std::vector<Point> intersect_convex_polygons(const std::vector<Point>& lhs,
                                                    const std::vector<Point>& rhs) {
    ScopedTimer timer(TimedOp::BoostPolygonIntersection);
    BgPolygon lhs_polygon = to_bg_polygon(lhs);
    BgPolygon rhs_polygon = to_bg_polygon(rhs);
    std::vector<BgPolygon> output;
    boost::geometry::intersection(lhs_polygon, rhs_polygon, output);

    if (output.empty()) {
        return {};
    }

    auto best = std::max_element(output.begin(), output.end(), [](const BgPolygon& a, const BgPolygon& b) {
        return boost::geometry::area(a) < boost::geometry::area(b);
    });

    return normalize_convex_polygon(from_bg_polygon(*best));
}

static void print_help() {
    std::cout << "Usage: simplify [options]\n"
              << "  --in <id>        Read input from data/taxi/<id>.txt (resolved absolutely)\n"
              << "  --out            Write output to data/taxi_simplified/<id>/original.txt & simplified.txt (resolved absolutely; requires --in <id>)\n"
              << "  --dist           After output, compute Frechet distance via ./frechet -n <id>\n"
              << "  --gui            Show GUI viewer \n"
              << "  --keep           Do not clear F/G/S polygons in GUI between steps\n"
              << "  --no-simp        Do not show the simplified curve in GUI\n"
              << "  --labels         Show vertex labels (indices) in GUI\n"
              << "  -d <delta>       Override DELTA (default " << DELTA << ")\n"
              << "  -e <epsilon>     Override EPSILON (default " << EPSILON << ")\n"
              << "  -F/-G/-S         Debug polygon display modes\n"
              << "  -h               Show this help and exit\n"
              << "\n"
              << "Shorthand: simplify <id> [flags] is equivalent to '--in <id> --out [flags]'\n";
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

// Print simple/orientation/area
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

bool point_in_convex(const Point &p, const std::vector<Point> &poly, bool ccw = true) {
    ScopedTimer timer(TimedOp::PointInConvex);
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
    ScopedTimer timer(TimedOp::WhichEdge);
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
    ScopedTimer timer(TimedOp::AppendRectPts);
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
    ScopedTimer timer(TimedOp::GetConvFromGrid);
    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    const double r = R_val();
    const double GRID = GRID_val();
    const double r2 = r * r;
    // treat p as (0, 0), find the topmost index of grid point that is contained
    // in the G_p
    int y_min = -(r / GRID) - 1;
    int y_max = -y_min;
    std::vector<Point> boundaries;
    boundaries.reserve(2 * y_max + 1);
    // TODO: improve this implementation
    for (int y = y_min; y <= y_max; y++) {
        const double y_actual = y * GRID;
        const int x_min = -sqrt(r2 - y_actual * y_actual) / GRID - 1;
        boundaries.emplace_back(px + x_min * GRID, py + y_actual);
    }
    for (int y = y_min; y <= y_max; y++) {
        const double y_actual = y * GRID;
        const int x_max = (sqrt(r2 - y_actual * y_actual) / GRID) + 1;
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

// refactor this function, too much duplication with above
std::vector<Point> get_points_from_grid(const Point &p) {
    ScopedTimer timer(TimedOp::GetPointsFromGrid);
    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    const double r = R_val();
    const double GRID = GRID_val();
    const double r2 = r * r;
    // treat p as (0, 0), find the topmost index of grid point that is contained
    // in the G_p
    int y_min = -(r / GRID) - 1;
    int y_max = -y_min;
    std::vector<Point> points;
    points.reserve(2 * y_max + 1);
    for (int y = y_min; y <= y_max; y++) {
        const double y_actual = y * GRID;
        const int x_min = -sqrt(r2 - y_actual * y_actual) / GRID - 1;
        const int x_max = -x_min;
        for (int x = x_min; x <= x_max; x++) {
            points.emplace_back(px + x * GRID, py + y_actual);
        }
    }
    return points;
}

// TODO: binary search
std::vector<int> find_tangent_idx(const Point &p,
                                 const std::vector<Point> &S) {
    ScopedTimer timer(TimedOp::FindTangentIdx);
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

    {
        ScopedTimer timer(TimedOp::CgalIntersection);
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
    }
    return std::nullopt;
}

std::vector<Point> find_F(const Point& p, const std::vector<Point>& S) {
    ScopedTimer timer(TimedOp::FindF);
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
    if (tangent.size() != 2) {
        std::cerr << "Unexpected tangent count: " << tangent.size() << "\n";
        auto F = current_bbox();
        return F;
    }

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
    Polygon polygon(S.begin(), S.end());
    CGAL_precondition(polygon.is_simple());
    CGAL_precondition(polygon.orientation() == CGAL::COUNTERCLOCKWISE);
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
    std::vector<Point> P = get_points_from_grid(p0);
    if (viewer) viewer->markP0(p0);
    std::array<Point, 2> buffer = {p0, p0};
    std::vector<std::vector<Point>> S(P.size(), std::vector<Point>{p0});
    int dead_cnt = 0;
    std::vector<int> dead(P.size());
    cur++;
    while (cur < int(stream.size())) {
        const Point& pi  = stream[cur];
        bool shown_debug = false;
        std::vector<std::vector<Point>> new_S(P.size());
        for (int i = 0; i < int(P.size()); i++) {
            if (dead[i]) {
                continue;
            }
            std::vector<Point> F = find_F(P[i], S[i]);
            std::vector<Point> Gi = get_conv_from_grid(pi);
            Polygon F_poly(F.begin(), F.end());
            Polygon Gi_poly(Gi.begin(), Gi.end());
            Polygon S_poly(S[i].begin(), S[i].end());

            if (!shown_debug) {
                const QColor stepColors[] = { 
                    Qt::red, 
                    QColor(255, 140, 0), // orange
                    Qt::blue, 
                    Qt::green, 
                    Qt::magenta, 
                    Qt::cyan 
                };
                QColor c = stepColors[cur % 6];

                if (showF && viewer) {
                    viewer->addPolygon(F_poly, c);
                } 
                if (showG && viewer) {
                    viewer->addPolygon(Gi_poly, c);
                }
                if (showS && viewer) {
                    viewer->addPolygon(S_poly, c);
                }
                shown_debug = true;
            }

            new_S[i] = intersect_convex_polygons(F, Gi);

            if (new_S[i].empty()) {
                dead[i] = true;
                dead_cnt++;
                continue;
            }
            buffer[0] = P[i];
            buffer[1] = new_S[i].front();
        }
        if (dead_cnt == int(P.size())) {
            break;
        }
        if (viewer) viewer->addOriginalPoint(pi);
        if (viewer) viewer->markPi(pi);
        for (int i = 0; i < int(P.size()); i++) {
            if (dead[i]) continue;
            S[i] = new_S[i];
        }
        cur++;
        if (viewer) viewer_process_events();
    }
    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    if (viewer && showSimp) viewer->addSimplifiedPoint(buffer[0]);
    if (viewer && showSimp) viewer->addSimplifiedPoint(buffer[1]);
    if (viewer && !keepPolygons) viewer->clearPolygons();
    if (viewer) viewer->clearMarkedP0();
    if (viewer) viewer->clearMarkedPi();
    if (viewer) viewer_process_events();
    return cur;
}

std::vector<Point> simplify(const std::vector<Point> &stream, MultiViewer* viewer = nullptr) {
    reset_timing_stats();
    std::vector<Point> simplified;
    std::cout << "Simplifying...\n";
    int cur = 0;
    while (cur != int(stream.size())) {
        std::cout << cur << '\n';
        cur = get_longest_stab(stream, cur, simplified, viewer);
    }
    std::cout << "Simplified!\n";
    print_timing_report();
    return simplified;
}

// TODO: clean up this mess of flags
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
        else if (strcmp(argv[i],"--keep") == 0) keepPolygons = true;
        else if (strcmp(argv[i],"--no-simp") == 0) showSimp = false;
        else if (strcmp(argv[i],"--labels") == 0) showLabels = true;
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

    // Shorthand: first positional numeric argument means '--in <id> --out'; allows extra flags (e.g., '--gui')
    if (test_case_no == -1 && argc >= 2 && argv[1][0] != '-') {
        try {
            test_case_no = std::stoi(argv[1]);
            out_flag = true;
        } catch (...) {
            // fall through to help if not parseable
        }
    }

    if (argc == 1) {
        print_help();
        return 0;
    }

    if (dist_flag) out_flag = true;

    std::vector<Point> stream;
    auto looks_like_repo_root = [](const std::filesystem::path& dir) {
        return std::filesystem::exists(dir / "CMakeLists.txt")
            && (std::filesystem::exists(dir / "simplify.cpp")
                || std::filesystem::exists(dir / "plot_curves.cpp")
                || std::filesystem::exists(dir / "algorithms"));
    };

    auto find_repo_root = [&](std::filesystem::path start) {
        std::error_code ec;
        if (start.empty()) return std::filesystem::path{};
        start = std::filesystem::weakly_canonical(start, ec);
        if (ec) start = start.lexically_normal();

        auto dir = std::filesystem::is_regular_file(start, ec) ? start.parent_path() : start;
        while (!dir.empty()) {
            if (looks_like_repo_root(dir)) return dir;
            auto parent = dir.parent_path();
            if (parent == dir) break;
            dir = parent;
        }
        return std::filesystem::path{};
    };

    // Determine repo_root once (used for input/output/frechet).
    std::filesystem::path repo_root;
    try {
        repo_root = find_repo_root(argv[0]);
    } catch (...) {}
    if (repo_root.empty()) {
        try {
            repo_root = find_repo_root(std::filesystem::current_path());
        } catch (...) {}
    }
    if (repo_root.empty()) {
        repo_root = std::filesystem::current_path();
    }
    if (test_case_no != -1) {
        // Prefer original in data/taxi_simplified/<id>/original.txt; fallback to data/taxi/<id>.txt
        auto simp_orig = repo_root / "data" / "taxi_simplified" / std::to_string(test_case_no) / "original.txt";
        auto raw_input = repo_root / "data" / "taxi" / (std::to_string(test_case_no) + ".txt");

        std::filesystem::path input_path;
        if (std::filesystem::exists(simp_orig)) {
            input_path = simp_orig;
        } else if (std::filesystem::exists(raw_input)) {
            input_path = raw_input;
        } else {
            std::cerr << "Cannot open " << simp_orig.string()
                      << " or " << raw_input.string()
                      << ".\nResolved repo root: " << repo_root.string()
                      << ".\nIf the dataset is not checked out yet, place inputs under data/taxi/<id>.txt.\n";
            return 1;
        }

        std::ifstream fin(input_path.string());
        if (!fin) {
            std::cerr << "Cannot open " << input_path.string() << "\n";
            return 1;
        }
        int N = 0;
        if (!(fin >> N)) { std::cerr << "Empty or invalid input in " << input_path.string() << "\n"; return 1; }
        stream.clear(); stream.reserve(N);
        for (int i = 0; i < N; ++i) {
            double x,y; if (!(fin >> x >> y)) { std::cerr << "Malformed pair at index " << i << " in " << input_path.string() << "\n"; return 1; }
            stream.emplace_back(x, y);
        }
    }

    QApplication app(argc, argv);
    MultiViewer viewer;
    MultiViewer* vptr = nullptr;
    if (gui_flag) {
        vptr = &viewer;
        viewer.setParameters(DELTA, EPSILON);
        viewer.setShowLabels(showLabels);
        viewer.show();
    }

    // Simplify
    const auto simplify_start = std::chrono::steady_clock::now();
    std::vector<Point> simplified = simplify(stream, vptr);
    const auto simplify_end = std::chrono::steady_clock::now();
    const auto simplify_ms = std::chrono::duration_cast<std::chrono::milliseconds>(simplify_end - simplify_start);
    std::cout << "Total elapsed: " << (static_cast<double>(simplify_ms.count()) / 1000.0) << "s\n";

    // Optional output
    if (out_flag) {
        if (test_case_no == -1) {
            std::cerr << "--out requires --in <id> to determine output location\n";
        } else {
            std::filesystem::path dir = repo_root / "data" / "taxi_simplified" / std::to_string(test_case_no);
            std::filesystem::create_directories(dir);
            
            // 1. our_simplified.txt
            std::ofstream simp(dir / "our_simplified.txt");
            std::size_t N = simplified.size();
            simp << N << '\n';
            for (const auto& p : simplified) {
                simp << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()) << '\n';
            }
            simp.close();

            // 2. <eps>_<delta>_our_simplified.txt
            std::string param_name = std::to_string(EPSILON) + "_" + std::to_string(DELTA) + "_our_simplified.txt";
            std::ofstream sad(dir / param_name);
            sad << N << '\n';
            for (const auto& p : simplified) {
                sad << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()) << '\n';
            }
            sad.close();
        }
        std::cout << "Output Written\n";
    }

    int gui_result = 0;
    if (gui_flag) {
        gui_result = app.exec();
    }

    // Run distance computation only after GUI is closed (or immediately if no GUI)
    if (dist_flag && test_case_no != -1) {
        std::filesystem::path frechet_path = repo_root / "frechet";
        // Quote the path to handle spaces (e.g., "Mobile Documents") and pass explicit args
        std::string cmd1 = std::string("\"") + frechet_path.string() + "\" --in " + std::to_string(test_case_no) + " --simplified";
        int rc = std::system(cmd1.c_str());
        (void)rc;
    }

    return gui_result;
}