#ifndef SIMPLIFY_GEOMETRY_H
#define SIMPLIFY_GEOMETRY_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point  = Kernel::Point_2;

using Polygon = CGAL::Polygon_2<Kernel>;

// ===========================================================================
//  Polygon intersection: double-precision Sutherland-Hodgman
// ===========================================================================
//
// Used for F(S,p) ∩ conv(G_i). Both are convex here, but SH only requires the
// cropping polygon to be convex. O(n*m) in doubles; n,m are small (~4..20).
namespace sh_double {

using Vec2 = std::array<double, 2>;

struct ClipEdge {
    Vec2 start;
    double dx, dy;
};

inline Vec2 sub(const Vec2& u, const Vec2& v) {
    return {u[0] - v[0], u[1] - v[1]};
}

inline double cross(const Vec2& u, const Vec2& v) {
    return u[0] * v[1] - u[1] * v[0];
}

inline Vec2 lerp(const Vec2& p, const Vec2& q, double t) {
    return {p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])};
}

// orient(a,b,p) = cross(b-a, p-a) = 2 * signed area of triangle (a,b,p).
// Positive => p is left of directed edge a→b.
inline double orient(const Vec2& a, const Vec2& b, const Vec2& p) {
    return cross(sub(b, a), sub(p, a));
}

// Point on segment p→q where orient along that segment is zero:
//   (1-t)*orient_p + t*orient_q = 0  =>  t = orient_p / (orient_p - orient_q).
inline Vec2 crossing_on_segment(const Vec2& p, const Vec2& q,
                                double orient_p, double orient_q) {
    return lerp(p, q, orient_p / (orient_p - orient_q));
}

inline void make_ccw(std::vector<Vec2>& polygon) {
    double twice_area = 0.0;
    const int n = static_cast<int>(polygon.size());
    for (int i = 0; i < n; ++i)
        twice_area += cross(polygon[i], polygon[(i + 1) % n]);
    if (twice_area < 0.0) std::reverse(polygon.begin(), polygon.end());
}

inline Vec2 to_vec2(const Point& p) {
    return {CGAL::to_double(p.x()), CGAL::to_double(p.y())};
}

inline void assign_ccw_doubles(const std::vector<Point>& src, std::vector<Vec2>& dst) {
    dst.clear();
    dst.reserve(src.size());
    for (const auto& p : src) dst.push_back(to_vec2(p));
    make_ccw(dst);
}

inline void assign_points(const std::vector<Vec2>& src, std::vector<Point>& dst) {
    dst.clear();
    dst.reserve(src.size());
    for (const auto& p : src) dst.emplace_back(p[0], p[1]);
}

inline const Vec2& next_ccw_vertex(const std::vector<Vec2>& polygon, int i) {
    const int n = static_cast<int>(polygon.size());
    return polygon[i + 1 == n ? 0 : i + 1];
}

// Keep the part of `polygon` that lies in the closed left half-plane of edge
// edge_start→edge_end. Writes the cropped ring into `cropped` (O(|polygon|)).
inline bool crop_to_left_of_edge(const std::vector<Vec2>& polygon,
                                 const Vec2& edge_start,
                                 const Vec2& edge_end,
                                 std::vector<Vec2>& cropped) {
    cropped.clear();
    const int n = static_cast<int>(polygon.size());
    if (n == 0) return false;

    Vec2 prev = polygon.back();
    double orient_prev = orient(edge_start, edge_end, prev);
    bool changed = false;
    for (int i = 0; i < n; ++i) {
        const Vec2& curr = polygon[i];
        const double orient_curr = orient(edge_start, edge_end, curr);
        const bool curr_inside = orient_curr >= 0.0;
        const bool prev_inside = orient_prev >= 0.0;

        if (!changed && (!curr_inside || !prev_inside)) {
            cropped.insert(cropped.end(), polygon.begin(), polygon.begin() + i);
            changed = true;
        }

        if (changed && curr_inside) {
            if (!prev_inside)
                cropped.push_back(
                    crossing_on_segment(prev, curr, orient_prev, orient_curr));
            cropped.push_back(curr);
        } else if (changed && prev_inside) {
            cropped.push_back(
                crossing_on_segment(prev, curr, orient_prev, orient_curr));
        }
        prev = curr;
        orient_prev = orient_curr;
    }
    return changed;
}

// The core's polygons are usually under 64 vertices. Keep both clip rings in
// reusable contiguous storage and write by index, while allowing larger rings.
struct FastClipBuffer {
    std::array<Vec2, 64> small;
    std::unique_ptr<Vec2[]> large;
    size_t capacity = small.size();
    size_t size = 0;

    Vec2* data() { return large ? large.get() : small.data(); }
    const Vec2* data() const { return large ? large.get() : small.data(); }

    void ensure(size_t needed) {
        if (needed <= capacity) return;
        const size_t next = std::max(needed, capacity * 2);
        auto replacement = std::make_unique<Vec2[]>(next);
        std::copy_n(data(), size, replacement.get());
        large = std::move(replacement);
        capacity = next;
    }
};

struct FastClipBuffers {
    FastClipBuffer first, second;
};

inline FastClipBuffers& fast_clip_buffers() {
    thread_local FastClipBuffers buffers;
    return buffers;
}

__attribute__((always_inline)) inline bool crop_to_left_of_edge_fast(
        const Vec2* polygon, size_t n, const ClipEdge& edge,
        Vec2* cropped, size_t& cropped_size) {
    cropped_size = 0;
    if (n == 0) return false;
    const double x0 = edge.start[0], y0 = edge.start[1];
    const double dx = edge.dx, dy = edge.dy;
    auto side = [&](const Vec2& p) {
        return dx * (p[1] - y0) - dy * (p[0] - x0);
    };
    double orient_prev = side(polygon[n - 1]);
    size_t i = 0;
    double orient_curr = 0.0;
    for (; i < n; ++i) {
        orient_curr = side(polygon[i]);
        if (orient_curr < 0.0 || orient_prev < 0.0) break;
        orient_prev = orient_curr;
    }
    if (i == n) return false;

    Vec2 prev = polygon[i == 0 ? n - 1 : i - 1];
    std::copy_n(polygon, i, cropped);
    cropped_size = i;
    auto emit = [&](const Vec2& curr, double orientation) {
        const bool curr_inside = orientation >= 0.0;
        const bool prev_inside = orient_prev >= 0.0;
        if (curr_inside) {
            if (!prev_inside)
                cropped[cropped_size++] =
                    crossing_on_segment(prev, curr, orient_prev, orientation);
            cropped[cropped_size++] = curr;
        } else if (prev_inside) {
            cropped[cropped_size++] =
                crossing_on_segment(prev, curr, orient_prev, orientation);
        }
        prev = curr;
        orient_prev = orientation;
    };
    emit(polygon[i], orient_curr);
    for (++i; i < n; ++i) {
        const Vec2& curr = polygon[i];
        emit(curr, side(curr));
    }
    return true;
}

// Retains vector capacity across clip() calls on this thread so steady-state
// intersection work does not heap-allocate.
struct ReusableClipBuffers {
    std::vector<Vec2> current_polygon;
    std::vector<Vec2> cropped_polygon;
    std::vector<Vec2> cropping_polygon;
    std::vector<Point> intersection;
};

inline ReusableClipBuffers& reusable_clip_buffers() {
    thread_local ReusableClipBuffers buffers;
    return buffers;
}

// Crop current_polygon by each left half-plane of CCW convex cropping_polygon.
// O(n·m). Returned reference is valid until the next clip() on this thread.
inline const std::vector<Point>& clip_prepared(const std::vector<Point>& current_polygon,
                                               const std::vector<Vec2>& cropping_polygon) {
    auto& buffers = reusable_clip_buffers();

    if (current_polygon.size() < 3 || cropping_polygon.size() < 3) {
        buffers.intersection.clear();
        return buffers.intersection;
    }

    assign_ccw_doubles(current_polygon, buffers.current_polygon);

    const int num_halfplanes = static_cast<int>(cropping_polygon.size());
    for (int e = 0; e < num_halfplanes && buffers.current_polygon.size() >= 3; ++e) {
        const Vec2& edge_start = cropping_polygon[e];
        const Vec2& edge_end   = next_ccw_vertex(cropping_polygon, e);
        if (crop_to_left_of_edge(buffers.current_polygon, edge_start, edge_end,
                                 buffers.cropped_polygon))
            buffers.current_polygon.swap(buffers.cropped_polygon);
    }

    if (buffers.current_polygon.size() < 3) {
        buffers.intersection.clear();
        return buffers.intersection;
    }
    assign_points(buffers.current_polygon, buffers.intersection);
    return buffers.intersection;
}

inline const std::vector<Point>& clip(const std::vector<Point>& current_polygon,
                                      const std::vector<Point>& cropping_polygon) {
    auto& buffers = reusable_clip_buffers();
    assign_ccw_doubles(cropping_polygon, buffers.cropping_polygon);
    return clip_prepared(current_polygon, buffers.cropping_polygon);
}

}  // namespace sh_double

// ===========================================================================
//  Intersection (public API)
// ===========================================================================

/**
 * @brief Remove consecutive near-duplicates (incl. wrap-around first/last).
 *
 * @param poly Input ring (must not alias @p out).
 * @param out  Cleared then filled; reused by callers to avoid allocation.
 *
 * Threshold EPS2 = 1e-12 (= (1e-6)^2): above ULP noise from find_F, below
 * real feature size. Prevents near-zero-length edges from reaching clip().
 */
inline void dedup_into(const std::vector<Point>& poly, std::vector<Point>& out) {
    out.clear();
    const int n = static_cast<int>(poly.size());
    if (n < 2) { out = poly; return; }
    constexpr double EPS2 = 1e-12;
    auto close = [](const Point& a, const Point& b) {
        const double dx = CGAL::to_double(a.x()) - CGAL::to_double(b.x());
        const double dy = CGAL::to_double(a.y()) - CGAL::to_double(b.y());
        return dx * dx + dy * dy <= EPS2;
    };
    out.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (!out.empty() && close(poly[i], out.back())) continue;
        out.push_back(poly[i]);
    }
    while (out.size() >= 2 && close(out.front(), out.back())) out.pop_back();
}

/** @brief Allocating wrapper around dedup_into(). */
inline std::vector<Point> dedup_consecutive(const std::vector<Point>& poly) {
    std::vector<Point> out;
    dedup_into(poly, out);
    return out;
}

/**
 * @brief Convex-convex intersection P ∩ Q via dedup + Sutherland-Hodgman.
 *
 * @param P_in   Subject polygon (typically F(S,p)).
 * @param Q_in   Clip polygon; must be convex (typically conv(G_i)).
 * @param result Cleared; on success holds CCW intersection (>=3 verts).
 * @return true iff result is a non-degenerate polygon.
 */
inline bool intersect(const std::vector<Point>& P_in,
                      const std::vector<Point>& Q_in,
                      std::vector<Point>& result) {
    thread_local std::vector<Point> P_verts, Q_verts;
    dedup_into(P_in, P_verts);
    dedup_into(Q_in, Q_verts);

    // clip's thread-local buffer is distinct from result, so this is safe.
    dedup_into(sh_double::clip(P_verts, Q_verts), result);
    if (result.size() < 3) {
        result.clear();
        return false;
    }
    return true;
}

// The grid hull is shared by every live anchor at a stream step. Prepare its
// deduplicated CCW vertices once instead of repeating that work in each clip.
struct PreparedClipPolygon {
    std::vector<sh_double::Vec2> vertices;
    std::vector<sh_double::ClipEdge> edges;
    double min_x, max_x, min_y, max_y;
};

struct AxisBounds {
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    void include(double x, double y) {
        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
    }

    bool contains(const Point& p) const {
        const double x = CGAL::to_double(p.x()), y = CGAL::to_double(p.y());
        return x >= min_x && x <= max_x && y >= min_y && y <= max_y;
    }
};

inline void prepare_clip_polygon(const std::vector<Point>& Q_in,
                                 PreparedClipPolygon& prepared) {
    thread_local std::vector<Point> unique;
    dedup_into(Q_in, unique);
    sh_double::assign_ccw_doubles(unique, prepared.vertices);
    prepared.edges.clear();
    prepared.edges.reserve(prepared.vertices.size());
    prepared.min_x = prepared.min_y = std::numeric_limits<double>::infinity();
    prepared.max_x = prepared.max_y = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < prepared.vertices.size(); ++i) {
        const auto& a = prepared.vertices[i];
        const auto& b = prepared.vertices[(i + 1) % prepared.vertices.size()];
        prepared.edges.push_back({a, b[0] - a[0], b[1] - a[1]});
        prepared.min_x = std::min(prepared.min_x, a[0]);
        prepared.max_x = std::max(prepared.max_x, a[0]);
        prepared.min_y = std::min(prepared.min_y, a[1]);
        prepared.max_y = std::max(prepared.max_y, a[1]);
    }
}

__attribute__((always_inline)) inline bool intersect_prepared(
        const std::vector<Point>& P_in, const PreparedClipPolygon& Q,
        std::vector<Point>& result, AxisBounds* result_bounds = nullptr,
        sh_double::FastClipBuffers* clip_buffers = nullptr) {
    if (P_in.size() < 3 || Q.vertices.size() < 3) {
        result.clear();
        return false;
    }
    auto& buffers = clip_buffers ? *clip_buffers : sh_double::fast_clip_buffers();
    sh_double::FastClipBuffer* subject = &buffers.first;
    sh_double::FastClipBuffer* scratch = &buffers.second;
    subject->size = 0;
    subject->ensure(P_in.size());
    constexpr double EPS2 = 1e-12;
    auto close = [](const sh_double::Vec2& a, const sh_double::Vec2& b) {
        const double dx = a[0] - b[0], dy = a[1] - b[1];
        return dx * dx + dy * dy <= EPS2;
    };
    // find_F supplies an ordered ring whose inherited S vertices were already
    // deduplicated when the previous stab polygon was produced.
    for (const Point& p : P_in) {
        const sh_double::Vec2 v = sh_double::to_vec2(p);
        subject->data()[subject->size++] = v;
    }

    // find_F emits CCW vertices, so the subject needs no area scan/reversal.
    for (const auto& edge : Q.edges) {
        if (subject->size < 3) break;
        scratch->size = 0;
        scratch->ensure(subject->size * 2 + 2);
        if (sh_double::crop_to_left_of_edge_fast(
                subject->data(), subject->size, edge, scratch->data(), scratch->size))
            std::swap(subject, scratch);
    }

    result.clear();
    if (result_bounds) *result_bounds = {};
    if (subject->size >= 3) {
        result.reserve(subject->size);
        sh_double::Vec2 first{}, last{};
        bool have_last = false;
        for (size_t i = 0; i < subject->size; ++i) {
            const auto& v = subject->data()[i];
            if (!have_last || !close(v, last)) {
                result.emplace_back(v[0], v[1]);
                if (result_bounds) result_bounds->include(v[0], v[1]);
                if (!have_last) first = v;
                last = v;
                have_last = true;
            }
        }
        while (result.size() >= 2 && close(first, last)) {
            result.pop_back();
            last = sh_double::to_vec2(result.back());
        }
    }
    if (result.size() < 3) {
        result.clear();
        return false;
    }
    return true;
}

// ===========================================================================
//  Trajectory-simplification geometry
// ===========================================================================
//
// The reachable-region construction from the streaming simplification
// algorithm: the bounding box, the delta-disk grid, and the free-space wedge
// F(S,p).  Shared by both the headless (simplify.cpp) and GUI
// (simplify_with_gui.cpp) front-ends, so it lives here as inline definitions.

// Axis-aligned working bounding box, sized per input by configure_bbox().
// inline (C++17) so both translation units share one definition.
inline double BMIN = -10000;
inline double BMAX = 10000;
inline constexpr double GEOM_TOL = 1e-6;

// Squared radius (in grid-offset units) of the farthest grid corner ever used;
// sqrt gives the algorithm's a-priori Frechet-distance guarantee.  Written by
// get_points_from_grid(), read by the front-ends after simplification.
inline double expected_frechet_squared = 0.0;

// Point-in-CCW-convex-polygon via filtered orientation: a double determinant
// guarded by its a-priori rounding bound resolves the clear majority; only
// genuinely ambiguous (near-boundary) corners defer to the exact predicate.
inline bool point_in_convex(const Point& p, const std::vector<Point>& poly, bool ccw = true) {
    const int n = static_cast<int>(poly.size());
    // A point cannot be inside a polygon with fewer than 3 vertices
    if (n < 3) return false;
    
    const int bad = ccw ? -1 : 1;   // sign that means "outside" (right for CCW, left for CW)
    const double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y());
    for (int i = 0; i < n; ++i) {
        const double ax = CGAL::to_double(poly[i].x()), ay = CGAL::to_double(poly[i].y());
        const double bx = CGAL::to_double(poly[(i + 1) % n].x()), by = CGAL::to_double(poly[(i + 1) % n].y());
        const double t1 = (bx - ax) * (py - ay);
        const double t2 = (by - ay) * (px - ax);
        const double det = t1 - t2;
        const double bound = 8.0 * std::numeric_limits<double>::epsilon() *
                             (std::abs(t1) + std::abs(t2));
        int s;
        if (det >  bound)      s =  1;
        else if (det < -bound) s = -1;
        else {
            switch (CGAL::orientation(poly[i], poly[(i + 1) % n], p)) {
                case CGAL::LEFT_TURN:  s =  1; break;
                case CGAL::RIGHT_TURN: s = -1; break;
                default:               s =  0; break;
            }
        }
        if (s == bad) return false;
    }
    return true;
}

// True when q lies on the closed segment ab (within GEOM_TOL).
inline bool point_on_segment(const Point& a, const Point& b, const Point& q) {
    const double ax = CGAL::to_double(a.x()), ay = CGAL::to_double(a.y());
    const double bx = CGAL::to_double(b.x()), by = CGAL::to_double(b.y());
    const double qx = CGAL::to_double(q.x()), qy = CGAL::to_double(q.y());
    const double t1 = (bx - ax) * (qy - ay);
    const double t2 = (by - ay) * (qx - ax);
    const double det = t1 - t2;
    const double bound = 8.0 * std::numeric_limits<double>::epsilon() *
                         (std::abs(t1) + std::abs(t2));
    if (std::abs(det) > bound) return false;
    if (CGAL::orientation(a, b, q) != CGAL::COLLINEAR) return false;
    return qx >= std::min(ax, bx) - GEOM_TOL && qx <= std::max(ax, bx) + GEOM_TOL &&
           qy >= std::min(ay, by) - GEOM_TOL && qy <= std::max(ay, by) + GEOM_TOL;
}

// Strict interior of a CCW convex polygon; boundary and exterior return false.
inline bool strictly_inside_convex(const Point& p, const std::vector<Point>& poly, bool ccw = true) {
    const int n = static_cast<int>(poly.size());
    if (n < 3) return false;

    const int bad = ccw ? -1 : 1;
    const double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y());
    for (int i = 0; i < n; ++i) {
        const double ax = CGAL::to_double(poly[i].x()), ay = CGAL::to_double(poly[i].y());
        const double bx = CGAL::to_double(poly[(i + 1) % n].x()), by = CGAL::to_double(poly[(i + 1) % n].y());
        const double t1 = (bx - ax) * (py - ay);
        const double t2 = (by - ay) * (px - ax);
        const double det = t1 - t2;
        const double bound = 8.0 * std::numeric_limits<double>::epsilon() *
                             (std::abs(t1) + std::abs(t2));
        int s;
        if (det >  bound)      s =  1;
        else if (det < -bound) s = -1;
        else {
            if (point_on_segment(poly[i], poly[(i + 1) % n], p)) return false;
            switch (CGAL::orientation(poly[i], poly[(i + 1) % n], p)) {
                case CGAL::LEFT_TURN:  s =  1; break;
                case CGAL::RIGHT_TURN: s = -1; break;
                default:               s =  0; break;
            }
        }
        if (s != 1) return false;
    }
    return true;
}

// First intersection of the ray p->dir with the working bbox, in doubles.
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
    if (std::abs(dx) > 1e-18)
        consider(((dx < 0.0 ? BMIN : BMAX) - px) / dx);
    if (std::abs(dy) > 1e-18)
        consider(((dy < 0.0 ? BMIN : BMAX) - py) / dy);
    if (!std::isfinite(best)) return std::nullopt;
    return Point(hx, hy);
}

// Corners (BL,BR,TR,TL) and edges of the bbox, indexed CCW for append_rect_pts.
enum class Bbox_edge {
    BL = 0, BOTTOM = 1, BR = 2, RIGHT = 3,
    TR = 4, TOP = 5,    TL = 6, LEFT = 7
};

inline std::array<Point, 4> current_bbox_corner() {
    return {
        Point(BMIN, BMIN), Point(BMAX, BMIN), Point(BMAX, BMAX), Point(BMIN, BMAX)
    };
}

inline std::vector<Point> current_bbox() {
    const auto corners = current_bbox_corner();
    return std::vector<Point>(corners.begin(), corners.end());
}

// Classify which bbox edge (or corner) the point s lies on.
inline std::optional<Bbox_edge> which_edge(const Point& s) {
    double x = CGAL::to_double(s.x()), y = CGAL::to_double(s.y());
    bool on_left   = std::abs(x - BMIN) < GEOM_TOL;
    bool on_right  = std::abs(x - BMAX) < GEOM_TOL;
    bool on_bottom = std::abs(y - BMIN) < GEOM_TOL;
    bool on_top    = std::abs(y - BMAX) < GEOM_TOL;

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

// Append the bbox corners strictly between edges `from` and `to`, walking CCW
// (or CW when ccw==false), when stitching the wedge boundary along the bbox.
inline void append_rect_pts(std::vector<Point>& out, Bbox_edge from, Bbox_edge to, bool ccw) {
    const auto corners = current_bbox_corner();
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

// Grid spacing and disk radius from the (epsilon, delta) parameters.
inline double GRID_val(double EPSILON, double DELTA, int multiplier = 1) {
    return EPSILON * DELTA / (2.0 * std::sqrt(2.0)) / multiplier;
}

inline double R_val(double EPSILON, double DELTA) {
    return (1.0 + EPSILON / 2.0) * DELTA;
}

// Size the working bbox to the input plus one grid-reach of padding, so every
// grid corner and wedge ray stays inside the box.
inline void configure_bbox(const std::vector<Point>& stream, double EPSILON, double DELTA) {
    double min_coord = std::numeric_limits<double>::infinity();
    double max_coord = -std::numeric_limits<double>::infinity();
    for (const Point& point : stream) {
        min_coord = std::min({min_coord, CGAL::to_double(point.x()), CGAL::to_double(point.y())});
        max_coord = std::max({max_coord, CGAL::to_double(point.x()), CGAL::to_double(point.y())});
    }
    const double grid_reach = R_val(EPSILON, DELTA) + GRID_val(EPSILON, DELTA) * std::sqrt(2.0);
    const double padding = grid_reach + std::max(1.0, grid_reach * 1e-6);
    BMIN = min_coord - padding;
    BMAX = max_coord + padding;
}

// All distinct grid corners within radius R of p (the delta-disk sample set).
// Updates expected_frechet_squared with the farthest corner offset seen.
inline std::vector<Point> get_points_from_grid(const Point& p, double EPSILON, double DELTA, int multiplier = 1) {
    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    const double GRID = GRID_val(EPSILON, DELTA, multiplier);
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
    std::vector<uint8_t> seen(size_t(corner_count) * corner_count, 0);
    std::vector<std::pair<int, int>> corner_coords;
    corner_coords.reserve(size_t(corner_count) * corner_count);

    auto add_corner = [&](int ji, int ki) {
        const int ix = ji - j_min, iy = ki - j_min;
        const size_t index = size_t(ix) * corner_count + iy;
        if (seen[index]) return;
        seen[index] = 1;
        const double corner_x = ji * GRID;
        const double corner_y = ki * GRID;
        expected_frechet_squared = std::max(
            expected_frechet_squared, corner_x * corner_x + corner_y * corner_y);
        corner_coords.push_back({ji, ki});
    };

    // Collect all corners within the disk
    for (int k = j_min; k <= j_max; ++k) {
        const double y0 = k * GRID, y1 = (k + 1) * GRID;
        for (int j = j_min; j <= j_max; ++j) {
            const double x0 = j * GRID, x1 = (j + 1) * GRID;
            const double nearest_x = x0 > 0.0 ? x0 : (x1 < 0.0 ? x1 : 0.0);
            const double nearest_y = y0 > 0.0 ? y0 : (y1 < 0.0 ? y1 : 0.0);
            if (nearest_x * nearest_x + nearest_y * nearest_y > r2) continue;
            add_corner(j,     k);
            add_corner(j + 1, k);
            add_corner(j + 1, k + 1);
            add_corner(j,     k + 1);
        }
    }

    // Sort corners in row-major order: top-to-bottom (descending y), then left-to-right (ascending x)
    std::sort(corner_coords.begin(), corner_coords.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second) return a.second > b.second; // y descending (top first)
        return a.first < b.first; // x ascending (left to right)
    });

    // Build the final points vector
    std::vector<Point> points;
    points.reserve(corner_coords.size());
    for (const auto& [ji, ki] : corner_coords) {
        const double corner_x = ji * GRID;
        const double corner_y = ki * GRID;
        points.emplace_back(px + corner_x, py + corner_y);
    }

    return points;
}

// Convex hull of the delta-disk grid samples around p (the region conv(G_i)).
//
// The grid-corner offsets depend only on (EPSILON, DELTA), never on p: every
// call to get_points_from_grid(p, ...) yields the same corner set merely
// translated by p.  Convex hull is translation-equivariant, so conv(G_i) has
// an identical shape for every point in the stream and only its position
// changes.  We therefore build the origin-centred hull once per (EPSILON,
// DELTA) and translate that cached template by p on each call, turning a
// per-step O(m log m) hull build into an O(h) copy-with-offset.
inline std::vector<Point> get_conv_from_grid(const Point& p, double EPSILON, double DELTA, int multiplier = 1) {
    thread_local double cached_eps   = std::numeric_limits<double>::quiet_NaN();
    thread_local double cached_delta = std::numeric_limits<double>::quiet_NaN();
    thread_local int cached_mult = 0;
    thread_local std::vector<std::array<double, 2>> hull_offsets;

    if (EPSILON != cached_eps || DELTA != cached_delta || multiplier != cached_mult) {
        std::vector<Point> points = get_points_from_grid(Point(0, 0), EPSILON, DELTA, multiplier);
        std::vector<Point> conv;
        CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(conv));
        hull_offsets.clear();
        hull_offsets.reserve(conv.size());
        for (const auto& q : conv)
            hull_offsets.push_back({CGAL::to_double(q.x()), CGAL::to_double(q.y())});
        cached_eps   = EPSILON;
        cached_delta = DELTA;
        cached_mult  = multiplier;
    }

    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    std::vector<Point> conv;
    conv.reserve(hull_offsets.size());
    for (const auto& off : hull_offsets)
        conv.emplace_back(px + off[0], py + off[1]);
    return conv;
}

// Boundary anchors for P: discrete convex outline of the grid samples -
// leftmost and rightmost on every y-row, plus every sample on the topmost
// and bottommost rows.
inline std::vector<Point> get_boundary_points_from_grid(const Point& p, double EPSILON, double DELTA, int multiplier = 1) {
    thread_local double cached_eps   = std::numeric_limits<double>::quiet_NaN();
    thread_local double cached_delta = std::numeric_limits<double>::quiet_NaN();
    thread_local int cached_mult = 0;
    thread_local std::vector<std::array<double, 2>> boundary_offsets;

    if (EPSILON != cached_eps || DELTA != cached_delta || multiplier != cached_mult) {
        // get_points_from_grid is row-major: y descending, then x ascending.
        const std::vector<Point> all = get_points_from_grid(Point(0, 0), EPSILON, DELTA, multiplier);
        std::vector<Point> boundary;
        boundary.reserve(all.size());
        if (!all.empty()) {
            const double y_top = CGAL::to_double(all.front().y());
            const double y_bot = CGAL::to_double(all.back().y());
            for (size_t i = 0; i < all.size(); ) {
                size_t j = i + 1;
                const double y = CGAL::to_double(all[i].y());
                while (j < all.size() && CGAL::to_double(all[j].y()) == y)
                    ++j;
                if (y == y_top || y == y_bot) {
                    boundary.insert(boundary.end(), all.begin() + static_cast<std::ptrdiff_t>(i),
                                    all.begin() + static_cast<std::ptrdiff_t>(j));
                } else {
                    boundary.push_back(all[i]);
                    if (j - 1 != i)
                        boundary.push_back(all[j - 1]);
                }
                i = j;
            }
        }

        boundary_offsets.clear();
        boundary_offsets.reserve(boundary.size());
        for (const auto& q : boundary)
            boundary_offsets.push_back({CGAL::to_double(q.x()), CGAL::to_double(q.y())});
        cached_eps   = EPSILON;
        cached_delta = DELTA;
        cached_mult  = multiplier;
    }

    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    std::vector<Point> out;
    out.reserve(boundary_offsets.size());
    for (const auto& off : boundary_offsets)
        out.emplace_back(px + off[0], py + off[1]);
    return out;
}

/**
 * @brief Indices of the two supporting (tangent) vertices from p to convex S.
 *
 * @pre p lies outside convex S. Then all of S sits in a <180° angular wedge
 *      from p, so bearing is a total order.
 * @param p External query point.
 * @param S Convex polygon vertices.
 * @return {min_idx, max_idx} of the most-CW and most-CCW vertices, or empty
 *         if both coincide (degenerate).
 *
 * Algorithm: O(n) scan; keep argmin / argmax of orientation(p, S[i], S[j])
 * in pure doubles. Collinear ties: either vertex is a valid support.
 */
inline std::array<int, 2> find_tangent_idx(const Point& p, const std::vector<Point>& S) {
    const int n = static_cast<int>(S.size());

    const double pdx = CGAL::to_double(p.x());
    const double pdy = CGAL::to_double(p.y());

    int rt = 0, lt = 0;   // most-clockwise, most-counterclockwise
    double rtx = CGAL::to_double(S[0].x()) - pdx;
    double rty = CGAL::to_double(S[0].y()) - pdy;
    double ltx = rtx, lty = rty;
    for (int j = 1; j < n; ++j) {
        const double wx = CGAL::to_double(S[j].x()) - pdx;
        const double wy = CGAL::to_double(S[j].y()) - pdy;
        if (rtx * wy - rty * wx < 0.0) {
            rt = j; rtx = wx; rty = wy;
        }
        if (ltx * wy - lty * wx > 0.0) {
            lt = j; ltx = wx; lty = wy;
        }
    }
    if (rt == lt) return {-1, -1};
    return {std::min(rt, lt), std::max(rt, lt)};
}

/**
 * @brief Conservative prune: true ⇒ Gi is provably disjoint from F(S,p).
 *
 * Lets the stab loop skip find_F + clip when separation is robust.
 *
 * Idea: F ⊆ cone(p; rays to tangent verts t0,t1) = H0 ∩ H1. If every vertex
 * of convex Gi lies strictly outside H0 (or outside H1), then Gi ∩ F = ∅.
 *
 * One-sided: uses raw-double tangents + a-priori rounding margin. Near-boundary
 * cases return false and fall through to the exact clip. Never drops a
 * candidate that exact intersection would keep.
 *
 * Early exits (return false = "not proven disjoint"):
 *   |S|<3, p ∈ S, no two tangents, or degenerate cone (orient == 0).
 */
inline bool wedge_gi_disjoint(const Point& p, const std::vector<Point>& S,
                              const std::vector<Point>& Gi) {
    const int sn = static_cast<int>(S.size());
    if (sn < 3) return false;                    // single point / degenerate: F = bbox
    if (point_in_convex(p, S)) return false;     // p inside S: F = whole bbox
    const auto tangent = find_tangent_idx(p, S);
    if (tangent[0] < 0) return false;             // find_F bails here anyway

    const double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y());
    const double t0x = CGAL::to_double(S[tangent[0]].x()) - px;
    const double t0y = CGAL::to_double(S[tangent[0]].y()) - py;
    const double t1x = CGAL::to_double(S[tangent[1]].x()) - px;
    const double t1y = CGAL::to_double(S[tangent[1]].y()) - py;

    // orient = cross(t0, t1); its sign tells which side of each ray is interior
    // (the side containing the other ray).
    const double orient = t0x * t1y - t0y * t1x;
    if (orient == 0.0) return false;             // degenerate cone
    const bool ccw = orient > 0.0;               // t1 is CCW of t0

    constexpr double K = 8.0 * std::numeric_limits<double>::epsilon();
    const int m = static_cast<int>(Gi.size());

    bool sep0 = true, sep1 = true;
    for (int k = 0; k < m && (sep0 || sep1); ++k) {
        const double qx = CGAL::to_double(Gi[k].x()) - px;
        const double qy = CGAL::to_double(Gi[k].y()) - py;
        if (sep0) {
            const double a = t0x * qy, b = t0y * qx;   // cross(t0, q) = a - b
            const double c0 = a - b;
            const double tol = K * (std::abs(a) + std::abs(b));
            // Interior of H0 has sign(c0) == sign(orient); exterior is the
            // opposite sign.  Require robustly-exterior to keep separation.
            if (ccw ? (c0 >= -tol) : (c0 <= tol)) sep0 = false;
        }
        if (sep1) {
            const double a = t1x * qy, b = t1y * qx;   // cross(t1, q) = a - b
            const double c1 = a - b;
            const double tol = K * (std::abs(a) + std::abs(b));
            // Interior of H1 has sign(c1) == -sign(orient); exterior opposite.
            if (ccw ? (c1 <= tol) : (c1 >= -tol)) sep1 = false;
        }
    }
    return sep0 || sep1;
}

/**
 * @brief Build free-space wedge F(S,p) into @p F (cleared first).
 *
 * F is the region of the working bbox reachable from convex stab region S
 * through external point p.
 *
 * Cases:
 *   1. |S|==1 or p ∈ S     → F = whole bbox.
 *   2. No two tangents     → F left empty (return).
 *   3. Ray miss / unclassified hit → F = whole bbox.
 *   4. Otherwise (p outside S):
 *        F = S-arc between tangents  ∪  bbox chain between ray hits
 *            ∪  the two ray-hit points.
 *      Arc choice: if (S[t0]-p)×(S[t1]-p) < 0 (right turn), copy t0..t1;
 *      else copy the wrap-around arc t1..end + begin..t0.
 *
 * @param p External (or interior) query point.
 * @param S Convex stab region; must not have exactly 2 vertices (assert).
 * @param F Output polygon, CCW.
 */
// Returns true when F is the full working bbox.
__attribute__((always_inline)) inline bool find_F(
        const Point& p, const std::vector<Point>& S, std::vector<Point>& F,
        const PreparedClipPolygon* clip_bounds = nullptr,
        bool* disjoint = nullptr,
        const AxisBounds* stab_bounds = nullptr,
        uint8_t* anchor_outside = nullptr) {
    F.clear();
    if (disjoint) *disjoint = false;
    assert(S.size() != 2);
    auto use_bbox = [&] {
        const auto corners = current_bbox_corner();
        F.assign(corners.begin(), corners.end());
        if (anchor_outside) *anchor_outside = false;
    };
    if (S.size() == 1) {
        use_bbox();
        return true;
    }
    if (!anchor_outside || !*anchor_outside) {
        if ((!stab_bounds || stab_bounds->contains(p)) && point_in_convex(p, S)) {
            use_bbox();
            return true;
        }
        // Each later stab polygon is contained in the continuation wedge.
        // The anchor lies outside that wedge once it lies outside S.
        if (anchor_outside) *anchor_outside = true;
    }

    const auto tangent = find_tangent_idx(p, S);
    if (tangent[0] < 0) {
        return false;
    }

    const double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y());
    const double ax = CGAL::to_double(S[tangent[0]].x()) - px;
    const double ay = CGAL::to_double(S[tangent[0]].y()) - py;
    const double bx = CGAL::to_double(S[tangent[1]].x()) - px;
    const double by = CGAL::to_double(S[tangent[1]].y()) - py;
    const double turn = ax * by - ay * bx;
    if (clip_bounds && disjoint && turn != 0.0) {
        // Every Gi vertex lies in this box. A tangent half-plane that misses
        // the whole box also misses Gi and therefore the free-space wedge.
        auto range = [&](double tx, double ty) {
            const double xmin = clip_bounds->min_x - px;
            const double xmax = clip_bounds->max_x - px;
            const double ymin = clip_bounds->min_y - py;
            const double ymax = clip_bounds->max_y - py;
            const double lo = tx * (tx >= 0 ? ymin : ymax) -
                              ty * (ty >= 0 ? xmax : xmin);
            const double hi = tx * (tx >= 0 ? ymax : ymin) -
                              ty * (ty >= 0 ? xmin : xmax);
            const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                     (std::abs(tx) + std::abs(ty)) *
                                     (std::abs(xmin) + std::abs(xmax) +
                                      std::abs(ymin) + std::abs(ymax) + 1.0);
            return std::array<double, 3>{lo, hi, tolerance};
        };
        const auto first = range(ax, ay);
        const auto second = range(bx, by);
        *disjoint = turn > 0.0
            ? (first[1] < -first[2] || second[0] > second[2])
            : (first[0] > first[2] || second[1] < -second[2]);
        if (*disjoint) return false;
    }

    auto hit1 = ray_hit_bbox(p, S[tangent[0]]);
    auto hit2 = ray_hit_bbox(p, S[tangent[1]]);
    if (!hit1 || !hit2) {
        use_bbox();
        return true;
    }

    auto e1 = which_edge(hit1.value());
    auto e2 = which_edge(hit2.value());
    if (!e1 || !e2) {
        use_bbox();
        return true;
    }

    int n = int(S.size());
    assert(n >= 3);
    assert(tangent[1] - tangent[0] - 1 >= 1 || tangent[0] + n - tangent[1] - 1 >= 1);

    F.reserve(n + 4);
    // Raw-double right_turn: sign((S[t1] - p) x (S[t2] - p))
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
    return false;
}

#endif // SIMPLIFY_GEOMETRY_H
