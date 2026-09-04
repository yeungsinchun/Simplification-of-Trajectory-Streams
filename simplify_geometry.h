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
#include <optional>
#include <vector>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point  = Kernel::Point_2;

using Polygon = CGAL::Polygon_2<Kernel>;

// ===========================================================================
//  Convex-convex intersection by double-precision Sutherland-Hodgman
// ===========================================================================
//
// Both inputs are convex (F(S,p) is a convex wedge clipped to the axis-aligned
// bbox; conv(G) is a convex grid-cell polygon), so their intersection is the
// subject polygon successively clipped by each half-plane of the clip polygon.
// Sutherland-Hodgman only ever *removes* subject area against each clip
// half-plane, so the output is always a subset of the subject and can never
// balloon; it needs no convexity guard, no containment guard, and no exact
// fallback.  Complexity is O(n*m) but both polygons are tiny (n,m ~ 4..20) and
// every operation is a plain double.
namespace sh_double {

// Reverse the polygon in place if it is clockwise, so callers can assume CCW.
inline void to_ccw(std::vector<std::array<double, 2>>& v) {
    double a = 0.0;
    const int n = static_cast<int>(v.size());
    for (int i = 0; i < n; ++i) {
        const auto& p = v[i];
        const auto& q = v[(i + 1) % n];
        a += p[0] * q[1] - q[0] * p[1];
    }
    if (a < 0.0) std::reverse(v.begin(), v.end());
}

// Thread-local workspace to avoid repeated allocations in clip().
// Each thread gets its own workspace, reused across all clip() calls.
struct ClipWorkspace {
    std::vector<std::array<double, 2>> subj;
    std::vector<std::array<double, 2>> clipv;
    std::vector<std::array<double, 2>> out;
    std::vector<std::array<double, 2>> in;
    std::vector<Point> result;
};

inline ClipWorkspace& get_clip_workspace() {
    thread_local ClipWorkspace ws;
    return ws;
}

// Clip convex subject P by convex clip Q; both given as Epick points.  Returns
// a const reference to the thread-local intersection polygon (CCW, no closing
// duplicate) in Epick points, or an empty vector when the intersection is
// degenerate (<3 vertices).  The returned reference is valid until the next
// clip() call on the same thread.
inline const std::vector<Point>& clip(const std::vector<Point>& P_verts,
                                      const std::vector<Point>& Q_verts) {
    // Reuse thread-local workspace to avoid allocations
    ClipWorkspace& ws = get_clip_workspace();

    const int n = static_cast<int>(P_verts.size());
    const int m = static_cast<int>(Q_verts.size());
    if (n < 3 || m < 3) { ws.result.clear(); return ws.result; }

    ws.subj.clear();
    ws.clipv.clear();
    ws.subj.reserve(n);
    ws.clipv.reserve(m);
    for (const auto& p : P_verts)
        ws.subj.push_back({CGAL::to_double(p.x()), CGAL::to_double(p.y())});
    for (const auto& q : Q_verts)
        ws.clipv.push_back({CGAL::to_double(q.x()), CGAL::to_double(q.y())});
    to_ccw(ws.subj);
    to_ccw(ws.clipv);

    // subj is not read after this; swap avoids copying the subject polygon.
    ws.out.swap(ws.subj);
    ws.in.clear();
    for (int e = 0; e < m && ws.out.size() >= 3; ++e) {
        const int e1 = (e + 1 == m) ? 0 : e + 1;   // avoid per-iteration modulo
        const double ax = ws.clipv[e][0],  ay = ws.clipv[e][1];
        const double bx = ws.clipv[e1][0], by = ws.clipv[e1][1];
        const double ex = bx - ax, ey = by - ay;
        // Signed area of triangle (a,b,pt); >=0 means pt is on/left of edge a->b
        // (inside, since the clip is CCW).
        auto side = [&](const std::array<double, 2>& pt) {
            return ex * (pt[1] - ay) - ey * (pt[0] - ax);
        };
        ws.in.swap(ws.out);
        ws.out.clear();
        const int sz = static_cast<int>(ws.in.size());
        int prev = sz - 1;                          // walk the previous index
        for (int i = 0; i < sz; ++i) {
            const auto& cur = ws.in[i];
            const auto& prv = ws.in[prev];
            prev = i;
            const double sc = side(cur), sp = side(prv);
            const bool cin = sc >= 0.0, pin = sp >= 0.0;
            if (cin) {
                if (!pin) {
                    const double t = sp / (sp - sc);
                    ws.out.push_back({prv[0] + t * (cur[0] - prv[0]),
                                      prv[1] + t * (cur[1] - prv[1])});
                }
                ws.out.push_back(cur);
            } else if (pin) {
                const double t = sp / (sp - sc);
                ws.out.push_back({prv[0] + t * (cur[0] - prv[0]),
                                  prv[1] + t * (cur[1] - prv[1])});
            }
        }
    }
    ws.result.clear();
    if (ws.out.size() < 3) return ws.result;

    ws.result.reserve(ws.out.size());
    for (const auto& p : ws.out) ws.result.emplace_back(p[0], p[1]);
    return ws.result;
}

}  // namespace sh_double

// ===========================================================================
//  Intersection (public API)
// ===========================================================================

// Collapse consecutive near-coincident vertices (including the wrap-around
// pair).  Near-duplicate points reach the intersector via find_F, producing
// edges only a few ULPs long; removing them up front avoids spurious
// near-zero-length clip edges.  EPS2 is a squared distance: 1e-12 == (1e-6)^2,
// far above the ~1e-13 duplicate noise yet far below any real feature.
// Core dedup, writing into a caller-provided buffer (reused across calls to
// avoid per-call allocation).  `out` must not alias `poly`.
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

inline std::vector<Point> dedup_consecutive(const std::vector<Point>& poly) {
    std::vector<Point> out;
    dedup_into(poly, out);
    return out;
}

// Convex-convex intersection of P_in and Q_in.  Deduplicates both inputs,
// runs the Sutherland-Hodgman clip, and returns true with the CCW result in
// `result` when the intersection is a non-degenerate polygon (>=3 vertices).
inline bool intersect(const std::vector<Point>& P_in,
                      const std::vector<Point>& Q_in,
                      std::vector<Point>& result) {
    // Thread-local scratch for the deduplicated inputs; reused across calls so
    // the steady state performs no heap allocation here.
    thread_local std::vector<Point> P_verts, Q_verts;
    dedup_into(P_in, P_verts);
    dedup_into(Q_in, Q_verts);

    // clip returns a reference to its own thread-local buffer, distinct from
    // `result`, so deduping straight into `result` is safe.
    dedup_into(sh_double::clip(P_verts, Q_verts), result);
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

// Boundary anchors for P: every grid sample on the convex-hull boundary.
// CGAL's hull keeps extreme vertices only; flat hull edges still contain
// collinear grid corners that must remain boundary anchors.
inline std::vector<Point> get_boundary_points_from_grid(const Point& p, double EPSILON, double DELTA, int multiplier = 1) {
    thread_local double cached_eps   = std::numeric_limits<double>::quiet_NaN();
    thread_local double cached_delta = std::numeric_limits<double>::quiet_NaN();
    thread_local int cached_mult = 0;
    thread_local std::vector<std::array<double, 2>> boundary_offsets;

    if (EPSILON != cached_eps || DELTA != cached_delta || multiplier != cached_mult) {
        const std::vector<Point> all = get_points_from_grid(Point(0, 0), EPSILON, DELTA, multiplier);
        std::vector<Point> hull;
        if (all.size() <= 2) {
            hull = all;
        } else {
            CGAL::convex_hull_2(all.begin(), all.end(), std::back_inserter(hull));
        }

        std::vector<Point> boundary;
        boundary.reserve(all.size());
        if (hull.size() < 3) {
            boundary = all;
        } else {
            for (const Point& q : all) {
                if (!strictly_inside_convex(q, hull)) boundary.push_back(q);
            }
        }

        std::sort(boundary.begin(), boundary.end(), [](const Point& a, const Point& b) {
            const double ay = CGAL::to_double(a.y()), by = CGAL::to_double(b.y());
            if (ay != by) return ay > by;
            return CGAL::to_double(a.x()) < CGAL::to_double(b.x());
        });

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

// The two tangent (supporting) vertices from external point p to convex
// polygon S.  Since p lies outside S, all of S falls within an angular wedge
// of <180 degrees where bearing is a total order, so a single linear scan
// keeping the most-clockwise and most-counterclockwise vertices finds both
// tangents in O(n).  The orientation sign is pure-double: on a collinear tie
// either vertex is an equally valid supporting vertex.
//
// Writes 0 or 2 indices into `out` (cleared first).
inline void find_tangent_idx(const Point& p, const std::vector<Point>& S,
                             std::vector<int>& out) {
    out.clear();
    const int n = static_cast<int>(S.size());
    if (n < 1) return;

    // Cache p and S coordinates as doubles once (S is scanned 2n times).
    thread_local std::vector<double> sx, sy;
    sx.resize(n);
    sy.resize(n);
    for (int i = 0; i < n; ++i) {
        sx[i] = CGAL::to_double(S[i].x());
        sy[i] = CGAL::to_double(S[i].y());
    }
    const double pdx = CGAL::to_double(p.x());
    const double pdy = CGAL::to_double(p.y());

    // Sign of orientation(p, S[i], S[j]) in pure doubles.
    // (+1 left, -1 right, 0 collinear.)
    auto orient = [&](int i, int j) -> int {
        const double ux = sx[i] - pdx, uy = sy[i] - pdy;
        const double wx = sx[j] - pdx, wy = sy[j] - pdy;
        const double det = ux * wy - uy * wx;
        if (det > 0.0) return  1;
        if (det < 0.0) return -1;
        return 0;
    };

    int rt = 0, lt = 0;   // most-clockwise, most-counterclockwise
    for (int j = 1; j < n; ++j) {
        if (orient(rt, j) < 0) rt = j;   // S[j] strictly right of p->S[rt]
        if (orient(lt, j) > 0) lt = j;   // S[j] strictly left  of p->S[lt]
    }
    if (rt != lt) {
        out.push_back(std::min(rt, lt));
        out.push_back(std::max(rt, lt));
    }
}

inline std::vector<int> find_tangent_idx(const Point& p, const std::vector<Point>& S) {
    std::vector<int> tangent;
    tangent.reserve(2);
    find_tangent_idx(p, S, tangent);
    return tangent;
}

// Conservative separating-axis prune for the stab loop.  Returns true only when
// the delta-disk region Gi is *provably* disjoint from the free-space wedge
// F(S,p), so the caller can skip the expensive clip and declare the candidate
// dead without computing F at all.
//
// When the prune does not fire and two supporting vertices were found,
// `tangent_out` (if non-null) holds those indices so find_F can reuse them
// instead of scanning S again.
//
// F(S,p) is always contained in the cone with apex p bounded by the two tangent
// rays p->S[t0] and p->S[t1] (all of S, its arc, and the bbox hits lie inside
// that angular wedge).  The cone is the intersection of two half-planes H0, H1
// whose boundary lines pass through p.  Hence if every vertex of the convex Gi
// lies strictly outside H0 (or strictly outside H1), then Gi cannot meet the
// cone and therefore cannot meet F, so F n Gi is empty.
//
// The test is one-sided: it uses the same raw-double tangents as find_F and an
// a-priori rounding margin, so it fires only when the separation is robust.
// Near-boundary cases fall through to the exact clip.  It can therefore never
// drop a candidate the exact intersection would have kept -- the Frechet
// guarantee is preserved (a wrongly-kept candidate is impossible; a wrongly
// dropped one cannot occur because we only prune on robust separation).
inline bool wedge_gi_disjoint(const Point& p, const std::vector<Point>& S,
                              const std::vector<Point>& Gi,
                              std::vector<int>* tangent_out = nullptr) {
    const int sn = static_cast<int>(S.size());
    if (sn < 3) return false;                    // single point / degenerate: F = bbox
    if (point_in_convex(p, S)) return false;     // p inside S: F = whole bbox

    // Prefer writing tangents straight into the caller's buffer to avoid an
    // extra allocation on the hot miss path.
    thread_local std::vector<int> local_tangent;
    std::vector<int>& tangent = tangent_out ? *tangent_out : local_tangent;
    find_tangent_idx(p, S, tangent);
    if (tangent.size() != 2) return false;       // find_F bails here anyway

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

// The free-space wedge F(S,p): the region of the bbox reachable from the
// convex stab region S through the external point p.  When p is inside S (or S
// is a single point) every bbox point is reachable, so F is the whole bbox;
// otherwise F is bounded by the two tangent rays from p to S, the arc of S
// between the tangent vertices, and the bbox boundary between the two ray hits.
//
// Optional `tangent_in`: when non-null and size==2, skips find_tangent_idx
// (caller already computed the supporting vertices, e.g. during a prune miss).
inline void find_F(const Point& p, const std::vector<Point>& S,
                   std::vector<Point>& F,
                   const std::vector<int>* tangent_in = nullptr) {
    F.clear();
    assert(S.size() != 2);

    // Precomputed tangents imply p is outside S (prune already established that).
    if (!(tangent_in && tangent_in->size() == 2)) {
        bool p_in_S = point_in_convex(p, S);
        if (S.size() == 1 || p_in_S) {
            F = current_bbox();
            return;
        }
    } else if (S.size() == 1) {
        F = current_bbox();
        return;
    }

    std::vector<int> tangent_storage;
    const std::vector<int>* tangent = tangent_in;
    if (!tangent || tangent->size() != 2) {
        find_tangent_idx(p, S, tangent_storage);
        tangent = &tangent_storage;
    }
    if (tangent->size() != 2) {
        return;
    }

    auto hit1 = ray_hit_bbox(p, S[(*tangent)[0]]);
    auto hit2 = ray_hit_bbox(p, S[(*tangent)[1]]);
    if (!hit1 || !hit2) {
        F = current_bbox();
        return;
    }

    auto e1 = which_edge(hit1.value());
    auto e2 = which_edge(hit2.value());
    if (!e1 || !e2) {
        F = current_bbox();
        return;
    }

    int n = int(S.size());
    assert(n >= 3);
    assert((*tangent)[1] - (*tangent)[0] - 1 >= 1 || (*tangent)[0] + n - (*tangent)[1] - 1 >= 1);

    F.reserve(n + 4);
    // Raw-double right_turn: sign((S[t1] - p) x (S[t2] - p))
    const double ax = CGAL::to_double(S[(*tangent)[0]].x() - p.x()),
                  ay = CGAL::to_double(S[(*tangent)[0]].y() - p.y());
    const double bx = CGAL::to_double(S[(*tangent)[1]].x() - p.x()),
                  by = CGAL::to_double(S[(*tangent)[1]].y() - p.y());
    const bool is_right_turn = (ax * by - ay * bx) < 0;
    if (is_right_turn) {
        std::copy(S.begin() + (*tangent)[0], S.begin() + (*tangent)[1] + 1, std::back_inserter(F));
        F.push_back(hit2.value());
        append_rect_pts(F, e2.value(), e1.value(), true);
        F.push_back(hit1.value());
    } else {
        std::copy(S.begin() + (*tangent)[1], S.end(), std::back_inserter(F));
        std::copy(S.begin(), S.begin() + (*tangent)[0] + 1, std::back_inserter(F));
        F.push_back(hit1.value());
        append_rect_pts(F, e1.value(), e2.value(), true);
        F.push_back(hit2.value());
    }
}

#endif // SIMPLIFY_GEOMETRY_H
