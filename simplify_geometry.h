#ifndef SIMPLIFY_GEOMETRY_H
#define SIMPLIFY_GEOMETRY_H

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <vector>

using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;

static inline CGAL::Cartesian_converter<Epick, Epeck> conv_to_exact;
static inline CGAL::Cartesian_converter<Epeck, Epick> conv_to_inexact;

using Kernel = Epick;
using Point  = Kernel::Point_2;
using Vector = Kernel::Vector_2;
using Segment = Kernel::Segment_2;
using Ray    = Kernel::Ray_2;
using Line   = Kernel::Line_2;
using Bbox   = CGAL::Bbox_2;
using Rect   = CGAL::Iso_rectangle_2<Kernel>;

using Polygon = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes = CGAL::Polygon_with_holes_2<Kernel>;

// ===========================================================================
//  Helpers
// ===========================================================================

inline double kernel_signed_area(const std::vector<Point>& P) {
    double s = 0.0;
    int sz = (int)P.size();
    for (int i = 0; i < sz; ++i) {
        const Point& a = P[i];
        const Point& b = P[(i + 1) % sz];
        s += CGAL::to_double(a.x() * b.y() - b.x() * a.y());
    }
    return s / 2.0;
}

// ===========================================================================
//  O'Rourke's ConvexIntersect (Ch. 7, "Computational Geometry in C", 2nd ed.)
// ===========================================================================
namespace orourke_cgal {

enum class InFlag { Pin, Qin, Unknown };

// Sign of the orientation determinant, IDENTICAL to CGAL::orientation but
// cheaper on the common case.  CGAL's Epick predicate runs an interval-
// arithmetic filter (Interval_nt: paired doubles with rounding-mode switches)
// before any exact stage; here we replace that filter with a single double
// determinant guarded by its a-priori forward-error bound (Shewchuk's orient2d
// bound is ~3*eps*(|t1|+|t2|); 8*eps is safely conservative).  When the double
// result clears the bound we return it directly; otherwise we defer to the
// exact predicate.  The result is bit-identical to CGAL::orientation, so no
// approximation is introduced -- only the redundant interval filter is skipped
// on the ~99% of calls that are numerically unambiguous.
inline int area_sign(const Point& a, const Point& b, const Point& c) {
    const double ax = CGAL::to_double(a.x()), ay = CGAL::to_double(a.y());
    const double bx = CGAL::to_double(b.x()), by = CGAL::to_double(b.y());
    const double cx = CGAL::to_double(c.x()), cy = CGAL::to_double(c.y());
    const double t1 = (bx - ax) * (cy - ay);
    const double t2 = (by - ay) * (cx - ax);
    const double det = t1 - t2;
    const double bound = 8.0 * std::numeric_limits<double>::epsilon() *
                         (std::abs(t1) + std::abs(t2));
    if (det >  bound) return  1;
    if (det < -bound) return -1;
    switch (CGAL::orientation(a, b, c)) {
        case CGAL::LEFT_TURN:  return  1;
        case CGAL::RIGHT_TURN: return -1;
        default:               return  0;
    }
}

// Sign of the 2x2 vector determinant, identical to CGAL::determinant's sign
// but with the same double-filter / exact-fallback treatment as area_sign.
inline int area_sign_vec(const Vector& a, const Vector& b) {
    const double ax = CGAL::to_double(a.x()), ay = CGAL::to_double(a.y());
    const double bx = CGAL::to_double(b.x()), by = CGAL::to_double(b.y());
    const double t1 = ax * by, t2 = ay * bx;
    const double det = t1 - t2;
    const double bound = 8.0 * std::numeric_limits<double>::epsilon() *
                         (std::abs(t1) + std::abs(t2));
    if (det >  bound) return  1;
    if (det < -bound) return -1;
    const auto c = CGAL::determinant(a, b);
    if (c > 0) return  1;
    if (c < 0) return -1;
    return 0;
}

inline bool between(const Point& a, const Point& b, const Point& c) {
    if (CGAL::orientation(a, b, c) != CGAL::COLLINEAR) return false;
    if (CGAL::abs(b.x() - a.x()) >= CGAL::abs(b.y() - a.y())) {
        return (a.x() <= c.x() && c.x() <= b.x()) ||
               (a.x() >= c.x() && c.x() >= b.x());
    }
    return (a.y() <= c.y() && c.y() <= b.y()) ||
           (a.y() >= c.y() && c.y() >= b.y());
}

inline int dot_sign(const Point& a, const Point& b,
                    const Point& c, const Point& d) {
    auto ux = b.x() - a.x(), uy = b.y() - a.y();
    auto vx = d.x() - c.x(), vy = d.y() - c.y();
    auto s  = ux * vx + uy * vy;
    if (s > 0) return  1;
    if (s < 0) return -1;
    return 0;
}

inline char parallel_int(const Point& a, const Point& b,
                         const Point& c, const Point& d,
                         Point& out_p, Point& out_q) {
    if (CGAL::orientation(a, b, c) != CGAL::COLLINEAR) return '0';

    if (between(a, b, c) && between(a, b, d)) { out_p = c; out_q = d; return 'e'; }
    if (between(c, d, a) && between(c, d, b)) { out_p = a; out_q = b; return 'e'; }
    if (between(a, b, c) && between(c, d, b)) { out_p = c; out_q = b; return 'e'; }
    if (between(a, b, c) && between(c, d, a)) { out_p = c; out_q = a; return 'e'; }
    if (between(a, b, d) && between(c, d, b)) { out_p = d; out_q = b; return 'e'; }
    if (between(a, b, d) && between(c, d, a)) { out_p = d; out_q = a; return 'e'; }
    return '0';
}

inline std::optional<Point> line_point(const Point& a, const Point& b,
                                       const Point& c, const Point& d) {
    const auto intersection = CGAL::intersection(Line(a, b), Line(c, d));
    if (!intersection) return std::nullopt;
    if (const Point* point = std::get_if<Point>(&*intersection)) return *point;
    return std::nullopt;
}

inline char seg_seg_int(const Point& a, const Point& b,
                        const Point& c, const Point& d,
                        Point& p, Point& q) {
    int d1 = area_sign(c, d, a);
    int d2 = area_sign(c, d, b);
    int d3 = area_sign(a, b, c);
    int d4 = area_sign(a, b, d);

    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
        ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
        const auto intersection = line_point(a, b, c, d);
        if (!intersection) return '0';
        p = *intersection;
        return '1';
    }

    if (d1 == 0 && d2 == 0 && d3 == 0 && d4 == 0) {
        return parallel_int(a, b, c, d, p, q);
    }

    if (d1 == 0 && between(c, d, a)) { p = a; return 'v'; }
    if (d2 == 0 && between(c, d, b)) { p = b; return 'v'; }
    if (d3 == 0 && between(a, b, c)) { p = c; return 'v'; }
    if (d4 == 0 && between(a, b, d)) { p = d; return 'v'; }
    return '0';
}

inline void append_unique(std::vector<Point>& out, const Point& point) {
    if (out.empty() || out.back() != point) out.push_back(point);
}

inline int advance(int a, int& aa, int n, bool inside,
                   const Point& v, std::vector<Point>& out) {
    if (inside) append_unique(out, v);
    aa++;
    return (a + 1) % n;
}

inline InFlag in_out(const Point& /*p*/, InFlag inflag, int aHB, int bHA) {
    if      (aHB > 0) return InFlag::Pin;
    else if (bHA > 0) return InFlag::Qin;
    return inflag;
}

inline bool point_in_convex_poly(const Point& p, const std::vector<Point>& poly) {
    int n = (int)poly.size();
    for (int i = 0; i < n; ++i) {
        const Point& a = poly[i];
        const Point& b = poly[(i + 1) % n];
        if (CGAL::orientation(a, b, p) == CGAL::RIGHT_TURN) return false;
    }
    return true;
}

inline bool all_points_in_convex_poly(const std::vector<Point>& inner,
                                      const std::vector<Point>& outer) {
    for (const auto& v : inner) {
        if (!point_in_convex_poly(v, outer)) return false;
    }
    return true;
}

// --- Tolerance-aware orientation for INEXACTLY-constructed points ----------
// O'Rourke emits intersection vertices via Epick line-line constructions, so a
// vertex that lies mathematically ON an input edge can land a few ULPs to
// either side.  Exact-sign predicates (CGAL::orientation) then misclassify it,
// reading a real edge as a spurious kink and discarding a valid intersection.
//
// We classify the turn by the PERPENDICULAR DISTANCE of c from the supporting
// line of a->b (= cross / |a->b|), compared against an ABSOLUTE band in
// coordinate units.  A relative band (REL*|cross|) was tried first but fails on
// the short edges of the result polygon: there |cross| ~ e*noise while the band
// ~ REL*e^2, so for e below ~100 units construction noise (~1e-7) exceeds the
// band and a boundary-resting vertex reads as a definite kink.  An absolute
// perpendicular band is scale-stable: genuine kinks from ill-conditioned
// constructions deviate >=1 unit, boundary noise is <1e-6, and TURN_TOL sits in
// the empty gap between them.
inline constexpr double TURN_TOL = 1e-3;

inline int tol_turn_sign(const Point& a, const Point& b, const Point& c) {
    const double abx = CGAL::to_double(b.x()) - CGAL::to_double(a.x());
    const double aby = CGAL::to_double(b.y()) - CGAL::to_double(a.y());
    const double acx = CGAL::to_double(c.x()) - CGAL::to_double(a.x());
    const double acy = CGAL::to_double(c.y()) - CGAL::to_double(a.y());
    const double cross = abx * acy - aby * acx;
    const double len = std::sqrt(abx * abx + aby * aby);
    if (len < 1e-9) return 0;               // degenerate edge: no turn info
    const double perp = cross / len;         // signed distance of c from line ab
    if (perp >  TURN_TOL) return  1;         // clear left turn
    if (perp < -TURN_TOL) return -1;         // clear right turn
    return 0;                                // (near-)collinear within band
}

// Absolute perpendicular-distance tolerance (in coordinate units) for the
// containment guard.  Empirically the O'Rourke result vertices split into two
// cleanly separated groups: correct ones sit <1e-6 outside either input (pure
// boundary rounding noise), while genuinely ill-conditioned line-line
// constructions poke >=1 unit out; the 1e-6..1e-3 band is empty.  A threshold
// in that gap rejects the ~40k false positives (which a relative cross-product
// band could not, because F's axis-aligned bbox edges make the relative band
// collapse to REL*|cross| and flag even a 1e-12 excursion) while still catching
// every real failure with a >=3-order margin on both sides.
inline constexpr double CONTAINMENT_TOL = 1e-3;

// Point-in-CCW-convex-polygon by absolute perpendicular distance: a point
// counts as outside only when it lies more than CONTAINMENT_TOL beyond some
// edge's supporting line.  Degenerate edges (near-zero length) carry no
// half-plane information and are skipped.
inline bool point_in_convex_poly_tol(const Point& p,
                                      const std::vector<Point>& poly) {
    const int n = (int)poly.size();
    const double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y());
    for (int i = 0; i < n; ++i) {
        const double ax = CGAL::to_double(poly[i].x()), ay = CGAL::to_double(poly[i].y());
        const double bx = CGAL::to_double(poly[(i + 1) % n].x()), by = CGAL::to_double(poly[(i + 1) % n].y());
        const double ex = bx - ax, ey = by - ay;
        const double len = std::sqrt(ex * ex + ey * ey);
        if (len < 1e-9) continue;  // degenerate edge: no constraint
        const double perp = (ex * (py - ay) - ey * (px - ax)) / len;
        if (perp < -CONTAINMENT_TOL) return false;  // clearly outside
    }
    return true;
}

inline bool all_points_in_convex_poly_tol(const std::vector<Point>& inner,
                                          const std::vector<Point>& outer) {
    for (const auto& v : inner) {
        if (!point_in_convex_poly_tol(v, outer)) return false;
    }
    return true;
}

// Collapse consecutive near-coincident vertices (including the wrap-around
// pair).  Near-duplicate points accumulate in the stab chain S across
// iterations and reach the intersector via find_F, producing edges only a few
// ULPs long (~1e-13 at these coordinate scales).  O'Rourke's ConvexIntersect
// derives its advance direction from the edge vector P[a]-P[a-1]; on such a
// near-zero edge that vector's exact cross-product sign is meaningless rounding
// noise, so the walk takes the wrong branch and silently drops the following
// vertex.  Removing these degenerate edges up front eliminates the failure at
// its source (Sutherland-Hodgman is immune because it never derives a
// direction from a subject edge, which is why it disagreed with O'Rourke).
// EPS2 is a squared distance: 1e-12 == (1e-6)^2, far above the ~1e-13 duplicate
// noise yet far below any real feature (adjacent distinct vertices are >=1
// apart, i.e. squared >=1).
inline std::vector<Point> dedup_consecutive(const std::vector<Point>& poly) {
    const int n = static_cast<int>(poly.size());
    if (n < 2) return poly;
    constexpr double EPS2 = 1e-12;
    auto close = [](const Point& a, const Point& b) {
        const double dx = CGAL::to_double(a.x()) - CGAL::to_double(b.x());
        const double dy = CGAL::to_double(a.y()) - CGAL::to_double(b.y());
        return dx * dx + dy * dy <= EPS2;
    };
    std::vector<Point> out;
    out.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (!out.empty() && close(poly[i], out.back())) continue;
        out.push_back(poly[i]);
    }
    while (out.size() >= 2 && close(out.front(), out.back())) out.pop_back();
    return out;
}

inline std::optional<Point> line_intersect(const Point& a, const Point& b,
                                           const Point& p1, const Point& p2) {
    // Guard against degenerate inputs.
    if (p1.x() == p2.x() && p1.y() == p2.y()) return p1;
    if (a.x() == b.x() && a.y() == b.y()) {
        return (CGAL::orientation(a, b, p1) != CGAL::RIGHT_TURN)
                   ? std::optional<Point>(p1) : std::nullopt;
    }

    // The clip edge defines an INFINITE line, not a finite segment (a, b)
    // — using Segment_2 here would drop vertices that cross the line
    // outside [a, b].  Epick's Line_2 intersection produces a
    // double-precision result; the rounding is bounded by 1-2 ulps and
    // doesn't flip the orientation predicate for non-degenerate inputs.
    Line line_ab(a, b);
    Line line_p(p1, p2);
    auto inter = CGAL::intersection(line_ab, line_p);
    if (!inter) return std::nullopt;
    if (const Point* ip = std::get_if<Point>(&*inter)) return *ip;
    return std::nullopt;
}

inline std::vector<Point> clip_halfplane(const std::vector<Point>& subject,
                                          const Point& a, const Point& b) {
    int n = (int)subject.size();
    if (n == 0) return {};

    // Epick orientation is exact in sign — half-plane classification
    // never misclassifies a vertex.  This avoids the per-call Epeck
    // conversions the previous version paid for.
    std::vector<Point> out;
    out.reserve(n + 1);

    auto o_prev = CGAL::orientation(a, b, subject[(n - 1) % n]);
    bool prev_in = (o_prev != CGAL::RIGHT_TURN);

    for (int i = 0; i < n; ++i) {
        const Point& cur  = subject[i];
        auto o_cur  = CGAL::orientation(a, b, cur);
        bool cur_in = (o_cur != CGAL::RIGHT_TURN);

        if (cur_in) {
            if (!prev_in) {
                if (auto ip = line_intersect(a, b, subject[(i + n - 1) % n], cur))
                    out.push_back(*ip);
            }
            out.push_back(cur);
        } else if (prev_in) {
            if (auto ip = line_intersect(a, b, subject[(i + n - 1) % n], cur))
                out.push_back(*ip);
        }
        prev_in = cur_in;
    }
    return out;
}

inline std::vector<Point> convex_intersect_fast(const std::vector<Point>& Pin,
                                                const std::vector<Point>& Qin) {
    if (Pin.size() < 3 || Qin.size() < 3) return {};

    auto kernel_signed_area = [](const std::vector<Point>& P) {
        double s = 0.0;
        int sz = (int)P.size();
        for (int i = 0; i < sz; ++i) {
            const Point& a = P[i];
            const Point& b = P[(i + 1) % sz];
            s += CGAL::to_double(a.x() * b.y() - b.x() * a.y());
        }
        return s / 2.0;
    };

    std::vector<Point> P, Q;
    P.reserve(Pin.size());
    Q.reserve(Qin.size());
    if (kernel_signed_area(Pin) < 0) { P.assign(Pin.rbegin(), Pin.rend()); }
    else                            { P = Pin; }
    if (kernel_signed_area(Qin) < 0) { Q.assign(Qin.rbegin(), Qin.rend()); }
    else                            { Q = Qin; }

    if (point_in_convex_poly(P[0], Q) && all_points_in_convex_poly(P, Q)) return P;
    if (point_in_convex_poly(Q[0], P) && all_points_in_convex_poly(Q, P)) return Q;

    std::vector<Point> running = P;
    int m = (int)Q.size();
    for (int i = 0; i < m; ++i) {
        const Point& a = Q[i];
        const Point& b = Q[(i + 1) % m];
        running = clip_halfplane(running, a, b);
        if (running.size() < 3) return {};
    }
    return running;
}

// Sutherland-Hodgman convex polygon clipping.  Clips subject against the
// half-planes defined by each edge of clip.  Both inputs must be CCW (which
// the callers — convex_hull_2 and find_F — guarantee).  Returns CCW result
// or empty vector if the polygons don't intersect.
//
// S-H is O(n*m) where n=subject.size() and m=clip.size(), the same asymptotic
// as O'Rourke, but each step is a single orientation test + possibly one
// line intersection (no segment-segment code, no inflection-point tracking,
// no loop guards).  Empirically ~3-4x faster than O'Rourke on small (~5-12
// vertex) polygons.
inline bool convex_intersect(const std::vector<Point>& subject,
                             const std::vector<Point>& clip,
                             std::vector<Point>& result,
                             bool& numerically_reliable) {
    result.clear();
    numerically_reliable = true;
    const int n = (int)subject.size();
    const int m = (int)clip.size();
    if (n < 3 || m < 3) return false;

    double subject_min_x = CGAL::to_double(subject[0].x());
    double subject_max_x = subject_min_x;
    double subject_min_y = CGAL::to_double(subject[0].y());
    double subject_max_y = subject_min_y;
    for (int i = 1; i < n; ++i) {
        const double x = CGAL::to_double(subject[i].x());
        const double y = CGAL::to_double(subject[i].y());
        subject_min_x = std::min(subject_min_x, x);
        subject_max_x = std::max(subject_max_x, x);
        subject_min_y = std::min(subject_min_y, y);
        subject_max_y = std::max(subject_max_y, y);
    }

    double clip_min_x = CGAL::to_double(clip[0].x());
    double clip_max_x = clip_min_x;
    double clip_min_y = CGAL::to_double(clip[0].y());
    double clip_max_y = clip_min_y;
    for (int i = 1; i < m; ++i) {
        const double x = CGAL::to_double(clip[i].x());
        const double y = CGAL::to_double(clip[i].y());
        clip_min_x = std::min(clip_min_x, x);
        clip_max_x = std::max(clip_max_x, x);
        clip_min_y = std::min(clip_min_y, y);
        clip_max_y = std::max(clip_max_y, y);
    }
    if (subject_max_x < clip_min_x || clip_max_x < subject_min_x ||
        subject_max_y < clip_min_y || clip_max_y < subject_min_y) {
        return false;
    }

    // The intersection of convex n- and m-gons has at most n+m vertices.
    // Keep the common small-polygon path entirely on the stack.  The fallback
    // is rare and favors correctness over speed for unusually large inputs.
    constexpr int STACK_CAPACITY = 128;
    if (n + m > STACK_CAPACITY) {
        result = subject;
        for (int e = 0; e < m && result.size() >= 3; ++e) {
            result = clip_halfplane(result, clip[e], clip[(e + 1) % m]);
        }
        return result.size() >= 3;
    }

    std::array<double, STACK_CAPACITY> cx, cy;
    std::array<double, STACK_CAPACITY> x0, y0, x1, y1;
    for (int i = 0; i < n; ++i) {
        x0[i] = CGAL::to_double(subject[i].x());
        y0[i] = CGAL::to_double(subject[i].y());
    }
    for (int i = 0; i < m; ++i) {
        cx[i] = CGAL::to_double(clip[i].x());
        cy[i] = CGAL::to_double(clip[i].y());
    }

    double* ix = x0.data();
    double* iy = y0.data();
    double* ox = x1.data();
    double* oy = y1.data();
    int in_size = n;

    for (int e = 0; e < m; ++e) {
        const int e1 = (e + 1) % m;
        const double ax = cx[e], ay = cy[e];
        const double bx = cx[e1], by = cy[e1];
        // Inside test: for CCW clip polygon, point is inside half-plane iff
        // (b - a) x (p - a) >= 0, i.e. (bx-ax)*(py-ay) - (by-ay)*(px-ax) >= 0.
        auto inside = [&](double px, double py) -> bool {
            const double ux = bx - ax;
            const double uy = by - ay;
            const double vx = px - ax;
            const double vy = py - ay;
            const double determinant = ux * vy - uy * vx;
            const double error_bound = 8.0 * std::numeric_limits<double>::epsilon() *
                                       (std::abs(ux * vy) + std::abs(uy * vx));
            if (std::abs(determinant) <= error_bound) {
                return CGAL::orientation(Point(ax, ay), Point(bx, by),
                                         Point(px, py)) != CGAL::RIGHT_TURN;
            }
            return determinant > 0.0;
        };
        // Line-line intersection of segment (S_{k-1} -> S_k) with the clip
        // edge (a -> b).  Both segments are finite; we resolve t in [0,1].
        // Returns 1 if both denominators are non-degenerate and t in (0, 1).
        auto intersect = [&](double p1x, double p1y,
                             double p2x, double p2y,
                             double& ox_out, double& oy_out) -> bool {
            const double rx = p2x - p1x, ry = p2y - p1y;
            const double sx2 = bx - ax, sy2 = by - ay;
            const double denom_lhs = rx * sy2;
            const double denom_rhs = ry * sx2;
            const double denom = denom_lhs - denom_rhs;
            const double error_bound = 8.0 * std::numeric_limits<double>::epsilon() *
                                       (std::abs(denom_lhs) + std::abs(denom_rhs));
            const double t_numerator = (ax - p1x) * sy2 - (ay - p1y) * sx2;
            if (std::abs(denom) > error_bound) {
                const double t = t_numerator / denom;
                const double t_error = 16.0 * std::numeric_limits<double>::epsilon();
                if (t > t_error && t < 1.0 - t_error) {
                    ox_out = p1x + t * rx;
                    oy_out = p1y + t * ry;
                    return true;
                }
            }

            const auto exact_intersection = CGAL::intersection(
                Line(Point(ax, ay), Point(bx, by)),
                Line(Point(p1x, p1y), Point(p2x, p2y)));
            if (!exact_intersection) return false;
            if (const Point* point = std::get_if<Point>(&*exact_intersection)) {
                ox_out = CGAL::to_double(point->x());
                oy_out = CGAL::to_double(point->y());
                return true;
            }
            return false;
        };

        int out_size = 0;
        if (in_size == 0) return false;
        auto append_unique = [&](double px, double py) {
            if (out_size > 0 && ox[out_size - 1] == px && oy[out_size - 1] == py) return;
            ox[out_size] = px;
            oy[out_size] = py;
            ++out_size;
        };

        double prev_x = ix[in_size - 1], prev_y = iy[in_size - 1];
        bool prev_in = inside(prev_x, prev_y);

        for (int k = 0; k < in_size; ++k) {
            const double cur_x = ix[k], cur_y = iy[k];
            const bool cur_in = inside(cur_x, cur_y);

            if (cur_in) {
                if (!prev_in) {
                    double ix_o, iy_o;
                    if (intersect(prev_x, prev_y, cur_x, cur_y, ix_o, iy_o)) {
                        append_unique(ix_o, iy_o);
                    }
                }
                append_unique(cur_x, cur_y);
            } else if (prev_in) {
                double ix_o, iy_o;
                if (intersect(prev_x, prev_y, cur_x, cur_y, ix_o, iy_o)) {
                    append_unique(ix_o, iy_o);
                }
            }
            prev_x = cur_x; prev_y = cur_y; prev_in = cur_in;
        }

        if (out_size < 3) return false;
        std::swap(ix, ox);
        std::swap(iy, oy);
        in_size = out_size;
    }

    result.reserve(in_size);
    for (int i = 0; i < in_size; ++i) {
        result.emplace_back(ix[i], iy[i]);
    }
    return true;
}

inline std::vector<Point> convex_intersect_robust(const std::vector<Point>& Pin,
                                                   const std::vector<Point>& Qin) {
    const int n = static_cast<int>(Pin.size());
    const int m = static_cast<int>(Qin.size());
    if (n < 3 || m < 3) return {};

    auto leftmost_index = [](const std::vector<Point>& polygon) {
        int best = 0;
        for (int i = 1; i < static_cast<int>(polygon.size()); ++i) {
            if (CGAL::compare_y(polygon[i], polygon[best]) == CGAL::SMALLER ||
                (CGAL::compare_y(polygon[i], polygon[best]) == CGAL::EQUAL &&
                 CGAL::compare_x(polygon[i], polygon[best]) == CGAL::SMALLER)) {
                best = i;
            }
        }
        return best;
    };

    std::vector<Point> P = Pin;
    std::vector<Point> Q = Qin;
    if (kernel_signed_area(P) < 0) std::reverse(P.begin(), P.end());
    if (kernel_signed_area(Q) < 0) std::reverse(Q.begin(), Q.end());

    int a = leftmost_index(P);
    int b = leftmost_index(Q);
    int aa = 0;
    int ba = 0;
    InFlag inflag = InFlag::Unknown;
    bool first_point = true;
    Point p0;
    std::vector<Point> out;

    do {
        const int a1 = (a + n - 1) % n;
        const int b1 = (b + m - 1) % m;
        const Vector A = P[a] - P[a1];
        const Vector B = Q[b] - Q[b1];
        const int cross = area_sign_vec(A, B);
        const int aHB = area_sign(Q[b1], Q[b], P[a]);
        const int bHA = area_sign(P[a1], P[a], Q[b]);

        Point p, q;
        const char code = seg_seg_int(P[a1], P[a], Q[b1], Q[b], p, q);
        if (code == '1' || code == 'v') {
            if (inflag == InFlag::Unknown && first_point) {
                aa = ba = 0;
                first_point = false;
                p0 = p;
            }
            append_unique(out, p);
            inflag = in_out(p, inflag, aHB, bHA);
        }

        if (code == 'e' && dot_sign(P[a1], P[a], Q[b1], Q[b]) < 0) {
            append_unique(out, p);
            append_unique(out, q);
            return out;
        }
        if (cross == 0 && aHB < 0 && bHA < 0) return {};

        if (cross == 0 && aHB == 0 && bHA == 0) {
            if (inflag == InFlag::Pin) {
                b = advance(b, ba, m, inflag == InFlag::Qin, Q[b], out);
            } else {
                a = advance(a, aa, n, inflag == InFlag::Pin, P[a], out);
            }
        } else if (cross >= 0) {
            if (bHA > 0) {
                a = advance(a, aa, n, inflag == InFlag::Pin, P[a], out);
            } else {
                b = advance(b, ba, m, inflag == InFlag::Qin, Q[b], out);
            }
        } else if (aHB > 0) {
            b = advance(b, ba, m, inflag == InFlag::Qin, Q[b], out);
        } else {
            a = advance(a, aa, n, inflag == InFlag::Pin, P[a], out);
        }
    } while (((aa < n) || (ba < m)) && aa < 2 * n && ba < 2 * m);

    if (first_point) {
        if (point_in_convex_poly(P[0], Q)) return P;
        if (point_in_convex_poly(Q[0], P)) return Q;
        return {};
    }

    append_unique(out, p0);
    if (out.size() > 1 && out.front() == out.back()) out.pop_back();
    return out;
}

}  // namespace orourke_cgal

// ===========================================================================
//  Exact convex-intersection fallback (EPECK)
// ===========================================================================
//
// The fast Epick O'Rourke path (below) is correct for well-conditioned inputs
// but its double-precision line-line constructions break down on near-
// degenerate configurations (e.g. a huge bbox-corner polygon clipped by a tiny
// grid cell): it can emit a vertex hundreds of units outside the true region.
// When the fast path's validation detects such output, intersect() re-runs the
// clip here with exact constructions.  Sutherland-Hodgman is used because,
// for convex inputs, it cannot produce the kinks/self-intersections O'Rourke
// can, and in EPECK every keep/drop decision and crossing point is exact.  The
// final vertices are converted back to Epick doubles (topology is exact; the
// residual coordinate error is ~1e-12, far below the ~1e-6 noise the fast path
// exhibits).
inline bool intersect_exact_convex(const std::vector<Point>& P_verts,
                                   const std::vector<Point>& Q_verts,
                                   std::vector<Point>& result) {
    result.clear();
    const int n = static_cast<int>(P_verts.size());
    const int m = static_cast<int>(Q_verts.size());
    if (n < 3 || m < 3) return false;

    using EPoint = Epeck::Point_2;
    using ELine  = Epeck::Line_2;

    std::vector<EPoint> subject, clip;
    subject.reserve(n);
    clip.reserve(m);
    for (const auto& p : P_verts) subject.push_back(conv_to_exact(p));
    for (const auto& q : Q_verts) clip.push_back(conv_to_exact(q));

    auto signed_area2 = [](const std::vector<EPoint>& v) {
        Epeck::FT s = 0;
        const int sz = static_cast<int>(v.size());
        for (int i = 0; i < sz; ++i) {
            const EPoint& a = v[i];
            const EPoint& b = v[(i + 1) % sz];
            s += a.x() * b.y() - b.x() * a.y();
        }
        return s;
    };
    if (signed_area2(subject) < 0) std::reverse(subject.begin(), subject.end());
    if (signed_area2(clip) < 0) std::reverse(clip.begin(), clip.end());

    std::vector<EPoint> out = subject;
    for (int e = 0; e < m && out.size() >= 3; ++e) {
        const EPoint& a = clip[e];
        const EPoint& b = clip[(e + 1) % m];
        const ELine edge(a, b);
        std::vector<EPoint> in;
        in.swap(out);
        const int sz = static_cast<int>(in.size());
        for (int i = 0; i < sz; ++i) {
            const EPoint& cur = in[i];
            const EPoint& prv = in[(i + sz - 1) % sz];
            const bool cur_in = CGAL::orientation(a, b, cur) != CGAL::RIGHT_TURN;
            const bool prv_in = CGAL::orientation(a, b, prv) != CGAL::RIGHT_TURN;
            if (cur_in) {
                if (!prv_in) {
                    const auto ip = CGAL::intersection(edge, ELine(prv, cur));
                    if (ip) if (const EPoint* pp = std::get_if<EPoint>(&*ip)) out.push_back(*pp);
                }
                out.push_back(cur);
            } else if (prv_in) {
                const auto ip = CGAL::intersection(edge, ELine(prv, cur));
                if (ip) if (const EPoint* pp = std::get_if<EPoint>(&*ip)) out.push_back(*pp);
            }
        }
    }
    if (out.size() < 3) { result.clear(); return false; }

    result.reserve(out.size());
    for (const auto& p : out) {
        result.emplace_back(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
    }
    if (result.size() > 1 && result.front() == result.back()) result.pop_back();
    if (result.size() < 3) { result.clear(); return false; }
    return true;
}

// ===========================================================================
//  Convex-convex intersection by double-precision Sutherland-Hodgman
// ===========================================================================
//
// Both inputs are convex (F(S,p) is a convex wedge clipped to the axis-aligned
// bbox; conv(G) is a convex grid-cell polygon), so their intersection is the
// subject polygon successively clipped by each half-plane of the clip polygon.
//
// Why this replaces O'Rourke's ConvexIntersect:  O'Rourke's dual-advance walk
// derives its next move from the *sign of a cross product of two edge
// directions* and, on the near-degenerate slivers F produces (a cap only a few
// units wide attached to two edges thousands of units long, reaching the bbox
// border), that sign is dominated by rounding and the walk takes the wrong
// branch — it then traces the *clip* polygon and emits a region far larger
// than the true intersection (observed: 18 vertices tracing Q where the exact
// answer is a 9-vertex cap 500+ units smaller).  Sutherland-Hodgman cannot
// exhibit this failure mode by construction: it only ever *removes* subject
// area against each clip half-plane, so the output is always a subset of the
// subject and can never balloon.  A near-collinear misclassification of the
// in/out side test can at worst displace a boundary vertex by construction
// noise (~1e-9 here), never change the topology.  Because it is structurally
// robust it needs no convexity guard, no containment guard, and no exact
// fallback — which is exactly where the previous version spent its GMP time.
//
// Complexity is O(n*m) but both polygons are tiny (n,m ~ 4..20), and every
// operation is a plain double, so in practice it is far faster than the
// O(n+m) O'Rourke walk once that walk's exact-predicate/fallback cost is
// included.
namespace sh_double {

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

// Clip convex subject P by convex clip Q; both given as Epick points.  Returns
// the intersection polygon (CCW, no closing duplicate) in Epick points, or an
// empty vector when the intersection is degenerate (<3 vertices).
inline std::vector<Point> clip(const std::vector<Point>& P_verts,
                               const std::vector<Point>& Q_verts) {
    const int n = static_cast<int>(P_verts.size());
    const int m = static_cast<int>(Q_verts.size());
    if (n < 3 || m < 3) return {};

    std::vector<std::array<double, 2>> subj, clipv;
    subj.reserve(n);
    clipv.reserve(m);
    for (const auto& p : P_verts)
        subj.push_back({CGAL::to_double(p.x()), CGAL::to_double(p.y())});
    for (const auto& q : Q_verts)
        clipv.push_back({CGAL::to_double(q.x()), CGAL::to_double(q.y())});
    to_ccw(subj);
    to_ccw(clipv);

    std::vector<std::array<double, 2>> out = subj;
    std::vector<std::array<double, 2>> in;
    for (int e = 0; e < m && out.size() >= 3; ++e) {
        const double ax = clipv[e][0],           ay = clipv[e][1];
        const double bx = clipv[(e + 1) % m][0], by = clipv[(e + 1) % m][1];
        const double ex = bx - ax, ey = by - ay;
        // Signed area of triangle (a,b,pt); >=0 means pt is on/left of edge a->b
        // (inside, since the clip is CCW).
        auto side = [&](const std::array<double, 2>& pt) {
            return ex * (pt[1] - ay) - ey * (pt[0] - ax);
        };
        in.swap(out);
        out.clear();
        const int sz = static_cast<int>(in.size());
        for (int i = 0; i < sz; ++i) {
            const auto& cur = in[i];
            const auto& prv = in[(i + sz - 1) % sz];
            const double sc = side(cur), sp = side(prv);
            const bool cin = sc >= 0.0, pin = sp >= 0.0;
            if (cin) {
                if (!pin) {
                    const double t = sp / (sp - sc);
                    out.push_back({prv[0] + t * (cur[0] - prv[0]),
                                   prv[1] + t * (cur[1] - prv[1])});
                }
                out.push_back(cur);
            } else if (pin) {
                const double t = sp / (sp - sc);
                out.push_back({prv[0] + t * (cur[0] - prv[0]),
                               prv[1] + t * (cur[1] - prv[1])});
            }
        }
    }
    if (out.size() < 3) return {};

    std::vector<Point> result;
    result.reserve(out.size());
    for (const auto& p : out) result.emplace_back(p[0], p[1]);
    return result;
}

}  // namespace sh_double

// ===========================================================================
//  Intersection (public API)
// ===========================================================================

#ifdef PROFILE_INTERSECT
inline long g_intersect_calls = 0;
inline long g_fallback_calls = 0;
#endif

inline bool intersect(const std::vector<Point>& P_in,
                      const std::vector<Point>& Q_in,
                      std::vector<Point>& result) {
#ifdef PROFILE_INTERSECT
    g_intersect_calls++;
#endif
    // Remove degenerate (near-zero-length) edges: collinear duplicate vertices
    // carry no half-plane information and only add spurious near-zero-length
    // clip edges.
    const std::vector<Point> P_verts = orourke_cgal::dedup_consecutive(P_in);
    const std::vector<Point> Q_verts = orourke_cgal::dedup_consecutive(Q_in);

    result = orourke_cgal::dedup_consecutive(sh_double::clip(P_verts, Q_verts));
    if (result.size() < 3) {
        result.clear();
        return false;
    }
    return true;
}

inline bool intersect(const Polygon& subject, const Polygon& clip,
                      std::vector<Point>& result) {
    std::vector<Point> P_verts(subject.vertices_begin(), subject.vertices_end());
    std::vector<Point> Q_verts(clip.vertices_begin(), clip.vertices_end());
    return intersect(P_verts, Q_verts, result);
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
// Bit-identical to CGAL::orientation but skips GMP on the near-collinear stab
// vertices that otherwise dominate this test.
inline bool point_in_convex(const Point& p, const std::vector<Point>& poly, bool ccw = true) {
    const int want = ccw ? 1 : -1;   // sign that means "outside"
    const double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y());
    const int n = static_cast<int>(poly.size());
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
        if (s == want) return false;
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
        const double corner_x = ji * GRID;
        const double corner_y = ki * GRID;
        expected_frechet_squared = std::max(
            expected_frechet_squared, corner_x * corner_x + corner_y * corner_y);
        points.emplace_back(px + corner_x, py + corner_y);
    };

    for (int j = j_min; j <= j_max; ++j) {
        const double x0 = j * GRID, x1 = (j + 1) * GRID;
        for (int k = j_min; k <= j_max; ++k) {
            const double y0 = k * GRID, y1 = (k + 1) * GRID;
            const double nearest_x = x0 > 0.0 ? x0 : (x1 < 0.0 ? x1 : 0.0);
            const double nearest_y = y0 > 0.0 ? y0 : (y1 < 0.0 ? y1 : 0.0);
            if (nearest_x * nearest_x + nearest_y * nearest_y > r2) continue;
            add_corner(j,     k);
            add_corner(j + 1, k);
            add_corner(j + 1, k + 1);
            add_corner(j,     k + 1);
        }
    }

    return points;
}

inline std::vector<Point> get_conv_from_grid(const Point& p, double EPSILON, double DELTA) {
    std::vector<Point> points = get_points_from_grid(p, EPSILON, DELTA);
    std::vector<Point> conv;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(conv));
    return conv;
}

// The two tangent (supporting) vertices from external point p to a convex
// polygon S.  The two supporting vertices are the vertices of EXTREME bearing
// as seen from p: since p lies outside the convex hull, all of S falls within
// an angular wedge of less than 180 degrees, inside which bearing is a total
// order, so a single linear scan keeping the most-clockwise and
// most-counterclockwise vertices finds both tangents in O(n).
//
// The orientation sign is computed in pure doubles with no exact fallback: for
// tangent finding an exact result on (near-)collinear ties is unnecessary
// because if two candidate vertices are collinear as seen from p, either is an
// equally valid supporting vertex (the wedge F's boundary ray passes through
// that collinear point regardless, so F is the same point set either way).
inline std::vector<int> find_tangent_idx(const Point& p, const std::vector<Point>& S) {
    std::vector<int> tangent;
    tangent.reserve(2);
    const int n = static_cast<int>(S.size());

    // Cache p and S coordinates as doubles once (S is scanned 2n times).
    static std::vector<double> sx, sy;
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
        tangent.push_back(std::min(rt, lt));
        tangent.push_back(std::max(rt, lt));
    }
    return tangent;
}

inline std::optional<Point> intersect_ray_with_rect(const Point& p, const Point& direction) {
    if (CGAL::to_double(CGAL::squared_distance(p, direction)) < GEOM_TOL * GEOM_TOL)
        return std::nullopt;

    Ray ray(p, direction);
    Bbox box(BMIN, BMIN, BMAX, BMAX);

    if (auto obj = CGAL::intersection(Rect(box), ray)) {
        if (const Point* ip = std::get_if<Point>(&*obj)) {
            if (CGAL::to_double(CGAL::squared_distance(*ip, p)) < GEOM_TOL * GEOM_TOL)
                return std::nullopt;
            return *ip;
        }
        if (const Segment* seg = std::get_if<Segment>(&*obj)) {
            const Point& a = seg->source();
            const Point& b = seg->target();
            double da = CGAL::to_double(CGAL::squared_distance(a, p));
            double db = CGAL::to_double(CGAL::squared_distance(b, p));
            if (da < GEOM_TOL * GEOM_TOL) return b;
            if (db < GEOM_TOL * GEOM_TOL) return a;
            return (da < db) ? a : b;
        }
    }
    return std::nullopt;
}

// The free-space wedge F(S,p): the region of the bbox reachable from the
// convex stab region S through the external point p.  When p is inside S (or S
// is a single point) every bbox point is reachable, so F is the whole bbox;
// otherwise F is bounded by the two tangent rays from p to S, the arc of S
// between the tangent vertices, and the bbox boundary between the two ray hits.
inline void find_F(const Point& p, const std::vector<Point>& S,
                   std::vector<Point>& F) {
    F.clear();
    assert(S.size() != 2);
    if (S.size() == 1 || point_in_convex(p, S)) {
        F = current_bbox();
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
        F = current_bbox();
        return;
    }

    auto e1 = which_edge(hit1.value());
    auto e2 = which_edge(hit2.value());
    if (!e1 || !e2) {
        std::cerr << "Cannot determine which Bbox edge the ray intersect with.\n";
        F = current_bbox();
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

#endif // SIMPLIFY_GEOMETRY_H
