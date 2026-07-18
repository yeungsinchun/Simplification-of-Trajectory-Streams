#ifndef SIMPLIFY_GEOMETRY_H
#define SIMPLIFY_GEOMETRY_H

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>

#include <array>
#include <cmath>
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

inline int area_sign(const Point& a, const Point& b, const Point& c) {
    switch (CGAL::orientation(a, b, c)) {
        case CGAL::LEFT_TURN:  return  1;
        case CGAL::RIGHT_TURN: return -1;
        case CGAL::COLLINEAR:  return  0;
    }
    return 0;
}

inline int area_sign_vec(const Vector& a, const Vector& b) {
    auto c = CGAL::determinant(a, b);
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

inline Point line_point(const Point& a, const Point& b,
                        const Point& c, const Point& d) {
    double ax = CGAL::to_double(a.x()), ay = CGAL::to_double(a.y());
    double bx = CGAL::to_double(b.x()), by = CGAL::to_double(b.y());
    double cx = CGAL::to_double(c.x()), cy = CGAL::to_double(c.y());
    double dx = CGAL::to_double(d.x()), dy = CGAL::to_double(d.y());
    double rx = bx - ax, ry = by - ay;
    double sx = dx - cx, sy = dy - cy;
    double denom = rx * sy - ry * sx;
    double t = ((cx - ax) * sy - (cy - ay) * sx) / denom;
    return Point(ax + t * rx, ay + t * ry);
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
        p = line_point(a, b, c, d);
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

inline int advance(int a, int& aa, int n, bool inside,
                   const Point& v, std::vector<Point>& out) {
    if (inside) out.push_back(v);
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
            out.push_back(p);
            inflag = in_out(p, inflag, aHB, bHA);
        }

        if (code == 'e' && dot_sign(P[a1], P[a], Q[b1], Q[b]) < 0) {
            out.push_back(p);
            out.push_back(q);
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

    out.push_back(p0);
    return out;
}

}  // namespace orourke_cgal

// ===========================================================================
//  Intersection (public API — calls into orourke_cgal::convex_intersect)
// ===========================================================================

inline bool intersect(const std::vector<Point>& P_verts,
                      const std::vector<Point>& Q_verts,
                      std::vector<Point>& result) {
    bool numerically_reliable = false;
    const bool fast_intersects = orourke_cgal::convex_intersect(
        P_verts, Q_verts, result, numerically_reliable);
    if (!numerically_reliable) {
        result = orourke_cgal::convex_intersect_robust(P_verts, Q_verts);
    } else if (!fast_intersects) {
        result.clear();
        return false;
    }

    if (!result.empty() && result.front() == result.back()) result.pop_back();
    if (result.size() < 3) {
        result.clear();
        return false;
    }

    Polygon polygon(result.begin(), result.end());
    if (!polygon.is_simple()) {
        result = orourke_cgal::convex_intersect_robust(P_verts, Q_verts);
        if (!result.empty() && result.front() == result.back()) result.pop_back();
        if (result.size() < 3) {
            result.clear();
            return false;
        }
        polygon = Polygon(result.begin(), result.end());
        if (!polygon.is_simple()) {
            result.clear();
            return false;
        }
    }
    if (polygon.is_clockwise_oriented()) std::reverse(result.begin(), result.end());
    return true;
}

inline bool intersect(const Polygon& subject, const Polygon& clip,
                      std::vector<Point>& result) {
    std::vector<Point> P_verts(subject.vertices_begin(), subject.vertices_end());
    std::vector<Point> Q_verts(clip.vertices_begin(), clip.vertices_end());
    return intersect(P_verts, Q_verts, result);
}

#endif // SIMPLIFY_GEOMETRY_H
