#ifndef SIMPLIFY_GEOMETRY_H
#define SIMPLIFY_GEOMETRY_H

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>

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

inline std::vector<Point> convex_intersect(const std::vector<Point>& Pin,
                                           const std::vector<Point>& Qin) {
    int n = (int)Pin.size();
    int m = (int)Qin.size();
    if (n < 3 || m < 3) return {};

    auto signed_area = [](const std::vector<Point>& P) {
        return kernel_signed_area(P);
    };
    auto leftmost_index = [](const std::vector<Point>& P) {
        int best = 0;
        for (int i = 1; i < (int)P.size(); ++i) {
            if (CGAL::compare_y(P[i], P[best]) == CGAL::SMALLER ||
                (CGAL::compare_y(P[i], P[best]) == CGAL::EQUAL &&
                 CGAL::compare_x(P[i], P[best]) == CGAL::SMALLER)) {
                best = i;
            }
        }
        return best;
    };
    std::vector<Point> Pr, Qr;
    int a, b;
    Pr = Pin;
    Qr = Qin;
    if (signed_area(Pr) < 0) std::reverse(Pr.begin(), Pr.end());
    if (signed_area(Qr) < 0) std::reverse(Qr.begin(), Qr.end());
    a = leftmost_index(Pr);
    b = leftmost_index(Qr);

    int aa = 0, ba = 0;
    InFlag inflag = InFlag::Unknown;
    bool first_point = true;
    Point p0;

    std::vector<Point> out;

    do {
        int a1 = (a + n - 1) % n;
        int b1 = (b + m - 1) % m;

        Vector A = Pr[a] - Pr[a1];
        Vector B = Qr[b] - Qr[b1];

        int cross  = area_sign_vec(A, B);
        int aHB    = area_sign(Qr[b1], Qr[b], Pr[a]);
        int bHA    = area_sign(Pr[a1], Pr[a], Qr[b]);

        Point p, q;
        char code = seg_seg_int(Pr[a1], Pr[a], Qr[b1], Qr[b], p, q);

        if (code == '1' || code == 'v') {
            if (inflag == InFlag::Unknown && first_point) {
                aa = ba = 0;
                first_point = false;
                p0 = p;
            }
            out.push_back(p);
            inflag = in_out(p, inflag, aHB, bHA);
        }

        if ((code == 'e') && (dot_sign(Pr[a1], Pr[a], Qr[b1], Qr[b]) < 0)) {
            out.push_back(p);
            out.push_back(q);
            return out;
        }

        if ((cross == 0) && (aHB < 0) && (bHA < 0)) {
            return {};
        }

        if ((cross == 0) && (aHB == 0) && (bHA == 0)) {
            if (inflag == InFlag::Pin)
                b = advance(b, ba, m, inflag == InFlag::Qin, Qr[b], out);
            else
                a = advance(a, aa, n, inflag == InFlag::Pin, Pr[a], out);
        } else if (cross >= 0) {
            if (bHA > 0)
                a = advance(a, aa, n, inflag == InFlag::Pin, Pr[a], out);
            else
                b = advance(b, ba, m, inflag == InFlag::Qin, Qr[b], out);
        } else {
            if (aHB > 0)
                b = advance(b, ba, m, inflag == InFlag::Qin, Qr[b], out);
            else
                a = advance(a, aa, n, inflag == InFlag::Pin, Pr[a], out);
        }
    } while (((aa < n) || (ba < m)) && (aa < 2 * n) && (ba < 2 * m));

    if (first_point) {
        if (point_in_convex_poly(Pr[0], Qr)) return Pr;
        if (point_in_convex_poly(Qr[0], Pr)) return Qr;
        return {};
    }

    out.push_back(p0);
    return out;
}

}  // namespace orourke_cgal

// ===========================================================================
//  Intersection (public API — calls into orourke_cgal)
// ===========================================================================

inline void intersect(const Polygon& subject, const Polygon& clip,
                      std::back_insert_iterator<std::vector<Polygon_with_holes>> result) {
    std::vector<Point> P_verts, Q_verts;
    P_verts.reserve(subject.size());
    Q_verts.reserve(clip.size());
    for (auto it = subject.vertices_begin(); it != subject.vertices_end(); ++it) {
        P_verts.push_back(*it);
    }
    for (auto it = clip.vertices_begin(); it != clip.vertices_end(); ++it) {
        Q_verts.push_back(*it);
    }

    std::vector<Point> verts = orourke_cgal::convex_intersect_fast(P_verts, Q_verts);
    if (verts.size() < 4) return;
    if (verts.front() == verts.back())
        verts.pop_back();

    Polygon inter_poly(verts.begin(), verts.end());
    if (!inter_poly.is_simple()) return;
    if (inter_poly.is_clockwise_oriented()) inter_poly.reverse_orientation();
    *result++ = Polygon_with_holes(inter_poly);
}

inline void intersect(const std::vector<Point>& P_verts, const std::vector<Point>& Q_verts,
                      std::back_insert_iterator<std::vector<Polygon_with_holes>> result) {
    std::vector<Point> verts = orourke_cgal::convex_intersect_fast(P_verts, Q_verts);
    if (verts.size() < 4) return;
    if (verts.front() == verts.back())
        verts.pop_back();

    Polygon inter_poly(verts.begin(), verts.end());
    if (!inter_poly.is_simple()) return;
    if (inter_poly.is_clockwise_oriented()) inter_poly.reverse_orientation();
    *result++ = Polygon_with_holes(inter_poly);
}

#endif // SIMPLIFY_GEOMETRY_H
