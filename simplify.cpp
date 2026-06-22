#include <CGAL/Polygon_2.h>
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
#include "simplify_geometry.h"
#include "timer.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

// Signed polygon area.  Uses CGAL::to_double on the cross product per edge pair.
// For Epick this is exact — the cross product of two doubles is exact, and there
// are no rounding risks from cancellation in this edge-pair summation.
static double kernel_signed_area(const std::vector<Point>& P) {
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

// Exact orientation sign (-1, 0, +1) of (a, b, c).  Uses Kernel's
// exact-predicate orientation, so the sign is correct even for nearly-
// collinear input.
static int area_sign(const Point& a, const Point& b, const Point& c) {
    switch (CGAL::orientation(a, b, c)) {
        case CGAL::LEFT_TURN:  return  1;
        case CGAL::RIGHT_TURN: return -1;
        case CGAL::COLLINEAR:  return  0;
    }
    return 0;
}

// Same as area_sign but for vector triples.  Used inside the main loop
// where A and B are the *edge vectors* of the two polygons, not points.
// Kernel provides Vector_2 for point - point.
static int area_sign_vec(const Vector& a, const Vector& b) {
    auto c = CGAL::determinant(a, b);
    if (c > 0) return  1;
    if (c < 0) return -1;
    return 0;
}

// Is c on the closed segment ab?  Exact collinearity test, then exact
// dominance test on whichever coordinate has the larger dynamic range.
static bool between(const Point& a, const Point& b, const Point& c) {
    if (CGAL::orientation(a, b, c) != CGAL::COLLINEAR) return false;
    if (CGAL::abs(b.x() - a.x()) >= CGAL::abs(b.y() - a.y())) {
        return (a.x() <= c.x() && c.x() <= b.x()) ||
               (a.x() >= c.x() && c.x() >= b.x());
    }
    return (a.y() <= c.y() && c.y() <= b.y()) ||
           (a.y() >= c.y() && c.y() >= b.y());
}

// Sign of the inner product (b - a) . (d - c).  Used by the shared-edge
// special case to detect "edges overlap with opposite orientation".
static int dot_sign(const Point& a, const Point& b,
                    const Point& c, const Point& d) {
    auto ux = b.x() - a.x(), uy = b.y() - a.y();
    auto vx = d.x() - c.x(), vy = d.y() - c.y();
    auto s  = ux * vx + uy * vy;
    if (s > 0) return  1;
    if (s < 0) return -1;
    return 0;
}

// Handle parallel/collinear segments ab and cd.  Returns one of:
//   '0' - no intersection
//   'e' - collinear overlap; out_p, out_q are the two overlap endpoints
static char parallel_int(const Point& a, const Point& b,
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

// Advance the index for polygon P or Q depending on `inside`, and append
// the appropriate vertex to `out` if `inside` is true.  Returns the new index.
static int advance(int a, int& aa, int n, bool inside,
                   const Point& v, std::vector<Point>& out) {
    if (inside) out.push_back(v);
    aa++;
    return (a + 1) % n;
}

// O'Rourke's InOut: state transition for the intersection boundary.
// Reference: lines of InOut() in O'Rourke's C code.
//   aHB > 0  ->  Pin
//   bHA > 0  ->  Qin
//   otherwise keep current inflag (status quo)
static InFlag in_out(const Point& /*p*/, InFlag inflag, int aHB, int bHA) {
    if      (aHB > 0) return InFlag::Pin;
    else if (bHA > 0) return InFlag::Qin;
    return inflag;
}

// Point-in-convex-polygon test, used only for the "boundaries don't cross"
// closure (one polygon entirely inside the other, or single touching point).
// Epick's `orientation` is an exact predicate (it's the construction that
// is inexact), so a vertex on a shared boundary edge is classified exactly
// the same as in Epeck.  Returning the wrong polygon here would silently
// corrupt the Fréchet bound, but predicates are still bit-exact.
static bool point_in_convex_poly(const Point& p, const std::vector<Point>& poly) {
    int n = (int)poly.size();
    for (int i = 0; i < n; ++i) {
        const Point& a = poly[i];
        const Point& b = poly[(i + 1) % n];
        if (CGAL::orientation(a, b, p) == CGAL::RIGHT_TURN) return false;
    }
    return true;
}

// Verify that every vertex of `inner` is inside `outer`.  Used to validate the
// containment shortcut: a single false positive (vertex on a shared boundary
// edge misclassified) would return the wrong polygon.
static bool all_points_in_convex_poly(const std::vector<Point>& inner,
                                      const std::vector<Point>& outer) {
    for (const auto& v : inner) {
        if (!point_in_convex_poly(v, outer)) return false;
    }
    return true;
}

// Line-line intersection for the half-plane clipper.  Given a directed edge
// (a, b) defining the "inside" half-plane (a point q is inside iff
// orientation(a, b, q) != RIGHT_TURN) and a directed input segment (p1, p2),
// return the intersection point of the two lines, or std::nullopt if the
// segments don't intersect (parallel/collinear case).
static std::optional<Point> line_intersect(const Point& a, const Point& b,
                                           const Point& p1, const Point& p2) {
    // Guard against degenerate segments.  clip_halfplane's orientation checks
    // filter out the common case, but a polygon vertex that was itself
    // produced by a prior intersection (and rounds to the same Epick coords
    // as its neighbour) can produce a zero-length segment here.  The old
    // floating-point version handled this via len2 > 0; CGAL::intersection
    // with Epeck segs crashes instead.
    if (p1.x() == p2.x() && p1.y() == p2.y()) return p1;
    if (a.x() == b.x() && a.y() == b.y()) {
        // Clip edge is degenerate; keep the segment endpoint inside the half-plane.
        return (CGAL::orientation(a, b, p1) != CGAL::RIGHT_TURN) ? std::optional<Point>(p1) : std::nullopt;
    }

    // Use Epeck for exact intersection construction.  Floating-point arithmetic
    // (Epick) can place the intersection slightly off the segment, and when the
    // next clipping edge uses that point as input, the error compounds — the
    // point no longer satisfies orientation(a, b, intersection) == COLLINEAR and
    // the polygon boundary shifts.  Epeck keeps every vertex exactly on the
    // intersection of the two lines, which is what Sutherland-Hodgman requires.
    Epeck::Point_2 ae = conv_to_exact(a);
    Epeck::Point_2 be = conv_to_exact(b);
    Epeck::Point_2 p1e = conv_to_exact(p1);
    Epeck::Point_2 p2e = conv_to_exact(p2);

    // Build segments inside the guard to catch any CGAL assertion/abort.
    Epeck::Segment_2 seg_ab(ae, be);
    Epeck::Segment_2 seg_p(p1e, p2e);
    auto inter = CGAL::intersection(seg_ab, seg_p);
    if (!inter) {
        // No intersection — parallel or collinear non-overlapping segments.
        return std::nullopt;
    }
    if (const Epeck::Point_2* ip = std::get_if<Epeck::Point_2>(&*inter)) {
        return conv_to_inexact(*ip);
    }
    // Collinear overlapping segments — treat as no intersection.
    return std::nullopt;
}

// Sutherland-Hodgman half-plane clip.  `subject` is a CCW convex polygon;
// the half-plane is defined by directed edge (a, b) with the inside being
// the LEFT side (i.e. the half-plane orientation(a, b, p) != RIGHT_TURN).
// Returns the clipped polygon, also CCW.
//
// Used by convex_intersect_fast: for each edge of the `clip` polygon we
// clip the running subject against that half-plane.  After |clip| edges
// the running subject is exactly `subject ∩ clip`.
static std::vector<Point> clip_halfplane(const std::vector<Point>& subject,
                                          const Point& a, const Point& b) {
    int n = (int)subject.size();
    if (n == 0) return {};

    // Convert clip edge to Epeck for exact orientation tests.  Epick (floating-
    // point) orientation can misclassify near-parallel edges due to rounding,
    // causing clip_halfplane to call line_intersect on segments that don't
    // actually cross the half-plane boundary.  line_intersect then returns null,
    // and we append a bogus vertex that corrupts the polygon geometry.
    Epeck::Point_2 ae = conv_to_exact(a);
    Epeck::Point_2 be = conv_to_exact(b);

    std::vector<Point> out;
    out.reserve(n + 1);

    for (int i = 0; i < n; ++i) {
        const Point& cur  = subject[i];
        const Point& prev = subject[(i + n - 1) % n];
        // Exact orientation test: only COLLINEAR is ambiguous.
        auto o_cur  = CGAL::orientation(ae, be, conv_to_exact(cur));
        auto o_prev = CGAL::orientation(ae, be, conv_to_exact(prev));
        bool cur_in  = (o_cur  != CGAL::RIGHT_TURN);
        bool prev_in = (o_prev != CGAL::RIGHT_TURN);

        if (cur_in) {
            if (!prev_in) {
                // crossing INTO the half-plane: emit intersection if it exists
                if (auto ip = line_intersect(a, b, prev, cur)) out.push_back(*ip);
            }
            out.push_back(cur);
        } else if (prev_in) {
            // crossing OUT of the half-plane: emit intersection if it exists
            if (auto ip = line_intersect(a, b, prev, cur)) out.push_back(*ip);
        }
    }
    return out;
}

// Specialised convex-convex intersection using Sutherland-Hodgman clipping.
// Assumes both inputs are CCW convex polygons.
//
// Faster than O'Rourke for small inputs (typical here: |P|~12, |Q|~19).
// Total work is O(|P| * |Q|) where each unit is one `orientation` test.
// O'Rourke's ConvexIntersect runs in O(|P| + |Q|) — each polygon edge is
// examined at most twice, and each iteration does O(1) work — so for the
// small polygons typical in this application the clipper's simpler per-step
// cost (one orientation) beats O'Rourke's (four predicates + bookkeeping).
//
// Returns the intersection polygon in CCW order, or {} if empty.
static std::vector<Point> convex_intersect_fast(const std::vector<Point>& Pin,
                                                 const std::vector<Point>& Qin) {
    if (Pin.size() < 3 || Qin.size() < 3) return {};

    // Normalise both to CCW using exact (Epeck) signed area.  The previous
    // Epick (floating-point) kernel_signed_area could misclassify nearly-zero-
    // area polygons, reversing one of them and corrupting the intersection.
    auto signed_area_sign = [](const std::vector<Point>& P) {
        Epeck::FT s = 0;
        int sz = (int)P.size();
        for (int i = 0; i < sz; ++i) {
            Epeck::Point_2 a = conv_to_exact(P[i]);
            Epeck::Point_2 b = conv_to_exact(P[(i + 1) % sz]);
            s += a.x() * b.y() - b.x() * a.y();
        }
        return CGAL::sign(s);
    };

    std::vector<Point> P, Q;
    P.reserve(Pin.size());
    Q.reserve(Qin.size());
    if (signed_area_sign(Pin) < 0) { P.assign(Pin.rbegin(), Pin.rend()); }
    else                            { P = Pin; }
    if (signed_area_sign(Qin) < 0) { Q.assign(Qin.rbegin(), Qin.rend()); }
    else                            { Q = Qin; }

    // Containment shortcuts: if one polygon is entirely inside the other,
    // return the inner one.  We test P[0] inside Q (and vice versa), then
    // validate that ALL vertices of the inner polygon are inside the outer one
    // — a single false positive (e.g. a shared boundary vertex misclassified)
    // would silently return the wrong polygon and corrupt the Fréchet bound.
    if (point_in_convex_poly(P[0], Q) && all_points_in_convex_poly(P, Q)) return P;
    if (point_in_convex_poly(Q[0], P) && all_points_in_convex_poly(Q, P)) return Q;

    // Otherwise clip P against each edge of Q.
    std::vector<Point> running = P;
    int m = (int)Q.size();
    for (int i = 0; i < m; ++i) {
        const Point& a = Q[i];
        const Point& b = Q[(i + 1) % m];
        running = clip_halfplane(running, a, b);
        // The intersection is empty if clipping against a single edge produces
        // fewer than 3 vertices.  Continuing would accumulate degenerate vertices
        // from line_intersect null returns, inflating the output polygon.
        if (running.size() < 3) return {};
    }
    return running;
}

// O'Rourke's ConvexIntersect, faithful transcription of the reference in
// "Computational Geometry in C" (2nd ed.), Chapter 7.  The reference
// requires:
//   * Both polygons are CCW
//   * The first vertex of each polygon is the bottom-most (leftmost on tie)
// Returns the intersection polygon as a CCW point list.  If P is inside Q
// (or vice versa) and they don't cross, the inner polygon is returned.
//
// Implementation kernel: Epick for fast common-case intersections, verified
// with exact orientation and falling back to Epeck only when the Epick
// point is off the segment.  Predicates (`orientation`, `compare_*`) remain
// bit-exact in Epick.  The fallback is rare (~< 1% of calls) but critical
// for correctness on near-degenerate inputs.
static std::vector<Point> convex_intersect(const std::vector<Point>& Pin,
                                           const std::vector<Point>& Qin) {
    int n = (int)Pin.size();
    int m = (int)Qin.size();
    if (n < 3 || m < 3) return {};

    // Normalize both to CCW and rotate to start at the bottommost (leftmost
    // on tie) vertex.  The reference requires both; we do it explicitly.
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
    std::vector<Point> Pr = Pin, Qr = Qin;
    if (signed_area(Pr) < 0) std::reverse(Pr.begin(), Pr.end());
    if (signed_area(Qr) < 0) std::reverse(Qr.begin(), Qr.end());
    int a = leftmost_index(Pr);
    int b = leftmost_index(Qr);

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

        // Intersection: fast path in Epick, verified with exact orientation.
        // Epick's orientation predicate is exact even when the construction
        // is slightly off, so the verification below catches all cases where
        // the double-intersection is too far from the true intersection.
        Point p, q;
        char code;
        {
            Epick::Segment_2 seg_a(Pr[a1], Pr[a]);
            Epick::Segment_2 seg_b(Qr[b1], Qr[b]);
            auto inter = CGAL::intersection(seg_a, seg_b);
            if (!inter) {
                code = '0';
            } else if (const Epick::Segment_2* seg =
                           std::get_if<Epick::Segment_2>(&*inter)) {
                p = seg->source();
                q = seg->target();
                code = 'e';
            } else if (const Epick::Point_2* ip =
                           std::get_if<Epick::Point_2>(&*inter)) {
                p = *ip;
                // Verify p lies on both segments using exact orientation.
                // orientation == COLLINEAR means p is exactly on the segment
                // in exact arithmetic — a slight double error would show as
                // LEFT_TURN or RIGHT_TURN, triggering the Epeck fallback.
                bool on_seg_a = (CGAL::orientation(Pr[a1], Pr[a], p) == CGAL::COLLINEAR);
                bool on_seg_b = (CGAL::orientation(Qr[b1], Qr[b], p) == CGAL::COLLINEAR);
                if (!on_seg_a || !on_seg_b) {
                    // Epick intersection is off; redo in Epeck for correctness.
                    Epeck::Segment_2 eseg_a(conv_to_exact(Pr[a1]), conv_to_exact(Pr[a]));
                    Epeck::Segment_2 eseg_b(conv_to_exact(Qr[b1]), conv_to_exact(Qr[b]));
                    auto einter = CGAL::intersection(eseg_a, eseg_b);
                    if (einter) {
                        if (const Epeck::Point_2* ep =
                                std::get_if<Epeck::Point_2>(&*einter))
                            p = conv_to_inexact(*ep);
                        else if (const Epeck::Segment_2* es =
                                     std::get_if<Epeck::Segment_2>(&*einter)) {
                            p = conv_to_inexact(es->source());
                            q = conv_to_inexact(es->target());
                        }
                    }
                }
                // 'v' = an endpoint of one segment lies on the other.
                if (p == Qr[b1] || p == Qr[b] || p == Pr[a1] || p == Pr[a])
                    code = 'v';
                else
                    code = '1';
            } else {
                code = '0';
            }
        }

        // If A & B intersect, update inflag and emit p.
        if (code == '1' || code == 'v') {
            if (inflag == InFlag::Unknown && first_point) {
                aa = ba = 0;
                first_point = false;
                p0 = p;
            }
            out.push_back(p);
            inflag = in_out(p, inflag, aHB, bHA);
        }

        // Special case: A & B overlap and oppositely oriented -> shared edge.
        if ((code == 'e') && (dot_sign(Pr[a1], Pr[a], Qr[b1], Qr[b]) < 0)) {
            out.push_back(p);
            out.push_back(q);
            return out;
        }

        // Special case: A & B parallel and separated -> disjoint.
        if ((cross == 0) && (aHB < 0) && (bHA < 0)) {
            return {};
        }

        // Special case: A & B collinear.
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
        } else {  // cross < 0
            if (aHB > 0)
                b = advance(b, ba, m, inflag == InFlag::Qin, Qr[b], out);
            else
                a = advance(a, aa, n, inflag == InFlag::Pin, Pr[a], out);
        }
    } while (((aa < n) || (ba < m)) && (aa < 2 * n) && (ba < 2 * m));

    if (first_point) {
        // No intersection found: one polygon inside the other, or single touch.
        if (point_in_convex_poly(Pr[0], Qr)) return Pr;
        if (point_in_convex_poly(Qr[0], Pr)) return Qr;
        return {};
    }

    // Close the polygon.
    out.push_back(p0);
    return out;
}

}  // namespace orourke_cgal

static void intersect(const Polygon& subject, const Polygon& clip,
                        std::back_insert_iterator<std::vector<Polygon_with_holes>> result) {
    std::vector<Point> P_verts, Q_verts;
    {
        TIMER("copy_verts");
        P_verts.reserve(subject.size());
        Q_verts.reserve(clip.size());
        for (auto it = subject.vertices_begin(); it != subject.vertices_end(); ++it) {
            P_verts.push_back(*it);
        }
        for (auto it = clip.vertices_begin(); it != clip.vertices_end(); ++it) {
            Q_verts.push_back(*it);
        }
    }

    std::vector<Point> verts;
    {
        TIMER("convex_intersect");
        verts = orourke_cgal::convex_intersect(P_verts, Q_verts);
    }
    // Drop the duplicated closing vertex (O'Rourke appends p0 at the end to
    // close the ring; CGAL::Polygon stores an open ring, so leaving the
    // duplicate in produces a polygon with a spurious zero-length edge that
    // later breaks find_F's tangent search).
    if (verts.size() < 4) return;            // < 3 real vertices -> nothing
    if (verts.front() == verts.back())
        verts.pop_back();

    // The intersection should be convex, hence simple.  But the O'Rourke
    // main loop can emit kinks / collinear vertices that confuse CGAL's
    // is_simple_2 check; ensure the output is well-formed before constructing
    // a Polygon (whose constructor asserts is_simple_2).
    Polygon inter_poly(verts.begin(), verts.end());
    if (!inter_poly.is_simple()) return;
    if (inter_poly.is_clockwise_oriented()) inter_poly.reverse_orientation();
    *result++ = Polygon_with_holes(inter_poly);
}

bool showF = false; // true if -F is passed
bool showG = false; // true if -G is passed
bool showS = false; // true if -S is passed
bool out_flag = false;      // write outputs if true
bool gui_flag = false;      // show viewer if true
bool dist_flag = false;     // compute frechet if true

constexpr double TOL = 1e-6;
constexpr double SQRT2 =
    1.41421356237; // sqrt in STL does not have constexpr version !?

// Bounding square [BMIN, BMAX]^2
extern const double BMIN = -10000;
extern const double BMAX = 10000;

// TODO: DELTA should not be too high that convex hull goes out of bounding square, which may cause the program to crash
double DELTA = 200;
double EPSILON = 0.6;

double GRID_val() { return EPSILON * DELTA / (2 * SQRT2); }

double R_val() { return (1.0 + EPSILON / 2.0) * DELTA; }

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

// Iterate every grid row y in [y_min, y_max] for which the disk
// of radius r = (1 + EPSILON/2) * DELTA has non-zero width.
// Calls f(px, py, y_actual, x_min, x_max) for each surviving row.
template <class F>
static void for_each_grid_row(double px, double py, F &&f) {
    const double r = R_val();
    const double GRID = GRID_val();
    const double r2 = r * r;
    const int y_min = -(r / GRID) - 1;
    const int y_max = -y_min;
    for (int y = y_min; y <= y_max; y++) {
        const double y_actual = y * GRID;
        const double dy2 = y_actual * y_actual;
        if (dy2 > r2) continue;  // row entirely outside the disk
        const double dx = std::sqrt(r2 - dy2);
        const int x_min = static_cast<int>(-dx / GRID - 1);
        const int x_max = static_cast<int>(+dx / GRID + 1);
        f(px, py, y_actual, x_min, x_max);
    }
}

std::vector<Point> get_conv_from_grid(const Point &p) {
    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    const double GRID = GRID_val();
    std::vector<Point> boundaries;
    for_each_grid_row(px, py, [&](double px_, double py_,
                                  double y_actual, int x_min, int x_max) {
        boundaries.emplace_back(px_ + x_min * GRID, py_ + y_actual);
        if (x_max != 0)  // avoid duplicating the rightmost of the previous row
            boundaries.emplace_back(px_ + x_max * GRID, py_ + y_actual);
    });
    std::vector<Point> conv;
    // O(n log n) hull; here h = O(n) so the asymptotic gain is modest.
    CGAL::convex_hull_2(boundaries.begin(), boundaries.end(),
                        std::back_inserter(conv));
    return conv;
}

double expected_frechet = 0;

std::vector<Point> get_points_from_grid(const Point &p) {
    const double px = CGAL::to_double(p.x());
    const double py = CGAL::to_double(p.y());
    const double GRID = GRID_val();
    if (DELTA == 0) {
        return std::vector<Point>{p};
    }
    std::vector<Point> points;
    for_each_grid_row(px, py, [&](double px_, double py_,
                                  double y_actual, int x_min, int x_max) {
        for (int x = x_min; x <= x_max; x++) {
            points.emplace_back(px_ + x * GRID, py_ + y_actual);
            expected_frechet = std::max(expected_frechet, CGAL::to_double(CGAL::squared_distance(p, points.back())));
        }
    });
    return points;
}

// Find the two tangent indices of p with respect to convex polygon S.
// A vertex S[i] is a "tangent" iff orientation(p, S[i-1], S[i]) !=
// orientation(p, S[i], S[i+1]) and the trailing edge is non-collinear.
//
// We walk the polygon once.  The original O(n) code recomputed both
// orientations every iteration, even though orientation(p, S[i], S[i+1])
// at iteration i is the same value as orientation(p, S[i-1], S[i]) at
// iteration i+1.  We cache that value across iterations, halving the
// number of exact predicates per call.  We also break out of the loop
// once both tangents are found — for a convex polygon there are at most
// two, so we never need to scan all the way around.
std::vector<int> find_tangent_idx(const Point &p,
                                 const std::vector<Point> &S) {
    int n = S.size();
    std::vector<int> tangent;
    {
        TIMER("find_tangent_idx");
        tangent.reserve(2);

        auto orient_at = [&](int i) {
            return CGAL::orientation(p, S[i], S[(i + 1) % n]);
        };

        int prev = orient_at(n - 1);  // orientation(p, S[n-1], S[0])
        for (int i = 0; i < n; ++i) {
            int cur = orient_at(i);
            if (cur != prev && cur != CGAL::COLLINEAR) {
                tangent.push_back(i);
                if (tangent.size() == 2) break;
            }
            prev = cur;
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
            if (CGAL::to_double(CGAL::squared_distance(*ip, p)) < TOL*TOL) return std::nullopt;
            return *ip;
        }
        if (const Segment* seg = std::get_if<Segment>(&*obj)) {
            const Point& a = seg->source();
            const Point& b = seg->target();
            double da = CGAL::to_double(CGAL::squared_distance(a, p));
            double db = CGAL::to_double(CGAL::squared_distance(b, p));
            if (da < TOL*TOL) return b;
            if (db < TOL*TOL) return a;
            return (da < db) ? a : b;
        }
    }
    return std::nullopt;
}

std::vector<Point> find_F(const Point& p, const std::vector<Point>& S) {
    assert(S.size() != 2); // wait why this check??
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

    std::optional<Point> hit1, hit2;
    {
        TIMER("ray_hits");
        hit1 = intersect_ray_with_rect(p, S[tangent[0]]);
        hit2 = intersect_ray_with_rect(p, S[tangent[1]]);
    }
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
    {
        TIMER("build_F_poly");
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
    TIMER("get_longest_stab");
    const Point& p0 = stream[cur];
    int p0cur = cur;
    std::vector<Point> P;
    {
        TIMER("get_points_from_grid");
        P = get_points_from_grid(p0);
        // std::cerr << "P.size() = " << P.size() << '\n';
    }
    if (viewer) {
        TIMER("gui_update");
        viewer->markP0(p0);
        viewer->addOriginalPoint(p0);
    }
    std::array<Point, 2> buffer = {p0, p0};
    std::vector<std::vector<Point>> S(P.size(), std::vector<Point>{p0});
    int dead_cnt = 0;
    std::vector<int> dead(P.size());
    cur++;
    while (cur < int(stream.size())) {
        const Point& pi  = stream[cur];
        std::vector<Point> Gi;
        {
            TIMER("get_conv_from_grid");
            Gi = get_conv_from_grid(pi);
        }
        std::vector<std::vector<Polygon_with_holes>> new_S(P.size());
        for (int i = 0; i < int(P.size()); i++) {
            if (dead[i]) {
                continue;
            }
            std::vector<Point> F;
            Polygon F_poly, Gi_poly, S_poly;
            {
                TIMER("find_F");
                F = find_F(P[i], S[i]);
            }
            {
                TIMER("polygon_construction");
                F_poly = Polygon(F.begin(), F.end());
                Gi_poly = Polygon(Gi.begin(), Gi.end());
                S_poly = Polygon(S[i].begin(), S[i].end());
            }

            if (i == 0) {
                const QColor stepColors[] = {
                    Qt::red,
                    Qt::blue,
                    Qt::green,
                    Qt::magenta,
                    Qt::cyan
                };
                QColor c = stepColors[cur % 6];
                if (viewer) {
                    TIMER("gui_update");
                    if (showF) viewer->addPolygon(F_poly, c);
                    if (showG) viewer->addPolygon(Gi_poly, c);
                    if (showS) viewer->addPolygon(S_poly, c);
                }
            }

            {
                TIMER("intersect");
                intersect(F_poly, Gi_poly, back_inserter(new_S[i]));
            }

            if (new_S[i].size() == 0) {
                dead[i] = true;
                dead_cnt++;
                continue;
            }
            assert(new_S[i].size() == 1);
            buffer[0] = P[i];
            buffer[1] = *new_S[i].begin()->outer_boundary().vertices_begin();
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
                S[i].clear();
                std::copy(new_S[i].begin()->outer_boundary().vertices_begin(),
                          new_S[i].begin()->outer_boundary().vertices_end(),
                          std::back_inserter(S[i]));
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

std::vector<Point> simplify(const std::vector<Point> &stream, MultiViewer* viewer = nullptr) {
    TIMER("simplify");
    std::vector<Point> simplified;
    std::cout << "Simplifying...\n";
    int cur = 0;
    while (cur != int(stream.size())) {
        // std::cerr << "New segment starts at G_" << cur << '\n';
        cur = get_longest_stab(stream, cur, simplified, viewer);
    }
    std::cout << "Expected Frechet distance: " << std::sqrt(expected_frechet) << '\n';
    std::cout << "The original stream of size " << stream.size() << " is simplified to " << simplified.size() << " points.\n";
    return simplified;
}

// TODO: clean up this mess of flags
int main(int argc, char** argv) {
    int test_case_no = -1;      // provided by --in <id>

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

    // Shorthand: first positional numeric argument means '--in <id> --out'; allows extra flags (e.g., '--gui')
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
    // Determine repo_root once (used for input/output/frechet). Try executable directory upward.
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
    // search from the executable
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


    // Simplify
    std::vector<Point> simplified;
    {
        TIMER("total");
        simplified = simplify(stream, vptr);
        stream = std::move(simplified);
    }
    print_timing_summary();

    if (out_flag) {
        std::filesystem::path dir = repo_root / "data" / std::to_string(test_case_no);
        std::filesystem::create_directories(dir);
        
        std::ofstream simp(dir / "simplify.txt");
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