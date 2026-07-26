#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "simplify_geometry.h"
#include "timer.h"

// ===========================================================================
//  Global parameters
// ===========================================================================

double DELTA = 200;
double EPSILON = 0.5;
double expected_frechet_squared = 0.0;

bool out_flag = false;
bool dist_flag = false;

// ===========================================================================
//  Trajectory-simplification geometry
// ===========================================================================

double BMIN = -10000;
double BMAX = 10000;
static constexpr double TOL = 1e-6;

inline bool point_in_convex(const Point& p, const std::vector<Point>& poly, bool ccw = true) {
    // Filtered orientation: a double determinant guarded by its a-priori
    // rounding bound resolves the clear majority; only genuinely ambiguous
    // (near-boundary) corners defer to the exact predicate.  Bit-identical to
    // CGAL::orientation but skips GMP on the near-collinear stab vertices that
    // otherwise dominate this test.
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
// polygon S.  The original formulation tested, for every vertex i, whether ALL
// other vertices lie weakly on one side of line (p, S[i]) — an O(n^2) sweep of
// exact orientation predicates, which the profiler showed dominating runtime
// (the near-collinear ray/hull configurations force CGAL's filter to fall
// through to exact GMP arithmetic).
//
// For a CONVEX polygon this is equivalent to a purely local O(n) test: the line
// through p and S[i] is a supporting line iff the two edges incident to S[i]
// leave S[i] toward the SAME closed side of that line, i.e. the neighbours
// S[i-1] and S[i+1] are not on strictly opposite sides.  Convexity guarantees
// that if both incident edges lie in a closed half-plane bounded by a line
// through their common vertex, the entire polygon does too — so the local test
// and the global sweep classify every vertex identically.  S is always convex
// here (it is produced by convex-polygon intersection), and find_F has already
// ensured p is not interior, so the tangents are well defined.  This is an
// exact equivalence, not an approximation: it uses the same exact orientation
// predicate, only O(n) times instead of O(n^2).
inline std::vector<int> find_tangent_idx(const Point& p, const std::vector<Point>& S) {
    std::vector<int> tangent;
    tangent.reserve(2);
    const int n = static_cast<int>(S.size());

    // Cache p and S coordinates as doubles once (to_double on an Epick point is
    // cheap but S is scanned 2n times, so hoisting the conversions out of the
    // inner loop is a measurable win).
    static std::vector<double> sx, sy;
    sx.resize(n);
    sy.resize(n);
    for (int i = 0; i < n; ++i) {
        sx[i] = CGAL::to_double(S[i].x());
        sy[i] = CGAL::to_double(S[i].y());
    }
    const double pdx = CGAL::to_double(p.x());
    const double pdy = CGAL::to_double(p.y());

    // Sign of orientation(p, S[i], S[j]) in pure doubles.  (+1 left, -1 right,
    // 0 collinear.)  No exact fallback: for tangent finding an exact result on
    // (near-)collinear ties is unnecessary.  The two tangents are the vertices
    // of extreme bearing from p; if two candidate vertices are collinear as
    // seen from p, either is an equally valid supporting vertex because the
    // wedge F's boundary ray passes through that collinear point regardless, so
    // F is the same point set either way.  Dropping the CGAL::orientation
    // fallback here removes essentially all remaining GMP arithmetic (the stab
    // chain S is riddled with near-collinear vertices that made the filtered
    // predicate defer to exact constantly).
    auto orient = [&](int i, int j) -> int {
        const double ux = sx[i] - pdx, uy = sy[i] - pdy;
        const double wx = sx[j] - pdx, wy = sy[j] - pdy;
        const double det = ux * wy - uy * wx;
        if (det > 0.0) return  1;
        if (det < 0.0) return -1;
        return 0;
    };

    // The two supporting (tangent) vertices from an external point p to a
    // convex polygon are the vertices of EXTREME bearing as seen from p.  Since
    // p lies outside the convex hull, all of S falls within an angular wedge of
    // less than 180 degrees, inside which bearing is a total order; so a single
    // linear scan keeping the most-clockwise vertex (every other vertex to its
    // left) finds one tangent, and keeping the most-counterclockwise vertex
    // finds the other.  This is 2n exact-orientation comparisons instead of the
    // O(n^2) all-pairs survey, and is exact: orient() returns the same sign as
    // CGAL::orientation.  On collinear extreme vertices it deterministically
    // keeps one valid tangent (the survey instead found both and bailed out),
    // which is strictly more robust.
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

// ===========================================================================
//  Help & debug printing
// ===========================================================================

static void print_help() {
    std::cout << "Usage: simplify [options]\n"
              << "  --in <id>        Read input from data/taxi/<id>.txt (resolved absolutely)\n"
              << "  --out            Write output to data/<id>/original.txt & simplify.txt (resolved absolutely; requires --in <id>)\n"
              << "  --dist           After output, compute Frechet distance by invoking ./frechet (Julia wrapper) with --in <id> --path <simplify.txt>\n"
              << "  -d <delta>       Override DELTA (default " << DELTA << ")\n"
              << "  -e <epsilon>     Override EPSILON (default " << EPSILON << ")\n"
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
//  Exact intersection (ground truth) — EPECK CGAL::intersection
// ===========================================================================

// Computes the intersection of two convex polygons using EPECK exact
// constructions. Mirrors the semantics of the fast intersect() in
// simplify_geometry.h: returns true and fills `result` (CCW, no duplicate
// closing vertex, >= 3 vertices) when the intersection has positive area;
// otherwise clears `result` and returns false.
inline bool intersect_exact(const std::vector<Point>& P_verts,
                            const std::vector<Point>& Q_verts,
                            std::vector<Point>& result) {
    using EPolygon = CGAL::Polygon_2<Epeck>;
    using EPwh     = CGAL::Polygon_with_holes_2<Epeck>;

    result.clear();
    if (P_verts.size() < 3 || Q_verts.size() < 3) return false;

    EPolygon P, Q;
    for (const auto& p : P_verts) P.push_back(conv_to_exact(p));
    for (const auto& q : Q_verts) Q.push_back(conv_to_exact(q));
    if (P.orientation() == CGAL::CLOCKWISE) P.reverse_orientation();
    if (Q.orientation() == CGAL::CLOCKWISE) Q.reverse_orientation();

    std::vector<EPwh> out;
    CGAL::intersection(P, Q, std::back_inserter(out));
    if (out.empty()) return false;

    const EPolygon& ob = out.front().outer_boundary();
    if (ob.size() < 3) return false;

    result.reserve(ob.size());
    for (auto it = ob.vertices_begin(); it != ob.vertices_end(); ++it) {
        result.emplace_back(CGAL::to_double(it->x()), CGAL::to_double(it->y()));
    }
    return true;
}

// ===========================================================================
//  Core algorithm
// ===========================================================================

int get_longest_stab(const std::vector<Point>& stream, int cur,
                     std::vector<Point>& simplified,
                     double EPSILON, double DELTA) {
    SIGNPOST_EVENT(stab.start, "cur=%d", cur);
    const Point& p0 = stream[cur];
    std::vector<Point> P;
    {
        SIGNPOST_BEGIN(get_points_from_grid);
        P = get_points_from_grid(p0, EPSILON, DELTA);
        SIGNPOST_END(get_points_from_grid);
        SIGNPOST_EVENT(stab.grid_points, "Pn=%lu", (unsigned long)P.size());
    }
    std::array<Point, 2> buffer = {p0, p0};
    const int Pn = (int)P.size();
    std::vector<std::vector<Point>> S(Pn);
    for (int i = 0; i < Pn; ++i) S[i] = {P[i]};
    int dead_cnt = 0;
    std::vector<int> dead(Pn);
    std::vector<std::vector<Point>> new_S(Pn);
    std::vector<std::vector<Point>> F(Pn);
    std::vector<Point> Gi;

    cur++;
    while (cur < int(stream.size())) {
        const Point& pi = stream[cur];
        SIGNPOST_BEGIN(stab.step);
        SIGNPOST_EVENT(stab.consume, "cur=%d", cur);
        Gi = get_conv_from_grid(pi, EPSILON, DELTA);

        for (int i = 0; i < Pn; ++i) {
            if (dead[i]) continue;
            {
                SIGNPOST_BEGIN(find_F);
                find_F(P[i], S[i], F[i]);
                SIGNPOST_END(find_F);
            }
            bool hit;
            {
                SIGNPOST_BEGIN(intersect);
                hit = intersect(F[i], Gi, new_S[i]);
                SIGNPOST_END(intersect);
            }
#ifdef DUMP_DIVERGENCE
            {
                // Does RAW O'Rourke output (pre-guard) actually differ from the
                // exact result, or only the guard verdict?  Compare the two
                // polygons cyclically within a coordinate tolerance.
                std::vector<Point> Fd = orourke_cgal::dedup_consecutive(F[i]);
                std::vector<Point> Gd = orourke_cgal::dedup_consecutive(Gi);
                std::vector<Point> raw = orourke_cgal::convex_intersect_robust(Fd, Gd);
                if (!raw.empty() && raw.front() == raw.back()) raw.pop_back();
                std::vector<Point> ex;
                bool ex_hit = intersect_exact_convex(F[i], Gi, ex);

                // Length-independent shape comparison: symmetric Hausdorff
                // distance between the RAW and EXACT vertex sets.  Immune to
                // degenerate edges, vertex count, and cyclic ordering, so it
                // measures GENUINE geometric divergence only.
                auto directed_hausdorff = [](const std::vector<Point>& A,
                                             const std::vector<Point>& B) -> double {
                    double worst = 0.0;
                    for (const auto& a : A) {
                        double best = 1e300;
                        for (const auto& b : B) {
                            double dx = CGAL::to_double(a.x()) - CGAL::to_double(b.x());
                            double dy = CGAL::to_double(a.y()) - CGAL::to_double(b.y());
                            best = std::min(best, dx*dx + dy*dy);
                        }
                        worst = std::max(worst, best);
                    }
                    return std::sqrt(worst);
                };

                static long n_match = 0;      // RAW ~= EXACT geometrically
                static long n_genuine = 0;     // RAW genuinely differs
                static long n_both_empty = 0;
                static double max_hd = 0.0;
                struct Reporter {
                    long *m, *g, *be; double* mh;
                    ~Reporter() {
                        std::ofstream f("hausdorff_tally.txt");
                        f << "match=" << *m << " genuine_diff=" << *g
                          << " both_empty=" << *be << " max_hausdorff=" << *mh << "\n";
                    }
                };
                static Reporter rep{&n_match, &n_genuine, &n_both_empty, &max_hd};

                bool raw_ok = raw.size() >= 3;
                if (!raw_ok && !ex_hit) { n_both_empty++; }
                else if (raw_ok && ex_hit) {
                    double hd = std::max(directed_hausdorff(raw, ex),
                                         directed_hausdorff(ex, raw));
                    max_hd = std::max(max_hd, hd);
                    if (hd < 1e-3) n_match++;
                    else {
                        n_genuine++;
                        static int dumped = 0;
                        if (dumped < 5) {
                            dumped++;
                            std::ofstream d("divergence_case_" + std::to_string(dumped) + ".txt");
                            d << std::setprecision(std::numeric_limits<double>::max_digits10);
                            d << "# hausdorff=" << hd << " raw=" << raw.size()
                              << " exact=" << ex.size() << "\n";
                            d << "F " << F[i].size() << "\n";
                            for (auto& p : F[i]) d << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
                            d << "Gi " << Gi.size() << "\n";
                            for (auto& p : Gi) d << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
                            d << "RAW " << raw.size() << "\n";
                            for (auto& p : raw) d << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
                            d << "EXACT " << ex.size() << "\n";
                            for (auto& p : ex) d << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
                            d.close();
                        }
                    }
                } else {
                    // one empty, one not — a genuine hit/miss disagreement
                    n_genuine++;
                }
            }
#endif
            if (!hit) {
                dead[i] = true;
                dead_cnt++;
            }
        }

        bool has_candidate = false;
        for (int i = Pn - 1; i >= 0 && !has_candidate; --i) {
            if (dead[i] || new_S[i].empty()) continue;
            buffer[0] = P[i];
            buffer[1] = new_S[i].front();
            has_candidate = true;
        }
        SIGNPOST_END(stab.step);
        if (!has_candidate || dead_cnt == Pn) break;

        for (int i = 0; i < Pn; ++i) {
            if (!dead[i]) S[i].swap(new_S[i]);
        }
        cur++;
    }
    SIGNPOST_EVENT(stab.emit, "buffer=(%{public}lf,%{public}lf)->(%{public}lf,%{public}lf)",
                   (double)CGAL::to_double(buffer[0].x()), (double)CGAL::to_double(buffer[0].y()),
                   (double)CGAL::to_double(buffer[1].x()), (double)CGAL::to_double(buffer[1].y()));

    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    return cur;
}

std::vector<Point> simplify(const std::vector<Point>& stream,
                            double EPSILON, double DELTA) {
    SIGNPOST_BEGIN(simplify.run);
    configure_bbox(stream, EPSILON, DELTA);
    std::vector<Point> simplified;
    std::cout << "Simplifying...\n";
    int cur = 0;
    int prefix = 0;
    while (cur != int(stream.size())) {
        SIGNPOST_EVENT(simplify.prefix, "from=%d size=%lu", cur, (unsigned long)simplified.size());
        cur = get_longest_stab(stream, cur, simplified, EPSILON, DELTA);
        prefix++;
    }
    SIGNPOST_END(simplify.run);
    std::cout << "Expected Frechet distance: " << std::sqrt(expected_frechet_squared) << '\n';
    std::cout << "The original stream of size " << stream.size() << " is simplified to " << simplified.size() << " points.\n";
    SIGNPOST_EVENT(simplify.done, "n=%lu m=%lu prefixes=%d",
                   (unsigned long)stream.size(), (unsigned long)simplified.size(), prefix);
    return simplified;
}

// ===========================================================================
//  main
// ===========================================================================

// TODO: clean up this mess of flags
int main(int argc, char** argv) {
    int test_case_no = -1;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i],"--out") == 0) out_flag = true;
        else if (strcmp(argv[i],"--dist") == 0) dist_flag = true;
        else if (strcmp(argv[i],"--gui") == 0 || strcmp(argv[i],"-F") == 0 ||
                 strcmp(argv[i],"-G") == 0 || strcmp(argv[i],"-S") == 0) {
            std::cerr << "GUI options require simplify_with_gui\n";
            return 1;
        }
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

    SIGNPOST_BEGIN(app.run);
    std::vector<Point> simplified = simplify(stream, EPSILON, DELTA);
    stream = std::move(simplified);
    SIGNPOST_END(app.run);

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

    if (dist_flag && test_case_no != -1) {
        std::filesystem::path frechet_path = repo_root / "scripts" / "frechet";
        std::string cmd1 = std::string("\"") + frechet_path.string() + "\" --in " + std::to_string(test_case_no) + " --path \"" + (repo_root / "data" / std::to_string(test_case_no) / "simplify.txt").string() + "\"";
        int rc = std::system(cmd1.c_str());
    }

    return 0;
}
