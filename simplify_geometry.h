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
//  Polygon intersection: double-precision Sutherland-Hodgman
// ===========================================================================
//
// Used for F(S,p) ∩ conv(G_i). Both are convex here, but SH only requires the
// cropping polygon to be convex. O(n*m) in doubles; n,m are small (~4..20).
namespace sh_double {

using Vec2 = std::array<double, 2>;

// Directed edge start -> start + (dx, dy) of a cropping polygon.
struct ClipEdge {
    Vec2 start;
    double dx, dy;
};

inline Vec2 to_vec2(const Point &p) {
    return {CGAL::to_double(p.x()), CGAL::to_double(p.y())};
}

// Point on segment p→q where orient along that segment is zero:
//   (1-t)*orient_p + t*orient_q = 0  =>  t = orient_p / (orient_p - orient_q).
inline Vec2 crossing_on_segment(const Vec2& p, const Vec2& q,
                                double orient_p, double orient_q) {
    const double t = orient_p / (orient_p - orient_q);
    return {p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])};
}

// Keep the part of polygon[0..n) in the closed left half-plane of `edge`,
// writing it to `cropped`. Returns false, leaving `cropped` unused, when the
// whole polygon is already inside (the common case).
__attribute__((always_inline)) inline bool
crop_to_left_of_edge(const Vec2 *polygon, size_t n, const ClipEdge &edge,
                     Vec2 *cropped, size_t &cropped_size) {
    cropped_size = 0;
    if (n == 0)
        return false;
    const double x0 = edge.start[0], y0 = edge.start[1];
    const double dx = edge.dx, dy = edge.dy;
    // Twice the signed area of (start, start + (dx, dy), p): positive left.
    auto side = [&](const Vec2 &p) __attribute__((always_inline)) {
        return dx * (p[1] - y0) - dy * (p[0] - x0);
    };
    double orient_prev = side(polygon[n - 1]);
    size_t i = 0;
    double orient_curr = 0.0;
    for (; i < n; ++i) {
        orient_curr = side(polygon[i]);
        if (orient_curr < 0.0 || orient_prev < 0.0)
            break;
        orient_prev = orient_curr;
    }
    if (i == n)
        return false;

    Vec2 prev = polygon[i == 0 ? n - 1 : i - 1];
    std::copy_n(polygon, i, cropped);
    cropped_size = i;
    auto emit_vertex = [&](const Vec2 &curr,
                           double orientation) __attribute__((always_inline)) {
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
    emit_vertex(polygon[i], orient_curr);
    for (++i; i < n; ++i) {
        const Vec2 &curr = polygon[i];
        emit_vertex(curr, side(curr));
    }
    return true;
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
 * real feature size. Prevents near-zero-length edges from reaching the clip.
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

// Reverse a clockwise ring so it runs counter-clockwise.
inline void make_ccw(std::vector<Point> &ring) {
    double twice_area = 0.0;
    const size_t n = ring.size();
    for (size_t i = 0; i < n; ++i) {
        const sh_double::Vec2 a = sh_double::to_vec2(ring[i]);
        const sh_double::Vec2 b = sh_double::to_vec2(ring[(i + 1) % n]);
        twice_area += a[0] * b[1] - a[1] * b[0];
    }
    if (twice_area < 0.0)
        std::reverse(ring.begin(), ring.end());
}

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

    bool contains(const Point &p) const {
        const double x = CGAL::to_double(p.x()), y = CGAL::to_double(p.y());
        return x >= min_x && x <= max_x && y >= min_y && y <= max_y;
    }
};

// A convex polygon in the form a clip against it consumes: the CCW edges of
// its deduplicated ring, in doubles, plus its bounding box. Build it once and
// clip any number of polygons against it.
class ConvexRegion {
  public:
    void assign(const std::vector<Point> &ring) {
        dedup_into(ring, ring_);
        make_ccw(ring_);
        edges_.clear();
        edges_.reserve(ring_.size());
        bounds_ = {};
        for (size_t i = 0; i < ring_.size(); ++i) {
            const sh_double::Vec2 a = sh_double::to_vec2(ring_[i]);
            const sh_double::Vec2 b =
                sh_double::to_vec2(ring_[(i + 1) % ring_.size()]);
            edges_.push_back({a, b[0] - a[0], b[1] - a[1]});
            bounds_.include(a[0], a[1]);
        }
    }

    const std::vector<sh_double::ClipEdge> &edges() const { return edges_; }
    const AxisBounds &bounds() const { return bounds_; }

  private:
    std::vector<Point> ring_;
    std::vector<sh_double::ClipEdge> edges_;
    AxisBounds bounds_;
};

// Sutherland-Hodgman clipping against a ConvexRegion. Owns the two rings the
// subject is cropped back and forth between, so repeated clips reuse them.
class ConvexClipper {
  public:
    /**
     * @brief result = subject ∩ region, with consecutive near-duplicate
     *        vertices removed.
     *
     * @param subject CCW ring without consecutive duplicates (find_F emits
     *        one), so it is cropped as given.
     * @param result_bounds Optional; receives the bounding box of result.
     * @return true iff result has at least 3 vertices; otherwise it is empty.
     */
    __attribute__((always_inline)) bool
    clip(const std::vector<Point> &subject, const ConvexRegion &region,
         std::vector<Point> &result, AxisBounds *result_bounds = nullptr) {
        if (subject.size() < 3 || region.edges().size() < 3) {
            result.clear();
            return false;
        }
        Ring *current = &first_;
        Ring *cropped = &second_;
        current->reserve(subject.size());
        current->size = 0;
        for (const Point &p : subject)
            current->data()[current->size++] = sh_double::to_vec2(p);

        for (const auto &edge : region.edges()) {
            if (current->size < 3)
                break;
            // Cropping by one edge adds at most one vertex per crossing.
            cropped->reserve(current->size * 2 + 2);
            if (sh_double::crop_to_left_of_edge(current->data(), current->size,
                                                edge, cropped->data(),
                                                cropped->size))
                std::swap(current, cropped);
        }

        result.clear();
        if (result_bounds)
            *result_bounds = {};
        if (current->size >= 3) {
            constexpr double EPS2 = 1e-12;
            auto close = [](const sh_double::Vec2 &a,
                            const sh_double::Vec2 &b) {
                const double dx = a[0] - b[0], dy = a[1] - b[1];
                return dx * dx + dy * dy <= EPS2;
            };
            result.reserve(current->size);
            sh_double::Vec2 first{}, last{};
            for (size_t i = 0; i < current->size; ++i) {
                const auto &v = current->data()[i];
                if (i != 0 && close(v, last))
                    continue;
                result.emplace_back(v[0], v[1]);
                if (result_bounds)
                    result_bounds->include(v[0], v[1]);
                if (i == 0)
                    first = v;
                last = v;
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

  private:
    struct Ring {
        std::vector<sh_double::Vec2> storage;
        size_t size = 0;

        sh_double::Vec2 *data() { return storage.data(); }
        void reserve(size_t n) {
            if (storage.size() < n)
                storage.resize(n);
        }
    };
    Ring first_, second_;
};

/**
 * @brief Convex-convex intersection P ∩ Q via dedup + Sutherland-Hodgman.
 *
 * @param P_in   Subject polygon (typically F(S,p)).
 * @param Q_in   Clip polygon; must be convex (typically conv(G_i)).
 * @param result Cleared; on success holds CCW intersection (>=3 verts).
 * @return true iff result is a non-degenerate polygon.
 */
inline bool intersect(const std::vector<Point> &P_in,
                      const std::vector<Point> &Q_in,
                      std::vector<Point> &result) {
    thread_local std::vector<Point> subject;
    thread_local ConvexRegion region;
    thread_local ConvexClipper clipper;
    dedup_into(P_in, subject);
    make_ccw(subject);
    region.assign(Q_in);
    return clipper.clip(subject, region, result);
}

// ===========================================================================
//  Trajectory-simplification geometry
// ===========================================================================
//
// The reachable-region construction from the streaming simplification
// algorithm: the bounding box, the delta-disk grid, and the free-space wedge
// F(S,p). Shared by the headless core (simplify_core.h), web trace
// (simplify.cpp), and GUI (simplify_with_gui.cpp).

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

// A shape derived from the grid samples around p, stored as offsets from p.
//
// The grid-corner offsets depend only on (EPSILON, DELTA), never on p: every
// call to get_points_from_grid(p, ...) yields the same corner set merely
// translated by p.  Constructions that commute with translation (convex hull,
// row outline) therefore give the same shape for every stream point, so the
// shape is built once around the origin per (EPSILON, DELTA) and each query
// is an O(h) copy-with-offset instead of a per-step O(m log m) rebuild.
class TranslatedGridShape {
  public:
    // `build` maps the grid samples around the origin to the shape there.
    template <class Build>
    std::vector<Point> at(const Point &p, double EPSILON, double DELTA,
                          int multiplier, Build &&build) {
        if (EPSILON != eps_ || DELTA != delta_ || multiplier != multiplier_) {
            const std::vector<Point> shape = build(
                get_points_from_grid(Point(0, 0), EPSILON, DELTA, multiplier));
            offsets_.clear();
            offsets_.reserve(shape.size());
            for (const auto &q : shape)
                offsets_.push_back(
                    {CGAL::to_double(q.x()), CGAL::to_double(q.y())});
            eps_ = EPSILON;
            delta_ = DELTA;
            multiplier_ = multiplier;
        }
        const double px = CGAL::to_double(p.x());
        const double py = CGAL::to_double(p.y());
        std::vector<Point> out;
        out.reserve(offsets_.size());
        for (const auto &off : offsets_)
            out.emplace_back(px + off[0], py + off[1]);
        return out;
    }

  private:
    double eps_ = std::numeric_limits<double>::quiet_NaN();
    double delta_ = std::numeric_limits<double>::quiet_NaN();
    int multiplier_ = 0;
    std::vector<std::array<double, 2>> offsets_;
};

// Convex hull of the delta-disk grid samples around p (the region conv(G_i)).
inline std::vector<Point> get_conv_from_grid(const Point &p, double EPSILON,
                                             double DELTA, int multiplier = 1) {
    thread_local TranslatedGridShape hull;
    return hull.at(p, EPSILON, DELTA, multiplier,
                   [](const std::vector<Point> &samples) {
                       std::vector<Point> conv;
                       CGAL::convex_hull_2(samples.begin(), samples.end(),
                                           std::back_inserter(conv));
                       return conv;
                   });
}

// Boundary anchors for P: discrete convex outline of the grid samples -
// leftmost and rightmost on every y-row, plus every sample on the topmost
// and bottommost rows.
inline std::vector<Point> get_boundary_points_from_grid(const Point& p, double EPSILON, double DELTA, int multiplier = 1) {
    thread_local TranslatedGridShape outline;
    return outline.at(
        p, EPSILON, DELTA, multiplier, [](const std::vector<Point> &all) {
            // get_points_from_grid is row-major: y descending, then x
            // ascending.
            std::vector<Point> boundary;
            boundary.reserve(all.size());
            if (all.empty())
                return boundary;
            const double y_top = CGAL::to_double(all.front().y());
            const double y_bot = CGAL::to_double(all.back().y());
            for (size_t i = 0; i < all.size();) {
                size_t j = i + 1;
                const double y = CGAL::to_double(all[i].y());
                while (j < all.size() && CGAL::to_double(all[j].y()) == y)
                    ++j;
                if (y == y_top || y == y_bot) {
                    boundary.insert(
                        boundary.end(),
                        all.begin() + static_cast<std::ptrdiff_t>(i),
                        all.begin() + static_cast<std::ptrdiff_t>(j));
                } else {
                    boundary.push_back(all[i]);
                    if (j - 1 != i)
                        boundary.push_back(all[j - 1]);
                }
                i = j;
            }
            return boundary;
        });
}

/**
 * @brief Indices of the two supporting (tangent) vertices from p to convex S.
 *
 * @pre p lies outside convex S. Then all of S sits in a <180° angular wedge
 *      from p, so bearing is a total order.
 * @param p External query point.
 * @param S Convex polygon vertices.
 * @return {min_idx, max_idx} of the most-CW and most-CCW vertices, or
 *         {-1, -1} if both coincide (degenerate).
 *
 * Algorithm: O(n) scan; keep argmin / argmax of orientation(p, S[i], S[j])
 * in pure doubles. Collinear ties: either vertex is a valid support.
 */
inline std::array<int, 2> find_tangent_idx(const Point &p,
                                           const std::vector<Point> &S) {
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
            rt = j;
            rtx = wx;
            rty = wy;
        }
        if (ltx * wy - lty * wx > 0.0) {
            lt = j;
            ltx = wx;
            lty = wy;
        }
    }
    if (rt == lt)
        return {-1, -1};
    return {std::min(rt, lt), std::max(rt, lt)};
}

// What find_F built.
enum class Wedge {
    whole_box,     // F is the whole working bbox
    cone,          // F is the part of the bbox beyond S as seen from p (empty
                   // when S offers no two distinct tangents)
    misses_target, // F provably misses the target box; F is left empty
};

/**
 * @brief Build free-space wedge F(S,p) into @p F (cleared first).
 *
 * F is the region of the working bbox reachable from convex stab region S
 * through external point p.
 *
 * Cases:
 *   1. |S|==1 or p ∈ S     → F = whole bbox.
 *   2. No two tangents     → F left empty (return).
 *   2b. Optional prune: if @p target is given and a tangent half-plane
 *       misses that box, return misses_target and leave F empty.
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
 * @param target Optional bbox of the region F is about to be intersected
 *        with, used only for the case-2b prune.
 * @param S_bounds Optional bbox of S; skips the p ∈ S test when p is outside.
 * @param p_outside_S Optional latch: once p lies outside S it stays outside
 *        every later S, so the p ∈ S test is skipped. Cleared on whole_box.
 */
__attribute__((always_inline)) inline Wedge
find_F(const Point &p, const std::vector<Point> &S, std::vector<Point> &F,
       const AxisBounds *target = nullptr, const AxisBounds *S_bounds = nullptr,
       bool *p_outside_S = nullptr) {
    F.clear();
    assert(S.size() != 2);
    auto use_bbox = [&]() __attribute__((always_inline)) {
        const auto corners = current_bbox_corner();
        F.assign(corners.begin(), corners.end());
        if (p_outside_S)
            *p_outside_S = false;
        return Wedge::whole_box;
    };
    if (S.size() == 1)
        return use_bbox();
    if (!p_outside_S || !*p_outside_S) {
        if ((!S_bounds || S_bounds->contains(p)) && point_in_convex(p, S))
            return use_bbox();
        // Each later stab polygon is contained in the continuation wedge.
        // The anchor lies outside that wedge once it lies outside S.
        if (p_outside_S)
            *p_outside_S = true;
    }

    const auto tangent = find_tangent_idx(p, S);
    if (tangent[0] < 0)
        return Wedge::cone;

    const double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y());
    const double ax = CGAL::to_double(S[tangent[0]].x()) - px;
    const double ay = CGAL::to_double(S[tangent[0]].y()) - py;
    const double bx = CGAL::to_double(S[tangent[1]].x()) - px;
    const double by = CGAL::to_double(S[tangent[1]].y()) - py;
    const double turn = ax * by - ay * bx;
    if (target && turn != 0.0) {
        // A tangent half-plane that misses the whole target box also misses
        // everything inside it, so F ∩ target is empty.
        auto range = [&](double tx, double ty) __attribute__((always_inline)) {
            const double xmin = target->min_x - px;
            const double xmax = target->max_x - px;
            const double ymin = target->min_y - py;
            const double ymax = target->max_y - py;
            const double lo =
                tx * (tx >= 0 ? ymin : ymax) - ty * (ty >= 0 ? xmax : xmin);
            const double hi =
                tx * (tx >= 0 ? ymax : ymin) - ty * (ty >= 0 ? xmin : xmax);
            const double tolerance = 64.0 *
                                     std::numeric_limits<double>::epsilon() *
                                     (std::abs(tx) + std::abs(ty)) *
                                     (std::abs(xmin) + std::abs(xmax) +
                                      std::abs(ymin) + std::abs(ymax) + 1.0);
            return std::array<double, 3>{lo, hi, tolerance};
        };
        const auto first = range(ax, ay);
        const auto second = range(bx, by);
        const bool misses =
            turn > 0.0 ? (first[1] < -first[2] || second[0] > second[2])
                       : (first[0] > first[2] || second[1] < -second[2]);
        if (misses)
            return Wedge::misses_target;
    }

    auto hit1 = ray_hit_bbox(p, S[tangent[0]]);
    auto hit2 = ray_hit_bbox(p, S[tangent[1]]);
    if (!hit1 || !hit2)
        return use_bbox();

    auto e1 = which_edge(hit1.value());
    auto e2 = which_edge(hit2.value());
    if (!e1 || !e2)
        return use_bbox();

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
    return Wedge::cone;
}

#endif // SIMPLIFY_GEOMETRY_H
