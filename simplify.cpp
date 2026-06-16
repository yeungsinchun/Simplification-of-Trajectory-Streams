#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
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
#include <chrono>
#include <iomanip>
#include <memory>
#include <QApplication>
#include "drawing.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>


// Timing instrumentation
static std::map<std::string, double> g_timing;
static std::map<std::string, long long> g_counters;
static std::map<std::string, double> g_child_time;
static std::vector<const char*> g_timer_stack;


struct Timer {
    std::chrono::high_resolution_clock::time_point start;
    const char* name;
    Timer(const char* n) : name(n) { 
        start = std::chrono::high_resolution_clock::now();
        g_timer_stack.push_back(name);
    }
    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration<double, std::milli>(end - start).count();
        
        // Add to this timer's total (wall time)
        g_timing[name] += dur;
        g_counters[name]++;
        
        // Add this duration to parent's child_time accumulator
        if (g_timer_stack.size() >= 2) {
            const char* parent = g_timer_stack[g_timer_stack.size() - 2];
            g_child_time[parent] += dur;
        }
        
        g_timer_stack.pop_back();
    }
};

#define TIMER(name) Timer _timer(name)

void print_timing_summary() {
    if (g_timing.empty()) return;
    
    // Calculate overhead for each timer: wall_time - sum_of_children_time
    std::map<std::string, double> overhead;
    for (const auto& kv : g_timing) {
        overhead[kv.first] = kv.second - g_child_time[kv.first];
    }
    
    // Define parent-child relationships
    std::map<std::string, std::vector<std::string>> children_of;
    children_of["total"] = {"simplify"};
    children_of["simplify"] = {"get_longest_stab"};
    children_of["get_longest_stab"] = {"intersect", "get_conv_from_grid", "find_F", "update_S", "polygon_construction", "get_points_from_grid", "gui_update"};
    children_of["intersect"] = {"intersection_total"};
    children_of["intersection_total"] = {"intersection_convert"};
    children_of["intersection_convert"] = {"intersection_convex_hull", "intersection_segment_intersect", "intersection_point_in_polygon_P", "intersection_point_in_polygon_Q"};
    children_of["find_F"] = {};  // find_F has no named children in our instrumentation
    
    // Build sorted list by wall time
    std::vector<std::pair<std::string, double>> sorted;
    for (const auto& kv : g_timing) {
        if (kv.second > 0.1) {
            sorted.emplace_back(kv.first, kv.second);
        }
    }
    std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });
    
    // Map timer to its wall time for quick lookup
    std::map<std::string, double> wall_time;
    for (const auto& kv : sorted) wall_time[kv.first] = kv.second;
    
    // Find root timers (those that are never children)
    std::set<std::string> is_child;
    for (const auto& kv : children_of) {
        for (const auto& child : kv.second) {
            is_child.insert(child);
        }
    }
    std::vector<std::string> roots;
    for (const auto& kv : wall_time) {
        if (is_child.find(kv.first) == is_child.end()) {
            roots.push_back(kv.first);
        }
    }
    std::sort(roots.begin(), roots.end(), [&](const std::string& a, const std::string& b) {
        return wall_time[a] > wall_time[b];
    });
    
    double total_wall = 0;
    // Sum all root-level timers (those that have no parent in our tree)
    std::set<std::string> parent_of;
    for (const auto& kv : children_of) {
        for (const auto& child : kv.second) {
            parent_of.insert(child);
        }
    }
    for (const auto& kv : wall_time) {
        if (parent_of.find(kv.first) == parent_of.end()) {
            total_wall += kv.second;
        }
    }
    
    auto visual_cols = [](const std::string& s) {
        int cols = 0;
        for (size_t i = 0; i < s.size(); ++i) {
            unsigned char c = s[i];
            if ((c & 0x80) == 0) {
                cols += 1; // ASCII
            } else if ((c & 0xF0) == 0xE0) {
                cols += 1; i += 2; // 3-byte UTF-8 (CJK, box drawing, etc.)
            } else if ((c & 0xE0) == 0xC0) {
                cols += 1; i += 1; // 2-byte UTF-8
            } else if ((c & 0xF8) == 0xF0) {
                cols += 2; i += 3; // 4-byte (emoji)
            } else {
                cols += 1; // continuation or unknown
            }
        }
        return cols;
    };
    
    std::cout << "\n========== TIMING SUMMARY ==========\n";
    std::cout << std::left << std::setw(50) << "Operation" 
              << std::right << std::setw(10) << "Wall (ms)"
              << std::setw(10) << "Self (ms)"
              << std::setw(10) << "Calls" 
              << std::setw(8) << "% Total" << "\n";
    std::cout << std::string(88, '-') << "\n";
    
    std::function<void(const std::string&, const std::string&, bool)> print_tree;
    print_tree = [&](const std::string& name, const std::string& prefix, bool is_last) {
        double wall = wall_time[name];
        double self = overhead[name];
        long long c = g_counters[name];
        double pct = 100.0 * wall / total_wall;
        
        std::string connector = is_last ? "└─ " : "├─ ";
        
        // Calculate visual columns so numbers always align at fixed position
        int prefix_visual = visual_cols(prefix);
        int name_visual = visual_cols(name);
        int connector_visual = 3; // └─ or ├─ are each 3 visual columns
        int total_visual = prefix_visual + connector_visual + name_visual;
        int padding = std::max(0, 48 - total_visual);
        
        std::cout << prefix << connector << name 
                  << std::string(padding, ' ')
                  << std::right << std::fixed << std::setprecision(2) << std::setw(10) << wall
                  << std::setw(10) << self
                  << std::setw(10) << c
                  << std::setw(7) << std::fixed << std::setprecision(1) << pct << "%\n";
        
        if (children_of.count(name)) {
            auto& kids = children_of[name];
            std::sort(kids.begin(), kids.end(), [&](const std::string& a, const std::string& b) {
                return wall_time[a] > wall_time[b];
            });
            for (size_t i = 0; i < kids.size(); i++) {
                bool last = (i == kids.size() - 1);
                std::string child_prefix = prefix + (is_last ? "  " : "│ ");
                print_tree(kids[i], child_prefix, last);
            }
        }
    };
    
    for (size_t i = 0; i < roots.size(); i++) {
        print_tree(roots[i], "", i == roots.size() - 1);
    }
    
    std::cout << std::string(88, '-') << "\n";
    std::cout << std::left << std::setw(50) << "TOTAL (wall time)" 
              << std::right << std::setw(10) << std::fixed << std::setprecision(2) << total_wall << "\n";
    std::cout << "====================================\n";
}

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = Kernel::Point_2;
using Segment = Kernel::Segment_2;
using Ray = Kernel::Ray_2;
using Bbox = CGAL::Bbox_2;
using Rect = CGAL::Iso_rectangle_2<Kernel>;
using Polygon = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes =  CGAL::Polygon_with_holes_2<Kernel>;

// Filtered kernel for intersection: exact predicates with double arithmetic
// Faster than Epeck for intersection computation, more robust than plain Simple_cartesian
using IntersectKernel = CGAL::Filtered_kernel<CGAL::Simple_cartesian<double>>;
using IntersectPoint = IntersectKernel::Point_2;
using IntersectSegment = IntersectKernel::Segment_2;
using IntersectPolygon = CGAL::Polygon_2<IntersectKernel>;

// Helper to convert from Kernel::Point_2 to IntersectPoint using double as intermediate
static IntersectPoint to_intersect(const Point& p) {
    return IntersectPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
}

// Helper to convert from IntersectPoint back to Point using double as intermediate
static Point from_intersect(const IntersectPoint& p) {
    return Point(static_cast<double>(p.x()), static_cast<double>(p.y()));
}

// Convex polygon intersection using boundary point collection + convex hull
// Uses Simple_cartesian<long double> for faster intersection computations
static std::vector<Point> get_vertices(const Polygon& poly) {
    std::vector<Point> verts;
    for (auto it = poly.vertices_begin(); it != poly.vertices_end(); ++it) {
        verts.push_back(*it);
    }
    return verts;
}

static bool point_in_polygon(const IntersectPoint& p, const std::vector<IntersectPoint>& poly) {
    int n = poly.size();
    for (int i = 0; i < n; i++) {
        const IntersectPoint& a = poly[i];
        const IntersectPoint& b = poly[(i + 1) % n];
        if (!CGAL::left_turn(a, b, p)) {
            return false;
        }
    }
    return true;
}

static std::optional<IntersectPoint> segment_intersection(const IntersectPoint& p1, const IntersectPoint& p2,
                                                         const IntersectPoint& q1, const IntersectPoint& q2) {
    auto inter = CGAL::intersection(IntersectSegment(p1, p2), IntersectSegment(q1, q2));
    if (inter && std::holds_alternative<IntersectPoint>(*inter)) {
        return std::get<IntersectPoint>(*inter);
    }
    return std::nullopt;
}

// Original: Epeck kernel + filtered intersection kernel + convex hull
static void convex_polygon_intersection(const Polygon& P, const Polygon& Q,
                                       std::vector<Polygon_with_holes>& result) {
    TIMER("intersection_total");
    {
        TIMER("intersection_convert");
        std::vector<IntersectPoint> P_verts, Q_verts;
        P_verts.reserve(P.size());
        Q_verts.reserve(Q.size());
        for (auto it = P.vertices_begin(); it != P.vertices_end(); ++it) {
            P_verts.push_back(to_intersect(*it));
        }
        for (auto it = Q.vertices_begin(); it != Q.vertices_end(); ++it) {
            Q_verts.push_back(to_intersect(*it));
        }

        if (P_verts.size() < 3 || Q_verts.size() < 3) return;

        int n = P_verts.size();
        int m = Q_verts.size();

        std::vector<IntersectPoint> candidates;
        std::set<std::pair<long double, long double>> seen;
        auto add_unique = [&](const IntersectPoint& p) {
            long double x = p.x();
            long double y = p.y();
            auto key = std::make_pair(x, y);
            if (seen.insert(key).second) {
                candidates.push_back(p);
            }
        };

        for (int i = 0; i < n; i++) {
            if (point_in_polygon(P_verts[i], Q_verts)) {
                add_unique(P_verts[i]);
            }
        }
        for (int j = 0; j < m; j++) {
            if (point_in_polygon(Q_verts[j], P_verts)) {
                add_unique(Q_verts[j]);
            }
        }

        for (int i = 0; i < n; i++) {
            const IntersectPoint& p1 = P_verts[i];
            const IntersectPoint& p2 = P_verts[(i + 1) % n];
            for (int j = 0; j < m; j++) {
                const IntersectPoint& q1 = Q_verts[j];
                const IntersectPoint& q2 = Q_verts[(j + 1) % m];
                auto inter = segment_intersection(p1, p2, q1, q2);
                if (inter) {
                    add_unique(*inter);
                }
            }
        }

        if (candidates.size() < 3) return;

        std::vector<IntersectPoint> hull;
        CGAL::convex_hull_2(candidates.begin(), candidates.end(), std::back_inserter(hull));

        if (hull.size() < 3) return;

        std::vector<Point> result_verts;
        result_verts.reserve(hull.size());
        for (const auto& hp : hull) {
            result_verts.push_back(from_intersect(hp));
        }

        Polygon inter_poly(result_verts.begin(), result_verts.end());

        if (!inter_poly.is_simple()) return;

        if (inter_poly.is_clockwise_oriented()) {
            inter_poly.reverse_orientation();
        }

        result.push_back(Polygon_with_holes(inter_poly));
    }
}

static void intersect(const Polygon& subject, const Polygon& clip,
                        std::back_insert_iterator<std::vector<Polygon_with_holes>> result) {
    std::vector<Polygon_with_holes> temp;
    convex_polygon_intersection(subject, clip, temp);
    for (auto& pwh : temp) {
        *result++ = std::move(pwh);
    }
}

// Convert CGAL Polygon to OpenCV Points
static std::vector<cv::Point2f> to_cv_points(const Polygon& poly) {
    std::vector<cv::Point2f> pts;
    pts.reserve(poly.size());
    for (auto it = poly.vertices_begin(); it != poly.vertices_end(); ++it) {
        pts.emplace_back(static_cast<float>(CGAL::to_double(it->x())),
                        static_cast<float>(CGAL::to_double(it->y())));
    }
    return pts;
}

// Convert OpenCV Points back to CGAL Polygon
static Polygon from_cv_points(const std::vector<cv::Point2f>& pts) {
    std::vector<Point> cg_pts;
    cg_pts.reserve(pts.size());
    for (const auto& p : pts) {
        cg_pts.emplace_back(p.x, p.y);
    }
    return Polygon(cg_pts.begin(), cg_pts.end());
}

// OpenCV-based polygon intersection using cv::intersectConvexConvex
static void intersect_opencv(const Polygon& A, const Polygon& B,
                            std::back_insert_iterator<std::vector<Polygon_with_holes>> result) {
    auto cvA = to_cv_points(A);
    auto cvB = to_cv_points(B);

    if (cvA.size() < 3 || cvB.size() < 3) {
        return;
    }

    std::vector<cv::Point2f> intersection;
    cv::intersectConvexConvex(cvA, cvB, intersection, true);

    if (intersection.empty() || intersection.size() < 3) {
        return;
    }

    // Convert back to CGAL
    Polygon inter_poly = from_cv_points(intersection);
    if (inter_poly.is_simple() && !inter_poly.is_clockwise_oriented()) {
        *result++ = Polygon_with_holes(inter_poly);
    }
}

bool showF = false; // true if -F is passed
bool showG = false; // true if -G is passed
bool showS = false; // true if -S is passed
bool showLabels = false; // true if --labels is passed
bool out_flag = false;      // write outputs if true
bool gui_flag = false;      // show viewer if true
bool dist_flag = false;     // compute frechet if true
bool dump_intersect_flag = false; // dump (F_poly, Gi_poly) pairs to disk

// Set by main() when --dump-intersect is passed.  Non-owning pointer to
// the open log file; nullptr disables dumping.  The dump helper at the
// intersect() call site reads this without changing any function signatures.
static std::ofstream* g_dump_intersect_stream = nullptr;

// Forward-declared: writes the (F_poly, Gi_poly) pair that is about to be
// passed to intersect() to the per-run log file (g_dump_intersect_stream).
// Each pair is appended; one record = "P <nv>\n<x y>\n... Q <nv>\n<x y>\n..."
// followed by a blank line.
static void dump_intersect_pair(const Polygon& F, const Polygon& G) {
    if (!g_dump_intersect_stream || !g_dump_intersect_stream->is_open()) return;
    auto& out = *g_dump_intersect_stream;
    auto write_one = [&](const char* tag, const Polygon& p) {
        out << tag << ' ' << p.size() << '\n';
        for (auto it = p.vertices_begin(); it != p.vertices_end(); ++it) {
            out << CGAL::to_double(it->x()) << ' '
                << CGAL::to_double(it->y()) << '\n';
        }
    };
    write_one("P", F);
    write_one("Q", G);
    out << '\n';
}

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
              << "  --pad <frac>     Data-bbox padding as a fraction of the larger extent (default 0.08)\n"
              << "  --vshift <px>    Shift canvas up by N pixels (default 0)\n"
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

// TODO: binary search
std::vector<int> find_tangent_idx(const Point &p,
                                 const std::vector<Point> &S) {
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
    assert(tangent.size() == 2);

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
    CGAL_precondition(Polygon(S.begin(), S.end()).is_counterclockwise_oriented());
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
    TIMER("get_longest_stab");
    const Point& p0 = stream[cur];
    int p0cur = cur;
    std::vector<Point> P;
    {
        TIMER("get_points_from_grid");
        P = get_points_from_grid(p0);
        // P = std::vector<Point>{p0};
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
        bool shown_debug = false;
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
            if (i == 0) {
                // std::cerr << "  [cur=" << cur << "] F.size()=" << F.size()
                //   << " Gi.size()=" << Gi.size() << '\n';
            }
            {
                TIMER("polygon_construction");
                F_poly = Polygon(F.begin(), F.end());
                Gi_poly = Polygon(Gi.begin(), Gi.end());
                S_poly = Polygon(S[i].begin(), S[i].end());
            }

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
                if (viewer) {
                    TIMER("gui_update");
                    if (showF) viewer->addPolygon(F_poly, c);
                    if (showG) viewer->addPolygon(Gi_poly, c);
                    if (showS) viewer->addPolygon(S_poly, c);
                }
                shown_debug = true;
            }

            {
                TIMER("intersect");
                if (dump_intersect_flag) dump_intersect_pair(F_poly, Gi_poly);
                intersect(F_poly, Gi_poly, back_inserter(new_S[i]));
                // intersect_opencv(F_poly, Gi_poly, back_inserter(new_S[i]));
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
    double pad_frac = 0.08;     // default data-bbox padding (8% of larger extent)
    int vshift = 0;             // pixels to shift the canvas up (positive = up)

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i],"-F") == 0) showF = true;
        else if (strcmp(argv[i],"-S") == 0) showS = true;
        else if (strcmp(argv[i],"-G") == 0) showG = true;
        else if (strcmp(argv[i],"--gui") == 0) gui_flag = true;
        else if (strcmp(argv[i],"--labels") == 0) showLabels = true;
        else if (strcmp(argv[i],"--out") == 0) out_flag = true;
        else if (strcmp(argv[i],"--dist") == 0) dist_flag = true;
        else if (strcmp(argv[i],"--dump-intersect") == 0) dump_intersect_flag = true;
        else if (strcmp(argv[i],"-d") == 0 && i+1 < argc) {
            try { DELTA = std::stod(argv[++i]); } catch(...) { std::cerr << "Invalid -d value\n"; return 1; }
        }
        else if (strcmp(argv[i],"-e") == 0 && i+1 < argc) {
            try { EPSILON = std::stod(argv[++i]); } catch(...) { std::cerr << "Invalid -e value\n"; return 1; }
        }
        else if (strcmp(argv[i],"--pad") == 0 && i+1 < argc) {
            try { pad_frac = std::stod(argv[++i]); }
            catch(...) { std::cerr << "Invalid --pad value\n"; return 1; }
        }
        else if (strcmp(argv[i],"--vshift") == 0 && i+1 < argc) {
            try { vshift = std::stoi(argv[++i]); }
            catch(...) { std::cerr << "Invalid --vshift value\n"; return 1; }
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
    try {
        repo_root = find_repo_root(argv[0], 5);
    } catch (const std::filesystem::filesystem_error&) {
        // argv[0] may be unresolvable; fall through to cwd.
    }
    if (repo_root.empty()) {
        repo_root = find_repo_root(std::filesystem::current_path(), 5);
        if (repo_root.empty()) {
            // Last resort: use cwd as-is even without a `data/` subdir.
            repo_root = std::filesystem::current_path();
        }
    }
    if (test_case_no != -1) {
        // Prefer original in data/<id>/original.txt; fallback to data/taxi/<id>.txt
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
        viewer.setShowLabels(showLabels);
        viewer.setPadFraction(pad_frac);
        viewer.setVShift(vshift);
        // Freeze the view to the full input stream extent so the canvas
        // does not re-zoom as points stream in during simplification.
        if (!stream.empty()) {
            double mn_x =  1e300, mn_y =  1e300;
            double mx_x = -1e300, mx_y = -1e300;
            for (const auto& p : stream) {
                double x = CGAL::to_double(p.x());
                double y = CGAL::to_double(p.y());
                if (x < mn_x) mn_x = x; if (x > mx_x) mx_x = x;
                if (y < mn_y) mn_y = y; if (y > mx_y) mx_y = y;
            }
            viewer.setDataBBox(mn_x, mn_y, mx_x, mx_y);
        }
        viewer.show();
    }

    // Optionally dump every (F_poly, Gi_poly) pair fed to intersect() to
    // data/<id>/intersect_pairs.txt.  This produces a stream of
    // real convex-convex intersection problems that can be replayed
    // against alternative intersection implementations.
    std::unique_ptr<std::ofstream> dump_stream;
    if (dump_intersect_flag) {
        if (test_case_no == -1) {
            std::cerr << "--dump-intersect requires --in <id> to determine "
                         "output location\n";
            return 1;
        }
        std::filesystem::path dir = repo_root / "data" / std::to_string(test_case_no);
        std::filesystem::create_directories(dir);
        std::filesystem::path dump_path = dir / "intersect_pairs.txt";
        dump_stream = std::make_unique<std::ofstream>(dump_path);
        if (!dump_stream->is_open()) {
            std::cerr << "Cannot open " << dump_path << " for writing\n";
            return 1;
        }
        g_dump_intersect_stream = dump_stream.get();
        std::cout << "Dumping (F_poly, Gi_poly) pairs to " << dump_path << '\n';
    }

    // Simplify
    std::vector<Point> simplified;
    {
        TIMER("total");
        simplified = simplify(stream, vptr);
        stream = std::move(simplified);
    }
    g_dump_intersect_stream = nullptr;
    if (dump_stream) {
        dump_stream->close();
        std::cout << "Intersection-pair dump complete\n";
    }
    print_timing_summary();

    // Optional output
    if (out_flag) {
        if (test_case_no == -1) {
            std::cerr << "--out requires --in <id> to determine output location\n";
        } else {
            std::filesystem::path dir = repo_root / "data" / std::to_string(test_case_no);
            std::filesystem::create_directories(dir);
            
            // simplify_*.txt only
            std::ofstream simp(dir / "simplify.txt");
            std::size_t N = stream.size();
            simp << N << '\n';
            for (const auto& p : stream) {
                simp << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y()) << '\n';
            }
            simp.close();
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
        std::string cmd1 = std::string("\"") + frechet_path.string() + "\" --in " + std::to_string(test_case_no) + " --path \"" + (repo_root / "data" / std::to_string(test_case_no) / "simplify.txt").string() + "\"";
        int rc = std::system(cmd1.c_str());
        (void)rc;
    }

    return gui_result;
}