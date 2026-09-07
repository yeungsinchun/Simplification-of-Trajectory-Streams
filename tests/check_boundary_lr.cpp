// Assert get_boundary_points_from_grid forms the discrete convex outline:
// leftmost/rightmost on every y-row, and every sample on topmost/bottommost rows.
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "simplify_geometry.h"

static void check_hull_outline(double epsilon, double delta) {
    const Point origin(0, 0);
    const auto all = get_points_from_grid(origin, epsilon, delta);
    const auto boundary = get_boundary_points_from_grid(origin, epsilon, delta);

    std::map<double, std::set<double>> xs_by_y;
    for (const auto& q : all) {
        xs_by_y[CGAL::to_double(q.y())].insert(CGAL::to_double(q.x()));
    }

    std::map<double, std::set<double>> boundary_xs_by_y;
    for (const auto& q : boundary) {
        boundary_xs_by_y[CGAL::to_double(q.y())].insert(CGAL::to_double(q.x()));
    }

    assert(!xs_by_y.empty());
    assert(boundary_xs_by_y.size() == xs_by_y.size());

    const double y_top = xs_by_y.rbegin()->first;  // highest y
    const double y_bot = xs_by_y.begin()->first;   // lowest y

    for (const auto& [y, xs] : xs_by_y) {
        const auto& bx = boundary_xs_by_y.at(y);
        assert(bx.count(*xs.begin()) == 1);
        assert(bx.count(*xs.rbegin()) == 1);
        if (y == y_top || y == y_bot) {
            assert(bx == xs);
        } else {
            // Middle rows: leftmost and rightmost only (one when singleton).
            assert(bx.size() == (xs.size() == 1 ? 1u : 2u));
            assert(bx.size() < xs.size() || xs.size() <= 2);
        }
    }
}

int main() {
    check_hull_outline(0.5, 1.0);
    check_hull_outline(1.0, 2.0);
    check_hull_outline(0.25, 0.5);
    std::cout << "check_boundary_lr: ok\n";
    return 0;
}
