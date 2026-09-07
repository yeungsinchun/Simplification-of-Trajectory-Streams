// Assert get_boundary_points_from_grid keeps only leftmost/rightmost sample
// on every y-row (no top/bottom full-row special case).
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "simplify_geometry.h"

static void check_lr_only(double epsilon, double delta) {
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

    assert(boundary_xs_by_y.size() == xs_by_y.size());
    for (const auto& [y, xs] : xs_by_y) {
        const auto& bx = boundary_xs_by_y.at(y);
        assert(bx.count(*xs.begin()) == 1);
        assert(bx.count(*xs.rbegin()) == 1);
        // At most two anchors per row (one when the row is a singleton).
        assert(bx.size() == (xs.size() == 1 ? 1u : 2u));
        assert(bx.size() < xs.size() || xs.size() <= 2);
    }
}

int main() {
    check_lr_only(0.5, 1.0);
    check_lr_only(1.0, 2.0);
    check_lr_only(0.25, 0.5);
    std::cout << "check_boundary_lr: ok\n";
    return 0;
}
