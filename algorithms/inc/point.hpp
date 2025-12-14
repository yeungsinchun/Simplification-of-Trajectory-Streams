#ifndef POINT_HPP
#define POINT_HPP
#include <cmath>
// x,y is 2d coordinate calculated by Mercator Projection
// https://en.wikipedia.org/wiki/Mercator_projection

class Point {
  public:
    double x, y;
    Point() {}
    Point(double x_, double y_) : x{x_}, y{y_} {}
    double distance(const Point &rhs) {
        return sqrt((x - rhs.x) * (x - rhs.x) + (y - rhs.y) * (y - rhs.y));
    }
};

#endif