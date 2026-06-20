#ifndef SIMPLIFY_GEOMETRY_H
#define SIMPLIFY_GEOMETRY_H

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/intersections.h>

using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;

// Convert Epick <-> Epeck for the intersection construction functions.
// Defined inline so each translation unit gets its own instance (avoids ODR).
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

#endif // SIMPLIFY_GEOMETRY_H
