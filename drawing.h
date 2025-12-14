#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <vector>
#include <QWidget>
#include <QColor>
#include <QString>
#include <limits>
#include <optional>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using Point   = Kernel::Point_2;
using Polygon = CGAL::Polygon_2<Kernel>;

class MultiViewer : public QWidget {
    Q_OBJECT
public:
    explicit MultiViewer(QWidget* parent = nullptr);

    void addPolygon(const Polygon& poly, const QColor& color = Qt::black);
    void addPoints(const std::vector<Point>& pts);
    // Special debug helpers: add polygons/points that should be highlighted
    // separately from regular content.
    void addSpecialPolygon(const Polygon& poly);
    void addSpecialPoint(const Point& p);
    void clearAll();
    void clearSpecials();
    // Clear all overlaid polygons (F/G/S debug shapes)
    void clearPolygons();

    void addOriginalPoint(const Point& p);
    void addSimplifiedPoint(const Point& p);
    void addOriginalPoints(const std::vector<Point>& v);
    void addSimplifiedPoints(const std::vector<Point>& v);

    void clearOriginal();
    void clearSimplified();

    // Generic colored curves with labels (for plotting multiple results)
    void addCurve(const std::vector<Point>& pts, const QColor& color, const QString& label);
    void clearCurves();

    const std::vector<Point>& original()   const { return original_; }
    const std::vector<Point>& simplified() const { return simplified_; }

    void setParameters(double delta, double epsilon);
    void setShowLabels(bool show);

    // mark a special reference point (p0) for debugging/visualization
    void markP0(const Point& p);
    void clearMarkedP0();
    // mark the currently processed point (pi) distinctly from p0
    void markPi(const Point& p);
    void clearMarkedPi();

protected:
    void paintEvent(QPaintEvent*) override;

private:
    struct ColoredPoly { Polygon poly; QColor color; };
    std::vector<ColoredPoly> polys_;
    std::vector<std::vector<Point>> pointSets_;

    // Trajectory streams
    std::vector<Point> original_;
    std::vector<Point> simplified_;

    // Display parameters
    double delta_   = std::numeric_limits<double>::quiet_NaN();
    double epsilon_ = std::numeric_limits<double>::quiet_NaN();
    bool showLabels_ = false;

    // Optional marked p0 for debugging
    std::optional<Point> marked_p0_;
    // Optional marked pi for debugging
    std::optional<Point> marked_pi_;
    // Special debug objects
    std::vector<Polygon> special_polys_;
    std::vector<Point>   special_points_;

    struct Curve { std::vector<Point> pts; QColor color; QString label; };
    std::vector<Curve> curves_;

    void compute_bbox(double& minx,dou