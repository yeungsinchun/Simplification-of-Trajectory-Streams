#pragma once

#include <CGAL/Polygon_2.h>
#include <vector>
#include <QWidget>
#include <QColor>
#include <QString>
#include <QRect>
#include <QMouseEvent>
#include <limits>
#include <optional>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include "simplify_geometry.h"



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

    void setOriginalVisible(bool v)   { original_visible_   = v; update(); }
    void setSimplifiedVisible(bool v) { simplified_visible_ = v; update(); }
    bool isOriginalVisible()   const { return original_visible_; }
    bool isSimplifiedVisible() const { return simplified_visible_; }

    // Generic colored curves with labels (for plotting multiple results)
    void addCurve(const std::vector<Point>& pts, const QColor& color, const QString& label);
    void clearCurves();

    const std::vector<Point>& original()   const { return original_; }
    const std::vector<Point>& simplified() const { return simplified_; }

    void setParameters(double delta, double epsilon);
    void setShowLabels(bool show);
    bool isViewFrozen() const { return view_frozen_; }

    // mark a special reference point (p0) for debugging/visualization
    void markP0(const Point& p);
    void clearMarkedP0();
    // mark the currently processed point (pi) distinctly from p0
    void markPi(const Point& p);
    void clearMarkedPi();

protected:
    void paintEvent(QPaintEvent*) override;
    void mousePressEvent(QMouseEvent*) override;

private:
    struct ColoredPoly { Polygon poly; QColor color; };
    std::vector<ColoredPoly> polys_;
    std::vector<std::vector<Point>> pointSets_;

    // Trajectory streams
    std::vector<Point> original_;
    std::vector<Point> simplified_;
    bool original_visible_   = true;
    bool simplified_visible_ = true;

    // Display parameters
    double delta_   = std::numeric_limits<double>::quiet_NaN();
    double epsilon_ = std::numeric_limits<double>::quiet_NaN();
    bool showLabels_ = false;

    // Frozen view extent (kept for API compatibility).
    bool   view_frozen_ = false;

    // Optional marked p0 for debugging
    std::optional<Point> marked_p0_;
    // Optional marked pi for debugging
    std::optional<Point> marked_pi_;
    // Special debug objects
    std::vector<Polygon> special_polys_;
    std::vector<Point>   special_points_;

    struct Curve { std::vector<Point> pts; QColor color; QString label; bool visible = true; };
    std::vector<Curve> curves_;

    // Legend hit-test rectangles filled during paintEvent; cleared on repaint.
    // Index into the legend entries in the same order they are drawn:
    //   [-2] original, [-1] simplified, [0..N-1] curves.
    std::vector<QRect> legend_rects_;
};

// process pending GUI events (call after adding points for incremental display)
void viewer_process_events();