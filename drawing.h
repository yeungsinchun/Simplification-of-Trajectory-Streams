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

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;
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
    // Set the data-bbox padding as a fraction of the larger data extent.
    // Default 0.08 (8%). Larger values push the data away from the canvas edges.
    void setPadFraction(double frac) { pad_frac_ = frac; update(); }
    double padFraction() const { return pad_frac_; }

    // Vertical shift of the canvas in pixels. Positive value moves the
    // data upward on screen (gives more room at the bottom). Default 0.
    void setVShift(int px) { vshift_ = px; update(); }
    int vShift() const { return vshift_; }

    // Freeze the view extent to the supplied data bounding box. When frozen,
    // subsequent add*() calls do not affect the visible region — useful when
    // the data streams in incrementally and we don't want the canvas to
    // re-zoom on every repaint. Pass NaNs to unfreeze (revert to dynamic).
    void setDataBBox(double minx, double miny, double maxx, double maxy);
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
    double pad_frac_ = 0.08;
    int    vshift_   = 0;

    // Frozen view extent. When view_frozen_ is true, paintEvent uses
    // these values instead of recomputing from the data each frame.
    bool   view_frozen_ = false;
    double frozen_minx_ = 0, frozen_miny_ = 0;
    double frozen_maxx_ = 0, frozen_maxy_ = 0;

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

    void compute_bbox(double& minx,double& miny,double& maxx,double& maxy) const;
};

// process pending GUI events (call after adding points for incremental display)
void viewer_process_events();