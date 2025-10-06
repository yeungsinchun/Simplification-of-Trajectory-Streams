#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <vector>
#include <QWidget>

using Kernel  = CGAL::Exact_predicates_exact_constructions_kernel;
using Point   = Kernel::Point_2;
using Polygon = CGAL::Polygon_2<Kernel>;

class MultiViewer : public QWidget {
    Q_OBJECT
public:
    explicit MultiViewer(QWidget* parent = nullptr);

    // Existing (polygons & arbitrary point sets)
    void addPolygon(const Polygon& poly);
    void addPoints(const std::vector<Point>& pts);
    void clearAll();

    // New: two incremental trajectory streams
    void addOriginalPoint(const Point& p);
    void addSimplifiedPoint(const Point& p);
    void addOriginalPoints(const std::vector<Point>& v);
    void addSimplifiedPoints(const std::vector<Point>& v);

    void clearOriginal();
    void clearSimplified();

    const std::vector<Point>& original()   const { return original_; }
    const std::vector<Point>& simplified() const { return simplified_; }

protected:
    void paintEvent(QPaintEvent*) override;

private:
    std::vector<Polygon> polys_;
    std::vector<std::vector<Point>> pointSets_;

    // Trajectory streams
    std::vector<Point> original_;
    std::vector<Point> simplified_;

    void compute_bbox(double& minx,double& miny,double& maxx,double& maxy) const;
};

// process pending GUI events (call after adding points for incremental display)
void viewer_process_events();