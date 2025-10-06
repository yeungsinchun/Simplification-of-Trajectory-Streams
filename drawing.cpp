#include "drawing.h"
#include <QPainter>
#include <QApplication>

MultiViewer::MultiViewer(QWidget* parent) : QWidget(parent) {
    setWindowTitle("Incremental Viewer");
    resize(1000, 800);
}

void MultiViewer::addPolygon(const Polygon& poly) {
    polys_.push_back(poly);
    update();
}

void MultiViewer::addPoints(const std::vector<Point>& pts) {
    pointSets_.push_back(pts);
    update();
}

void MultiViewer::clearAll() {
    polys_.clear();
    pointSets_.clear();
    original_.clear();
    simplified_.clear();
    update();
}

void MultiViewer::addOriginalPoint(const Point& p) {
    original_.push_back(p);
    update();
}

void MultiViewer::addSimplifiedPoint(const Point& p) {
    simplified_.push_back(p);
    update();
}

void MultiViewer::addOriginalPoints(const std::vector<Point>& v) {
    original_.insert(original_.end(), v.begin(), v.end());
    update();
}

void MultiViewer::addSimplifiedPoints(const std::vector<Point>& v) {
    simplified_.insert(simplified_.end(), v.begin(), v.end());
    update();
}

void MultiViewer::clearOriginal() {
    original_.clear();
    update();
}

void MultiViewer::clearSimplified() {
    simplified_.clear();
    update();
}

void MultiViewer::compute_bbox(double& minx,double& miny,double& maxx,double& maxy) const {
    minx =  1e300; miny =  1e300;
    maxx = -1e300; maxy = -1e300;
    auto acc = [&](const Point& p){
        double x = CGAL::to_double(p.x());
        double y = CGAL::to_double(p.y());
        minx = std::min(minx,x); miny = std::min(miny,y);
        maxx = std::max(maxx,x); maxy = std::max(maxy,y);
    };
    for (auto& poly : polys_)
        for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
            acc(*v);
    for (auto& set : pointSets_)
        for (auto& p : set) acc(p);
    for (auto& p : original_)   acc(p);
    for (auto& p : simplified_) acc(p);
    if (!(minx <= maxx)) { minx = -1; maxx = 1; miny = -1; maxy = 1; }
}

void MultiViewer::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing, true);

    double minx,miny,maxx,maxy;
    compute_bbox(minx,miny,maxx,maxy);
    double dx = maxx - minx; if (dx == 0) dx = 1;
    double dy = maxy - miny; if (dy == 0) dy = 1;
    double margin = 30;
    int W = width(), H = height();
    double scale = std::min((W - 2*margin)/dx, (H - 2*margin)/dy);
    auto mapX = [&](double x){ return margin + (x - minx)*scale; };
    auto mapY = [&](double y){ return H - (margin + (y - miny)*scale); };

    p.fillRect(rect(), Qt::white);

    // arbitrary point sets
    p.setPen(Qt::NoPen);
    p.setBrush(Qt::black);
    for (auto& set : pointSets_) {
        for (auto& pt : set)
            p.drawEllipse(QPointF(mapX(CGAL::to_double(pt.x())),
                                  mapY(CGAL::to_double(pt.y()))), 3, 3);
    }

    // polygons
    const std::vector<QColor> colors { Qt::blue, Qt::red, Qt::darkGreen,
                                       Qt::magenta, Qt::darkCyan, Qt::darkYellow };
    int ci = 0;
    p.setBrush(Qt::NoBrush);
    for (auto& poly : polys_) {
        QPolygonF qp;
        for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
            qp << QPointF(mapX(CGAL::to_double(v->x())),
                          mapY(CGAL::to_double(v->y())));
        p.setPen(QPen(colors[ci % colors.size()], 2));
        p.drawPolygon(qp);
        ++ci;
    }

    // original trajectory (polyline + points)
    if (!original_.empty()) {
        p.setPen(QPen(Qt::darkGray, 2, Qt::SolidLine));
        for (size_t i = 1; i < original_.size(); ++i) {
            p.drawLine(QPointF(mapX(CGAL::to_double(original_[i-1].x())),
                               mapY(CGAL::to_double(original_[i-1].y()))),
                       QPointF(mapX(CGAL::to_double(original_[i].x())),
                               mapY(CGAL::to_double(original_[i].y()))));
        }
        p.setPen(Qt::NoPen);
        p.setBrush(Qt::darkGray);
        for (auto& pt : original_)
            p.drawEllipse(QPointF(mapX(CGAL::to_double(pt.x())),
                                  mapY(CGAL::to_double(pt.y()))), 2.5, 2.5);
    }

    // simplified trajectory (overlay, highlight)
    if (!simplified_.empty()) {
        p.setPen(QPen(Qt::red, 3, Qt::SolidLine));
        for (size_t i = 1; i < simplified_.size(); ++i) {
            p.drawLine(QPointF(mapX(CGAL::to_double(simplified_[i-1].x())),
                               mapY(CGAL::to_double(simplified_[i-1].y()))),
                       QPointF(mapX(CGAL::to_double(simplified_[i].x())),
                               mapY(CGAL::to_double(simplified_[i].y()))));
        }
        p.setPen(Qt::NoPen);
        p.setBrush(Qt::red);
        for (auto& pt : simplified_)
            p.drawEllipse(QPointF(mapX(CGAL::to_double(pt.x())),
                                  mapY(CGAL::to_double(pt.y()))), 3.5, 3.5);
    }
}

void viewer_process_events() {
    if (qApp) qApp->processEvents();
}