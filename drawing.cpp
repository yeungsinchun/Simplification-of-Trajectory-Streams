#include "drawing.h"
#include <QPainter>
#include <QApplication>
#include <QFontMetrics>
#include <QString>

// Bounding box declared/defined in simplify.cpp
extern const double BMIN;
extern const double BMAX;

MultiViewer::MultiViewer(QWidget* parent) : QWidget(parent) {
    setWindowTitle("Incremental Viewer");
    resize(1000, 800);
}

void MultiViewer::addPolygon(const Polygon& poly) {
    polys_.push_back(poly);
    update();
}

void MultiViewer::addSpecialPolygon(const Polygon& poly) {
    special_polys_.push_back(poly);
    update();
}

void MultiViewer::addPoints(const std::vector<Point>& pts) {
    pointSets_.push_back(pts);
    update();
}

void MultiViewer::addSpecialPoint(const Point& p) {
    special_points_.push_back(p);
    update();
}

void MultiViewer::clearAll() {
    polys_.clear();
    pointSets_.clear();
    original_.clear();
    simplified_.clear();
    update();
}

void MultiViewer::clearPolygons() {
    polys_.clear();
    update();
}

void MultiViewer::clearSpecials() {
    special_polys_.clear();
    special_points_.clear();
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

void MultiViewer::addCurve(const std::vector<Point>& v, const QColor& color, const QString& label) {
    curves_.push_back(Curve{v, color, label});
    update();
}

void MultiViewer::clearCurves() {
    curves_.clear();
    update();
}

void MultiViewer::setParameters(double delta, double epsilon) {
    delta_ = delta;
    epsilon_ = epsilon;
    update();
}

void MultiViewer::markP0(const Point& p) {
    marked_p0_ = p;
    update();
}

void MultiViewer::clearMarkedP0() {
    marked_p0_.reset();
    update();
}

void MultiViewer::markPi(const Point& p) {
    marked_pi_ = p;
    update();
}

void MultiViewer::clearMarkedPi() {
    marked_pi_.reset();
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

    // Map using the global bounding square [BMIN, BMAX]^2 so the
    // bounding box is always visible and drawn at the beginning.
    double minx = BMIN, miny = BMIN, maxx = BMAX, maxy = BMAX;
    double dx = maxx - minx; if (dx == 0) dx = 1;
    double dy = maxy - miny; if (dy == 0) dy = 1;
    double margin = 30;
    int W = width(), H = height();
    double scale = std::min((W - 2*margin)/dx, (H - 2*margin)/dy);
    auto mapX = [&](double x){ return margin + (x - minx)*scale; };
    auto mapY = [&](double y){ return H - (margin + (y - miny)*scale); };

    p.fillRect(rect(), Qt::white);

    // Draw the global bounding square outline first (no fill)
    {
        double x0 = mapX(BMIN);
        double x1 = mapX(BMAX);
        double y0 = mapY(BMIN);
        double y1 = mapY(BMAX);
        double left = std::min(x0, x1);
        double top  = std::min(y0, y1);
        double w    = std::abs(x1 - x0);
        double h    = std::abs(y1 - y0);
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(QColor(80,80,80), 2, Qt::SolidLine));
        p.drawRect(QRectF(left, top, w, h));
    }

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

    // special debug polygons (distinct style: dashed orange with translucent fill)
    if (!special_polys_.empty()) {
        QColor edge(255, 140, 0);         // dark orange
        QColor fill(255, 165, 0, 60);     // light orange, translucent
        p.setBrush(fill);
        p.setPen(QPen(edge, 2, Qt::DashLine));
        for (auto& poly : special_polys_) {
            QPolygonF qp;
            for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
                qp << QPointF(mapX(CGAL::to_double(v->x())),
                              mapY(CGAL::to_double(v->y())));
            p.drawPolygon(qp);
        }
        p.setBrush(Qt::NoBrush);
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

    // additional colored curves with labels
    for (const auto& c : curves_) {
        if (c.pts.empty()) continue;
        p.setPen(QPen(c.color, 2, Qt::SolidLine));
        for (size_t i = 1; i < c.pts.size(); ++i) {
            p.drawLine(QPointF(mapX(CGAL::to_double(c.pts[i-1].x())),
                               mapY(CGAL::to_double(c.pts[i-1].y()))),
                       QPointF(mapX(CGAL::to_double(c.pts[i].x())),
                               mapY(CGAL::to_double(c.pts[i].y()))));
        }
        p.setPen(Qt::NoPen);
        p.setBrush(c.color);
        for (auto& pt : c.pts)
            p.drawEllipse(QPointF(mapX(CGAL::to_double(pt.x())),
                                  mapY(CGAL::to_double(pt.y()))), 3.0, 3.0);
    }

    // marked p0 (draw late so it appears on top)
    if (marked_p0_) {
        const auto& mp = *marked_p0_;
        double mx = mapX(CGAL::to_double(mp.x()));
        double my = mapY(CGAL::to_double(mp.y()));
        p.setPen(Qt::NoPen);
        p.setBrush(QColor(50, 200, 50)); // bright green
        p.drawEllipse(QPointF(mx, my), 3.5, 3.5);
        p.setPen(Qt::black);
    }

    // marked pi (current point) - use a different color from p0
    if (marked_pi_) {
        const auto& mp = *marked_pi_;
        double mx = mapX(CGAL::to_double(mp.x()));
        double my = mapY(CGAL::to_double(mp.y()));
        p.setPen(Qt::NoPen);
        p.setBrush(QColor(240, 140, 0)); // orange
        p.drawEllipse(QPointF(mx, my), 3.5, 3.5);
        p.setPen(Qt::black);
    }

    // special debug points (distinct style: cyan circles with black outline)
    if (!special_points_.empty()) {
        p.setPen(QPen(Qt::black, 1));
        p.setBrush(QColor(0, 200, 255)); // cyan
        for (const auto& sp : special_points_) {
            p.drawEllipse(QPointF(mapX(CGAL::to_double(sp.x())),
                                  mapY(CGAL::to_double(sp.y()))), 4.0, 4.0);
        }
        p.setPen(Qt::black);
        p.setBrush(Qt::NoBrush);
    }

    // Heads-up display (top-right): counts, ratio, and parameters
    {
        const int W = width();
        const int pad = 8;
        const int marginTR = 10; // margin from top-right edges

        const size_t n_orig = original_.size();
        const size_t n_simp = simplified_.size();
        double ratio = (n_orig == 0) ? 0.0 : (100.0 * double(n_simp) / double(n_orig));

        QString line1 = QString("orig: %1  simp: %2  ratio: %3%")
                            .arg(qulonglong(n_orig))
                            .arg(qulonglong(n_simp))
                            .arg(ratio, 0, 'f', 1);
        QString line2 = QString("delta=%1  epsilon=%2")
                            .arg(std::isfinite(delta_) ? QString::number(delta_, 'g', 4) : QString("-"))
                            .arg(std::isfinite(epsilon_) ? QString::number(epsilon_, 'g', 4) : QString("-"));

        QFontMetrics fm(p.font());
        int textW = std::max(fm.horizontalAdvance(line1), fm.horizontalAdvance(line2));
        int textH = fm.height() * 2 + 4; // 2 lines + small spacing

        QRect rectBg(W - marginTR - (textW + 2*pad), marginTR, textW + 2*pad, textH + 2*pad);
        p.setPen(Qt::NoPen);
        p.setBrush(QColor(255, 255, 255, 210));
        p.drawRect(rectBg);

        p.setPen(Qt::black);
        int tx = rectBg.left() + pad;
        int ty = rectBg.top() + pad + fm.ascent();
        p.drawText(tx, ty, line1);
        ty += fm.height();
        p.drawText(tx, ty, line2);
    }

    // Legend (top-left): labels and colors for curves (with point counts)
    {
        const int marginTL = 50;
        const int swatch = 10;
        const int pad = 6;
        int x0 = marginTL, y0 = marginTL;

        QFontMetrics fm(p.font());
    // Compose legend entries: original, simplified, then added curves, each with point count
        struct Entry { QColor color; QString text; };
        std::vector<Entry> entries;
        if (!original_.empty()) entries.push_back({Qt::darkGray, QString("original (%1)").arg(qulonglong(original_.size()))});
        if (!simplified_.empty()) entries.push_back({Qt::red, QString("simplified (%1)").arg(qulonglong(simplified_.size()))});
        for (const auto& c : curves_) entries.push_back({c.color, QString("%1 (%2)").arg(c.label).arg(qulonglong(c.pts.size()))});
    if (!special_polys_.empty()) entries.push_back({QColor(255, 140, 0), QString("special polys (%1)").arg(qulonglong(special_polys_.size()))});
    if (!special_points_.empty()) entries.push_back({QColor(0, 200, 255), QString("special points (%1)").arg(qulonglong(special_points_.size()))});

        for (const auto& e : entries) {
            // draw color box
            p.setPen(Qt::NoPen);
            p.setBrush(e.color);
            p.drawRect(x0, y0, swatch, swatch);
            // draw label
            p.setPen(Qt::black);
            p.drawText(x0 + swatch + pad, y0 + fm.ascent(), e.text);
            y0 += std::max(swatch, fm.height()) + 4;
        }
    }
}

void viewer_process_events() {
    if (qApp) qApp->processEvents();
}