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

void MultiViewer::addPolygon(const Polygon& poly, const QColor& color) {
    polys_.push_back({poly, color});
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

void MultiViewer::setShowLabels(bool show) {
    showLabels_ = show;
    update();
}

void MultiViewer::setDataBBox(double minx, double miny, double maxx, double maxy) {
    // NaN coordinates mean "unfreeze" (revert to dynamic per-paint view).
    if (!std::isfinite(minx) || !std::isfinite(miny) ||
        !std::isfinite(maxx) || !std::isfinite(maxy)) {
        view_frozen_ = false;
        update();
        return;
    }
    frozen_minx_ = minx;
    frozen_miny_ = miny;
    frozen_maxx_ = maxx;
    frozen_maxy_ = maxy;
    view_frozen_ = true;
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
    for (auto& cp : polys_)
        for (auto v = cp.poly.vertices_begin(); v != cp.poly.vertices_end(); ++v)
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

    // Effective viewport: anchored on the data bounding box with a small
    // world-space padding. The algorithmic grid [BMIN, BMAX]^2 is just a
    // reference frame drawn on top; it should not constrain the view.
    //
    // If the view was frozen via setDataBBox(), use that bbox and don't
    // re-fit on every paint (so the canvas doesn't re-zoom as more points
    // stream in). Otherwise compute the bbox from the current data.
    double dminx, dminy, dmaxx, dmaxy;
    if (view_frozen_) {
        dminx = frozen_minx_; dminy = frozen_miny_;
        dmaxx = frozen_maxx_; dmaxy = frozen_maxy_;
    } else {
        compute_bbox(dminx, dminy, dmaxx, dmaxy);
    }
    double dx_data = std::max(dmaxx - dminx, 1.0);
    double dy_data = std::max(dmaxy - dminy, 1.0);
    // The BBOX rectangle drawn on screen is the *unpadded* data extent
    // [dminx, dmaxx] x [dminy, dmaxy] (so its edges sit on the actual
    // extreme data points). The world view we project to the canvas is
    // a square that contains this data extent plus 8% padding on every
    // side, so the display is always a square regardless of the data
    // aspect ratio.
    double pad = std::max(dx_data, dy_data) * pad_frac_;
    double side = std::max(dx_data, dy_data) + 2.0 * pad;
    double cx = 0.5 * (dminx + dmaxx);
    double cy = 0.5 * (dminy + dmaxy);
    double minx = cx - 0.5 * side;
    double maxx = cx + 0.5 * side;
    double miny = cy - 0.5 * side;
    double maxy = cy + 0.5 * side;
    double dx = maxx - minx; if (dx == 0) dx = 1;
    double dy = maxy - miny; if (dy == 0) dy = 1;
    double margin = 30;
    int W = width(), H = height();
    double availW = W - 2*margin;
    double availH = H - 2*margin;
    double scale  = std::min(availW/dx, availH/dy);
    // Center the data within the canvas: any slack on the looser axis is
    // split equally so the data is centered both horizontally and vertically.
    // vshift_ (pixels) shifts the data upward on screen — positive vshift_
    // moves the canvas up, giving more room at the bottom. Clamp vshift_
    // so the data extent can never be pushed off the top of the window:
    // we need mapY(dminy) >= 0, i.e. y_off <= H, i.e. vshift_ <= (H+margin)/2 - dy*scale/2.
    double vshift_eff = vshift_;
    {
        double y_off_no_shift = margin + (availH - dy * scale) * 0.5;
        double max_vshift = H - y_off_no_shift; // mapY(dminy) = H - y_off = H - y_off_no_shift - vshift; require >= 0
        if (vshift_eff > max_vshift) vshift_eff = max_vshift;
        if (vshift_eff < -y_off_no_shift) vshift_eff = -y_off_no_shift; // mapY(dmaxy) = H - y_off - dy*scale; require <= H
    }
    double x_off = margin + (availW - dx * scale) * 0.5;
    double y_off = margin + (availH - dy * scale) * 0.5 + vshift_eff;
    auto mapX = [&](double x){ return x_off + (x - minx) * scale; };
    auto mapY = [&](double y){ return H - (y_off + (y - miny) * scale); };

    p.fillRect(rect(), Qt::white);

    // Draw the data bounding box (no fill) first. The box outlines the
    // *unpadded* data extent [dminx, dmaxx] x [dminy, dmaxy], with a
    // tiny visual inset (1% of the data side) so points and the index
    // labels we draw on top of them don't sit directly on the box
    // stroke. The world view used for projection is a slightly larger
    // square (centered on the data center, side = max(dx,dy) + 2*pad),
    // which keeps the display square and gives labels room to render
    // just outside the box.
    double box_left, box_right, box_top, box_bot;
    {
        double inset = std::max(dx_data, dy_data) * 0.01;
        double bx0 = dminx - inset, bx1 = dmaxx + inset;
        double by0 = dminy - inset, by1 = dmaxy + inset;
        double x_lo = mapX(bx0);
        double x_hi = mapX(bx1);
        double y_lo = mapY(by0);
        double y_hi = mapY(by1);
        box_left = std::min(x_lo, x_hi);
        box_right = std::max(x_lo, x_hi);
        box_top  = std::min(y_lo, y_hi);
        box_bot  = std::max(y_lo, y_hi);
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(QColor(80,80,80), 2, Qt::SolidLine));
        p.drawRect(QRectF(box_left, box_top, box_right - box_left, box_bot - box_top));
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
    p.setBrush(Qt::NoBrush);
    for (auto& cp : polys_) {
        QPolygonF qp;
        for (auto v = cp.poly.vertices_begin(); v != cp.poly.vertices_end(); ++v)
            qp << QPointF(mapX(CGAL::to_double(v->x())),
                          mapY(CGAL::to_double(v->y())));
        p.setPen(QPen(cp.color, 2));
        p.drawPolygon(qp);
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
    if (!original_.empty() && original_visible_) {
        p.setPen(QPen(Qt::darkGray, 2, Qt::SolidLine));
        for (size_t i = 1; i < original_.size(); ++i) {
            p.drawLine(QPointF(mapX(CGAL::to_double(original_[i-1].x())),
                               mapY(CGAL::to_double(original_[i-1].y()))),
                       QPointF(mapX(CGAL::to_double(original_[i].x())),
                               mapY(CGAL::to_double(original_[i].y()))));
        }
        p.setPen(Qt::NoPen);
        p.setBrush(Qt::darkGray);
        for (size_t i = 0; i < original_.size(); ++i) {
            const auto& pt = original_[i];
            QPointF pf(mapX(CGAL::to_double(pt.x())), mapY(CGAL::to_double(pt.y())));
            p.drawEllipse(pf, 2.5, 2.5);
            if (showLabels_) {
                p.setPen(Qt::black);
                QString txt = QString::number(i);
                QFontMetrics fm(p.font());
                int tw = fm.horizontalAdvance(txt);
                int th = fm.height();
                // Default offset: 5px right and 5px above the point.
                double lx = pf.x() + 5;
                double ly = pf.y() - 5;
                // If the label would cross the top edge of the bbox, push
                // it below the point instead.
                if (ly - th < box_top) ly = pf.y() + 5 + th;
                // If the label would cross the right edge, push it left of
                // the point instead.
                if (lx + tw > box_right) lx = pf.x() - 5 - tw;
                // If the label would cross the left edge, push it right
                // beyond the point (with a little extra room).
                if (lx < box_left) lx = pf.x() + 5;
                // If the label would cross the bottom edge, push it above
                // the point (with a little extra room).
                if (ly > box_bot) ly = pf.y() - 5 - th;
                p.drawText(QPointF(lx, ly), txt);
                p.setPen(Qt::NoPen);
            }
        }
    }

    // simplified trajectory (overlay, highlight)
    if (!simplified_.empty() && simplified_visible_) {
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
        if (!c.visible || c.pts.empty()) continue;
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

    // Heads-up display (top-right): single line: ratio, e, d
    {
        const int W = width();
        const int pad = 8;
        const int marginTR = 10; // margin from top-right edges

        const size_t n_orig = original_.size();
        const size_t n_simp = simplified_.size();
        QString ratioStr = (n_orig == 0) ? QString("ratio=-")
                                         : QString("ratio=%1% ").arg(100.0 * double(n_simp) / double(n_orig), 0, 'f', 1);
        QString eStr = QString("e=%1 ").arg(std::isfinite(epsilon_) ? QString::number(epsilon_, 'g', 4) : QString("-"));
        QString dStr = QString("d=%1").arg(std::isfinite(delta_) ? QString::number(delta_, 'g', 4) : QString("-"));
        QString line = ratioStr + eStr + dStr;

        QFontMetrics fm(p.font());
        int textW = fm.horizontalAdvance(line);
        int textH = fm.height();

        QRect rectBg(W - marginTR - (textW + 2*pad), marginTR, textW + 2*pad, textH + 2*pad);
        p.setPen(Qt::NoPen);
        p.setBrush(QColor(255, 255, 255, 210));
        p.drawRect(rectBg);

        p.setPen(Qt::black);
        int tx = rectBg.left() + pad;
        int ty = rectBg.top() + pad + fm.ascent();
        p.drawText(tx, ty, line);
    }

    // Legend (top-left): show only number of points for each entry.
    // Each entry records a hit-test rectangle used by mousePressEvent to
    // toggle visibility on click. The first two slots (indices 0 and 1)
    // refer to the original and simplified trajectories respectively;
    // the rest (2..N+1) index into curves_.
    {
        const int marginTL = 50;
        const int swatch = 10;
        const int pad = 6;
        int x0 = marginTL, y0 = marginTL;

        QFontMetrics fm(p.font());
    // Compose legend entries: original, simplified, then added curves, each with label and point count
        struct Entry { QColor color; QString text; bool hidden = false; };
        std::vector<Entry> entries;
        if (!original_.empty()) {
            entries.push_back({Qt::darkGray, QString("original (%1)").arg(qulonglong(original_.size())), !original_visible_});
        }
        if (!simplified_.empty()) {
            entries.push_back({Qt::red, QString("simplified (%1)").arg(qulonglong(simplified_.size())), !simplified_visible_});
        }
        for (const auto& c : curves_) {
            entries.push_back({c.color, QString("%1 (%2)").arg(c.label).arg(qulonglong(c.pts.size())), !c.visible});
        }
        if (!special_polys_.empty()) entries.push_back({QColor(255, 140, 0), QString("special polys (%1)").arg(qulonglong(special_polys_.size()))});
        if (!special_points_.empty()) entries.push_back({QColor(0, 200, 255), QString("special points (%1)").arg(qulonglong(special_points_.size()))});

        legend_rects_.clear();
        legend_rects_.reserve(entries.size());

        for (const auto& e : entries) {
            int rowH = std::max(swatch, fm.height());
            int textX = x0 + swatch + pad;
            int textW = fm.horizontalAdvance(e.text);
            QRect rowRect(x0, y0, (textX - x0) + textW, rowH);
            legend_rects_.push_back(rowRect);

            // dim the entry when its underlying series is hidden
            QColor swatchColor = e.color;
            QColor textColor   = Qt::black;
            if (e.hidden) {
                swatchColor.setAlpha(90);
                textColor = QColor(150, 150, 150);
            }

            // draw color box
            p.setPen(Qt::NoPen);
            p.setBrush(swatchColor);
            p.drawRect(x0, y0, swatch, swatch);
            // draw label
            p.setPen(textColor);
            p.drawText(textX, y0 + fm.ascent(), e.text);
            y0 += rowH + 4;
        }
    }
}

void MultiViewer::mousePressEvent(QMouseEvent* ev) {
    if (ev->button() != Qt::LeftButton) {
        QWidget::mousePressEvent(ev);
        return;
    }
    const QPoint pos = ev->pos();
    if (legend_rects_.empty()) return;
    // Walk entries in reverse (top-of-stack first) so overlapping rows resolve to the topmost.
    for (int i = static_cast<int>(legend_rects_.size()) - 1; i >= 0; --i) {
        if (!legend_rects_[i].contains(pos)) continue;
        // Map legend index back to underlying object.
        // First two slots are original/simplified (toggle visibility).
        size_t idx = 0;
        if (!original_.empty()) {
            if (i == static_cast<int>(idx)) {
                original_visible_ = !original_visible_;
                update();
                return;
            }
            ++idx;
        }
        if (!simplified_.empty()) {
            if (i == static_cast<int>(idx)) {
                simplified_visible_ = !simplified_visible_;
                update();
                return;
            }
            ++idx;
        }
        if (idx < curves_.size() &&
            i >= static_cast<int>(idx) &&
            (i - static_cast<int>(idx)) < static_cast<int>(curves_.size())) {
            size_t ci = static_cast<size_t>(i) - idx;
            if (ci < curves_.size()) {
                curves_[ci].visible = !curves_[ci].visible;
                update();
            }
            return;
        }
        // special polys / points entries are non-toggleable
        return;
    }
}

void viewer_process_events() {
    if (qApp) qApp->processEvents();
}