#include "drawing.h"
#include <QPainter>
#include <QApplication>
#include <QFontMetrics>
#include <QString>

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

void MultiViewer::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing, true);

    // Fixed world-space: always [-10000, 10000]^2
    const double WORLD_MIN = -10000.0;
    const double WORLD_MAX =  10000.0;
    const double WORLD_SIDE = WORLD_MAX - WORLD_MIN; // 20000
    const double margin = 30.0;

    int W = width(), H = height();
    double availW = W - 2.0 * margin;
    double availH = H - 2.0 * margin;
    double scale = std::min(availW / WORLD_SIDE, availH / WORLD_SIDE);

    // Center horizontally and vertically: any slack is split equally
    double x_off = margin + (availW - WORLD_SIDE * scale) * 0.5;
    double y_off = margin + (availH - WORLD_SIDE * scale) * 0.5;

    auto mapX = [&](double x) { return x_off + (x - WORLD_MIN) * scale; };
    auto mapY = [&](double y) { return H - (y_off + (y - WORLD_MIN) * scale); };

    // Viewport rect in screen pixels
    double vp_left   = mapX(WORLD_MIN);
    double vp_right  = mapX(WORLD_MAX);
    double vp_top    = mapY(WORLD_MAX);
    double vp_bottom = mapY(WORLD_MIN);

    p.fillRect(rect(), Qt::white);

    // Draw the fixed bounding box for [-10000, 10000]^2
    {
        double bx_lo = mapX(WORLD_MIN);
        double bx_hi = mapX(WORLD_MAX);
        double by_lo = mapY(WORLD_MIN);
        double by_hi = mapY(WORLD_MAX);
        double box_left  = std::min(bx_lo, bx_hi);
        double box_right = std::max(bx_lo, bx_hi);
        double box_top   = std::min(by_lo, by_hi);
        double box_bot   = std::max(by_lo, by_hi);
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
                // Clamp to the actual viewport so labels never go off-screen.
                if (ly - th < vp_top)   ly = pf.y() + 5 + th;
                if (lx + tw > vp_right) lx = pf.x() - 5 - tw;
                if (lx < vp_left)       lx = vp_left;
                if (ly > vp_bottom)     ly = pf.y() - 5 - th;
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

    // Heads-up display (top-right, inside viewport): ratio, e, d
    {
        const int pad = 6;
        const int marginTR = 10; // margin from viewport right/top edges

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

        // Position inside the viewport (top-right corner)
        double hud_x = vp_right - marginTR - (textW + 2*pad);
        double hud_y = vp_top + marginTR;
        QRect rectBg(static_cast<int>(hud_x), static_cast<int>(hud_y),
                     textW + 2*pad, textH + 2*pad);
        p.setPen(Qt::NoPen);
        p.setBrush(QColor(255, 255, 255, 210));
        p.drawRect(rectBg);

        p.setPen(Qt::black);
        int tx = rectBg.left() + pad;
        int ty = rectBg.top() + pad + fm.ascent();
        p.drawText(tx, ty, line);
    }

    // Legend (top-left, inside the viewport): show number of points for each entry.
    // Each entry records a hit-test rectangle used by mousePressEvent to
    // toggle visibility on click. The first two slots (indices 0 and 1)
    // refer to the original and simplified trajectories respectively;
    // the rest (2..N+1) index into curves_.
    {
        const int swatch = 10;
        const int pad = 6;
        // Position the legend just inside the top-left corner of the viewport
        double x0 = vp_left + swatch;
        double y0 = vp_top + swatch;

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
            int ix0 = static_cast<int>(x0);
            int iy0 = static_cast<int>(y0);
            int textX = ix0 + swatch + pad;
            int textW = fm.horizontalAdvance(e.text);
            QRect rowRect(ix0, iy0, (textX - ix0) + textW, rowH);
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
            p.drawRect(ix0, iy0, swatch, swatch);
            // draw label
            p.setPen(textColor);
            p.drawText(textX, iy0 + fm.ascent(), e.text);
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