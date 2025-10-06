#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/number_utils.h>

#include <QApplication>
#include <QPainter>
#include <QPolygonF>
#include <QWidget>

using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;

std::vector<Point_2> generate_stream(int N, double min_coord = -1e6,
                                     double max_coord = 1e6) {
  std::vector<Point_2> out;
  out.reserve(std::max(0, N));

  std::mt19937_64 rng(std::random_device{}());
  std::uniform_real_distribution<double> pos_dist(min_coord, max_coord);
  std::uniform_real_distribution<double> heading_dist(-M_PI, M_PI);
  std::normal_distribution<double> step_dist(5.0 * 1e4, 1.5 * 1e4);
  std::normal_distribution<double> steer_noise(0.0, 0.02);
  std::bernoulli_distribution turn_event(0.2);
  std::uniform_real_distribution<double> turn_large(-M_PI / 4, M_PI / 4);
  std::normal_distribution<double> gps_noise(0.0, 3.0);

  double x = pos_dist(rng);
  double y = pos_dist(rng);
  double heading = heading_dist(rng);

  for (int i = 0; i < N; ++i) {
    double step = step_dist(rng);
    if (turn_event(rng))
      heading += turn_large(rng);
    heading += steer_noise(rng);

    x += step * std::cos(heading);
    y += step * std::sin(heading);

    if (x < min_coord) {
      x = min_coord + (min_coord - x);
      heading += 0.5;
    }
    if (x > max_coord) {
      x = max_coord - (x - max_coord);
      heading -= 0.5;
    }
    if (y < min_coord) {
      y = min_coord + (min_coord - y);
      heading += 0.5;
    }
    if (y > max_coord) {
      y = max_coord - (y - max_coord);
      heading -= 0.5;
    }

    out.emplace_back(x + gps_noise(rng), y + gps_noise(rng));
  }

  return out;
}

class StreamWidget : public QWidget {
public:
  StreamWidget(std::vector<Point_2> pts, QWidget *parent = nullptr)
      : QWidget(parent), pts_(std::move(pts)) {
    setWindowTitle("Generated GPS-like stream");
    resize(900, 700);
    compute_bbox();
  }

protected:
  void paintEvent(QPaintEvent *) override {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing, true);

    const int W = width();
    const int H = height();
    const double margin = 20.0;

    double dx = maxx_ - minx_;
    double dy = maxy_ - miny_;
    if (dx == 0)
      dx = 1.0;
    if (dy == 0)
      dy = 1.0;
    double scale = std::min((W - 2 * margin) / dx, (H - 2 * margin) / dy);

    auto mapx = [&](double x) { return margin + (x - minx_) * scale; };
    auto mapy = [&](double y) { return H - (margin + (y - miny_) * scale); };

    // draw path
    p.setPen(QPen(Qt::blue, 1));
    QPolygonF poly;
    poly.reserve(pts_.size());
    for (const auto &pt : pts_) {
      double x = CGAL::to_double(pt.x());
      double y = CGAL::to_double(pt.y());
      poly << QPointF(mapx(x), mapy(y));
    }
    if (!poly.isEmpty())
      p.drawPolyline(poly);

    // draw points
    p.setBrush(Qt::black);
    p.setPen(Qt::NoPen);
    for (const auto &pt : pts_) {
      double x = CGAL::to_double(pt.x());
      double y = CGAL::to_double(pt.y());
      QPointF center(mapx(x), mapy(y));
      p.drawEllipse(center, 3.0, 3.0);
    }
    // draw points as filled dots with thin outline (more visible)
    const double base_radius = 4.0;
    for (const auto &pt : pts_) {
      double x = CGAL::to_double(pt.x());
      double y = CGAL::to_double(pt.y());
      QPointF center(mapx(x), mapy(y));

      // filled dot
      p.setPen(Qt::NoPen);
      p.setBrush(Qt::red);
      p.drawEllipse(center, base_radius, base_radius);

      // thin black outline
      p.setBrush(Qt::NoBrush);
      p.setPen(QPen(Qt::black, 0.8));
      p.drawEllipse(center, base_radius, base_radius);
    }
  }

private:
  void compute_bbox() {
    if (pts_.empty()) {
      minx_ = miny_ = -1.0;
      maxx_ = maxy_ = 1.0;
      return;
    }
    minx_ = miny_ = 1e300;
    maxx_ = maxy_ = -1e300;
    for (auto &pt : pts_) {
      double x = CGAL::to_double(pt.x());
      double y = CGAL::to_double(pt.y());
      minx_ = std::min(minx_, x);
      miny_ = std::min(miny_, y);
      maxx_ = std::max(maxx_, x);
      maxy_ = std::max(maxy_, y);
    }
  }

  std::vector<Point_2> pts_;
  double minx_, miny_, maxx_, maxy_;
};

int main(int argc, char **argv) {
  int N;
  if (!(std::cin >> N))
    return 0;

  std::vector<Point_2> stream;
  if (N > 0) {
    stream = generate_stream(N);
  } else {
    stream = generate_stream(1000); // default
  }

  // write stream to data/input.txt in the format:
  // <N>
  // x y
  // x y
  // ...
  try {
    std::ofstream ofs("../data/input2.txt");
    ofs << std::setprecision(15);
    ofs << stream.size() << "\n";
    for (const auto &pt : stream) {
      std::cout << CGAL::to_double(pt.x()) << ' ' << CGAL::to_double(pt.y())
                << '\n';
      ofs << CGAL::to_double(pt.x()) << " " << CGAL::to_double(pt.y()) << "\n";
    }
    std::cout << "Stream has been written to ../data/input2.txt\n";
  } catch (...) {
    std::cout << "Not ok";
  }

  QApplication app(argc, argv);
  StreamWidget w(std::move(stream));
  w.show();
  return app.exec();
}