#include "drawing.h"
#include <QApplication>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <cctype>

// Provide the global bounding square expected by drawing.cpp
extern const double BMIN = -10000.0;
extern const double BMAX =  10000.0;

using P = Point; // CGAL kernel point

static bool parse_two_line_xy(const std::string& path, std::vector<P>& out) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line1, line2;
    if (!std::getline(in, line1)) return false;
    // Skip possible empty/comment lines
    while (line1.size() && std::all_of(line1.begin(), line1.end(), [](unsigned char c){return std::isspace(c);} )) {
        if (!std::getline(in, line1)) return false;
    }
    if (!std::getline(in, line2)) return false;
    std::istringstream xs(line1), ys(line2);
    std::vector<double> vx, vy;
    for (double v; xs >> v; ) vx.push_back(v);
    for (double v; ys >> v; ) vy.push_back(v);
    if (vx.empty() || vx.size() != vy.size()) return false;
    out.clear(); out.reserve(vx.size());
    for (size_t i = 0; i < vx.size(); ++i) out.emplace_back(vx[i], vy[i]);
    return true;
}

static bool parse_n_pairs(const std::string& path, std::vector<P>& out) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    // Attempt to detect an optional leading N
    size_t expected = 0;
    if (std::getline(in, line)) {
        std::istringstream ss(line);
        size_t ntmp; double tx, ty;
        if ((ss >> ntmp) && !(ss >> tx)) {
            expected = ntmp;
        } else {
            // The first line already contains coordinates; process it below
            in.clear();
            in.seekg(0);
        }
    }
    out.clear();
    double x, y;
    while (in >> x >> y) out.emplace_back(x, y);
    if (expected && out.size() != expected) {
        // Not a hard error; accept what was parsed
    }
    return !out.empty();
}

static bool load_curve(const std::string& path, std::vector<P>& out) {
    // Try two-line format first, then fallback
    if (parse_two_line_xy(path, out)) return true;
    out.clear();
    if (parse_n_pairs(path, out)) return true;
    return false;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <original.txt> <simplified.txt>\n";
        std::cerr << "Example: " << argv[0] << " ../data/taxi_simplified/10/original.txt ../data/taxi_simplified/10/simplified.txt\n";
        return 1;
    }
    std::string origPath = argv[1];
    std::string simpPath = argv[2];

    std::vector<P> orig, simp;
    if (!load_curve(origPath, orig)) {
        std::cerr << "Failed to read curve from " << origPath << "\n";
        return 2;
    }
    if (!load_curve(simpPath, simp)) {
        std::cerr << "Failed to read curve from " << simpPath << "\n";
        return 3;
    }

    QApplication app(argc, argv);
    MultiViewer viewer;
    viewer.setWindowTitle("Curve Viewer: original (gray) vs simplified (red)");
    viewer.addOriginalPoints(orig);
    viewer.addSimplifiedPoints(simp);
    viewer.resize(1000, 800);
    viewer.show();
    return app.exec();
}
