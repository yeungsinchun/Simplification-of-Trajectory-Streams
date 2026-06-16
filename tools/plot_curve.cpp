// plot_curve: visualize all curves in data/<id>/ side by side.
//
// Usage: ./plot_curve <id>
//
// Loads data/<id>/original.txt (if present) and every other *.txt file in
// data/<id> as a separately-colored curve. Click on a legend entry to
// toggle that series on/off; hidden entries are dimmed in the legend.
//
// Paths are resolved from the executable location when possible, falling
// back to the current working directory. Repository root is identified by
// the presence of a `data/` directory.

#include "drawing.h"
#include <QApplication>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mach-o/dyld.h>
#include <string>
#include <vector>

// Bounding square expected by drawing.cpp (matches simplify.cpp).
extern const double BMIN = -10000.0;
extern const double BMAX =  10000.0;

using P = Point;

static bool parse_n_pairs(const std::string& path, std::vector<P>& out) {
    std::ifstream in(path);
    if (!in) return false;
    size_t N = 0;
    if (!(in >> N)) return false;
    out.clear();
    out.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        double x, y;
        if (!(in >> x >> y)) return false;
        out.emplace_back(x, y);
    }
    return true;
}

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

static QString label_from_filename(const std::string& filename) {
    std::string base = filename;
    auto pos = base.rfind('/');
    if (pos != std::string::npos) base = base.substr(pos + 1);
    pos = base.rfind('\\');
    if (pos != std::string::npos) base = base.substr(pos + 1);
    base = to_lower(base);

    if (base == "original.txt") return "Original";
    if (base == "simplified.txt" || base == "our_simplified.txt" || base == "simplify_sh.txt")
        return "Simplified";

    // Strip a trailing ".txt" before applying the more specific rules so
    // `simplify_against_dp_0.5_355.2885.txt` becomes
    // `simplify_against_dp_0.5_355.2885` (no extension).
    auto strip_ext = [](std::string& s) {
        auto dot = s.rfind('.');
        if (dot != std::string::npos) s = s.substr(0, dot);
    };
    strip_ext(base);

    if (base == "dp" || base == "dp_simplified")    return "DP";
    if (base == "dots" || base == "dots_simplified") return "DOTS";
    if (base == "squish" || base == "squish_simplified") return "SQUISH";

    // Generic *_simplified fallback
    const std::string suffix = "_simplified";
    if (base.size() > suffix.size() &&
        base.substr(base.size() - suffix.size()) == suffix) {
        std::string algo = base.substr(0, base.size() - suffix.size());
        if (algo == "operab") return "OPERBA";
        std::string label;
        for (char c : algo) label += std::toupper(c);
        return QString::fromStdString(label);
    }

    return QString::fromStdString(base);
}

namespace fs = std::filesystem;

static fs::path find_repo_root() {
    fs::path repo_root;
    try {
        char self_path[4096] = {0};
        uint32_t size = sizeof(self_path);
        if (_NSGetExecutablePath(self_path, &size) == 0) {
            fs::path exe = fs::canonical(self_path);
            fs::path dir = exe.parent_path();
            for (int i = 0; i < 6 && !dir.empty(); ++i) {
                if (fs::exists(dir / "data") && fs::is_directory(dir / "data")) {
                    return dir;
                }
                dir = dir.parent_path();
            }
        }
    } catch (...) {}
    fs::path cwd = fs::current_path();
    if (fs::exists(cwd / "data") && fs::is_directory(cwd / "data")) {
        return cwd;
    }
    return cwd;
}

int main(int argc, char** argv) {
    if (argc < 2 || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
        std::cerr << "Usage: " << argv[0] << " <id>\n"
                  << "Plots every .txt file in data/<id>/. Click a legend\n"
                  << "entry to toggle that series on/off.\n";
        return argc < 2 ? 1 : 0;
    }

    int id = 0;
    try {
        id = std::stoi(argv[1]);
    } catch (...) {
        std::cerr << "Invalid id: " << argv[1] << "\n";
        return 1;
    }

    fs::path repo_root = find_repo_root();
    fs::path data_dir = repo_root / "data" / std::to_string(id);

    if (!fs::exists(data_dir) || !fs::is_directory(data_dir)) {
        std::cerr << "Directory not found: " << data_dir << "\n";
        return 1;
    }

    // Gather curves from every .txt file in the directory.
    std::vector<P> orig;
    std::vector<std::pair<QString, std::vector<P>>> curves;

    bool have_original = false;
    for (const auto& entry : fs::directory_iterator(data_dir)) {
        if (!entry.is_regular_file()) continue;
        const auto& path = entry.path();
        if (path.extension() != ".txt") continue;
        std::string name = path.filename().string();
        std::vector<P> pts;
        if (!parse_n_pairs(path.string(), pts)) continue;
        QString label = label_from_filename(name);
        if (to_lower(name) == "original.txt") {
            orig = std::move(pts);
            have_original = true;
        } else {
            curves.emplace_back(label, std::move(pts));
        }
    }

    if (!have_original) {
        std::cerr << "Note: no original.txt in " << data_dir << "\n";
    }
    if (curves.empty() && !have_original) {
        std::cerr << "No .txt curves found in " << data_dir << "\n";
        return 1;
    }

    // Stable, predictable legend order: special ones first, then alphabetical.
    auto rank = [](const QString& s) -> int {
        if (s == "Original")   return 0;
        if (s == "Simplified") return 1;
        if (s == "DP")         return 2;
        if (s == "DOTS")       return 3;
        if (s == "SQUISH")     return 4;
        if (s == "OPERB")      return 5;
        if (s == "OPERBA")     return 6;
        if (s == "FBQS")       return 7;
        return 100;
    };
    std::sort(curves.begin(), curves.end(),
              [&](auto& a, auto& b) {
                  int ra = rank(a.first), rb = rank(b.first);
                  if (ra != rb) return ra < rb;
                  return a.first < b.first;
              });

    QApplication app(argc, argv);
    MultiViewer viewer;
    viewer.setWindowTitle(QString("Trajectory curves for id %1").arg(id));
    if (have_original) viewer.addOriginalPoints(orig);

    auto color_for = [](const QString& label) {
        if (label == "Simplified") return QColor(255, 0, 0);
        if (label == "DP")         return QColor(255, 127, 0);
        if (label == "DOTS")       return QColor(166, 86, 40);
        if (label == "SQUISH")     return QColor(147, 112, 219);
        if (label == "OPERB")      return QColor(55, 126, 184);
        if (label == "OPERBA")     return QColor(77, 175, 74);
        if (label == "FBQS")       return QColor(152, 78, 163);
        return QColor(120, 120, 120);
    };

    for (auto& [label, pts] : curves) {
        if (label == "Simplified") {
            viewer.addSimplifiedPoints(pts);
        } else {
            viewer.addCurve(pts, color_for(label), label);
        }
    }

    viewer.resize(1000, 800);
    viewer.show();
    return app.exec();
}
