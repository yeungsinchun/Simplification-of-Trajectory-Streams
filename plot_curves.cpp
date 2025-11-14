#include "drawing.h"
#include <QApplication>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <cctype>
#include <filesystem>
#include <algorithm>

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

static std::string to_lower(std::string s) { std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); }); return s; }

static QString label_from_filename(const std::string& filename) {
    // filename like: dp_simplified.txt, operb_simplified.txt, simplified.txt, original.txt
    std::string base = filename;
    auto pos = base.rfind('/');
    if (pos != std::string::npos) base = base.substr(pos+1);
    pos = base.rfind('\\');
    if (pos != std::string::npos) base = base.substr(pos+1);
    base = to_lower(base);
    if (base == "original.txt") return "Original";
    if (base == "simplified.txt") return "Simplified";
    if (base.find("dp_") == 0) return "DP";
    // Accept both spellings: operba_* and operab_* map to OPERBA
    if (base.find("operba_") == 0 || base.find("operab_") == 0) return "OPERBA";
    if (base.find("operb_") == 0) return "OPERB";
    if (base.find("fbqs_") == 0) return "FBQS";
    // Fallback: strip extension
    auto dot = base.rfind('.');
    if (dot != std::string::npos) base = base.substr(0, dot);
    return QString::fromStdString(base);
}

int main(int argc, char** argv) {
    // Expect: ./plot_curves <id>
    if (argc < 2 || (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))) {
        std::cerr << "Usage: " << argv[0] << " <id> [--orig] [--all | -dp -operb -operba -fbqs -simplify]\n";
        std::cerr << "Examples:\n  " << argv[0] << " 1 --all\n  " << argv[0] << " 1 -dp -operb\n  " << argv[0] << " 47 --orig\n";
        std::cerr << "Loads data/taxi/<id>.txt as original (if needed) and curves from data/taxi_simplified/<id>/*.txt (paths resolved from executable or cwd).\n";
        std::cerr << "When --orig is provided, only the original curve is displayed.\n";
        return argc < 2 ? 1 : 0;
    }

    int id = 0;
    try { id = std::stoi(argv[1]); } catch (...) {
        std::cerr << "Invalid id: " << argv[1] << "\n"; return 1;
    }

    // Parse selection flags (optional)
    bool selAny = false, selAll = false;
    bool selDP = false, selOPERB = false, selOPERBA = false, selFBQS = false, selSimplified = false;
    bool onlyOrig = false;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--orig") { onlyOrig = true; }
        else if (a == "--all") { selAll = true; selAny = true; }
        else if (a == "--dp") { selDP = true; selAny = true; }
        else if (a == "--operb") { selOPERB = true; selAny = true; }
        else if (a == "--operba" || a == "--operab") { selOPERBA = true; selAny = true; }
        else if (a == "--fbqs") { selFBQS = true; selAny = true; }
        else if (a == "--simplify" || a == "--simplified") { selSimplified = true; selAny = true; }
        else {
            std::cerr << "Warning: unknown flag '" << a << "' ignored.\n";
        }
    }

    // Resolve repo root by searching from executable directory upwards for a folder containing `data`
    std::filesystem::path repo_root;
    try {
        auto exe = std::filesystem::canonical(argv[0]);
        auto dir = exe.parent_path();
        for (int i = 0; i < 5 && !dir.empty(); ++i) {
            if (std::filesystem::exists(dir / "data")) { repo_root = dir; break; }
            dir = dir.parent_path();
        }
    } catch (...) {}
    if (repo_root.empty()) {
        auto cwd = std::filesystem::current_path();
        if (std::filesystem::exists(cwd / "data")) repo_root = cwd; else repo_root = cwd;
    }

    // Paths
    std::filesystem::path simplified_dir = repo_root / "data" / "taxi_simplified" / std::to_string(id);
    std::filesystem::path orig_data_path  = repo_root / "data" / "taxi" / (std::to_string(id) + ".txt");

    // Gather curves
    std::vector<P> orig;
    std::vector<std::pair<QString, std::vector<P>>> curves; // label -> points

    // Prefer original.txt in the simplified dir; else read from data/taxi/<id>.txt
    std::vector<P> tmp;
    bool have_original = false;
    {
        std::filesystem::path p = simplified_dir / "original.txt";
        if (std::filesystem::exists(p) && load_curve(p.string(), tmp)) {
            orig = tmp;
            have_original = true;
        } else if (std::filesystem::exists(orig_data_path) && parse_n_pairs(orig_data_path.string(), tmp)) {
            orig = tmp;
            have_original = true;
        }
    }

    // Read simplified curves from directory unless only the original is requested
    if (!onlyOrig) {
        if (!std::filesystem::exists(simplified_dir) || !std::filesystem::is_directory(simplified_dir)) {
            std::cerr << "Directory not found: " << simplified_dir << "\n";
            return 1;
        }

        for (auto& entry : std::filesystem::directory_iterator(simplified_dir)) {
            if (!entry.is_regular_file()) continue;
            auto path = entry.path();
            if (path.extension() != ".txt") continue;

            auto name = path.filename().string();
            // We'll treat original/simplified specially for consistent colors
            if (to_lower(name) == "original.txt") continue; // already handled above

            std::vector<P> pts;
            if (!load_curve(path.string(), pts)) continue;
            QString label = label_from_filename(name);
            curves.emplace_back(label, std::move(pts));
        }
    }

    if (!have_original) {
        std::cerr << "Warning: original curve not found at " << (simplified_dir / "original.txt")
                  << " or " << orig_data_path << ". Proceeding without original.\n";
    }

    // Sort curves for deterministic legend order: Original, Simplified, DP, OPERB, OPERBA, FBQS, then others
    if (!onlyOrig) {
        auto rank = [](const QString& s)->int{
            if (s == "Simplified") return 1;
            if (s == "DP") return 2;
            if (s == "OPERB") return 3;
            if (s == "OPERBA") return 4;
            if (s == "FBQS") return 5;
            return 100;
        };
        std::sort(curves.begin(), curves.end(), [&](auto& a, auto& b){ return rank(a.first) < rank(b.first); });
    }

    // Filter curves by selection if any flags were provided (ignored when --orig)
    if (!onlyOrig && selAny && !selAll) {
        std::vector<std::pair<QString, std::vector<P>>> filtered;
        auto wanted = [&](const QString& lbl){
            if (lbl == "DP") return selDP;
            if (lbl == "OPERB") return selOPERB;
            if (lbl == "OPERBA") return selOPERBA;
            if (lbl == "FBQS") return selFBQS;
            if (lbl == "Simplified") return selSimplified;
            return false; // unknown labels hidden when filtering
        };
        for (auto& it : curves) if (wanted(it.first)) filtered.push_back(std::move(it));
        curves.swap(filtered);
    }

    QApplication app(argc, argv);
    MultiViewer viewer;
    viewer.setWindowTitle("Simplified curves with different algorithms");

    // Add original if present
    if (have_original) viewer.addOriginalPoints(orig);

    // Add curves with colors unless only original is requested
    if (!onlyOrig) {
        auto color_for = [&](const QString& label){
            if (label == "Simplified") return QColor(255, 0, 0);        // red (match viewer simplified)
            if (label == "DP")         return QColor(255, 127, 0);      // orange (distinct from Simplified)
            if (label == "OPERB")      return QColor(55, 126, 184);     // blue   rgba(55, 126, 184, 1)
            if (label == "OPERBA")     return QColor(77, 175, 74);      // green  rgba(77, 175, 74, 1)
            if (label == "FBQS")       return QColor(152, 78, 163);     // purple #984ea3
            return QColor(120, 120, 120);                               // gray
        };

        for (auto& [label, pts] : curves) {
            if (label == "Simplified") {
                viewer.addSimplifiedPoints(pts);
            } else {
                viewer.addCurve(pts, color_for(label), label);
            }
        }
    }

    viewer.resize(1000, 800);
    viewer.show();
    return app.exec();
}
