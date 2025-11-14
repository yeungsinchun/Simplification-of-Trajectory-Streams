#include <algorithm>
#include <cctype>
#include <cerrno>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static std::string trim(const std::string& s) {
    size_t i = 0, j = s.size();
    while (i < j && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    while (j > i && std::isspace(static_cast<unsigned char>(s[j-1]))) --j;
    return s.substr(i, j - i);
}

static bool parse_last_two_csv_fields(const std::string& line, double& x, double& y) {
    // Find last two comma positions efficiently
    int last = -1, second_last = -1;
    for (int i = static_cast<int>(line.size()) - 1; i >= 0; --i) {
        if (line[i] == ',') {
            if (last == -1) last = i;
            else { second_last = i; break; }
        }
    }
    if (last == -1) return false; // not enough fields
    std::string sx, sy;
    sy = trim(line.substr(last + 1));
    if (second_last == -1) {
        // Only one comma, take from start to last as x
        sx = trim(line.substr(0, last));
    } else {
        sx = trim(line.substr(second_last + 1, last - (second_last + 1)));
    }
    try {
        x = std::stod(sx);
        y = std::stod(sy);
    } catch (...) {
        return false;
    }
    return true;
}

struct Points { std::vector<double> xs, ys; };

static bool read_raw_xy(const fs::path& in, Points& pts) {
    std::ifstream fin(in);
    if (!fin) return false;
    pts.xs.clear();
    pts.ys.clear();
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        double x, y;
        if (parse_last_two_csv_fields(line, x, y)) {
            pts.xs.push_back(x);
            pts.ys.push_back(y);
        }
    }
    return !pts.xs.empty();
}

static void normalize_minx_to_range(Points& pts, double targetMinX, double targetMaxX) {
    const size_t n = pts.xs.size();
    if (n == 0) return;
    double minx = pts.xs[0], maxx = pts.xs[0], miny = pts.ys[0], maxy = pts.ys[0];
    double sumy = 0.0;
    for (size_t i = 0; i < n; ++i) {
        minx = std::min(minx, pts.xs[i]);
        maxx = std::max(maxx, pts.xs[i]);
        miny = std::min(miny, pts.ys[i]);
        maxy = std::max(maxy, pts.ys[i]);
        sumy += pts.ys[i];
    }
    const double xrange = maxx - minx;
    const double ymean = sumy / static_cast<double>(n);

    // Compute the maximum absolute deviation from the mean for Y
    double maxAbsDevY = 0.0;
    for (size_t i = 0; i < n; ++i) {
        maxAbsDevY = std::max(maxAbsDevY, std::abs(pts.ys[i] - ymean));
    }

    // Scale for X to fit into [targetMinX, targetMaxX]
    const double targetWidth = (targetMaxX - targetMinX);
    const double sx = (xrange > 0.0) ? (targetWidth / xrange) : 1.0;

    // Scale for Y so that centered-at-mean values stay within [-8000, 8000]
    const double yLimit = 8000.0;
    const double sy = (maxAbsDevY > 0.0) ? (yLimit / maxAbsDevY) : 1.0;

    // Use the smaller scale so both axes remain within range
    const double scale = std::min(sx, sy);

    for (size_t i = 0; i < n; ++i) {
        // Map x linearly so that minx -> targetMinX (range may be inset if scale < sx)
        pts.xs[i] = targetMinX + (pts.xs[i] - minx) * scale;
        // Scale y linearly around its mean to preserve shape/aspect and stay within bounds
        pts.ys[i] = (pts.ys[i] - ymean) * scale; // centered around 0
    }
}

static bool write_cleaned(const fs::path& out, const Points& pts) {
    fs::create_directories(out.parent_path());
    std::ofstream fout(out);
    if (!fout) return false;
    const size_t n = pts.xs.size();
    fout << n << '\n';
    for (size_t i = 0; i < n; ++i) {
        // 15 significant digits similar to awk %.15g
        fout.setf(std::ios::fmtflags(0), std::ios::floatfield);
        fout << std::setprecision(15) << pts.xs[i] << ' ' << pts.ys[i] << '\n';
    }
    return true;
}

static void usage(const char* prog, const fs::path& src, const fs::path& out) {
    std::cerr << "Usage: " << prog << " [--all] | [--in ID]\n\n"
                 "Options:\n"
                 "  --all        Process all .txt files in " << src << "\n"
                 "  --in ID      Process single file ID (e.g. --in 16 -> " << (src/"16.txt") << ")\n"
                 "  ID           Shorthand positional (e.g. '" << prog << " 16')\n"
                 "  -h, --help   Show this message\n\n"
                 "This tool reads the last two comma-separated fields as (x,y),\n"
                 "maps min x -> -8000 and max x -> 8000, and scales y with the same\n"
                 "uniform factor around its mean (preserving shape/aspect). Output\n"
                 "is written to " << out << " with first line N, followed by N lines 'x y'.\n";
}

int main(int argc, char** argv) {
    // Locate repo root by walking up until we find the source folder
    fs::path repo_root = fs::current_path();
    fs::path found_root;
    for (int up = 0; up < 5; ++up) {
        if (fs::exists(repo_root / "taxi_log_2008_by_id")) { found_root = repo_root; break; }
        if (repo_root.has_parent_path()) repo_root = repo_root.parent_path();
        else break;
    }
    if (!found_root.empty()) repo_root = found_root;

    const fs::path source_dir = repo_root / "taxi_log_2008_by_id";
    const fs::path out_dir    = repo_root / "data" / "taxi";

    if (argc < 2) {
        usage(argv[0], source_dir, out_dir);
        return 1;
    }

    bool do_all = false; int single_id = -1;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--all") { do_all = true; }
        else if (a == "--in" && i + 1 < argc) {
            single_id = std::stoi(argv[++i]);
        } else if (a == "-h" || a == "--help") {
            usage(argv[0], source_dir, out_dir);
            return 0;
        } else if (!a.empty() && a[0] != '-' && single_id == -1) {
            // Shorthand positional ID: normalize <id>
            try { single_id = std::stoi(a); }
            catch (...) {
                std::cerr << "Invalid ID: " << a << "\n";
                return 1;
            }
        } else {
            std::cerr << "Unknown argument: " << a << "\n";
            usage(argv[0], source_dir, out_dir);
            return 1;
        }
    }

    auto process_one = [&](int id) {
        fs::path in = source_dir / (std::to_string(id) + ".txt");
        if (!fs::exists(in)) {
            std::cerr << "File not found: " << in << "\n";
            return;
        }
        Points pts;
        if (!read_raw_xy(in, pts)) {
            std::cerr << "No points parsed from: " << in << "\n";
            return;
        }
        normalize_minx_to_range(pts, -8000.0, 8000.0);
        fs::path simp_dir = repo_root / "data" / "taxi_simplified" / std::to_string(id);
        fs::create_directories(simp_dir);
        fs::path orig_path = simp_dir / "original.txt";
        if (!write_cleaned(orig_path, pts)) {
            std::cerr << "Failed to write: " << orig_path << "\n";
        } else {
            std::cout << "Wrote: " << orig_path << " (" << pts.xs.size() << " points)\n";
        }
    };

    if (do_all) {
        if (!fs::exists(source_dir)) {
            std::cerr << "Source directory not found: " << source_dir << "\n";
            return 1;
        }
        // collect ids by scanning .txt files
        std::vector<int> ids;
        for (auto& entry : fs::directory_iterator(source_dir)) {
            if (!entry.is_regular_file()) continue;
            auto p = entry.path();
            if (p.extension() != ".txt") continue;
            auto name = p.stem().string();
            try { ids.push_back(std::stoi(name)); } catch (...) { /* skip */ }
        }
        std::sort(ids.begin(), ids.end());
        for (int id : ids) process_one(id);
    } else if (single_id >= 0) {
        process_one(single_id);
    } else {
        usage(argv[0], source_dir, out_dir);
        return 1;
    }

    return 0;
}
