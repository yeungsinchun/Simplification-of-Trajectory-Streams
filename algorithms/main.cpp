#include "./inc/dp.hpp"
#include "./inc/fbqs.hpp"
#include "./inc/operb.hpp"
#include "./inc/trajectory.hpp"
#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

#pragma comment(linker, "/STACK:1024000000,1024000000")

std::string algorithm_type[4] = {"dp", "operb", "operba", "fbqs"};

static std::filesystem::path resolve_repo_root(const char* argv0) {
    try {
        auto exe = std::filesystem::canonical(argv0);
        auto dir = exe.parent_path();
        // Walk up a few levels to find a directory containing `data`
        for (int i = 0; i < 4 && !dir.empty(); ++i) {
            if (std::filesystem::exists(dir / "data")) return dir;
            dir = dir.parent_path();
        }
    } catch (...) {
        // fall through to cwd
    }
    auto cwd = std::filesystem::current_path();
    if (std::filesystem::exists(cwd / "data")) return cwd;
    return cwd; // fallback
}

static void print_help(const char* prog) {
    std::cerr << "Usage: " << prog << " <id> [error_bound]\n";
    std::cerr << "  <id>           Required trajectory id (integer)\n";
    std::cerr << "  [error_bound]  Optional error bound (default: 200)\n";
    std::cerr << "Example: " << prog << " 1 100\n";
    std::cerr << "Runs all algorithms (dp, operb, operba, fbqs) on taxi/<id>.txt and writes outputs to data/taxi_simplified/<id>/.\n";
}

int main(int argc, char *argv[]) {
    // Expect: ./main <id> [error_bound]
    if (argc < 2 || (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))) {
        print_help(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    int id = 0;
    double error_bound = 200.0; // default
    try {
        id = std::stoi(argv[1]);
        if (argc >= 3) {
            error_bound = std::stod(argv[2]);
        } else {
            std::cout << "No error_bound provided. Using default error_bound=200\n";
        }
    } catch (const std::exception&) {
        std::cerr << "Invalid arguments. id must be integer; error_bound must be a number (if provided).\n";
        print_help(argv[0]);
        return 1;
    }

    double tx, ty;

    double start_time = clock();

    // Resolve repository root and build absolute paths
    auto repo_root = resolve_repo_root(argv[0]);
    Trajectory<Point> *traj = new Trajectory<Point>;
    // Prefer data/taxi_simplified/<id>/original.txt (N+pairs), fallback to data/taxi/<id>.txt
    std::filesystem::path simp_orig = repo_root / "data" / "taxi_simplified" / std::to_string(id) / "original.txt";
    std::filesystem::path file_path = simp_orig;
    std::string file_name = file_path.string();

    // Parse strictly as N then N lines/pairs of x y
    std::ifstream fin(file_name);
    if (!fin) {
        std::cerr << "Failed to open " << file_name << "\n";
        delete traj; return 1;
    }
    int N = 0;
    if (!(fin >> N)) {
        std::cerr << "Invalid header in " << file_name << " (expected N)\n";
        delete traj; return 1;
    }
    for (int i = 0; i < N; ++i) {
        if (!(fin >> tx >> ty)) {
            std::cerr << "Unexpected end of file while reading coordinates in " << file_name << "\n";
            delete traj; return 1;
        }
        traj->push(Point{tx, ty});
    }
    fin.close();

    std::cout << "Running on trajectory id=" << id << " (points=" << traj->size() << ") with error_bound=" << error_bound << "\n";

    // Ensure output directory exists once
    std::filesystem::path out_dir = repo_root / "data" / "taxi_simplified" / std::to_string(id);
    try {
        std::filesystem::create_directories(out_dir);
    } catch (const std::exception& e) {
        std::cerr << "Failed to create output directory: " << out_dir << ": " << e.what() << "\n";
        delete traj;
        return 1;
    }

    // Helper lambda to run one algorithm and write output
    auto run_and_write = [&](const std::string& alg_name) {
        double t0 = clock();
        Trajectory<Line>* result = nullptr;
        if (alg_name == algorithm_type[0]) {
            DP algo{error_bound};
            result = algo.compress(traj);
        } else if (alg_name == algorithm_type[1]) {
            OPERB algo{error_bound};
            result = algo.compress(traj);
        } else if (alg_name == algorithm_type[2]) {
            OPERBA algo{error_bound};
            result = algo.compress(traj);
        } else if (alg_name == algorithm_type[3]) {
            FBQS algo{error_bound};
            result = algo.compress(traj);
        } else {
            std::cerr << "Unknown algorithm: " << alg_name << "\n";
            return;
        }
        double t1 = clock();

        // Build simplified vertex sequence from segments
        std::vector<Point> simplified_pts;
        if (result && result->size() > 0) {
            simplified_pts.reserve(result->size() + 1);
            auto first_seg = (*result)[0];
            simplified_pts.push_back(first_seg.start_point());
            for (std::size_t k = 0; k < result->size(); ++k) {
                auto seg = (*result)[k];
                simplified_pts.push_back(seg.end_point());
            }
        } else {
            simplified_pts.reserve(traj->size());
            for (std::size_t k = 0; k < traj->size(); ++k) simplified_pts.push_back((*traj)[k]);
        }

        // Write to data/taxi_simplified/<id>/<alg>_simplified.txt in N + pairs format
        try {
            std::filesystem::path out_path = out_dir / (alg_name + std::string("_simplified.txt"));
            std::ofstream fout(out_path);
            if (!fout) {
                std::cerr << "Failed to open output file: " << out_path << "\n";
            } else {
                std::size_t Nsim = simplified_pts.size();
                fout << Nsim << '\n';
                for (std::size_t i = 0; i < Nsim; ++i) {
                    fout << simplified_pts[i].x << ' ' << simplified_pts[i].y << '\n';
                }
                std::cout << "[" << alg_name << "] points=" << simplified_pts.size() << ", time="
                          << (double)(t1 - t0) / CLOCKS_PER_SEC << "s -> " << out_path << "\n";
            }
        } catch (const std::exception& e) {
            std::cerr << "Error writing output for " << alg_name << ": " << e.what() << "\n";
        }

        delete result;
    };

    // Run all algorithms
    for (const auto& name : algorithm_type) {
        run_and_write(name);
    }

    delete traj;

    double end_time = clock();
    std::cout << "Total elapsed: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s\n";
    return 0;
}