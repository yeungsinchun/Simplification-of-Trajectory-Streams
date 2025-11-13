#include "./inc/algorithm.hpp"
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
int main(int argc, char *argv[]) {

    if (argc < 4) {
        return 0;
    }

    double error_bound = std::stod(argv[1]);
    int size = std::stoi(argv[2]);
    Algorithm *pta = nullptr;
    std::string alg_name = argv[3];

    if (argv[3] == algorithm_type[0]) {
        pta = new DP{error_bound};
    } else if (argv[3] == algorithm_type[1]) {
        pta = new OPERB{error_bound};
    } else if (argv[3] == algorithm_type[2]) {
        pta = new OPERBA{error_bound};
    } else if (argv[3] == algorithm_type[3]) {
        pta = new FBQS{error_bound};
    }

    double tx, ty;

    double start_time = clock();

    // Read single trajectory id = <size>
    Trajectory<Point> *traj = new Trajectory<Point>;
    int id = size; // interpret argv[2] as the file id (e.g., 1)
    std::string file_name = "../../data/taxi/" + std::to_string(id) + ".txt";

    // Use C++ streams to read: first an integer N, then N pairs (x y)
    std::ifstream fin(file_name);
    if (!fin) {
        std::cerr << "Failed to open " << file_name << "\n";
        delete traj;
        return 1;
    }
    int N = 0;
    if (!(fin >> N)) {
        std::cerr << "Invalid header in " << file_name << " (expected N)\n";
        delete traj;
        return 1;
    }
    for (int j = 0; j < N; ++j) {
        if (!(fin >> tx >> ty)) {
            std::cerr << "Unexpected end of file while reading coordinates in " << file_name << "\n";
            fin.close();

            std::cout << "Running on trajectory id=" << id << " (points=" << traj->size() << ")\n";

            // Run compression
            Trajectory<Line>* result = pta ? pta->compress(traj) : nullptr;

            // Build simplified vertex sequence from segments
            std::vector<Point> simplified_pts;
            if (result && result->size() > 0) {
                simplified_pts.reserve(result->size() + 1);
                // start point of first segment
                auto first_seg = (*result)[0];
                simplified_pts.push_back(first_seg.start_point());
                for (std::size_t k = 0; k < result->size(); ++k) {
                    auto seg = (*result)[k];
                    simplified_pts.push_back(seg.end_point());
                }
            } else {
                // Fallback: no segments produced; use the original points
                simplified_pts.reserve(traj->size());
                for (std::size_t k = 0; k < traj->size(); ++k) simplified_pts.push_back((*traj)[k]);
            }

            // Write to ../../data/taxi_simplified/<id>/<alg>_simplified.txt in two-line x/y format
            try {
                std::filesystem::path out_dir = std::filesystem::path("../../data/taxi_simplified") / std::to_string(id);
                std::filesystem::create_directories(out_dir);
                std::filesystem::path out_path = out_dir / (alg_name + std::string("_simplified.txt"));
                std::ofstream fout(out_path);
                if (!fout) {
                    std::cerr << "Failed to open output file: " << out_path << "\n";
                } else {
                    for (std::size_t i = 0; i < simplified_pts.size(); ++i) {
                        if (i) fout << ' ';
                        fout << simplified_pts[i].x;
                    }
                    fout << "\n";
                    for (std::size_t i = 0; i < simplified_pts.size(); ++i) {
                        if (i) fout << ' ';
                        fout << simplified_pts[i].y;
                    }
                    fout << "\n";
                    std::cout << "Wrote simplified to " << out_path << " (points=" << simplified_pts.size() << ")\n";
                }
            } catch (const std::exception& e) {
                std::cerr << "Error writing output: " << e.what() << "\n";
            }

            delete result;
            delete traj;

            double end_time = clock();
            std::cout << "Elapsed: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s\n";
            return 0;
        }
    }
}