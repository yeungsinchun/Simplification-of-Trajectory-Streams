// Long-lived, single-threaded workload for Instruments / perf. Not a timing
// benchmark: repetitions warm the thread-local grid caches after the first run.
#include "simplify_core.h"

#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cerr << "Usage: profile_sots original.txt epsilon delta repetitions\n";
        return 1;
    }
    std::ifstream input(argv[1]);
    size_t count;
    if (!(input >> count) || count == 0)
        return 1;
    std::vector<Point> stream;
    for (size_t i = 0; i < count; ++i) {
        double x, y;
        if (!(input >> x >> y))
            return 1;
        stream.emplace_back(x, y);
    }
    const double epsilon = std::stod(argv[2]), delta = std::stod(argv[3]);
    const int repetitions = std::stoi(argv[4]);
    if (!(epsilon > 0) || !(delta > 0) || repetitions < 1)
        return 1;
    configure_bbox(stream, epsilon, delta);
    double checksum = 0;
    for (int i = 0; i < repetitions; ++i) {
        const auto result = Simplifier(epsilon, delta).simplify(stream);
        checksum += CGAL::to_double(result.back().x()) + result.size();
    }
    std::cout << checksum << '\n';
}
