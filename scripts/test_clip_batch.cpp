// Compare the optimized edge crop with a straightforward Sutherland-Hodgman
// oracle, including non-convex rings, wraparound crossings and exact
// boundaries.
#include "simplify_geometry.h"

#include <bit>
#include <iostream>
#include <random>
#include <stdexcept>

using sh_double::ClipEdge;
using sh_double::Vec2;

static std::vector<Vec2> scalar_crop(const std::vector<Vec2> &ring,
                                     const ClipEdge &edge) {
    std::vector<Vec2> out;
    if (ring.empty())
        return out;
    auto side = [&](const Vec2 &p) {
        return edge.dx * (p[1] - edge.start[1]) -
               edge.dy * (p[0] - edge.start[0]);
    };
    Vec2 prev = ring.back();
    double a = side(prev);
    for (const auto &curr : ring) {
        const double b = side(curr);
        if ((a >= 0) != (b >= 0)) {
            const double t = a / (a - b);
            out.push_back({prev[0] + t * (curr[0] - prev[0]),
                           prev[1] + t * (curr[1] - prev[1])});
        }
        if (b >= 0)
            out.push_back(curr);
        prev = curr;
        a = b;
    }
    return out;
}

static void check(const std::vector<Vec2> &ring, const ClipEdge &edge) {
    const auto expected = scalar_crop(ring, edge);
    std::vector<Vec2> actual(ring.size() * 2 + 2);
    size_t size = 0;
    if (sh_double::crop_to_left_of_edge(ring.data(), ring.size(), edge,
                                        actual.data(), size))
        actual.resize(size);
    else
        actual = ring;
    if (actual.size() != expected.size())
        throw std::runtime_error("size differs");
    for (size_t i = 0; i < actual.size(); ++i)
        for (int axis = 0; axis < 2; ++axis)
            if (std::bit_cast<uint64_t>(actual[i][axis]) !=
                std::bit_cast<uint64_t>(expected[i][axis]))
                throw std::runtime_error("coordinate bits differ");
}

int main() {
    // Every combination of inside/outside/on-edge vertices around block sizes.
    for (size_t n = 0; n <= 10; ++n) {
        size_t combinations = 1;
        for (size_t i = 0; i < n; ++i)
            combinations *= 3;
        for (size_t pattern = 0; pattern < combinations; ++pattern) {
            size_t code = pattern;
            std::vector<Vec2> ring;
            for (size_t i = 0; i < n; ++i, code /= 3)
                ring.push_back({double(i), double(int(code % 3) - 1)});
            check(ring, {{0, 0}, 1, 0});
        }
    }
    std::mt19937_64 random(20260928);
    std::uniform_real_distribution<double> coord(-1e6, 1e6);
    for (size_t trial = 0; trial < 20000; ++trial) {
        std::vector<Vec2> ring(trial % 129);
        for (auto &point : ring)
            point = {coord(random), coord(random)};
        check(ring,
              {{coord(random), coord(random)}, coord(random), coord(random)});
    }
    std::cout << "Edge cropping matches scalar oracle bit for bit\n";
}
