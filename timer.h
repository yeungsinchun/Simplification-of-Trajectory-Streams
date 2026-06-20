#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <functional>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdio>

namespace timer_detail {
    inline std::map<std::string, double>& timing() {
        static std::map<std::string, double> t;
        return t;
    }
    inline std::map<std::string, long long>& counters() {
        static std::map<std::string, long long> c;
        return c;
    }
    inline std::map<std::string, double>& child_time() {
        static std::map<std::string, double> c;
        return c;
    }
    inline std::map<std::string, std::vector<std::string>>& children() {
        static std::map<std::string, std::vector<std::string>> ch;
        return ch;
    }
    inline std::vector<const char*>& stack() {
        static std::vector<const char*> s;
        return s;
    }
}

struct Timer {
    std::chrono::high_resolution_clock::time_point start;
    const char* name;
    Timer(const char* n) : name(n) {
        start = std::chrono::high_resolution_clock::now();
        timer_detail::stack().push_back(name);
    }
    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration<double, std::milli>(end - start).count();
        timer_detail::timing()[name] += dur;
        timer_detail::counters()[name]++;
        auto& stack = timer_detail::stack();
        if (stack.size() >= 2) {
            const char* parent = stack[stack.size() - 2];
            timer_detail::child_time()[parent] += dur;
            timer_detail::children()[parent].push_back(name);
        }
        stack.pop_back();
    }
};

#define TIMER(name) Timer _timer(name)

inline void print_timing_summary() {
    auto& t = timer_detail::timing();
    if (t.empty()) return;
    double grand = t.count("total") ? t["total"] : 0.0;

    fprintf(stderr, "\n========== TIMING SUMMARY ==========\n");
    // Column widths must match between header and data.  Numeric columns
    // are right-aligned at fixed positions; the name column always has
    // the same total width (indent + name = 28) so the numbers stay
    // anchored regardless of depth.  When the name doesn't fit in the
    // remaining space it is truncated.
    const int NAME_TOTAL_WIDTH = 28;
    std::string header_left = "Operation";
    header_left.append(NAME_TOTAL_WIDTH - header_left.size(), ' ');
    fprintf(stderr, "%s %12s %12s %10s %9s %9s\n",
            header_left.c_str(),
            "Wall (ms)", "Self (ms)", "Calls", "% Total", "% Parent");
    fprintf(stderr, "--------------------------------------------------------------------------------\n");

    auto report_one = [&](const std::string& name, int depth,
                          double parent_wall) {
        auto it = t.find(name);
        if (it == t.end()) return;
        double wall = it->second;
        long long calls = timer_detail::counters()[name];
        double self = wall - timer_detail::child_time()[name];
        double pct_total = (grand > 0.0 ? wall / grand * 100.0 : 0.0);
        double pct_parent = (parent_wall > 0.0 ? wall / parent_wall * 100.0 : 0.0);
        int indent_cols = depth * 2;
        int name_cols = NAME_TOTAL_WIDTH - indent_cols;
        if (name_cols < 1) name_cols = 1;
        std::string shown = name;
        if ((int)shown.size() > name_cols) shown.resize(name_cols);
        // Build the left field and pad/truncate to exactly NAME_TOTAL_WIDTH.
        std::string left = std::string(indent_cols, ' ') + shown;
        if ((int)left.size() < NAME_TOTAL_WIDTH) left.append(NAME_TOTAL_WIDTH - left.size(), ' ');
        else if ((int)left.size() > NAME_TOTAL_WIDTH) left.resize(NAME_TOTAL_WIDTH);
        fprintf(stderr, "%s %12.2f %12.2f %10lld %8.1f%% %8.1f%%\n",
                left.c_str(), wall, self, calls, pct_total, pct_parent);
    };

    // Recursive tree print starting at `total`.  `simplify` is the
    // outermost function-level timer; it is structurally identical to
    // `total` so we skip it and reparent its children to `total`.
    std::function<void(const std::string&, int, double)> print_recursive;
    print_recursive = [&](const std::string& name, int depth, double parent_wall) {
        if (!t.count(name)) return;
        if (name == "simplify") {
            // Re-emit children of simplify at depth-1 (under total).
            std::set<std::string> printed;
            for (const auto& child : timer_detail::children()[name]) {
                if (printed.insert(child).second && t.count(child)) {
                    print_recursive(child, depth, parent_wall);
                }
            }
            return;
        }
        report_one(name, depth, parent_wall);
        double my_wall = t[name];
        std::set<std::string> printed;
        for (const auto& child : timer_detail::children()[name]) {
            if (printed.insert(child).second && t.count(child)) {
                print_recursive(child, depth + 1, my_wall);
            }
        }
    };

    print_recursive("total", 0, grand);
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    fprintf(stderr, "%-28s %12.2f\n", "TOTAL (wall time)", grand);
    fprintf(stderr, "====================================\n");
}

#endif // TIMER_H
