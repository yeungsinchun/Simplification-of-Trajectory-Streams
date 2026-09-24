#ifndef SIMPLIFY_CORE_H
#define SIMPLIFY_CORE_H

#include <array>
#include <cstdint>
#include <vector>

#include "simplify_geometry.h"
#include "timer.h"

// ===========================================================================
//  Core algorithm (opt-in TIMER sites; active only when --time is set)
// ===========================================================================
//
// One stab loop shared by the headless (simplify.cpp) and GUI
// (simplify_with_gui.cpp) front-ends, so the viewer shows exactly the anchors,
// wedges and stab regions the headless run computes.

struct StabScratch {
    std::vector<Point> P, Gi, F;
    std::vector<std::vector<Point>> S;
    std::vector<int> active;
    std::vector<AxisBounds> stab_bounds;
    std::vector<uint8_t> anchor_outside;
    PreparedClipPolygon prepared_Gi;
    sh_double::FastClipBuffers clip_buffers;
};

// Hooks into get_longest_stab. The default no-op observer compiles away, so
// the headless build runs the bare loop.
struct NoStabObserver {
    // p0 starts a stab.
    void stab_begin(const Point & /*p0*/) {}
    // Anchor `anchor` reached intersect with stab region S (before this step)
    // and wedge F = F(S, P[anchor]). Anchors dropped by find_F's disjoint
    // prune never reach intersect and are not reported.
    void anchor_wedge(int /*anchor*/, const std::vector<Point> & /*S*/,
                      const std::vector<Point> & /*F*/) {}
    // The anchor last passed to anchor_wedge intersected Gi and stays live.
    void anchor_survived(int /*anchor*/) {}
    // stream[cur] was consumed with hull Gi and at least one live anchor.
    void step_end(int /*cur*/, const Point & /*pi*/,
                  const std::vector<Point> & /*Gi*/) {}
    // The stab emitted segment {anchor, S.front()}.
    void stab_end(const std::array<Point, 2> & /*segment*/) {}
};

template <class Observer = NoStabObserver>
int get_longest_stab(const std::vector<Point> &stream, int cur,
                     std::vector<Point> &simplified, double EPSILON,
                     double DELTA, StabScratch &scratch,
                     Observer &&observer = {}) {
    TIMER("get_longest_stab");
    const Point &p0 = stream[cur];
    auto &P = scratch.P;
    auto &Gi = scratch.Gi;
    auto &S = scratch.S;
    auto &F = scratch.F;
    auto &active = scratch.active;
    auto &stab_bounds = scratch.stab_bounds;
    auto &anchor_outside = scratch.anchor_outside;
    observer.stab_begin(p0);
    {
        TIMER("boundary_P");
        P = get_boundary_points_from_grid(p0, EPSILON, DELTA);
    }
    std::array<Point, 2> buffer = {p0, p0};
    const int Pn = (int)P.size();
    S.resize(Pn);
    stab_bounds.resize(Pn);
    anchor_outside.assign(Pn, 0);
    active.clear();
    active.reserve(Pn);
    for (int i = 0; i < Pn; ++i) {
        S[i].clear();
        S[i].push_back(P[i]);
        active.push_back(i);
    }

    cur++;
    while (cur < int(stream.size())) {
        {
            TIMER("hull_Gi");
            Gi = get_conv_from_grid(stream[cur], EPSILON, DELTA);
        }
        prepare_clip_polygon(Gi, scratch.prepared_Gi);
        std::vector<Point> bbox_result;
        AxisBounds bbox_bounds;
        bool bbox_cached = false, bbox_hit = false;
        size_t surviving = 0;
        for (int i : active) {
            bool full_bbox, disjoint;
            {
                TIMER("find_F");
                full_bbox =
                    find_F(P[i], S[i], F, &scratch.prepared_Gi, &disjoint,
                           &stab_bounds[i], &anchor_outside[i]);
            }
            if (disjoint)
                continue;
            observer.anchor_wedge(i, S[i], F);
            bool hit;
            {
                TIMER("intersect");
                if (full_bbox && bbox_cached) {
                    hit = bbox_hit;
                    S[i] = bbox_result;
                    stab_bounds[i] = bbox_bounds;
                } else {
                    hit = intersect_prepared(F, scratch.prepared_Gi, S[i],
                                             anchor_outside[i] && !full_bbox
                                                 ? nullptr
                                                 : &stab_bounds[i],
                                             scratch.clip_buffers);
                    if (full_bbox) {
                        bbox_result = S[i];
                        bbox_bounds = stab_bounds[i];
                        bbox_hit = hit;
                        bbox_cached = true;
                    }
                }
            }
            if (!hit)
                continue;
            observer.anchor_survived(i);
            active[surviving++] = i;
        }
        active.resize(surviving);
        if (active.empty())
            break;
        const int chosen = active.back();
        buffer[0] = P[chosen];
        buffer[1] = S[chosen].front();
        observer.step_end(cur, stream[cur], Gi);
        cur++;
    }
    simplified.emplace_back(buffer[0]);
    simplified.emplace_back(buffer[1]);
    observer.stab_end(buffer);
    return cur;
}

#endif // SIMPLIFY_CORE_H
