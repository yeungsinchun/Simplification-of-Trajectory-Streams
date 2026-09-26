#ifndef SIMPLIFY_CORE_H
#define SIMPLIFY_CORE_H

#include <array>
#include <vector>

#include "simplify_geometry.h"
#include "timer.h"

// ===========================================================================
//  Core algorithm (opt-in TIMER sites; active only when --time is set)
// ===========================================================================
//
// Each stab starts at p0 = stream[cur] with the anchors P, the boundary grid
// samples around p0. Anchor p keeps a stab region S, the endpoints q such
// that segment pq still passes near every point consumed so far; initially
// S = {p}. Consuming p_i updates every live anchor to
//
//     S ← F(S, p) ∩ conv(G_i)
//
// and drops the anchors whose S becomes empty. The stab ends when no anchor
// survives (or the stream ends) and emits the segment from the last anchor,
// in P order, alive after the final covered point to the first vertex of
// its S.

// Reusable working memory for advancing anchors through one stream step. The
// whole-box clip is shared by every anchor with that wedge at the step.
class StepWorkspace {
  public:
    // Called before the first anchor of each step.
    void begin_step() { box_clip_ready_ = false; }

  private:
    friend class Anchor;

    std::vector<Point> F_; // wedge of the anchor being advanced
    ConvexClipper clipper_;
    std::vector<Point> box_clip_; // F = whole box: box ∩ conv(G_i)
    AxisBounds box_clip_bounds_;
    bool box_clip_nonempty_ = false;
    bool box_clip_ready_ = false;
};

// A candidate start point p of the current stab and its stab region S.
class Anchor {
  public:
    void reset(const Point &p) {
        p_ = p;
        S_.assign(1, p);
        p_outside_S_ = false;
    }

    // S ← F(S, p) ∩ G. Returns false when the intersection is empty, after
    // which the anchor is dead for the rest of the stab.
    bool advance(const ConvexRegion &G, StepWorkspace &ws) {
        Wedge wedge;
        {
            TIMER("find_F");
            wedge =
                find_F(p_, S_, ws.F_, &G.bounds(), &S_bounds_, &p_outside_S_);
        }
        if (wedge == Wedge::misses_target)
            return false;
        TIMER("intersect");
        // Every anchor whose wedge is the whole box gets the same S.
        const bool whole_box = wedge == Wedge::whole_box;
        if (whole_box && ws.box_clip_ready_) {
            S_ = ws.box_clip_;
            S_bounds_ = ws.box_clip_bounds_;
            return ws.box_clip_nonempty_;
        }
        // S_bounds_ only serves the p ∈ S test, which the latch now skips.
        const bool nonempty =
            ws.clipper_.clip(ws.F_, G, S_, p_outside_S_ ? nullptr : &S_bounds_);
        if (whole_box) {
            ws.box_clip_ = S_;
            ws.box_clip_bounds_ = S_bounds_;
            ws.box_clip_nonempty_ = nonempty;
            ws.box_clip_ready_ = true;
        }
        return nonempty;
    }

    // Segment from the anchor into its stab region.
    std::array<Point, 2> segment() const { return {p_, S_.front()}; }

  private:
    Point p_;
    std::vector<Point> S_;
    // Shortcuts for the p ∈ S test in find_F: a p outside S's bbox is
    // outside S, and once p is outside S it stays outside every later S.
    AxisBounds S_bounds_;
    bool p_outside_S_ = false;
};

// The anchors of the current stab, in P order, and which of them are live.
// Storage is kept across stabs, so steady-state stabs do not allocate.
class AnchorSet {
  public:
    void reset(const std::vector<Point> &P) {
        anchors_.resize(P.size());
        live_.resize(P.size());
        for (size_t i = 0; i < P.size(); ++i) {
            anchors_[i].reset(P[i]);
            live_[i] = int(i);
        }
    }

    // Keeps the live anchors for which keep(anchor) is true, in order.
    template <class Keep> void retain(Keep &&keep) {
        size_t kept = 0;
        for (const int i : live_)
            if (keep(anchors_[i]))
                live_[kept++] = i;
        live_.resize(kept);
    }

    bool empty() const { return live_.empty(); }
    const Anchor &last() const { return anchors_[live_.back()]; }

  private:
    std::vector<Anchor> anchors_;
    std::vector<int> live_;
};

// Streaming simplification by repeated longest stabs.
//
// CRTP base: Derived may shadow advance_anchors() to change how one step
// updates the live anchors; calls are resolved at compile time, so the
// default sequential loop pays nothing for the customization point.
template <class Derived> class StreamSimplifier {
  public:
    StreamSimplifier(double epsilon, double delta)
        : epsilon_(epsilon), delta_(delta) {}

    std::vector<Point> simplify(const std::vector<Point> &stream) {
        std::vector<Point> simplified;
        int cur = 0;
        while (cur != int(stream.size()))
            cur = longest_stab(stream, cur, simplified);
        return simplified;
    }

  protected:
    // S ← F(S, p) ∩ G for every live anchor, dropping the emptied ones.
    void advance_anchors(const ConvexRegion &G) {
        workspace_.begin_step();
        anchors_.retain(
            [&](Anchor &anchor) { return anchor.advance(G, workspace_); });
    }

    AnchorSet anchors_;
    StepWorkspace workspace_;

  private:
    Derived &derived() { return static_cast<Derived &>(*this); }

    // Consumes stream[cur..] while some anchor survives, appends the stab's
    // segment to `simplified`, and returns the index of the first point the
    // stab could not cover.
    int longest_stab(const std::vector<Point> &stream, int cur,
                     std::vector<Point> &simplified) {
        TIMER("get_longest_stab");
        const Point &p0 = stream[cur];
        {
            TIMER("boundary_P");
            anchors_.reset(get_boundary_points_from_grid(p0, epsilon_, delta_));
        }
        std::array<Point, 2> segment = {p0, p0};
        for (++cur; cur < int(stream.size()); ++cur) {
            {
                TIMER("hull_Gi");
                Gi_.assign(get_conv_from_grid(stream[cur], epsilon_, delta_));
            }
            derived().advance_anchors(Gi_);
            if (anchors_.empty())
                break;
            segment = anchors_.last().segment();
        }
        simplified.push_back(segment[0]);
        simplified.push_back(segment[1]);
        return cur;
    }

    double epsilon_, delta_;
    ConvexRegion Gi_;
};

// The sequential simplifier.
class Simplifier final : public StreamSimplifier<Simplifier> {
  public:
    using StreamSimplifier::StreamSimplifier;
};

#endif // SIMPLIFY_CORE_H
