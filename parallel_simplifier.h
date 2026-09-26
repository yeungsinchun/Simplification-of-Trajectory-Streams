#ifndef PARALLEL_SIMPLIFIER_H
#define PARALLEL_SIMPLIFIER_H

#include <cstdint>
#include <optional>
#include <vector>

#include "simplify_core.h"
#include "timer.h"
#include "worker_pool.h"

// StreamSimplifier that advances the live anchors of a step on a WorkerPool.
//
// Within one step, Anchor::advance reads only the shared hull G and writes
// only its own anchor, so anchors can be advanced in any order on any worker
// (each with its own StepWorkspace). Survivors are then compacted in P order,
// so the output is identical to the sequential Simplifier for any number of
// workers.
class ParallelSimplifier final : public StreamSimplifier<ParallelSimplifier> {
  public:
    ParallelSimplifier(double epsilon, double delta, int workers)
        : StreamSimplifier(epsilon, delta), workers_(workers) {}

  private:
    friend StreamSimplifier;

    // A parallel step pays a hand-off and a join, so only steps with at least
    // this much work (live anchors x hull edges) use the pool; small steps
    // (coarse epsilon, few live anchors) stay on the calling thread.
    static constexpr size_t kMinParallelWork = 256;

    void advance_anchors(const ConvexRegion &G) {
        const size_t count = anchors_.live_count();
        // Timers keep global state, so --time runs stay sequential.
        if (timer_detail::enabled() ||
            count * G.edges().size() < kMinParallelWork) {
            StreamSimplifier::advance_anchors(G);
            return;
        }
        alive_.resize(count);
        for (auto &worker : workers_)
            worker.workspace.begin_step();
        auto advance = [&](int begin, int end, int worker) {
            for (int k = begin; k < end; ++k)
                alive_[k] =
                    anchors_.live(k).advance(G, workers_[worker].workspace);
        };
        // Started on the first parallel step: runs that never have one
        // (coarse epsilon) do not pay for starting threads.
        if (!pool_)
            pool_.emplace(int(workers_.size()));
        pool_->run(int(count), advance);
        anchors_.retain_flagged(alive_);
    }

    // Workers write their workspace on every anchor; a cache line each keeps
    // them from invalidating one another's.
    struct alignas(128) Worker {
        StepWorkspace workspace;
    };
    std::vector<Worker> workers_;
    std::vector<uint8_t> alive_;
    std::optional<WorkerPool> pool_; // last: its threads stop first
};

#endif // PARALLEL_SIMPLIFIER_H
