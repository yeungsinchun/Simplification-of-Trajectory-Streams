#ifndef WORKER_POOL_H
#define WORKER_POOL_H

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <thread>
#include <vector>

// Persistent helper threads for short, frequent data-parallel loops.
//
// run(count, fn) calls fn(begin, end, worker) on chunks covering [0, count)
// across all workers, the calling thread being worker 0, and returns once
// every call has finished. `worker` identifies the calling participant, so
// callers can keep per-worker state without synchronisation.
//
// Work is claimed in small chunks through one ticket holding (generation,
// next index); a claim is valid only while the generation still matches, and
// run() waits for finished items rather than for every helper, so a helper
// the OS has descheduled never stalls a run. Idle helpers spin for up to
// kIdleSpin and then sleep on a condition variable.
class WorkerPool {
  public:
    explicit WorkerPool(int workers) : workers_(workers) {
        for (int w = 1; w < workers; ++w)
            threads_.emplace_back([this, w] { helper_loop(w); });
    }
    ~WorkerPool() {
        stop_.store(true);
        publish(0, nullptr, nullptr);
        for (auto &thread : threads_)
            thread.join();
    }
    WorkerPool(const WorkerPool &) = delete;
    WorkerPool &operator=(const WorkerPool &) = delete;

    int workers() const { return workers_; }

    template <class Fn> void run(int count, Fn &fn) {
        const uint64_t generation =
            publish(count, &fn, [](void *job, int begin, int end, int worker) {
                (*static_cast<Fn *>(job))(begin, end, worker);
            });
        drain(0, generation);
        while (done_.load(std::memory_order_acquire) != count)
            cpu_relax();
    }

  private:
    using Call = void (*)(void *, int, int, int);
    struct Job {
        std::atomic<void *> ctx{nullptr};
        std::atomic<Call> call{nullptr};
        std::atomic<int> count{0};
    };
    static constexpr int kChunk = 4;

    // Runs arrive every few microseconds while a caller is busy; helpers idle
    // this long (e.g. after the caller finished) go to sleep.
    static constexpr auto kIdleSpin = std::chrono::milliseconds(1);

    static void cpu_relax() {
#if defined(__aarch64__)
        asm volatile("yield");
#elif defined(__x86_64__) || defined(__i386__)
        __builtin_ia32_pause();
#endif
    }

    // Jobs alternate between two slots; a slot is rewritten only after the
    // generation using it has finished and the next one has been published.
    uint64_t publish(int count, void *ctx, Call call) {
        const uint64_t generation = (ticket_.load() >> 32) + 1;
        Job &job = jobs_[generation & 1];
        job.ctx.store(ctx, std::memory_order_relaxed);
        job.call.store(call, std::memory_order_relaxed);
        job.count.store(count, std::memory_order_relaxed);
        done_.store(0, std::memory_order_relaxed);
        ticket_.store(generation << 32);
        if (sleepers_.load() > 0) {
            std::lock_guard<std::mutex> lock(mutex_);
            wake_.notify_all();
        }
        return generation;
    }

    void drain(int worker, uint64_t generation) {
        const Job &job = jobs_[generation & 1];
        uint64_t ticket = ticket_.load(std::memory_order_acquire);
        for (;;) {
            if ((ticket >> 32) != generation)
                return;
            const int begin = int(uint32_t(ticket));
            const int count = job.count.load(std::memory_order_relaxed);
            if (begin >= count)
                return;
            if (!ticket_.compare_exchange_weak(ticket, ticket + kChunk,
                                               std::memory_order_acq_rel,
                                               std::memory_order_acquire))
                continue;
            const int end = std::min(begin + kChunk, count);
            void *ctx = job.ctx.load(std::memory_order_relaxed);
            const Call call = job.call.load(std::memory_order_relaxed);
            call(ctx, begin, end, worker);
            done_.fetch_add(end - begin, std::memory_order_release);
            ticket = ticket_.load(std::memory_order_acquire);
        }
    }

    void helper_loop(int worker) {
        uint64_t seen = 0;
        for (;;) {
            uint64_t generation;
            auto idle_since = std::chrono::steady_clock::now();
            for (int spins = 0; (generation = ticket_.load() >> 32) == seen;
                 ++spins) {
                cpu_relax();
                if (spins % 64 != 63 ||
                    std::chrono::steady_clock::now() - idle_since < kIdleSpin)
                    continue;
                std::unique_lock<std::mutex> lock(mutex_);
                sleepers_.fetch_add(1);
                wake_.wait(lock,
                           [&] { return (ticket_.load() >> 32) != seen; });
                sleepers_.fetch_sub(1);
            }
            seen = generation;
            if (stop_.load())
                return;
            drain(worker, generation);
        }
    }

    const int workers_;
    std::vector<std::thread> threads_;
    Job jobs_[2];
    std::atomic<uint64_t> ticket_{0};
    std::atomic<int> done_{0}, sleepers_{0};
    std::atomic<bool> stop_{false};
    std::mutex mutex_;
    std::condition_variable wake_;
};

#endif // WORKER_POOL_H
