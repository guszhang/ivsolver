/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * Minimal portable parallel-for built on std::thread. Work is split into fixed-size
 * chunks that are handed out dynamically, so results that are reduced per chunk are
 * independent of the thread count (deterministic output).
 *
 * In single-threaded builds (e.g. WebAssembly without pthreads) everything runs
 * on the calling thread.
 */
#pragma once

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <vector>

#if !defined(__EMSCRIPTEN__) || defined(__EMSCRIPTEN_PTHREADS__)
#include <thread>
#define IVSOLVER_HAS_THREADS 1
#endif

namespace ivsolver
{
    inline int resolve_thread_count(int requested)
    {
#ifdef IVSOLVER_HAS_THREADS
        if (requested > 0)
        {
            return requested;
        }
        unsigned hw = std::thread::hardware_concurrency();
        return hw == 0 ? 1 : static_cast<int>(hw);
#else
        (void)requested;
        return 1;
#endif
    }

    // Calls body(chunk_begin, chunk_end, chunk_index) for consecutive chunks of
    // [0, count). Exceptions thrown inside a worker are rethrown on the caller.
    template <typename Body>
    void parallel_chunks(size_t count, size_t chunk_size, int threads, Body body)
    {
        if (count == 0)
        {
            return;
        }
        chunk_size = std::max<size_t>(1, chunk_size);
        const size_t chunks = (count + chunk_size - 1) / chunk_size;
        int workers = std::min<int>(resolve_thread_count(threads), static_cast<int>(chunks));
        std::atomic<size_t> next{0};
        std::exception_ptr failure;
        std::mutex failure_mutex;
        auto run = [&]()
        {
            try
            {
                for (size_t chunk = next++; chunk < chunks; chunk = next++)
                {
                    size_t begin = chunk * chunk_size;
                    size_t end = std::min(count, begin + chunk_size);
                    body(begin, end, chunk);
                }
            }
            catch (...)
            {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure)
                {
                    failure = std::current_exception();
                }
                next = chunks;
            }
        };
#ifdef IVSOLVER_HAS_THREADS
        if (workers > 1)
        {
            std::vector<std::thread> pool;
            pool.reserve(workers - 1);
            for (int i = 0; i < workers - 1; i++)
            {
                pool.emplace_back(run);
            }
            run();
            for (auto &t : pool)
            {
                t.join();
            }
        }
        else
        {
            run();
        }
#else
        (void)workers;
        run();
#endif
        if (failure)
        {
            std::rethrow_exception(failure);
        }
    }
}
