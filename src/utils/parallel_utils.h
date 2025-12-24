#ifndef PARALLEL_UTILS_H
#define PARALLEL_UTILS_H

#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>

namespace parallel {

// Get optimal number of threads
inline int getNumThreads(int requested = 0) {
#ifdef _OPENMP
    if (requested <= 0) {
        return omp_get_max_threads();
    }
    return std::min(requested, omp_get_max_threads());
#else
    return 1;
#endif
}

// Set number of threads
inline void setNumThreads(int n) {
#ifdef _OPENMP
    omp_set_num_threads(n);
#endif
}

} // namespace parallel

#endif // PARALLEL_UTILS_H

