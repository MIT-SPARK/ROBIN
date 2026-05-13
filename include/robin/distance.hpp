// Copyright (c) 2020, Massachusetts Institute of Technology,
// Cambridge, MA 02139
// All Rights Reserved
// Authors: Jingnan Shi, et al. (see THANKS for the full author list)
// See LICENSE for the license information

#pragma once

#include <cmath>
#include <cstddef>

#include <omp.h>

#include <robin/math.hpp>

namespace robin {

/**
 * @brief Precompute all N*(N-1)/2 pairwise Euclidean distances for a set of points.
 * Points are stored column-major in a (dim x N) matrix as raw double*.
 * Output is float* for SIMD throughput (8 floats vs 4 doubles per 256-bit register).
 */
inline void precompute_pairwise_euclidean_f(const double* points, size_t dim, size_t N, float* out) {
  size_t num_pairs = N * (N - 1) / 2;
#pragma omp parallel for schedule(static)
  for (size_t k = 0; k < num_pairs; ++k) {
    auto [i, j] = pair_from_linear_index(k, N);
    const double* pi = points + i * dim;
    const double* pj = points + j * dim;
    double sum_sq = 0.0;
    for (size_t d = 0; d < dim; ++d) {
      double diff = pi[d] - pj[d];
      sum_sq += diff * diff;
    }
    out[k] = static_cast<float>(std::sqrt(sum_sq));
  }
}

/**
 * @brief Precompute all N*(N-1)/2 pairwise cosine similarities.
 * Vectors are stored column-major in a (dim x N) matrix as raw double*.
 */
inline void precompute_pairwise_cosine_f(const double* vectors, size_t dim, size_t N, float* out) {
  size_t num_pairs = N * (N - 1) / 2;
#pragma omp parallel for schedule(static)
  for (size_t k = 0; k < num_pairs; ++k) {
    auto [i, j] = pair_from_linear_index(k, N);
    const double* vi = vectors + i * dim;
    const double* vj = vectors + j * dim;
    double dot = 0.0, norm_i = 0.0, norm_j = 0.0;
    for (size_t d = 0; d < dim; ++d) {
      dot += vi[d] * vj[d];
      norm_i += vi[d] * vi[d];
      norm_j += vj[d] * vj[d];
    }
    double cosine = dot / (std::sqrt(norm_i) * std::sqrt(norm_j));
    if (cosine > 1.0)
      cosine = 1.0;
    else if (cosine < -1.0)
      cosine = -1.0;
    out[k] = static_cast<float>(cosine);
  }
}

} // namespace robin
