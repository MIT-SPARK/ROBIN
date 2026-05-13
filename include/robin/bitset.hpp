// Copyright (c) 2020, Massachusetts Institute of Technology,
// Cambridge, MA 02139
// All Rights Reserved
// Authors: Jingnan Shi, et al. (see THANKS for the full author list)
// See LICENSE for the license information

#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>

#include <omp.h>

#include <robin/math.hpp>

namespace robin {

inline void bitset_set(uint64_t* bits, size_t k) { bits[k / 64] |= (1ULL << (k % 64)); }

inline bool bitset_test(const uint64_t* bits, size_t k) {
  return (bits[k / 64] >> (k % 64)) & 1ULL;
}

/**
 * @brief Count set bits in a bitset using parallel popcount
 */
inline size_t bitset_popcount(const uint64_t* bits, size_t num_words) {
  size_t count = 0;
#pragma omp parallel for reduction(+ : count)
  for (size_t w = 0; w < num_words; ++w) {
    count += __builtin_popcountll(bits[w]);
  }
  return count;
}

/**
 * @brief Build CSR graph from edge bitset.
 * The bitset has one bit per pair (i,j) in linear order.
 * Builds undirected graph: each edge stored in both directions.
 * @param bits Edge bitset
 * @param N Number of vertices
 * @param num_pairs N*(N-1)/2
 * @param offsets_out Output: CSR offsets array (size N+1), caller must delete[]
 * @param edges_out Output: CSR edges array (size 2*edge_count), caller must delete[]
 * @param edge_count_out Output: number of undirected edges
 */
inline void bitset_to_csr(const uint64_t* bits, size_t N, size_t num_pairs, size_t** offsets_out,
                           size_t** edges_out, size_t* edge_count_out) {
  // Step 1: Compute per-vertex degree
  auto* degrees = new size_t[N]();
  for (size_t k = 0; k < num_pairs; ++k) {
    if (bitset_test(bits, k)) {
      auto [i, j] = pair_from_linear_index(k, N);
      degrees[i]++;
      degrees[j]++;
    }
  }

  // Step 2: Prefix sum to get offsets
  auto* offsets = new size_t[N + 1];
  offsets[0] = 0;
  for (size_t v = 0; v < N; ++v) {
    offsets[v + 1] = offsets[v] + degrees[v];
  }
  size_t total_entries = offsets[N]; // = 2 * edge_count
  *edge_count_out = total_entries / 2;

  // Step 3: Fill edges
  auto* edges = new size_t[total_entries];
  auto* pos = new size_t[N]();
  for (size_t k = 0; k < num_pairs; ++k) {
    if (bitset_test(bits, k)) {
      auto [i, j] = pair_from_linear_index(k, N);
      edges[offsets[i] + pos[i]++] = j;
      edges[offsets[j] + pos[j]++] = i;
    }
  }

  *offsets_out = offsets;
  *edges_out = edges;
  delete[] degrees;
  delete[] pos;
}

} // namespace robin
