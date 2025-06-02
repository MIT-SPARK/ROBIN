// Copyright (c) 2020, Massachusetts Institute of Technology,
// Cambridge, MA 02139
// All Rights Reserved
// Authors: Jingnan Shi, et al. (see THANKS for the full author list)
// See LICENSE for the license information

#include <robin/pkc.hpp>

#include <iostream>
#include <vector>
#include <atomic>
#include <omp.h>
#include <algorithm>
#include <cmath>

int free_graph(robin::pkc::graph_t* g) {
  if (g->num_edges != nullptr) {
    free(g->num_edges);
  }

  if (g->adj != nullptr) {
    free(g->adj);
  }

  return 0;
}

/**
 * @brief K-core decomposition with parallel PKC. This particular function is suitable for large
 * graph, with optimization detailed in the original paper.
 *
 * @param g
 * @param deg
 */
static void PKC_optimized(const robin::IGraph& g, std::atomic<size_t>* deg) {
  double frac = 0.98;
  size_t MAX_NUM_THREADS = omp_get_max_threads();
  long n = g.VertexCount();
  long reduceN = 0;
  std::atomic<long> Size(0);
  std::atomic<long> visited(0);

  std::atomic<size_t>* newDeg = NULL;
  size_t* mapIndexToVtx = NULL;
  robin::pkc::vid_t* vertexToIndex = (robin::pkc::vid_t*)malloc(n * sizeof(robin::pkc::vid_t));
  std::vector<size_t> cumNumEdges(MAX_NUM_THREADS);
  long part = 0;
  robin::pkc::graph_t g_small;
  g_small.num_edges = NULL;
  g_small.adj = NULL;

#pragma omp parallel default(none) shared(std::cout, std::cerr, reduceN, Size, frac, part, cumNumEdges, deg, newDeg, visited, g, g_small, vertexToIndex, mapIndexToVtx) firstprivate(n)
  {
    int NUM_THREADS = omp_get_num_threads();
    int level = 0;
    int tid = omp_get_thread_num();

    std::vector<robin::pkc::vid_t> buff_vec;
    //vid_t* buff = (vid_t*)malloc((n * sizeof(vid_t)) / NUM_THREADS);
    //assert(buff != NULL);

    size_t start = 0, end = 0;
    int useSmallQ = 0;

#pragma omp for schedule(static)
    for (long i = 0; i < n; i++) {
      deg[i] = (g.GetVertexDegree(i));
    }

    while (visited.load(std::memory_order_relaxed) < n) {

      if ((useSmallQ == 0) && (visited.load(std::memory_order_relaxed) >= (size_t)( static_cast<double>(n) * frac))) {

        useSmallQ = 1;
        if (tid == 0) {
          reduceN = n - visited.load(std::memory_order_relaxed);
          newDeg = new std::atomic<size_t> [reduceN * sizeof(size_t)];
          mapIndexToVtx = (size_t*)malloc(reduceN * sizeof(size_t));
          g_small.n = reduceN;
          g_small.num_edges = (robin::pkc::vid_t*)malloc((reduceN + 1) * sizeof(robin::pkc::vid_t));
          g_small.num_edges[0] = 0;

          part = reduceN / NUM_THREADS;
        }
#pragma omp barrier

#pragma omp for schedule(static)
        for (long i = 0; i < n; i++) {
          if (deg[i] >= level) {
            buff_vec.push_back(i);
            //buff[end] = i;
            end++;
          }
        }

        // Now add them atomically
        size_t begin = Size.fetch_add(end, std::memory_order_relaxed);

        for (size_t i = 0; i < end; i++) {
          //newDeg[begin + i] = deg[buff[i]];
          //mapIndexToVtx[begin + i] = buff[i];
          //vertexToIndex[buff[i]] = begin + i;
          newDeg[begin + i].store(deg[buff_vec[i]].load(std::memory_order_relaxed), std::memory_order_relaxed);
          mapIndexToVtx[begin + i] = buff_vec[i];
          vertexToIndex[buff_vec[i]] = begin + i;
        }

        end = 0;
        buff_vec.clear();
#pragma omp barrier

        // Make a graph with reduceN vertices
        size_t edgeCount = 0;

        if (tid == (NUM_THREADS - 1)) {
          for (size_t i = tid * part; i < Size.load(std::memory_order_relaxed); i++) {
            size_t v = mapIndexToVtx[i];
            size_t prevEdgeCount = edgeCount;
            size_t v_edge_count = g.GetVertexDegree(v);
            for (size_t e_id = 0; e_id < v_edge_count; ++e_id) {
              auto u = g.GetVertexEdge(v, e_id);
              if (deg[u] >= level) {
                edgeCount++;
              }
            }
            //for (const auto& u : g.GetVertexEdges(v)) {
            //  if (deg[u] >= level) {
            //    edgeCount++;
            //  }
            //}
            // for(eid_t j = g->num_edges[v]; j < g->num_edges[v+1]; j++) {
            //  if( deg[ g->adj[j] ] >= level ) {
            //    edgeCount ++;
            //  }
            //}
            g_small.num_edges[i] = prevEdgeCount;
          }
          cumNumEdges[tid] = edgeCount;
        } else {

          for (size_t i = tid * part; i < (tid + 1) * part; i++) {
            size_t v = mapIndexToVtx[i];
            size_t prevEdgeCount = edgeCount;
            size_t v_edge_count = g.GetVertexDegree(v);
            for (size_t e_id = 0; e_id < v_edge_count; ++e_id) {
              auto u = g.GetVertexEdge(v, e_id);
              if (deg[u] >= level) {
                edgeCount++;
              }
            }
            //for (const auto& u : g.GetVertexEdges(v)) {
            //  if (deg[u] >= level) {
            //    edgeCount++;
            //  }
            //}
            // for(eid_t j = g->num_edges[v]; j < g->num_edges[v+1]; j++) {
            //  if( deg[ g->adj[j] ] >= level ) {
            //    edgeCount ++;
            //  }
            //}
            g_small.num_edges[i] = prevEdgeCount;
          }
          cumNumEdges[tid] = edgeCount;
        }

#pragma omp barrier
        if (tid == 0) {
          size_t start = cumNumEdges[0];
          for (int i = 1; i < NUM_THREADS; i++) {
            size_t prevEdgeCount = start;
            start = start + cumNumEdges[i];
            cumNumEdges[i] = prevEdgeCount;
          }
          g_small.m = start;
          g_small.num_edges[Size.load(std::memory_order_relaxed)] = start;

          g_small.adj = (robin::pkc::eid_t*)malloc(g_small.m * sizeof(robin::pkc::eid_t));
          cumNumEdges[0] = 0;
        }
#pragma omp barrier
        if (tid == (NUM_THREADS - 1)) {

          for (size_t i = tid * part; i < Size.load(std::memory_order_relaxed); i++) {
            g_small.num_edges[i] = g_small.num_edges[i] + cumNumEdges[tid];

            size_t v = mapIndexToVtx[i];
            size_t v_edge_count = g.GetVertexDegree(v);
            for (size_t e_id = 0; e_id < v_edge_count; ++e_id) {
              auto u = g.GetVertexEdge(v, e_id);
              if (deg[u] >= level) {
                g_small.adj[g_small.num_edges[i]] = vertexToIndex[u];
                g_small.num_edges[i]++;
              }
            }
            //for (const auto& u : g.GetVertexEdges(v)) {
            //  if (deg[u] >= level) {
            //    g_small.adj[g_small.num_edges[i]] = vertexToIndex[u];
            //    g_small.num_edges[i]++;
            //  }
            //}
            // for(eid_t j = g->num_edges[v]; j < g->num_edges[v+1]; j++) {
            //  if( deg[ g->adj[j] ] >= level ) {
            //    g_small.adj[ g_small.num_edges[i] ] = vertexToIndex[ g->adj[j] ];
            //    g_small.num_edges[i]++;
            //  }
            //}
          }
        } else {
          for (size_t i = tid * part; i < (tid + 1) * part; i++) {
            g_small.num_edges[i] = g_small.num_edges[i] + cumNumEdges[tid];

            size_t v = mapIndexToVtx[i];
            size_t v_edge_count = g.GetVertexDegree(v);
            for (size_t e_id = 0; e_id < v_edge_count; ++e_id) {
              auto u = g.GetVertexEdge(v, e_id);
              if (deg[u] >= level) {
                g_small.adj[g_small.num_edges[i]] = vertexToIndex[u];
                g_small.num_edges[i]++;
              }
            }
            //for (const auto& u : g.GetVertexEdges(v)) {
            //  if (deg[u] >= level) {
            //    g_small.adj[g_small.num_edges[i]] = vertexToIndex[u];
            //    g_small.num_edges[i]++;
            //  }
            //}
            // for(eid_t j = g->num_edges[v]; j < g->num_edges[v+1]; j++) {
            //  if( deg[ g->adj[j] ] >= level ) {
            //    g_small.adj[ g_small.num_edges[i] ] = vertexToIndex[ g->adj[j] ];
            //    g_small.num_edges[i]++;
            //  }
            //}
          }
        }
#pragma omp barrier
        // Now fix num_edges array
        if (tid == NUM_THREADS - 1) {

          for (long i = (Size.load(std::memory_order_relaxed) - 1); i >= (tid * part + 1); i--) {
            g_small.num_edges[i] = g_small.num_edges[i - 1];
          }
          g_small.num_edges[tid * part] = cumNumEdges[tid];

        } else {

          for (long i = ((tid + 1) * part - 1); i >= (tid * part + 1); i--) {
            g_small.num_edges[i] = g_small.num_edges[i - 1];
          }
          g_small.num_edges[tid * part] = cumNumEdges[tid];
        }
#pragma omp barrier
      }

      if (useSmallQ == 0) {

#pragma omp for schedule(static)
        for (long i = 0; i < n; i++) {
          if (deg[i] == level) {
            //buff[end] = i;
            buff_vec.push_back(i);
            end++;
          }
        }

        // Get work from curr queue and also add work after the current size
        while (start < end) {

          //vid_t v = buff[start];
          robin::pkc::vid_t v = buff_vec[start];
          start++;

          size_t v_edge_count = g.GetVertexDegree(v);
          for (size_t e_id = 0; e_id < v_edge_count; ++e_id) {
            auto u = g.GetVertexEdge(v, e_id);
            size_t deg_u = deg[u];

            if (deg_u > level) {
              size_t du = deg[u].fetch_sub(1, std::memory_order_relaxed);

              if (du == (level + 1)) {
                //buff[end] = u;
                buff_vec.push_back(u);
                end++;
              }

              if (du <= level) {
                deg[u].fetch_add(1, std::memory_order_relaxed);
              }
            } // deg_u > level

          } // visit adjacencies
        }   // end of while loop
      } else {

#pragma omp for schedule(static)
        for (long i = 0; i < Size.load(std::memory_order_relaxed); i++) {
          if (newDeg[i] == level) {
            //buff[end] = i;
            buff_vec.push_back(i);
            end++;
          }
        }

        // Get work from curr queue and also add work after the current size
        while (start < end) {

          //vid_t v = buff[start];
          robin::pkc::vid_t v = buff_vec[start];
          start++;

          for (robin::pkc::eid_t j = g_small.num_edges[v]; j < g_small.num_edges[v + 1]; j++) {
            robin::pkc::vid_t u = g_small.adj[j];

            size_t deg_u = newDeg[u];

            if (deg_u > level) {
              size_t du = newDeg[u].fetch_sub(1, std::memory_order_relaxed);

              if (du == (level + 1)) {
                //buff[end] = u;
                buff_vec.push_back(u);
                end++;
              }

              if (du <= level) {
                newDeg[u].fetch_add(1, std::memory_order_relaxed);
              }
            } // deg_u > level
          }   // visit adjacencies
        }     // end of while loop
      }

      visited.fetch_add(end, std::memory_order_relaxed);

#pragma omp barrier
      start = 0;
      end = 0;
      level = level + 1;
      buff_vec.clear();

    } // end of #visited < n

    //free(buff);

    // copy core values from newDeg to deg
#pragma omp for schedule(static)
    for (long i = 0; i < Size.load(std::memory_order_relaxed); i++) {
      deg[mapIndexToVtx[i]].store(newDeg[i].load(std::memory_order_relaxed), std::memory_order_relaxed);
    }

  } //#end of parallel region

  delete[] newDeg;
  free(mapIndexToVtx);
  free(vertexToIndex);
  free_graph(&g_small);
}

/**
 * @brief K-core decomposition with parallel PKC. This is the non-optimized version.
 * Each thread has its own queue -- parallel PKC_org for k-core decomposition
 *
 * @param g
 * @param deg
 */
static void PKC_original(const robin::IGraph& g, std::atomic<size_t>* deg) {
  long n = g.VertexCount();
  std::atomic<size_t> visited(0);

#pragma omp parallel default(none) shared(std::cout, g, visited, deg) firstprivate(n)
  {
    int level = 0;

    // thread-local buffer
    std::vector<robin::pkc::vid_t> buff_vec;
    int NUM_THREADS = omp_get_num_threads();

    long start = 0, end = 0;

#pragma omp for schedule(static)
    for (long i = 0; i < n; i++) {
      deg[i] = g.GetVertexDegree(i);
    }

    while (visited.load(std::memory_order_relaxed) < n) {

#pragma omp for schedule(static)
      for (long i = 0; i < n; i++) {
        if (deg[i] == level) {
          buff_vec.push_back(i);
          end++;
        }
      }

      // Get work from curr queue and also add work after the current size
      while (start < end) {

        robin::pkc::vid_t v = buff_vec[start];
        start++;

        size_t v_edge_count = g.GetVertexDegree(v);
        for (size_t e_id = 0; e_id < v_edge_count; ++e_id) {
          auto u = g.GetVertexEdge(v, e_id);
          int deg_u = deg[u].load(std::memory_order_relaxed);

          if (deg_u > level) {
            // TODO: Consider https://en.cppreference.com/w/cpp/atomic/atomic/fetch_sub
            size_t du = deg[u].fetch_sub(1, std::memory_order_relaxed);

            if (du == (level + 1)) {
              buff_vec.push_back(u);
              end++;
            }
            if (du <= level) {
              deg[u].fetch_add(1, std::memory_order_relaxed);
            }
          } // deg_u > level

        } // visit adjacencies
      }   // end of while loop

      visited.fetch_add(end, std::memory_order_relaxed);

#pragma omp barrier
      start = 0;
      end = 0;
      buff_vec.clear();
      level = level + 1;

    } // end of #visited < n

  } //#end of parallel region
}

void robin::pkc::BZ_kCores(const IGraph& g, std::vector<size_t>* deg) {
  long n = g.VertexCount();

  // Two arrays of size n
  unsigned int* vert = (unsigned int*)malloc(n * sizeof(unsigned int));
  assert(vert != NULL);

  unsigned int* pos = (unsigned int*)malloc(n * sizeof(unsigned int));
  assert(pos != NULL);

  // Maximum degree
  int maxDeg = 0;

  // deg[i] -- contains degree of vertex i
  for (size_t i = 0; i < n; i++) {
    (*deg)[i] = g.GetVertexDegree(i);

    if ((*deg)[i] > maxDeg)
      maxDeg = (*deg)[i];
  }

  // Used for bin-sort -- is of size maxDeg + 1
  unsigned int* bin = (unsigned int*)calloc(maxDeg + 1, sizeof(unsigned int));

  // Sort the vertices by increasing degree using bin-sort
  // Count number of vertices with each degree in 0...maxDeg
  for (long i = 0; i < n; i++) {
    bin[(*deg)[i]]++;
  }

  unsigned int start = 0;
  for (int d = 0; d < maxDeg + 1; d++) {
    unsigned int num = bin[d];
    bin[d] = start; // Start index of vertex in vert with degree d
    start = start + num;
  }

  // Do bin-sort of the vertices
  // vert -- contains the vertices in sorted order of degree
  // pos -- contains the positon of a vertex in vert array
  for (long i = 0; i < n; i++) {
    pos[i] = bin[(*deg)[i]];
    vert[pos[i]] = i;
    bin[(*deg)[i]]++;
  }

  for (int d = maxDeg; d >= 1; d--) {
    bin[d] = bin[d - 1];
  }
  bin[0] = 0;

  // kcores computation
  for (long i = 0; i < n; i++) {
    // Process the vertices in increasing order of degree
    robin::pkc::vid_t v = vert[i];
    size_t v_edge_count = g.GetVertexDegree(v);
    for (size_t e_id = 0; e_id < v_edge_count; ++e_id) {
      auto u = g.GetVertexEdge(v, e_id);
      if ((*deg)[u] > (*deg)[v]) {

        // Swap u with the first vertex in bin[du]
        unsigned int du = (*deg)[u];
        unsigned int pu = pos[u];
        unsigned int pw = bin[du];
        unsigned int w = vert[pw];

        if (u != w) {
          pos[u] = pw;
          vert[pu] = w;
          pos[w] = pu;
          vert[pw] = u;
        }

        // Increase the starting index of bin du
        bin[du]++;

        // Decrease degree of u -- so u is in previous bin now
        (*deg)[u]--;
      }
    }
  }

  free(vert);
  free(pos);
  free(bin);
}

void robin::pkc::PKC_original_serial(const IGraph& g, std::vector<size_t>* deg) {
  long n = g.VertexCount();
  long visited = 0;

  int level = 0;

  robin::pkc::vid_t* buff = (robin::pkc::vid_t*)malloc(n * sizeof(robin::pkc::vid_t));
  assert(buff != NULL);

  long start = 0, end = 0;
  for (long i = 0; i < n; i++) {
    (*deg)[i] = (g.GetVertexDegree(i));
  }

  while (visited < n) {

    for (long i = 0; i < n; i++) {

      if ((*deg)[i] == level) {
        buff[end] = i;
        end++;
      }
    }

    while (start < end) {

      robin::pkc::vid_t v = buff[start];
      start++;

      // Check the adj list of vertex v
      size_t v_edge_count = g.GetVertexDegree(v);
      for (size_t e_id = 0; e_id < v_edge_count; ++e_id) {
        auto u = g.GetVertexEdge(v, e_id);
        if ((*deg)[u] > level) {
          (*deg)[u] = (*deg)[u] - 1;

          if ((*deg)[u] == level) {

            buff[end] = u;
            end++;
          }

        } // deg_u > level

      } // visit adjacencies
    }   // end of while loop

    visited = visited + end;

    start = 0;
    end = 0;
    level = level + 1;

  } // end of #visited < n

  free(buff);
}

void robin::pkc::PKC_parallel(const IGraph& g, std::vector<size_t>* core,
                              bool use_optimized) {
  assert(core != NULL);
  std::atomic<size_t>* temp_core = new std::atomic<size_t>[g.VertexCount()];
  assert(temp_core != NULL);
  if (use_optimized) {
    PKC_optimized(g, temp_core);
  } else {
    PKC_original(g, temp_core);
  }
  core->resize(g.VertexCount());
  for (size_t i = 0; i < core->size(); ++i) {
    (*core)[i] = temp_core[i].load(std::memory_order_relaxed);;
  }
  delete[] temp_core;
}