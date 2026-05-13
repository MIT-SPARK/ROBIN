// Copyright (c) 2020, Massachusetts Institute of Technology,
// Cambridge, MA 02139
// All Rights Reserved
// Authors: Jingnan Shi, et al. (see THANKS for the full author list)
// See LICENSE for the license information

#include <cmath>

#include <robin/bitset.hpp>
#include <robin/distance.hpp>
#include <robin/robin.hpp>

std::vector<size_t> robin::FindInlierStructure(const IGraph* g,
                                               InlierGraphStructure graph_structure) {
  // identify inlier structures
  switch (graph_structure) {
  case InlierGraphStructure::MAX_CORE: {
    KCoreDecompositionSolver k_core_decomposition_solver(
        robin::KCoreDecompositionSolver::KCORE_SOLVER_MODE::BZ_SERIAL);
    k_core_decomposition_solver.Solve(*g);
    return k_core_decomposition_solver.GetMaxKCore();
  }
  case InlierGraphStructure::MAX_CLIQUE: {
    robin::MaxCliqueSolver::Params clique_params;
    clique_params.solver_mode = robin::MaxCliqueSolver::CLIQUE_SOLVER_MODE::PMC_EXACT;
    robin::MaxCliqueSolver clique_solver(clique_params);
    return clique_solver.FindMaxClique(*g);
  }
  }
}

robin::IGraph* robin::MakeVecAvgInvGraph(const Eigen::MatrixXd& measurements, double noise_bound,
                                         GraphsStorageType graph_storage_type) {
  // Measurements set
  VectorY y_set(measurements);

  // Compatibility check function
  robin::SvaCompCheck sva_comp_check(noise_bound);

  // Compatibility graph constructor
  robin::CompGraphConstructor<VectorY, robin::SvaCompCheck, 2> graph_constructor;
  graph_constructor.SetCompCheckFunction(&sva_comp_check);
  graph_constructor.SetMeasurements(&y_set);

  // InvGraphConstructor
  auto* g = graph_constructor.BuildCompGraph(graph_storage_type);
  return g;
}

robin::IGraph* robin::MakeRotAvgInvGraph(const std::vector<Eigen::Matrix3d>& measurements,
                                         So3Distance dist_type, double noise_bound,
                                         GraphsStorageType graph_storage_type) {
  // Measurements sets
  So3Y y_set(measurements);

  // InvFunc
  switch (dist_type) {
  case So3Distance::CHORDAL_DISTANCE: {
    // Distance function
    So3ChordalDist inv_func;

    // Compatibility check function
    robin::SraCompCheck<So3ChordalDist> comp_check(&inv_func, noise_bound);

    // InvGraphConstructor
    robin::CompGraphConstructor<So3Y, robin::SraCompCheck<So3ChordalDist>, 2> graph_constructor;
    graph_constructor.SetCompCheckFunction(&comp_check);
    graph_constructor.SetMeasurements(&y_set);

    // Build graph
    auto* g = graph_constructor.BuildCompGraph(graph_storage_type);

    return g;
  }
  case robin::So3Distance::GEODESIC_DISTANCE: {
    // InvFunc
    So3GeodesicDist inv_func;

    // Compatibility check function
    robin::SraCompCheck<So3GeodesicDist> comp_check(&inv_func, noise_bound);

    // InvGraphConstructor
    robin::CompGraphConstructor<So3Y, robin::SraCompCheck<So3GeodesicDist>, 2> graph_constructor;
    graph_constructor.SetCompCheckFunction(&comp_check);
    graph_constructor.SetMeasurements(&y_set);

    // Build graph
    auto* g = graph_constructor.BuildCompGraph(graph_storage_type);

    return g;
  }
  }
}

robin::IGraph* robin::Make3dRegInvGraph(const Eigen::Matrix3Xd& src_3d_points,
                                        const Eigen::Matrix3Xd& dst_3d_points, double noise_bound,
                                        GraphsStorageType graph_storage_type) {
  // Try optimized precomputed path
  size_t N = dst_3d_points.cols();
  auto* g_precomputed =
      BuildCSRCompGraph_3dReg_precomputed(dst_3d_points.data(), src_3d_points.data(), N, noise_bound);
  if (g_precomputed != nullptr) {
    return g_precomputed;
  }

  // Fall back to existing path
  // Measurement 3D points set
  Points3d measurements(dst_3d_points);

  // Model set
  // Note: for our purpose measurements and model points can be flipped
  Points3d model(src_3d_points);

  // Compatibility function
  Points3dRegCompCheck comp_check(&model, noise_bound);

  // InvGraphConstructor
  robin::CompGraphConstructor<Points3d, robin::Points3dRegCompCheck, 2> graph_constructor;
  graph_constructor.SetCompCheckFunction(&comp_check);
  graph_constructor.SetMeasurements(&measurements);

  // Build graph
  auto* g = graph_constructor.BuildCompGraph(graph_storage_type);

  return g;
}

robin::IGraph* robin::Make3dNormalRegInvGraph(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& src_3d_points_with_normals,
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& dst_3d_points_with_normals,
    const Eigen::Vector2d& noise_bound, GraphsStorageType graph_storage_type) {
  assert(noise_bound.size() == 2);
  assert(src_3d_points_with_normals.cols() == dst_3d_points_with_normals.cols());

  size_t N = dst_3d_points_with_normals.cols();

  // Transform the normal noise bound the same way the existing path does
  double normal_noise_bound = std::cos(noise_bound(1));

  // Extract contiguous point and normal data from the 6xN matrices.
  // The 6xN matrix is column-major, so rows are interleaved. We need separate 3xN buffers.
  Eigen::Matrix3Xd dst_points = dst_3d_points_with_normals.topRows<3>();
  Eigen::Matrix3Xd dst_normals = dst_3d_points_with_normals.bottomRows<3>();
  Eigen::Matrix3Xd src_points = src_3d_points_with_normals.topRows<3>();
  Eigen::Matrix3Xd src_normals = src_3d_points_with_normals.bottomRows<3>();

  // Try optimized precomputed path
  auto* g_precomputed = BuildCSRCompGraph_3dNormalReg_precomputed(
      dst_points.data(), dst_normals.data(), src_points.data(), src_normals.data(), N,
      noise_bound[0], normal_noise_bound);
  if (g_precomputed != nullptr) {
    return g_precomputed;
  }

  // Fall back to existing path
  // Measurements
  Points3dWithNormals measurements(dst_3d_points_with_normals);

  // Models
  Points3dWithNormals models(src_3d_points_with_normals);

  // Compatibility check function
  PointsNormals3dRegCompCheck comp_check(&models, noise_bound[0], normal_noise_bound);

  // InvGraphConstructor
  robin::CompGraphConstructor<Points3dWithNormals, robin::PointsNormals3dRegCompCheck, 2>
      graph_constructor;
  graph_constructor.SetCompCheckFunction(&comp_check);
  graph_constructor.SetMeasurements(&measurements);

  // Build graph
  auto* g = graph_constructor.BuildCompGraph(graph_storage_type);

  return g;
}

//
// Precomputed distance pipeline implementations
//

robin::CSRGraph* robin::BuildCSRCompGraph_3dReg_precomputed(const double* src_points,
                                                            const double* model_points, size_t N,
                                                            double noise_bound) {
  size_t num_pairs = N * (N - 1) / 2;

  // Memory budget check: if > 512MB for the two distance arrays, fall back
  size_t mem_needed = num_pairs * sizeof(float) * 2;
  if (mem_needed > 512ULL * 1024 * 1024) {
    return nullptr;
  }

  // Precompute distances
  auto* src_dists = new float[num_pairs];
  auto* model_dists = new float[num_pairs];
  precompute_pairwise_euclidean_f(src_points, 3, N, src_dists);
  precompute_pairwise_euclidean_f(model_points, 3, N, model_dists);

  // Bitset scan
  size_t num_words = (num_pairs + 63) / 64;
  auto* edge_bits = new uint64_t[num_words]();
  float threshold_f = static_cast<float>(noise_bound);

#pragma omp parallel for schedule(static)
  for (size_t k = 0; k < num_pairs; ++k) {
    if (std::fabsf(src_dists[k] - model_dists[k]) <= threshold_f) {
      size_t word = k / 64;
      size_t bit = k % 64;
      __atomic_or_fetch(&edge_bits[word], 1ULL << bit, __ATOMIC_RELAXED);
    }
  }

  delete[] src_dists;
  delete[] model_dists;

  // Build CSR from bitset
  size_t* offsets = nullptr;
  size_t* edges = nullptr;
  size_t edge_count = 0;
  bitset_to_csr(edge_bits, N, num_pairs, &offsets, &edges, &edge_count);
  delete[] edge_bits;

  auto* g = new CSRGraph();
  g->SetEdges(edge_count, edges);
  g->SetOffsets(N, offsets);
  return g;
}

robin::CSRGraph* robin::BuildCSRCompGraph_3dNormalReg_precomputed(
    const double* src_points, const double* src_normals, const double* model_points,
    const double* model_normals, size_t N, double point_noise_bound, double normal_noise_bound) {
  size_t num_pairs = N * (N - 1) / 2;

  // Memory budget check: 4 float arrays (src_dists, model_dists, src_cosines, model_cosines)
  size_t mem_needed = num_pairs * sizeof(float) * 4;
  if (mem_needed > 512ULL * 1024 * 1024) {
    return nullptr;
  }

  // Precompute point distances
  auto* src_dists = new float[num_pairs];
  auto* model_dists = new float[num_pairs];
  precompute_pairwise_euclidean_f(src_points, 3, N, src_dists);
  precompute_pairwise_euclidean_f(model_points, 3, N, model_dists);

  // Precompute normal cosine similarities
  auto* src_cosines = new float[num_pairs];
  auto* model_cosines = new float[num_pairs];
  precompute_pairwise_cosine_f(src_normals, 3, N, src_cosines);
  precompute_pairwise_cosine_f(model_normals, 3, N, model_cosines);

  // Bitset scan with short-circuit: point check first, then normal check
  size_t num_words = (num_pairs + 63) / 64;
  auto* edge_bits = new uint64_t[num_words]();
  float pt_threshold_f = static_cast<float>(point_noise_bound);
  float normal_threshold_f = static_cast<float>(normal_noise_bound);

#pragma omp parallel for schedule(static)
  for (size_t k = 0; k < num_pairs; ++k) {
    // Point compatibility check
    if (std::fabsf(src_dists[k] - model_dists[k]) > pt_threshold_f) {
      continue;
    }
    // Normal compatibility check (same formula as PointsNormals3dRegCompCheck)
    float ca = src_cosines[k];
    float cb = model_cosines[k];
    float lhs = ca * cb + std::sqrtf((1.0f - ca * ca) * (1.0f - cb * cb));
    if (lhs < normal_threshold_f) {
      continue;
    }
    size_t word = k / 64;
    size_t bit = k % 64;
    __atomic_or_fetch(&edge_bits[word], 1ULL << bit, __ATOMIC_RELAXED);
  }

  delete[] src_dists;
  delete[] model_dists;
  delete[] src_cosines;
  delete[] model_cosines;

  // Build CSR from bitset
  size_t* offsets = nullptr;
  size_t* edges = nullptr;
  size_t edge_count = 0;
  bitset_to_csr(edge_bits, N, num_pairs, &offsets, &edges, &edge_count);
  delete[] edge_bits;

  auto* g = new CSRGraph();
  g->SetEdges(edge_count, edges);
  g->SetOffsets(N, offsets);
  return g;
}
