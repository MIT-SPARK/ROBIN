// Copyright (c) 2020, Massachusetts Institute of Technology,
// Cambridge, MA 02139
// All Rights Reserved
// Authors: Jingnan Shi, et al. (see THANKS for the full author list)
// See LICENSE for the license information

#include <robin/robin.hpp>

#include <nanoflann.hpp>

namespace {

/**
 * @brief Convert a directed half-edge adjacency list (where only edges i->j with j>i are stored)
 *        into a full undirected graph of the requested storage type.
 */
robin::IGraph* BuildGraphFromHalfEdges(std::vector<std::vector<size_t>>& half_adj_list,
                                       robin::GraphsStorageType graph_storage_type) {
  const size_t N = half_adj_list.size();

  size_t total_edges = 0;
  for (size_t i = 0; i < N; ++i) {
    total_edges += half_adj_list[i].size();
  }

  if (graph_storage_type == robin::GraphsStorageType::ADJ_LIST) {
    std::vector<std::vector<size_t>> full_adj_list(N);
    for (size_t i = 0; i < N; ++i) {
      for (size_t j : half_adj_list[i]) {
        full_adj_list[i].push_back(j);
        full_adj_list[j].push_back(i);
      }
    }
    return new robin::AdjListGraph(std::move(full_adj_list), total_edges);
  }

  // CSR and ATOMIC_CSR use edge-list constructors
  std::vector<std::pair<size_t, size_t>> edge_list;
  edge_list.reserve(total_edges);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j : half_adj_list[i]) {
      edge_list.emplace_back(i, j);
    }
  }

  if (graph_storage_type == robin::GraphsStorageType::CSR) {
    return new robin::CSRGraph(edge_list);
  }
  return new robin::AtomicCSRGraph(edge_list);
}

} // anonymous namespace

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
  const size_t N = measurements.cols();
  const size_t dim = measurements.rows();

  // Transpose to row-major (N x dim) for nanoflann
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> data_rowmajor =
      measurements.transpose();

  // Build KD-tree using nanoflann's Eigen matrix adaptor
  using kd_tree_t = nanoflann::KDTreeEigenMatrixAdaptor<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>,
      -1 /* dimensionality at runtime */,
      nanoflann::metric_L2>;

  kd_tree_t tree(dim, std::cref(data_rowmajor), {10 /* max leaf size */});

  // nanoflann radiusSearch uses squared L2 distance
  const double search_radius_sq = noise_bound * noise_bound;

  std::vector<std::vector<size_t>> half_adj_list(N);

  nanoflann::SearchParameters params;
  params.sorted = false;

#pragma omp parallel for schedule(dynamic, 64)
  for (size_t i = 0; i < N; ++i) {
    std::vector<nanoflann::ResultItem<size_t, double>> matches;
    tree.index->radiusSearch(data_rowmajor.row(i).data(), search_radius_sq, matches, params);

    for (const auto& match : matches) {
      if (match.first > i) {
        half_adj_list[i].push_back(match.first);
      }
    }
  }

  return BuildGraphFromHalfEdges(half_adj_list, graph_storage_type);
}

robin::IGraph* robin::MakeRotAvgInvGraph(const std::vector<Eigen::Matrix3d>& measurements,
                                         So3Distance dist_type, double noise_bound,
                                         GraphsStorageType graph_storage_type) {
  const size_t N = measurements.size();

  // Embed rotations as 9D vectors (flatten each 3x3 matrix)
  Eigen::Matrix<double, Eigen::Dynamic, 9, Eigen::RowMajor> embedded(N, 9);
  for (size_t i = 0; i < N; ++i) {
    Eigen::Map<const Eigen::Matrix<double, 1, 9, Eigen::RowMajor>> flat(measurements[i].data());
    embedded.row(i) = flat;
  }

  // Build KD-tree in R^9
  using kd_tree_t = nanoflann::KDTreeEigenMatrixAdaptor<
      Eigen::Matrix<double, Eigen::Dynamic, 9, Eigen::RowMajor>,
      9,
      nanoflann::metric_L2>;

  kd_tree_t tree(9, std::cref(embedded), {10 /* max leaf size */});

  // L2 in R^9 equals the Frobenius (chordal) distance.
  // For geodesic mode, use chordal upper bound as a pre-filter:
  // chordal = 2*sqrt(2)*sin(geodesic/2), then refine with exact geodesic check.
  double search_radius_sq;
  if (dist_type == So3Distance::GEODESIC_DISTANCE) {
    double chordal_ub = 2.0 * std::sqrt(2.0) * std::sin(noise_bound / 2.0);
    search_radius_sq = chordal_ub * chordal_ub;
  } else {
    search_radius_sq = noise_bound * noise_bound;
  }

  std::vector<std::vector<size_t>> half_adj_list(N);
  So3GeodesicDist geodesic_fn;

  nanoflann::SearchParameters params;
  params.sorted = false;

#pragma omp parallel for schedule(dynamic, 64)
  for (size_t i = 0; i < N; ++i) {
    std::vector<nanoflann::ResultItem<size_t, double>> matches;
    tree.index->radiusSearch(embedded.row(i).data(), search_radius_sq, matches, params);

    for (const auto& match : matches) {
      if (match.first > i) {
        bool keep = true;
        if (dist_type == So3Distance::GEODESIC_DISTANCE) {
          keep = (geodesic_fn(measurements[i], measurements[match.first]) <= noise_bound);
        }
        if (keep) {
          half_adj_list[i].push_back(match.first);
        }
      }
    }
  }

  return BuildGraphFromHalfEdges(half_adj_list, graph_storage_type);
}

robin::IGraph* robin::Make3dRegInvGraph(const Eigen::Matrix3Xd& src_3d_points,
                                        const Eigen::Matrix3Xd& dst_3d_points, double noise_bound,
                                        GraphsStorageType graph_storage_type) {
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

  // Measurements
  Points3dWithNormals measurements(dst_3d_points_with_normals);

  // Models
  Points3dWithNormals models(src_3d_points_with_normals);

  // Compatibility check function
  Eigen::Vector2d updated_noise_bound = noise_bound;
  updated_noise_bound(1) = std::cos(updated_noise_bound(1));
  PointsNormals3dRegCompCheck comp_check(&models, updated_noise_bound[0], updated_noise_bound[1]);

  // InvGraphConstructor
  robin::CompGraphConstructor<Points3dWithNormals, robin::PointsNormals3dRegCompCheck, 2>
      graph_constructor;
  graph_constructor.SetCompCheckFunction(&comp_check);
  graph_constructor.SetMeasurements(&measurements);

  // Build graph
  auto* g = graph_constructor.BuildCompGraph(graph_storage_type);

  return g;
}
