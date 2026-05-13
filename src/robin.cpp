// Copyright (c) 2020, Massachusetts Institute of Technology,
// Cambridge, MA 02139
// All Rights Reserved
// Authors: Jingnan Shi, et al. (see THANKS for the full author list)
// See LICENSE for the license information

#include <robin/robin.hpp>

namespace {
std::vector<size_t> solve_on_graph(const robin::IGraph* g,
                                   robin::InlierGraphStructure graph_structure) {
  switch (graph_structure) {
  case robin::InlierGraphStructure::MAX_CORE: {
    robin::KCoreDecompositionSolver k_core_decomposition_solver(
        robin::KCoreDecompositionSolver::KCORE_SOLVER_MODE::BZ_SERIAL);
    k_core_decomposition_solver.Solve(*g);
    return k_core_decomposition_solver.GetMaxKCore();
  }
  case robin::InlierGraphStructure::MAX_CLIQUE: {
    robin::MaxCliqueSolver::Params clique_params;
    clique_params.solver_mode = robin::MaxCliqueSolver::CLIQUE_SOLVER_MODE::PMC_EXACT;
    robin::MaxCliqueSolver clique_solver(clique_params);
    return clique_solver.FindMaxClique(*g);
  }
  }
  return {};
}
} // anonymous namespace

std::vector<size_t> robin::FindInlierStructure(const IGraph* g,
                                               InlierGraphStructure graph_structure) {
  size_t N = g->VertexCount();

  // Identify non-isolated vertices
  std::vector<size_t> new_to_old;
  std::vector<size_t> old_to_new(N, SIZE_MAX);

  for (size_t i = 0; i < N; ++i) {
    if (g->GetVertexDegree(i) > 0) {
      old_to_new[i] = new_to_old.size();
      new_to_old.push_back(i);
    }
  }

  size_t N_compact = new_to_old.size();

  // If no vertices were filtered, run solver on original graph
  if (N_compact == N) {
    return solve_on_graph(g, graph_structure);
  }

  // All vertices isolated — no inliers
  if (N_compact == 0) {
    return {};
  }

  // Build compacted adjacency list with remapped vertex IDs
  std::vector<std::vector<size_t>> compact_adj(N_compact);
  size_t compact_edges = 0;
  for (size_t new_i = 0; new_i < N_compact; ++new_i) {
    size_t old_i = new_to_old[new_i];
    size_t deg = g->GetVertexDegree(old_i);
    compact_adj[new_i].reserve(deg);
    for (size_t e = 0; e < deg; ++e) {
      size_t old_j = g->GetVertexEdge(old_i, e);
      size_t new_j = old_to_new[old_j];
      compact_adj[new_i].push_back(new_j);
    }
    compact_edges += deg;
  }
  compact_edges /= 2; // undirected: each edge counted twice

  AdjListGraph compact_graph(std::move(compact_adj), compact_edges);

  // Run solver on compact graph
  auto compact_result = solve_on_graph(&compact_graph, graph_structure);

  // Remap back to original vertex IDs
  std::vector<size_t> result;
  result.reserve(compact_result.size());
  for (size_t new_id : compact_result) {
    result.push_back(new_to_old[new_id]);
  }
  return result;
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
