// Copyright (c) 2020, Massachusetts Institute of Technology,
// Cambridge, MA 02139
// All Rights Reserved
// Authors: Jingnan Shi, et al. (see THANKS for the full author list)
// See LICENSE for the license information

#include "catch.hpp"
#include "test_utils.hpp"

#include <algorithm>
#include <chrono>
#include <iterator>
#include <omp.h>
#include <set>

#include <Eigen/Core>

#include <robin/core.hpp>
#include <robin/graph.hpp>
#include <robin/problems.hpp>
#include <robin/robin.hpp>
#include <robin/utils.hpp>

#include <xenium/ramalhete_queue.hpp>
#include <xenium/reclamation/stamp_it.hpp>

TEST_CASE("concurrent queues") {
  SECTION("simple large queue") {
    using EdgeType = std::pair<size_t, size_t>;
    xenium::ramalhete_queue<std::unique_ptr<EdgeType>,
                            xenium::policy::reclaimer<xenium::reclamation::stamp_it>>
        edge_queue;

    size_t N = 500000;
    for (size_t i = 0; i < N; ++i) {
      edge_queue.push(std::make_unique<EdgeType>(0, 0));
    }

    std::unique_ptr<EdgeType> edge;
    while (edge_queue.try_pop(edge)) {
      REQUIRE(edge->first == 0);
      REQUIRE(edge->second == 0);
      N--;
    }
    REQUIRE(N == 0);
  }
}

TEST_CASE("large compatibility graph csr") {
  // an arbitrary transformation matrix
  Eigen::Matrix4d T;
  // clang-format off
  T << 9.96926560e-01,  6.68735757e-02, -4.06664421e-02, -1.15576939e-01,
      -6.61289946e-02, 9.97617877e-01,  1.94008687e-02, -3.87705398e-02,
      4.18675510e-02, -1.66517807e-02,  9.98977765e-01, 1.14874890e-01,
      0,              0,                0,              1;
  // clang-format on

  size_t N = 200;

  // random 3d points
  Eigen::Matrix<double, 3, Eigen::Dynamic> src_points =
      Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, N);
  Eigen::Matrix<double, 4, Eigen::Dynamic> src_h;
  src_h.resize(4, src_points.cols());
  src_h.topRows(3) = src_points;
  src_h.bottomRows(1) = Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(N);

  // apply transformation
  Eigen::Matrix<double, 4, Eigen::Dynamic> tgt_h = T * src_h;
  Eigen::Matrix<double, 3, Eigen::Dynamic> tgt_points = tgt_h.topRows(3);

  // noise bound
  double noise_bound = 0.1;

  // add noise
  Eigen::Matrix<double, 3, Eigen::Dynamic> noise =
      Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, N) * noise_bound;
  tgt_points = tgt_points + noise;

  //
  // Measurement 3D points set
  robin::Points3d measurements(tgt_points);

  // Model set
  // Note: for our purpose measurements and model points can be flipped
  robin::Points3d model(src_points);

  // Compatibility function
  robin::Points3dRegCompCheck comp_check(&model, noise_bound);

  // InvGraphConstructor
  robin::CompGraphConstructor<robin::Points3d, robin::Points3dRegCompCheck, 2> graph_constructor;
  graph_constructor.SetCompCheckFunction(&comp_check);
  graph_constructor.SetMeasurements(&measurements);

  // Build graph
  size_t trials = 10;

  double graph_construction_time = 0;
  double clique_finding_time = 0;
  double core_finding_time = 0;
  for (size_t i = 0; i < trials; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
    auto* g = graph_constructor.BuildCompGraph(robin::GraphsStorageType::CSR);
    auto t1 = std::chrono::high_resolution_clock::now();

    auto actual_max_clique_indices =
        robin::FindInlierStructure(g, robin::InlierGraphStructure::MAX_CLIQUE);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto actual_max_core_indices =
        robin::FindInlierStructure(g, robin::InlierGraphStructure::MAX_CORE);
    auto t3 = std::chrono::high_resolution_clock::now();

    graph_construction_time +=
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - start).count();
    clique_finding_time += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    core_finding_time += std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    delete g;
  }
  std::cout << "Graph time  (ms): " << graph_construction_time / trials << std::endl;
  std::cout << "Clique time (ms): " << clique_finding_time / trials << std::endl;
  std::cout << "Core time   (ms): " << core_finding_time / trials << std::endl;
}

TEST_CASE("large compatibility graph atomic csr") {
  // an arbitrary transformation matrix
  Eigen::Matrix4d T;
  // clang-format off
  T << 9.96926560e-01,  6.68735757e-02, -4.06664421e-02, -1.15576939e-01,
      -6.61289946e-02, 9.97617877e-01,  1.94008687e-02, -3.87705398e-02,
      4.18675510e-02, -1.66517807e-02,  9.98977765e-01, 1.14874890e-01,
      0,              0,                0,              1;
  // clang-format on

  size_t N = 200;

  // random 3d points
  Eigen::Matrix<double, 3, Eigen::Dynamic> src_points =
      Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, N);
  Eigen::Matrix<double, 4, Eigen::Dynamic> src_h;
  src_h.resize(4, src_points.cols());
  src_h.topRows(3) = src_points;
  src_h.bottomRows(1) = Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(N);

  // apply transformation
  Eigen::Matrix<double, 4, Eigen::Dynamic> tgt_h = T * src_h;
  Eigen::Matrix<double, 3, Eigen::Dynamic> tgt_points = tgt_h.topRows(3);

  // noise bound
  double noise_bound = 0.1;

  // add noise
  Eigen::Matrix<double, 3, Eigen::Dynamic> noise =
      Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, N) * noise_bound;
  tgt_points = tgt_points + noise;

  //
  // Measurement 3D points set
  robin::Points3d measurements(tgt_points);

  // Model set
  // Note: for our purpose measurements and model points can be flipped
  robin::Points3d model(src_points);

  // Compatibility function
  robin::Points3dRegCompCheck comp_check(&model, noise_bound);

  // InvGraphConstructor
  robin::CompGraphConstructor<robin::Points3d, robin::Points3dRegCompCheck, 2> graph_constructor;
  graph_constructor.SetCompCheckFunction(&comp_check);
  graph_constructor.SetMeasurements(&measurements);

  // Build graph
  size_t trials = 10;

  double graph_construction_time = 0;
  double clique_finding_time = 0;
  double core_finding_time = 0;
  for (size_t i = 0; i < trials; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
    auto* g = graph_constructor.BuildCompGraph(robin::GraphsStorageType::ATOMIC_CSR);
    auto t1 = std::chrono::high_resolution_clock::now();

    auto actual_max_clique_indices =
        robin::FindInlierStructure(g, robin::InlierGraphStructure::MAX_CLIQUE);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto actual_max_core_indices =
        robin::FindInlierStructure(g, robin::InlierGraphStructure::MAX_CORE);
    auto t3 = std::chrono::high_resolution_clock::now();

    graph_construction_time +=
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - start).count();
    clique_finding_time += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    core_finding_time += std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    delete g;
  }
  std::cout << "Graph time  (ms): " << graph_construction_time / trials << std::endl;
  std::cout << "Clique time (ms): " << clique_finding_time / trials << std::endl;
  std::cout << "Core time   (ms): " << core_finding_time / trials << std::endl;
}

TEST_CASE("large compatibility graph adj list") {
  // an arbitrary transformation matrix
  Eigen::Matrix4d T;
  // clang-format off
  T << 9.96926560e-01,  6.68735757e-02, -4.06664421e-02, -1.15576939e-01,
      -6.61289946e-02, 9.97617877e-01,  1.94008687e-02, -3.87705398e-02,
      4.18675510e-02, -1.66517807e-02,  9.98977765e-01, 1.14874890e-01,
      0,              0,                0,              1;
  // clang-format on

  size_t N = 200;

  // random 3d points
  Eigen::Matrix<double, 3, Eigen::Dynamic> src_points =
      Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, N);
  Eigen::Matrix<double, 4, Eigen::Dynamic> src_h;
  src_h.resize(4, src_points.cols());
  src_h.topRows(3) = src_points;
  src_h.bottomRows(1) = Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(N);

  // apply transformation
  Eigen::Matrix<double, 4, Eigen::Dynamic> tgt_h = T * src_h;
  Eigen::Matrix<double, 3, Eigen::Dynamic> tgt_points = tgt_h.topRows(3);

  // noise bound
  double noise_bound = 0.1;

  // add noise
  Eigen::Matrix<double, 3, Eigen::Dynamic> noise =
      Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, N) * noise_bound;
  tgt_points = tgt_points + noise;

  //
  // Measurement 3D points set
  robin::Points3d measurements(tgt_points);

  // Model set
  // Note: for our purpose measurements and model points can be flipped
  robin::Points3d model(src_points);

  // Compatibility function
  robin::Points3dRegCompCheck comp_check(&model, noise_bound);

  // InvGraphConstructor
  robin::CompGraphConstructor<robin::Points3d, robin::Points3dRegCompCheck, 2> graph_constructor;
  graph_constructor.SetCompCheckFunction(&comp_check);
  graph_constructor.SetMeasurements(&measurements);

  // Build graph
  size_t trials = 10;

  double graph_construction_time = 0;
  double clique_finding_time = 0;
  double core_finding_time = 0;
  for (size_t i = 0; i < trials; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
    auto* g = graph_constructor.BuildCompGraph(robin::GraphsStorageType::ADJ_LIST);
    auto t1 = std::chrono::high_resolution_clock::now();

    auto actual_max_clique_indices =
        robin::FindInlierStructure(g, robin::InlierGraphStructure::MAX_CLIQUE);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto actual_max_core_indices =
        robin::FindInlierStructure(g, robin::InlierGraphStructure::MAX_CORE);
    auto t3 = std::chrono::high_resolution_clock::now();

    graph_construction_time +=
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - start).count();
    clique_finding_time += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    core_finding_time += std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();

    delete g;
  }
  std::cout << "Graph time  (ms): " << graph_construction_time / trials << std::endl;
  std::cout << "Clique time (ms): " << clique_finding_time / trials << std::endl;
  std::cout << "Core time   (ms): " << core_finding_time / trials << std::endl;
}

TEST_CASE("graph construction methods produce same edge sets") {
  // Validates that stack-allocation refactor doesn't break any construction path.
  // Uses a small deterministic problem so results are reproducible.
  robin::test::VecAveragingProblemFixture fixture;
  Eigen::Vector2d exp_vector;
  exp_vector << 2, 2;
  double noise_bound = 0.1;
  size_t N = 5;
  auto measurements = fixture.GenerateMeasurements(exp_vector, N, 0, noise_bound);
  robin::VectorY measurements_set(measurements);
  robin::SvaCompCheck sva_comp_check(2 * noise_bound);

  robin::CompGraphConstructor<robin::VectorY, robin::SvaCompCheck, 2> graph_constructor;
  graph_constructor.SetCompCheckFunction(&sva_comp_check);
  graph_constructor.SetMeasurements(&measurements_set);

  // Helper: extract sorted edge set from an AdjListGraph
  auto extract_edges = [](robin::AdjListGraph& g) {
    std::set<std::pair<size_t, size_t>> edges;
    for (size_t i = 0; i < g.VertexCount(); ++i) {
      for (size_t k = 0; k < g.GetVertexDegree(i); ++k) {
        size_t j = g.GetVertexEdge(i, k);
        // store each edge with smaller index first
        edges.insert(i < j ? std::make_pair(i, j) : std::make_pair(j, i));
      }
    }
    return edges;
  };

  SECTION("serial vs parallel_edge_buffer vs parallel_vertex") {
    robin::AdjListGraph g_serial;
    graph_constructor.BuildCompGraph_serial(&g_serial);
    auto edges_serial = extract_edges(g_serial);

    robin::AdjListGraph g_edge_buf;
    graph_constructor.BuildCompGraph_parallel_edge_buffer(&g_edge_buf);
    auto edges_edge_buf = extract_edges(g_edge_buf);

    robin::AdjListGraph g_vertex;
    graph_constructor.BuildCompGraph_parallel_vertex(&g_vertex);
    auto edges_vertex = extract_edges(g_vertex);

    // All methods must produce the same vertex and edge counts
    REQUIRE(g_serial.VertexCount() == N);
    REQUIRE(g_edge_buf.VertexCount() == N);
    REQUIRE(g_vertex.VertexCount() == N);

    REQUIRE(edges_serial.size() == edges_edge_buf.size());
    REQUIRE(edges_serial.size() == edges_vertex.size());

    // All methods must produce identical edge sets
    REQUIRE(edges_serial == edges_edge_buf);
    REQUIRE(edges_serial == edges_vertex);

    // For N=5 with no outliers and tight noise bound, expect complete graph
    REQUIRE(edges_serial.size() == N * (N - 1) / 2);
  }
}

TEST_CASE("sample comp graph construction") {
  SECTION("2d vector averaging") {
    //
    // Prepare test data
    //
    size_t N = 10;
    Eigen::Vector2d exp_vector;
    exp_vector << 2, 2;

    // noises between -0.1,0.1
    Eigen::Matrix2Xd random_noises = Eigen::Matrix2Xd::Random(2, N);
    random_noises.colwise().normalize();
    random_noises /= 10;
    double noise_bound = 0.1;

    // manually create outliers
    random_noises.col(0) *= 10;
    random_noises.col(9) *= 20;
    std::vector<size_t> expected_inliers_sizet = {1, 2, 3, 4, 5, 6, 7, 8};

    // Measurements set
    class H : public robin::SetBase<Eigen::Vector2d> {
    public:
      Eigen::Vector2d operator[](size_t i) const override { return h_set.col(i); }
      size_t size() const override { return h_set.cols(); }
      Eigen::Matrix2Xd h_set;
    };
    H H_set;

    // actual measurements
    H_set.h_set = Eigen::Matrix2Xd::Zero(2, N);
    for (size_t i = 0; i < random_noises.cols(); ++i) {
      H_set.h_set.col(i) = exp_vector + random_noises.col(i);
    }

    // Compatibility functor
    struct Vec2dComp {
      Vec2dComp(double thres) : noise_threshold(thres) {}
      bool operator()(H* measurements, size_t* subset_indices) {
        return ((*measurements)[subset_indices[0]] - (*measurements)[subset_indices[1]]).norm() <=
               noise_threshold;
      }
      double noise_threshold = 0;
    };
    Vec2dComp compatibility_func(2 * noise_bound);

    // Construct the invariant graph
    robin::CompGraphConstructor<H, Vec2dComp, 2> comp_graph_constructor;
    comp_graph_constructor.SetMeasurements(&H_set);
    comp_graph_constructor.SetCompCheckFunction(&compatibility_func);
    auto* g = comp_graph_constructor.BuildCompGraph(robin::GraphsStorageType::ADJ_LIST);

    // find k-core
    robin::KCoreDecompositionSolver k_core_decomposition_solver(
        robin::KCoreDecompositionSolver::KCORE_SOLVER_MODE::PKC_PARALLEL);
    k_core_decomposition_solver.Solve(*g);
    auto max_core = k_core_decomposition_solver.GetMaxKCore();
    std::sort(max_core.begin(), max_core.end());
    REQUIRE_THAT(max_core, Catch::Equals(expected_inliers_sizet));

    // find clique
    robin::MaxCliqueSolver::Params clique_params;
    clique_params.solver_mode = robin::MaxCliqueSolver::CLIQUE_SOLVER_MODE::PMC_EXACT;
    robin::MaxCliqueSolver clique_solver(clique_params);
    auto clique = clique_solver.FindMaxClique(*g);
    std::sort(clique.begin(), clique.end());
    REQUIRE_THAT(clique, Catch::Equals(expected_inliers_sizet));

    delete g;
  }
}

TEST_CASE("vertex parallel comp graph construction") {
  // test constructing comp graph in using atomic csr graph data structure
  SECTION("2d vector averaging no outliers") {
    robin::test::VecAveragingProblemFixture fixture;
    Eigen::Vector2d exp_vector;
    exp_vector << 2, 2;
    double noise_bound = 0.1;
    size_t N = 5;
    auto measurements = fixture.GenerateMeasurements(exp_vector, N, 0, noise_bound);
    robin::VectorY measurements_set(measurements);
    robin::SvaCompCheck sva_comp_check(2*noise_bound);

    // prepare comp graph constructor
    robin::CompGraphConstructor<robin::VectorY, robin::SvaCompCheck, 2> graph_constructor;
    graph_constructor.SetCompCheckFunction(&sva_comp_check);
    graph_constructor.SetMeasurements(&measurements_set);

    // construct the csr graph
    robin::AdjListGraph g;
    graph_constructor.BuildCompGraph_parallel_vertex(&g);

    // check the generated comp graph to be a complete graph
    REQUIRE(g.VertexCount() == N);
    REQUIRE(g.EdgeCount() == N*(N-1)/2);
    for (size_t i = 0; i < N; ++i) {
      auto c_deg = g.GetVertexDegree(i);
      std::set<size_t> c_edge_set;
      for (size_t k = 0; k < c_deg; ++k) {
        auto e = g.GetVertexEdge(i, k);
        c_edge_set.insert(e);
      }
      REQUIRE(c_edge_set.size() == (N-1));
    }
  }

  SECTION("2d vector averaging with outliers") {
    robin::test::VecAveragingProblemFixture fixture;
    Eigen::Vector2d exp_vector;
    exp_vector << 2, 2;
    double noise_bound = 0.1;
    size_t N = 5;
    auto measurements = fixture.GenerateMeasurements(exp_vector, N, 2, noise_bound);
    auto expected_outliers = fixture.ExpectedOutliers;
    robin::VectorY measurements_set(measurements);
    robin::SvaCompCheck sva_comp_check(2*noise_bound);

    // prepare comp graph constructor
    robin::CompGraphConstructor<robin::VectorY, robin::SvaCompCheck, 2> graph_constructor;
    graph_constructor.SetCompCheckFunction(&sva_comp_check);
    graph_constructor.SetMeasurements(&measurements_set);

    // construct the csr graph
    robin::AdjListGraph g;
    graph_constructor.BuildCompGraph_parallel_vertex(&g);

    // check the generated comp graph to be a complete graph
    REQUIRE(g.VertexCount() == N);
    for (size_t i = 0; i < N; ++i) {
      auto c_deg = g.GetVertexDegree(i);
      std::set<size_t> c_edge_set;
      for (size_t e_id = 0; e_id < c_deg; ++e_id) {
        auto e = g.GetVertexEdge(i, e_id);
        c_edge_set.insert(e);
      }
      for (const auto& m : expected_outliers) {
        REQUIRE(c_edge_set.count(m) == 0);
      }
    }
  }
}

/**
 * @brief Test constructing comp graph in using atomic csr graph data structure
 */
TEST_CASE("atomic csr comp graph construction") {
  SECTION("2d vector averaging no outliers") {
    robin::test::VecAveragingProblemFixture fixture;
    Eigen::Vector2d exp_vector;
    exp_vector << 2, 2;
    double noise_bound = 0.1;
    size_t N = 5;
    auto measurements = fixture.GenerateMeasurements(exp_vector, N, 0, noise_bound);
    robin::VectorY measurements_set(measurements);
    robin::SvaCompCheck sva_comp_check(2*noise_bound);

    // prepare comp graph constructor
    robin::CompGraphConstructor<robin::VectorY, robin::SvaCompCheck, 2> graph_constructor;
    graph_constructor.SetCompCheckFunction(&sva_comp_check);
    graph_constructor.SetMeasurements(&measurements_set);

    // construct the csr graph
    robin::AtomicCSRGraph g;
    graph_constructor.BuildCSRCompGraph_parallel_two_passes_atomic(&g);

    // check the generated comp graph to be a complete graph
    REQUIRE(g.VertexCount() == N);
    REQUIRE(g.EdgeCount() == N*(N-1)/2);
    for (size_t i = 0; i < N; ++i) {
      auto c_deg = g.GetVertexDegree(i);
      std::set<size_t> c_edge_set;
      for (size_t e_id = 0; e_id < c_deg; ++e_id) {
        auto e = g.GetVertexEdge(i, e_id);
        c_edge_set.insert(e);
      }
      REQUIRE(c_edge_set.size() == (N-1));
    }
  }

  SECTION("2d vector averaging with outliers") {
    robin::test::VecAveragingProblemFixture fixture;
    Eigen::Vector2d exp_vector;
    exp_vector << 2, 2;
    double noise_bound = 0.1;
    size_t N = 5;
    auto measurements = fixture.GenerateMeasurements(exp_vector, N, 2, noise_bound);
    auto expected_outliers = fixture.ExpectedOutliers;
    std::set<size_t> expected_outliers_set(expected_outliers.begin(), expected_outliers.end());

    robin::VectorY measurements_set(measurements);
    robin::SvaCompCheck sva_comp_check(2*noise_bound);

    // prepare comp graph constructor
    robin::CompGraphConstructor<robin::VectorY, robin::SvaCompCheck, 2> graph_constructor;
    graph_constructor.SetCompCheckFunction(&sva_comp_check);
    graph_constructor.SetMeasurements(&measurements_set);

    // construct the csr graph
    robin::AtomicCSRGraph g;
    graph_constructor.BuildCSRCompGraph_parallel_two_passes_atomic(&g);

    // check the generated comp graph to be a complete graph
    REQUIRE(g.VertexCount() == N);
    // 2 outliers, 3 inliers --> 3 clique with 3 edges
    REQUIRE(g.EdgeCount() == 3);
    for (size_t i = 0; i < N; ++i) {
      auto c_deg = g.GetVertexDegree(i);
      if (expected_outliers_set.count(i)) {
        REQUIRE(c_deg == 0);
      }
      std::set<size_t> c_edge_set;
      for (size_t e_id = 0; e_id < c_deg; ++e_id) {
        auto e = g.GetVertexEdge(i, e_id);
        c_edge_set.insert(e);
      }
      for (const auto& m : expected_outliers) {
        REQUIRE(c_edge_set.count(m) == 0);
      }
    }
  }
}