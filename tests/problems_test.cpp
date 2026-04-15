// Copyright (c) 2020, Massachusetts Institute of Technology,
// Cambridge, MA 02139
// All Rights Reserved
// Authors: Jingnan Shi, et al. (see THANKS for the full author list)
// See LICENSE for the license information

#include "catch.hpp"
#include "test_utils.hpp"

#include <algorithm>
#include <iterator>
#include <omp.h>
#include <set>

#include <Eigen/Core>

#include <robin/core.hpp>
#include <robin/graph.hpp>
#include <robin/problems.hpp>
#include <robin/utils.hpp>

TEST_CASE("So3Y ref() returns const reference") {
  std::vector<Eigen::Matrix3d> rotations = {Eigen::Matrix3d::Identity(),
                                             Eigen::Matrix3d::Identity() * 2};
  robin::So3Y y_set(rotations);
  const Eigen::Matrix3d& ref1 = y_set.ref(0);
  const Eigen::Matrix3d& ref2 = y_set.ref(0);
  REQUIRE(&ref1 == &ref2); // Same address -- true reference, not temporary
}

TEST_CASE("So3Y operator[] still works by value") {
  std::vector<Eigen::Matrix3d> rotations = {Eigen::Matrix3d::Identity()};
  robin::So3Y y_set(rotations);
  Eigen::Matrix3d val = y_set[0];
  REQUIRE(val.isApprox(Eigen::Matrix3d::Identity()));
}

