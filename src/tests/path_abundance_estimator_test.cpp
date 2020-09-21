
#include "catch.hpp"

#include "../path_abundance_estimator.hpp"
#include "../utils.hpp"


TEST_CASE("Weighted minimum path cover can be found") {
    
    auto path_abundance_estimator = MinimumPathAbundanceEstimator(1, 1, 1);

    Eigen::ColMatrixXb read_path_cover(4, 3);
	read_path_cover << 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1;

	Eigen::RowVectorXui read_counts(1, 4);
	read_counts << 1, 3, 1, 5;

	Eigen::RowVectorXd path_weights(1, 3);
	path_weights << 1, 1, 1;

	REQUIRE(path_abundance_estimator.weightedMinimumPathCover(read_path_cover, read_counts, path_weights) == vector<uint32_t>({0, 1}));

	SECTION("Weights influence minimum path cover") {

		path_weights(2) = 0.01;
		REQUIRE(path_abundance_estimator.weightedMinimumPathCover(read_path_cover, read_counts, path_weights) == vector<uint32_t>({0, 1, 2}));
	}
}

