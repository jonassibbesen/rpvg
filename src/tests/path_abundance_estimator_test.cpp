
#include "catch.hpp"

#include "../path_abundance_estimator.hpp"
#include "../utils.hpp"


TEST_CASE("Minimum path cover can be sampled") {
    
    auto path_abundance_estimator = MinimumPathAbundanceEstimator(1, 1, 1, 1);

    Eigen::ColMatrixXd read_path_noise_log_probs(4, 3);
	read_path_noise_log_probs << -8, 0, -1.5, 0, -0.01, 0, -2, -2, 0, 0, -7, -0.1;

	Eigen::RowVectorXui read_counts(1, 4);
	read_counts << 1, 3, 1, 5;

	REQUIRE(path_abundance_estimator.sampleMinimumPathCover(read_path_noise_log_probs, read_counts) == vector<uint32_t>({2, 1}));

	SECTION("Minimum path cover is random") {

    	auto path_abundance_estimator_2 = MinimumPathAbundanceEstimator(1, 1, 1, 10);
		REQUIRE(path_abundance_estimator_2.sampleMinimumPathCover(read_path_noise_log_probs, read_counts) == vector<uint32_t>({2, 0, 1}));
	}
}

