
#include "catch.hpp"

#include "gssw.h"

#include "../read_path_probabilities.hpp"
#include "../utils.hpp"

   
const double score_log_base = gssw_dna_recover_log_base(1, 4, 0.5, double_precision);

TEST_CASE("ReadPathProbabilities can be equal") {
    
	ReadPathProbabilities read_path_probs_1(2, score_log_base);
	ReadPathProbabilities read_path_probs_2(2, score_log_base);

    REQUIRE(read_path_probs_1 == read_path_probs_2);	

    read_path_probs_1.noise_prob = 0.01;    
    REQUIRE(read_path_probs_1 != read_path_probs_2);

    read_path_probs_2.noise_prob = 0.01;    
    REQUIRE(read_path_probs_1 == read_path_probs_2);

    read_path_probs_1.read_path_probs.push_back(0.1);
    REQUIRE(read_path_probs_1 != read_path_probs_2);

    read_path_probs_2.read_path_probs.push_back(0.2);
    REQUIRE(read_path_probs_1 != read_path_probs_2);
    
    read_path_probs_2.read_path_probs.back() = 0.1;
    REQUIRE(read_path_probs_1 == read_path_probs_2);
}

TEST_CASE("Multiple ReadPathProbabilities can be sorted") {
    
	vector<ReadPathProbabilities> multiple_read_path_probs(2, ReadPathProbabilities(2, score_log_base));

	multiple_read_path_probs.front().read_path_probs.back() = 1;
	sort(multiple_read_path_probs.begin(), multiple_read_path_probs.end());

	REQUIRE(doubleCompare(multiple_read_path_probs.front().read_path_probs.back(), 0));
	REQUIRE(doubleCompare(multiple_read_path_probs.back().read_path_probs.back(), 1));

	multiple_read_path_probs.back().noise_prob = 0;
	sort(multiple_read_path_probs.begin(), multiple_read_path_probs.end());

	REQUIRE(doubleCompare(multiple_read_path_probs.front().read_path_probs.back(), 1));
	REQUIRE(doubleCompare(multiple_read_path_probs.back().read_path_probs.back(), 0));
}

TEST_CASE("Read path probabilities can be calculated from alignment paths") {
    
	vector<AlignmentPath> alignment_paths(1, AlignmentPath(10, 10, 3, vector<gbwt::size_type>({100, 200})));

	unordered_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}});
	FragmentLengthDist fragment_length_dist(10, 2);

	ReadPathProbabilities read_path_probs(2, score_log_base);
	read_path_probs.calcReadPathProbabilities(alignment_paths, clustered_path_index, fragment_length_dist, false);

	REQUIRE(doubleCompare(read_path_probs.noise_prob, 0.1));
	REQUIRE(read_path_probs.read_path_probs.size() == 2);
	REQUIRE(doubleCompare(read_path_probs.read_path_probs.front(), 0.45));
	REQUIRE(doubleCompare(read_path_probs.read_path_probs.back(), 0.45));

    SECTION("Extremely improbable alignment path returns finite probabilities") {

    	alignment_paths.front().seq_length = 10000;

		ReadPathProbabilities read_path_probs_2(2, score_log_base);
		read_path_probs_2.calcReadPathProbabilities(alignment_paths, clustered_path_index, fragment_length_dist, false);

		REQUIRE(read_path_probs == read_path_probs_2);
	}

    SECTION("Probabilities are calculated from multiple alignment paths") {

		alignment_paths.emplace_back(AlignmentPath(15, 10, 5, vector<gbwt::size_type>({50})));
		
		clustered_path_index.emplace(10, 2);
		clustered_path_index.emplace(50, 3);

		ReadPathProbabilities read_path_probs_3(4, score_log_base);
		read_path_probs_3.calcReadPathProbabilities(alignment_paths, clustered_path_index, fragment_length_dist, false);

		REQUIRE(doubleCompare(read_path_probs_3.noise_prob, 0.1));
		REQUIRE(read_path_probs_3.read_path_probs.size() == 4);
		REQUIRE(doubleCompare(read_path_probs_3.read_path_probs.at(0), 0.3334779864688257));
		REQUIRE(doubleCompare(read_path_probs_3.read_path_probs.at(1), 0.3334779864688257));
		REQUIRE(doubleCompare(read_path_probs_3.read_path_probs.at(2), 0));
		REQUIRE(doubleCompare(read_path_probs_3.read_path_probs.at(3), 0.2330440270623483));
	}

    SECTION("Positional probabilities are calculated from path lengths") {

		read_path_probs.addPositionalProbabilities(vector<double>({3, 2}));

		REQUIRE(doubleCompare(read_path_probs.noise_prob, 0.1));
		REQUIRE(read_path_probs.read_path_probs.size() == 2);
		REQUIRE(doubleCompare(read_path_probs.read_path_probs.front(), 0.3599999999999999));
		REQUIRE(doubleCompare(read_path_probs.read_path_probs.back(), 0.5400000000000000));
	}	
}

