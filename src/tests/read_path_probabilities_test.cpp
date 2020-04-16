
#include "catch.hpp"

#include "gssw.h"

#include "../read_path_probabilities.hpp"
#include "../utils.hpp"

   
const double score_log_base = gssw_dna_recover_log_base(1, 4, 0.5, double_precision);

TEST_CASE("Read path probabilities can be calculated from alignment paths") {
    
	vector<AlignmentPath> alignment_paths(1, AlignmentPath(10, 10, 3, vector<gbwt::size_type>({100, 200})));

	unordered_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}});
	FragmentLengthDist fragment_length_dist(10, 2);

	vector<Path> paths(2);
	paths.front().effective_length = 3;
	paths.back().effective_length = 3;

	ReadPathProbabilities read_path_probs(1, 2, score_log_base, fragment_length_dist);
	read_path_probs.calcReadPathProbabilities(alignment_paths, clustered_path_index, paths, false);

	REQUIRE(doubleCompare(read_path_probs.noiseProbability(), 0.1));
	REQUIRE(read_path_probs.probabilities().size() == 2);
	REQUIRE(doubleCompare(read_path_probs.probabilities().front(), 0.45));
	REQUIRE(doubleCompare(read_path_probs.probabilities().back(), 0.45));

    SECTION("Extremely improbable alignment path returns finite probabilities") {

    	alignment_paths.front().seq_length = 10000;

		ReadPathProbabilities read_path_probs_2(1, 2, score_log_base, fragment_length_dist);
		read_path_probs_2.calcReadPathProbabilities(alignment_paths, clustered_path_index, paths, false);

		REQUIRE(read_path_probs == read_path_probs_2);
	}

    SECTION("Probabilities are calculated from multiple alignment paths") {

		alignment_paths.emplace_back(AlignmentPath(15, 10, 5, vector<gbwt::size_type>({50})));
		
		clustered_path_index.emplace(10, 2);
		clustered_path_index.emplace(50, 3);

		paths.emplace_back(Path());
		paths.back().effective_length = 3;

		paths.emplace_back(Path());
		paths.back().effective_length = 3;

		ReadPathProbabilities read_path_probs_3(1, 4, score_log_base, fragment_length_dist);
		read_path_probs_3.calcReadPathProbabilities(alignment_paths, clustered_path_index, paths, false);

		REQUIRE(doubleCompare(read_path_probs_3.noiseProbability(), 0.1));
		REQUIRE(read_path_probs_3.probabilities().size() == 4);
		REQUIRE(doubleCompare(read_path_probs_3.probabilities().at(0), 0.3334779864688257));
		REQUIRE(doubleCompare(read_path_probs_3.probabilities().at(1), 0.3334779864688257));
		REQUIRE(doubleCompare(read_path_probs_3.probabilities().at(2), 0));
		REQUIRE(doubleCompare(read_path_probs_3.probabilities().at(3), 0.2330440270623483));
	}

    SECTION("Effective path lengths affect probabilities") {

		paths.back().effective_length = 2;

		ReadPathProbabilities read_path_probs_4(1, 2, score_log_base, fragment_length_dist);
		read_path_probs_4.calcReadPathProbabilities(alignment_paths, clustered_path_index, paths, false);

		REQUIRE(doubleCompare(read_path_probs_4.noiseProbability(), 0.1));
		REQUIRE(read_path_probs_4.probabilities().size() == 2);
		REQUIRE(doubleCompare(read_path_probs_4.probabilities().front(), 0.3599999999999999));
		REQUIRE(doubleCompare(read_path_probs_4.probabilities().back(), 0.5400000000000000));
	}
}

