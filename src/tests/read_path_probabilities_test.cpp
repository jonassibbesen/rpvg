
#include "catch.hpp"

#include "gssw.h"

#include "../read_path_probabilities.hpp"
#include "../utils.hpp"

   
const double score_log_base = gssw_dna_recover_log_base(1, 4, 0.5, double_precision);

TEST_CASE("Read path probabilities can be calculated from alignment paths") {
    
	unordered_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}});
	FragmentLengthDist fragment_length_dist(10, 2);

	vector<AlignmentPath> alignment_paths(1, AlignmentPath(10, 10, 3, false, gbwt::SearchState()));
	auto alignment_path_ids = vector<vector<gbwt::size_type> >(1, vector<gbwt::size_type>({100, 200}));

	vector<PathInfo> paths(2);
	paths.front().effective_length = 3;
	paths.back().effective_length = 3;

	ReadPathProbabilities read_path_probs(1, 2, score_log_base, fragment_length_dist);
	read_path_probs.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, false);

	REQUIRE(read_path_probs.readCount() == 1);
	REQUIRE(doubleCompare(read_path_probs.noiseProbability(), 0.1));
	REQUIRE(read_path_probs.probabilities().size() == 2);
	REQUIRE(doubleCompare(read_path_probs.probabilities().front(), 0.45));
	REQUIRE(doubleCompare(read_path_probs.probabilities().back(), 0.45));

    SECTION("Improbable alignment path returns finite probabilities") {

    	alignment_paths.front().seq_length = 100000;

		ReadPathProbabilities read_path_probs_2(1, 2, score_log_base, fragment_length_dist);
		read_path_probs_2.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, false);

		REQUIRE(doubleCompare(read_path_probs_2.noiseProbability(), 0.1));
		REQUIRE(read_path_probs_2.probabilities().size() == 2);
		REQUIRE(read_path_probs.probabilities().front() - read_path_probs_2.probabilities().front() < pow(10, -8));
		REQUIRE(read_path_probs.probabilities().back() - read_path_probs_2.probabilities().back() < pow(10, -8));
	}

    SECTION("Probabilities are calculated from multiple alignment paths") {

		alignment_paths.emplace_back(AlignmentPath(15, 10, 5, false, gbwt::SearchState()));
		alignment_path_ids.emplace_back(vector<gbwt::size_type>({50}));
		
		clustered_path_index.emplace(10, 2);
		clustered_path_index.emplace(50, 3);

		paths.emplace_back(PathInfo());
		paths.back().effective_length = 3;

		paths.emplace_back(PathInfo());
		paths.back().effective_length = 3;

		ReadPathProbabilities read_path_probs_3(1, 4, score_log_base, fragment_length_dist);
		read_path_probs_3.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, false);

		REQUIRE(read_path_probs_3.readCount() == 1);
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
		read_path_probs_4.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, false);

		REQUIRE(read_path_probs_4.readCount() == 1);
		REQUIRE(doubleCompare(read_path_probs_4.noiseProbability(), 0.1));
		REQUIRE(read_path_probs_4.probabilities().size() == 2);
		REQUIRE(doubleCompare(read_path_probs_4.probabilities().front(), 0.3599999999999999));
		REQUIRE(doubleCompare(read_path_probs_4.probabilities().back(), 0.5400000000000000));
	}
}

TEST_CASE("Identical read path probabilities can be merged") {

	unordered_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}});
	FragmentLengthDist fragment_length_dist(10, 2);

	vector<AlignmentPath> alignment_paths(1, AlignmentPath(10, 10, 3, false, gbwt::SearchState()));
	auto alignment_path_ids = vector<vector<gbwt::size_type> >(1, vector<gbwt::size_type>({100, 200}));

	vector<PathInfo> paths(2);
	paths.front().effective_length = 3;
	paths.back().effective_length = 3;

	ReadPathProbabilities read_path_probs(1, 2, score_log_base, fragment_length_dist);
	read_path_probs.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, false);

	REQUIRE(read_path_probs.mergeIdenticalReadPathProbabilities(read_path_probs, 0.001));

	REQUIRE(read_path_probs.readCount() == 2);
	REQUIRE(doubleCompare(read_path_probs.noiseProbability(), 0.1));
	REQUIRE(read_path_probs.probabilities().size() == 2);
	REQUIRE(doubleCompare(read_path_probs.probabilities().front(), 0.45));
	REQUIRE(doubleCompare(read_path_probs.probabilities().back(), 0.45));

	SECTION("Probability precision affect merge") {

		vector<PathInfo> paths(2);
		paths.front().effective_length = 2;
		paths.back().effective_length = 3;

		ReadPathProbabilities read_path_probs_2(3, 2, score_log_base, fragment_length_dist);
		read_path_probs_2.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, false);

		REQUIRE(read_path_probs.mergeIdenticalReadPathProbabilities(read_path_probs_2, 0.1));
		REQUIRE(read_path_probs.readCount() == 5);

		REQUIRE(!read_path_probs.mergeIdenticalReadPathProbabilities(read_path_probs_2, 0.01));
		REQUIRE(read_path_probs.readCount() == 5);
	}
}

TEST_CASE("Read path probabilities can be collapsed") {

	
	unordered_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}, {10, 2}, {50, 3}});
	FragmentLengthDist fragment_length_dist(10, 2);

	vector<AlignmentPath> alignment_paths(1, AlignmentPath(10, 10, 3, false, gbwt::SearchState()));
	alignment_paths.emplace_back(AlignmentPath(15, 10, 5, false, gbwt::SearchState()));

	auto alignment_path_ids = vector<vector<gbwt::size_type> >(1, vector<gbwt::size_type>({100, 200}));
	alignment_path_ids.emplace_back(vector<gbwt::size_type>({50}));

	vector<PathInfo> paths(4);
	paths.at(0).effective_length = 3;
	paths.at(1).effective_length = 3;
	paths.at(2).effective_length = 3;
	paths.at(3).effective_length = 3;

	ReadPathProbabilities read_path_probs(1, 4, score_log_base, fragment_length_dist);
	read_path_probs.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, false);

	auto collapsed_probs = read_path_probs.collapsedProbabilities(0.01);
	REQUIRE(collapsed_probs.size() == 3);

	REQUIRE(doubleCompare(collapsed_probs.at(0).first, 0));
	REQUIRE(doubleCompare(collapsed_probs.at(1).first, 0.233044027062348));
	REQUIRE(doubleCompare(collapsed_probs.at(2).first, 0.3334779864688257));

	REQUIRE(collapsed_probs.at(0).second == vector<uint32_t>({2}));
	REQUIRE(collapsed_probs.at(1).second == vector<uint32_t>({3}));
	REQUIRE(collapsed_probs.at(2).second == vector<uint32_t>({0, 1}));

	SECTION("Probability precision affect collapse") {

		auto collapsed_probs = read_path_probs.collapsedProbabilities(0.2);
		REQUIRE(collapsed_probs.size() == 2);

		REQUIRE(doubleCompare(collapsed_probs.at(0).first, 0));
		REQUIRE(doubleCompare(collapsed_probs.at(1).first, 0.3334779864688257));

		REQUIRE(collapsed_probs.at(0).second == vector<uint32_t>({2}));
		REQUIRE(collapsed_probs.at(1).second == vector<uint32_t>({0, 1, 3}));
	}
}

