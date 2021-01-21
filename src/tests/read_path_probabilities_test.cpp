
#include "catch.hpp"
#include "sparsepp/spp.h"

#include "../read_path_probabilities.hpp"
#include "../utils.hpp"

   
TEST_CASE("Read path probabilities can be calculated from alignment paths") {
    
	spp::sparse_hash_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}});
	FragmentLengthDist fragment_length_dist(10, 2);

	vector<AlignmentPath> alignment_paths(1, AlignmentPath(make_pair(gbwt::SearchState(), 0), false, 10, 10, 3));
	auto alignment_path_ids = vector<vector<gbwt::size_type> >(1, vector<gbwt::size_type>({100, 200}));

	vector<PathInfo> paths(2, PathInfo(""));
	paths.front().effective_length = 3;
	paths.back().effective_length = 3;

	ReadPathProbabilities read_path_probs(1, pow(10, -8));
	read_path_probs.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

	REQUIRE(read_path_probs.readCount() == 1);
	REQUIRE(Utils::doubleCompare(read_path_probs.noiseProbability(), 0.1));
	REQUIRE(read_path_probs.probabilities().size() == 2);

	REQUIRE(read_path_probs.probabilities().front().first == 0);
	REQUIRE(Utils::doubleCompare(read_path_probs.probabilities().front().second, 0.45));
	REQUIRE(read_path_probs.probabilities().back().first == 1);
	REQUIRE(Utils::doubleCompare(read_path_probs.probabilities().back().second, 0.45));

    SECTION("Improbable alignment path returns finite probabilities") {

    	alignment_paths.front().frag_length = 100000;

		ReadPathProbabilities read_path_probs_2(1, pow(10, -8));
		read_path_probs_2.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProbability(), 0.1));
		REQUIRE(read_path_probs_2.probabilities().size() == 2);

		REQUIRE(read_path_probs.probabilities().front().first == read_path_probs_2.probabilities().front().first);
		REQUIRE(read_path_probs.probabilities().front().second - read_path_probs_2.probabilities().front().second < pow(10, -8));
		REQUIRE(read_path_probs.probabilities().back().first == read_path_probs_2.probabilities().back().first);
		REQUIRE(read_path_probs.probabilities().back().second - read_path_probs_2.probabilities().back().second < pow(10, -8));
	}

    SECTION("Probabilities are calculated from multiple alignment paths") {

		alignment_paths.emplace_back(AlignmentPath(make_pair(gbwt::SearchState(), 0), false, 15, 10, 5));
		alignment_path_ids.emplace_back(vector<gbwt::size_type>({50}));
		
		clustered_path_index.emplace(10, 2);
		clustered_path_index.emplace(50, 3);

		paths.emplace_back(PathInfo(""));
		paths.back().effective_length = 3;

		paths.emplace_back(PathInfo(""));
		paths.back().effective_length = 3;

		ReadPathProbabilities read_path_probs_3(1, pow(10, -8));
		read_path_probs_3.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_3.readCount() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_3.noiseProbability(), 0.1));
		REQUIRE(read_path_probs_3.probabilities().size() == 3);

		REQUIRE(read_path_probs_3.probabilities().at(0).first == 0);
		REQUIRE(Utils::doubleCompare(read_path_probs_3.probabilities().at(0).second, 0.333477986468937));
		REQUIRE(read_path_probs_3.probabilities().at(1).first == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_3.probabilities().at(1).second, 0.333477986468937));
		REQUIRE(read_path_probs_3.probabilities().at(2).first == 3);
		REQUIRE(Utils::doubleCompare(read_path_probs_3.probabilities().at(2).second, 0.233044027062125));
	}

    SECTION("Effective path lengths affect probabilities") {

		paths.back().effective_length = 2;

		ReadPathProbabilities read_path_probs_4(1, pow(10, -8));
		read_path_probs_4.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_4.readCount() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_4.noiseProbability(), 0.1));
		REQUIRE(read_path_probs_4.probabilities().size() == 2);

		REQUIRE(read_path_probs.probabilities().front().first == 0);
		REQUIRE(Utils::doubleCompare(read_path_probs_4.probabilities().front().second, 0.359999999999999));
		REQUIRE(read_path_probs.probabilities().back().first == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_4.probabilities().back().second, 0.540000000000000));
	}
}

TEST_CASE("Identical read path probabilities can be merged") {

	spp::sparse_hash_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}});
	FragmentLengthDist fragment_length_dist(10, 2);

	vector<AlignmentPath> alignment_paths(1, AlignmentPath(make_pair(gbwt::SearchState(), 0), false, 10, 10, 3));
	auto alignment_path_ids = vector<vector<gbwt::size_type> >(1, vector<gbwt::size_type>({100, 200}));

	vector<PathInfo> paths(2, PathInfo(""));
	paths.front().effective_length = 3;
	paths.back().effective_length = 3;

	ReadPathProbabilities read_path_probs(1, pow(10, -8));
	read_path_probs.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

	REQUIRE(read_path_probs.mergeIdenticalReadPathProbabilities(read_path_probs));

	REQUIRE(read_path_probs.readCount() == 2);
	REQUIRE(Utils::doubleCompare(read_path_probs.noiseProbability(), 0.1));
	REQUIRE(read_path_probs.probabilities().size() == 2);

	REQUIRE(read_path_probs.probabilities().front().first == 0);
	REQUIRE(Utils::doubleCompare(read_path_probs.probabilities().front().second, 0.45));
	REQUIRE(read_path_probs.probabilities().back().first == 1);
	REQUIRE(Utils::doubleCompare(read_path_probs.probabilities().back().second, 0.45));

	SECTION("Probability precision affect merge") {

		vector<PathInfo> paths(2, PathInfo(""));
		paths.front().effective_length = 2;
		paths.back().effective_length = 3;

		ReadPathProbabilities read_path_probs_2(3, 0.1);
		read_path_probs_2.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_2.mergeIdenticalReadPathProbabilities(read_path_probs));
		REQUIRE(read_path_probs_2.readCount() == 5);
	}
}

TEST_CASE("Read path probabilities can be collapsed") {

	
	spp::sparse_hash_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}, {10, 2}, {50, 3}});
	FragmentLengthDist fragment_length_dist(10, 2);

	vector<AlignmentPath> alignment_paths(1, AlignmentPath(make_pair(gbwt::SearchState(), 0), false, 10, 10, 3));
	alignment_paths.emplace_back(AlignmentPath(make_pair(gbwt::SearchState(), 0), false, 15, 10, 5));

	auto alignment_path_ids = vector<vector<gbwt::size_type> >(1, vector<gbwt::size_type>({100, 200}));
	alignment_path_ids.emplace_back(vector<gbwt::size_type>({50}));

	vector<PathInfo> paths(4, PathInfo(""));
	paths.at(0).effective_length = 3;
	paths.at(1).effective_length = 3;
	paths.at(2).effective_length = 3;
	paths.at(3).effective_length = 3;

	ReadPathProbabilities read_path_probs(1, 0.01);
	read_path_probs.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

	auto collapsed_probs = read_path_probs.collapsedProbabilities();
	REQUIRE(collapsed_probs.size() == 2);

	REQUIRE(Utils::doubleCompare(collapsed_probs.at(0).first, 0.233044027062125));
	REQUIRE(Utils::doubleCompare(collapsed_probs.at(1).first, 0.333477986468937));

	REQUIRE(collapsed_probs.at(0).second == vector<uint32_t>({3}));
	REQUIRE(collapsed_probs.at(1).second == vector<uint32_t>({0, 1}));

	SECTION("Probability precision affect collapse") {

		ReadPathProbabilities read_path_probs(1, 0.2);
		read_path_probs.calcReadPathProbabilities(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		auto collapsed_probs = read_path_probs.collapsedProbabilities();
		REQUIRE(collapsed_probs.size() == 1);

		REQUIRE(Utils::doubleCompare(collapsed_probs.at(0).first, 0.333477986468937));
		REQUIRE(collapsed_probs.at(0).second == vector<uint32_t>({0, 1, 3}));
	}
}

