
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
	read_path_probs.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

	REQUIRE(read_path_probs.readCount() == 1);
	REQUIRE(Utils::doubleCompare(read_path_probs.noiseProb(), 0.1));

	REQUIRE(read_path_probs.pathProbs().size() == 1);
	REQUIRE(Utils::doubleCompare(read_path_probs.pathProbs().front().first, 0.45));
	REQUIRE(read_path_probs.pathProbs().front().second == vector<uint32_t>({0, 1}));

    SECTION("Improbable alignment path returns finite path probabilities") {

    	alignment_paths.front().frag_length = 100000;

		ReadPathProbabilities read_path_probs_2(1, pow(10, -8));
		read_path_probs_2.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_2.readCount() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), 0.1));

		REQUIRE(read_path_probs_2.pathProbs().size() == 1);
		REQUIRE(abs(read_path_probs_2.pathProbs().front().first - read_path_probs.pathProbs().front().first) < pow(10, -8));
		REQUIRE(read_path_probs_2.pathProbs().front().second == read_path_probs.pathProbs().front().second);
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

		ReadPathProbabilities read_path_probs_2(1, pow(10, -8));
		read_path_probs_2.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_2.readCount() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), 0.1));
		
		REQUIRE(read_path_probs_2.pathProbs().size() == 2);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().front().first, 0.233044027062125));
		REQUIRE(read_path_probs_2.pathProbs().front().second == vector<uint32_t>({3}));
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().back().first, 0.333477986468937));
		REQUIRE(read_path_probs_2.pathProbs().back().second == vector<uint32_t>({0, 1}));

		SECTION("Probability precision affect number of unique path probabilities") {

			paths.back().effective_length = 2;

			ReadPathProbabilities read_path_probs_3(1, 0.1);
			read_path_probs_3.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

			REQUIRE(read_path_probs_3.readCount() == 1);
			REQUIRE(Utils::doubleCompare(read_path_probs_3.noiseProb(), 0.1));

			REQUIRE(read_path_probs_3.pathProbs().size() == 1);
			REQUIRE(Utils::doubleCompare(read_path_probs_3.pathProbs().front().first, 0.3));
			REQUIRE(read_path_probs_3.pathProbs().front().second == vector<uint32_t>({0, 1, 3}));
		}
	}

    SECTION("Effective path lengths affect path probabilities") {

		paths.back().effective_length = 2;

		ReadPathProbabilities read_path_probs_2(1, pow(10, -8));
		read_path_probs_2.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_2.readCount() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), 0.1));

		REQUIRE(read_path_probs_2.pathProbs().size() == 2);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().front().first, 0.36));
		REQUIRE(read_path_probs_2.pathProbs().front().second == vector<uint32_t>({0}));
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().back().first, 0.54));
		REQUIRE(read_path_probs_2.pathProbs().back().second == vector<uint32_t>({1}));
	}

    SECTION("Base noise probability affect path probabilities") {

		ReadPathProbabilities read_path_probs_2(1, pow(10, -8));
		read_path_probs_2.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0.3);

		REQUIRE(read_path_probs_2.readCount() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), 0.3));

		REQUIRE(read_path_probs_2.pathProbs().size() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().front().first, 0.35));
		REQUIRE(read_path_probs_2.pathProbs().front().second == read_path_probs.pathProbs().front().second);
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
	read_path_probs.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

	REQUIRE(read_path_probs.quickMergeIdentical(read_path_probs));

	REQUIRE(read_path_probs.readCount() == 2);
	REQUIRE(Utils::doubleCompare(read_path_probs.noiseProb(), 0.1));

	REQUIRE(read_path_probs.pathProbs().size() == 1);
	REQUIRE(Utils::doubleCompare(read_path_probs.pathProbs().front().first, 0.45));
	REQUIRE(read_path_probs.pathProbs().front().second == vector<uint32_t>({0, 1}));
}

