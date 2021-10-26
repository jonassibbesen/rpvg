
#include "catch.hpp"
#include "sparsepp/spp.h"

#include "../read_path_probabilities.hpp"
#include "../utils.hpp"

   
TEST_CASE("Read path probabilities can be calculated from alignment paths") {
    
	spp::sparse_hash_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}});
	FragmentLengthDist fragment_length_dist(10, 2);

	vector<AlignmentPath> alignment_paths;
	alignment_paths.emplace_back(make_pair(gbwt::SearchState(), 0), true, 10, 3, 5, 10);
	alignment_paths.emplace_back(make_pair(gbwt::SearchState(), 0), true, 10, numeric_limits<int32_t>::lowest(), 0, 0);
	
	vector<vector<gbwt::size_type> > alignment_path_ids;
	alignment_path_ids.emplace_back(vector<gbwt::size_type>({100, 200}));
	alignment_path_ids.emplace_back(vector<gbwt::size_type>());

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

		REQUIRE(read_path_probs_2.readCount() == read_path_probs.readCount());
		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), read_path_probs.noiseProb()));

		REQUIRE(read_path_probs_2.pathProbs().size() == 1);
		REQUIRE(abs(read_path_probs_2.pathProbs().front().first - read_path_probs.pathProbs().front().first) < pow(10, -8));
		REQUIRE(read_path_probs_2.pathProbs().front().second == read_path_probs.pathProbs().front().second);
	}

    SECTION("Probabilities are calculated from multiple alignment paths") {

		alignment_paths.at(1) = AlignmentPath(make_pair(gbwt::SearchState(), 0), true, 10, 5, 8, 15);
		alignment_paths.emplace_back(make_pair(gbwt::SearchState(), 0), true, 10, numeric_limits<int32_t>::lowest(), 0, 0);

		alignment_path_ids.at(1) = vector<gbwt::size_type>({50});
		alignment_path_ids.emplace_back(vector<gbwt::size_type>());
		
		clustered_path_index.emplace(10, 2);
		clustered_path_index.emplace(50, 3);

		paths.emplace_back(PathInfo(""));
		paths.back().effective_length = 3;

		paths.emplace_back(PathInfo(""));
		paths.back().effective_length = 3;

		ReadPathProbabilities read_path_probs_2(1, pow(10, -8));
		read_path_probs_2.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_2.readCount() == read_path_probs.readCount());
		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), read_path_probs.noiseProb()));
		
		REQUIRE(read_path_probs_2.pathProbs().size() == 2);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().front().first, 0.233044027062125));
		REQUIRE(read_path_probs_2.pathProbs().front().second == vector<uint32_t>({3}));
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().back().first, 0.333477986468937));
		REQUIRE(read_path_probs_2.pathProbs().back().second == vector<uint32_t>({0, 1}));

		SECTION("Probability precision affect number of unique path probabilities") {

			paths.back().effective_length = 2;

			ReadPathProbabilities read_path_probs_3(1, 0.1);
			read_path_probs_3.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

			REQUIRE(read_path_probs_2.readCount() == read_path_probs.readCount());
			REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), read_path_probs.noiseProb()));

			REQUIRE(read_path_probs_3.pathProbs().size() == 1);
			REQUIRE(Utils::doubleCompare(read_path_probs_3.pathProbs().front().first, 0.3));
			REQUIRE(read_path_probs_3.pathProbs().front().second == vector<uint32_t>({0, 1, 3}));
		}

		SECTION("Longest alignment path is always chosen") {

			alignment_paths.at(2) = AlignmentPath(make_pair(gbwt::SearchState(), 0), true, 10, 3, 10, 10);
			alignment_paths.emplace_back(make_pair(gbwt::SearchState(), 0), true, 10, numeric_limits<int32_t>::lowest(), 0, 0);

			alignment_path_ids.at(2) = vector<gbwt::size_type>({50});
			alignment_path_ids.emplace_back(vector<gbwt::size_type>());

			ReadPathProbabilities read_path_probs_3(1, 0.1);
			read_path_probs_3.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

			REQUIRE(read_path_probs_3.readCount() == read_path_probs.readCount());
			REQUIRE(Utils::doubleCompare(read_path_probs_3.noiseProb(), read_path_probs.noiseProb()));

			REQUIRE(read_path_probs_3.pathProbs().size() == 1);
			REQUIRE(Utils::doubleCompare(read_path_probs_3.pathProbs().front().first, 0.3));
			REQUIRE(read_path_probs_3.pathProbs().front().second == vector<uint32_t>({0, 1, 3}));
		}

		SECTION("Highest scoring alignment path is chosen if identical") {

			alignment_paths.at(2) = AlignmentPath(make_pair(gbwt::SearchState(), 0), true, 10, 3, 8, 15);
			alignment_paths.emplace_back(make_pair(gbwt::SearchState(), 0), true, 10, numeric_limits<int32_t>::lowest(), 0, 0);

			alignment_path_ids.at(2) = vector<gbwt::size_type>({50});
			alignment_path_ids.emplace_back(vector<gbwt::size_type>());

			ReadPathProbabilities read_path_probs_3(1, 0.1);
			read_path_probs_3.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

			REQUIRE(read_path_probs_3.readCount() == read_path_probs.readCount());
			REQUIRE(Utils::doubleCompare(read_path_probs_3.noiseProb(), read_path_probs.noiseProb()));

			REQUIRE(read_path_probs_3.pathProbs().size() == 2);
			REQUIRE(abs(read_path_probs_3.pathProbs().front().first - read_path_probs_2.pathProbs().front().first) < pow(10, -8));
			REQUIRE(read_path_probs_3.pathProbs().front().second == read_path_probs_2.pathProbs().front().second);
			REQUIRE(abs(read_path_probs_3.pathProbs().back().first - read_path_probs_2.pathProbs().back().first) < pow(10, -8));
			REQUIRE(read_path_probs_3.pathProbs().back().second == read_path_probs_2.pathProbs().back().second);
		}
	}

    SECTION("Noise alignment path affect noise probability affect") {

   		alignment_paths.back().score_sum = -2.302585 / Utils::noise_score_log_base;

		ReadPathProbabilities read_path_probs_2(1, pow(10, -8));
		read_path_probs_2.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_2.readCount() == read_path_probs.readCount());
		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), 0.190000008369464));

		REQUIRE(read_path_probs_2.pathProbs().size() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().front().first, 0.404999995815267));
		REQUIRE(read_path_probs_2.pathProbs().front().second == vector<uint32_t>({0, 1}));

		alignment_paths.back().score_sum = 0;

		ReadPathProbabilities read_path_probs_3(1, pow(10, -8));
		read_path_probs_3.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_3.readCount() == read_path_probs.readCount());
		REQUIRE(Utils::doubleCompare(read_path_probs_3.noiseProb(), 1));

		REQUIRE(read_path_probs_3.pathProbs().empty());
	}

    SECTION("Effective path lengths affect path probabilities") {

		paths.back().effective_length = 2;

		ReadPathProbabilities read_path_probs_2(1, pow(10, -8));
		read_path_probs_2.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0);

		REQUIRE(read_path_probs_2.readCount() == read_path_probs.readCount());
		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), read_path_probs.noiseProb()));

		REQUIRE(read_path_probs_2.pathProbs().size() == 2);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().front().first, 0.36));
		REQUIRE(read_path_probs_2.pathProbs().front().second == vector<uint32_t>({0}));
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().back().first, 0.54));
		REQUIRE(read_path_probs_2.pathProbs().back().second == vector<uint32_t>({1}));
	}

    SECTION("Base noise probability affect path probabilities") {

   		alignment_paths.back().score_sum = -5.0 / Utils::noise_score_log_base;

		ReadPathProbabilities read_path_probs_2(1, pow(10, -8));
		read_path_probs_2.calcAlignPathProbs(alignment_paths, alignment_path_ids, clustered_path_index, paths, fragment_length_dist, false, 0.3);
		
		REQUIRE(read_path_probs_2.readCount() == read_path_probs.readCount());
		REQUIRE(Utils::doubleCompare(read_path_probs_2.noiseProb(), 0.304716562899359));

		REQUIRE(read_path_probs_2.pathProbs().size() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs_2.pathProbs().front().first, 0.347641718550320));
		REQUIRE(read_path_probs_2.pathProbs().front().second == read_path_probs.pathProbs().front().second);
	}


	SECTION("Identical read path probabilities can be merged") {

		REQUIRE(read_path_probs.quickMergeIdentical(read_path_probs));

		REQUIRE(read_path_probs.readCount() == 2);
		REQUIRE(Utils::doubleCompare(read_path_probs.noiseProb(), 0.1));

		REQUIRE(read_path_probs.pathProbs().size() == 1);
		REQUIRE(Utils::doubleCompare(read_path_probs.pathProbs().front().first, 0.45));
		REQUIRE(read_path_probs.pathProbs().front().second == vector<uint32_t>({0, 1}));
	}
}

