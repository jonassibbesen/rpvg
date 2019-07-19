
#include "catch2/catch.hpp"

#include "../read_path_probs.hpp"
#include "../utils.hpp"


TEST_CASE("ReadPathProbs can be equal") {
    
	ReadPathProbs read_path_probs_1(2);
	ReadPathProbs read_path_probs_2(2);

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

TEST_CASE("Multiple ReadPathProbs can be sorted") {
    
	vector<ReadPathProbs> multiple_read_path_probs(2, ReadPathProbs(2));

	multiple_read_path_probs.front().read_path_probs.back() = 1;
	sort(multiple_read_path_probs.begin(), multiple_read_path_probs.end());

	REQUIRE(doubleCompare(multiple_read_path_probs.front().read_path_probs.back(), 0));
	REQUIRE(doubleCompare(multiple_read_path_probs.back().read_path_probs.back(), 1));

	multiple_read_path_probs.back().noise_prob = 0;
	sort(multiple_read_path_probs.begin(), multiple_read_path_probs.end());

	REQUIRE(doubleCompare(multiple_read_path_probs.front().read_path_probs.back(), 1));
	REQUIRE(doubleCompare(multiple_read_path_probs.back().read_path_probs.back(), 0));
}

TEST_CASE("Relative probabilities can be calculated from alignment scores") {
    
	vector<AlignmentPath> alignment_paths(1, AlignmentPath());
	alignment_paths.front().scores.push_back(1);

	ReadPathProbs read_path_probs(0);
	auto align_paths_log_probs = read_path_probs.calcRelativeAlignmentScoreLogProbs(alignment_paths);

	REQUIRE(align_paths_log_probs.size() == 1);
	REQUIRE(doubleCompare(align_paths_log_probs.front(), 0));	

	alignment_paths.emplace_back(AlignmentPath());
	alignment_paths.back().scores.push_back(2);
	align_paths_log_probs = read_path_probs.calcRelativeAlignmentScoreLogProbs(alignment_paths);

	REQUIRE(align_paths_log_probs.size() == 2);
	REQUIRE(doubleCompare(align_paths_log_probs.front(), -1.607063344188675));			
	REQUIRE(doubleCompare(align_paths_log_probs.back(), -0.223738075450029));			
}

TEST_CASE("Read path probabilities can be calculated from alignment paths") {
    
	vector<AlignmentPath> alignment_paths(1, AlignmentPath());
	alignment_paths.front().path_ids.push_back(100);
	alignment_paths.front().path_ids.push_back(200);
	alignment_paths.front().seq_length = 10;
	alignment_paths.front().mapqs.push_back(10);
	alignment_paths.front().mapqs.push_back(20);
	alignment_paths.front().scores.push_back(1);
	alignment_paths.front().scores.push_back(2);

	unordered_map<uint32_t, uint32_t> clustered_path_index({{100, 0}, {200, 1}});
	FragmentLengthDist fragment_length_dist(10, 2);

	ReadPathProbs read_path_probs(2);
	read_path_probs.calcReadPathProbs(alignment_paths, clustered_path_index, fragment_length_dist);

	REQUIRE(doubleCompare(read_path_probs.noise_prob, 0.109));
	REQUIRE(read_path_probs.read_path_probs.size() == 2);
	REQUIRE(doubleCompare(read_path_probs.read_path_probs.front(), 0.5));
	REQUIRE(doubleCompare(read_path_probs.read_path_probs.back(), 0.5));

    SECTION("Extremely improbable alignment path returns finite probabilities") {

    	alignment_paths.front().seq_length = 10000; 

		ReadPathProbs read_path_probs_2(2);
		read_path_probs_2.calcReadPathProbs(alignment_paths, clustered_path_index, fragment_length_dist);

		REQUIRE(read_path_probs == read_path_probs_2);
	}

    SECTION("Probabilities can be calculated from multiple alignment paths") {

		alignment_paths.emplace_back(AlignmentPath());
		alignment_paths.back().path_ids.push_back(50);
		alignment_paths.back().seq_length = 15;
		alignment_paths.back().mapqs.push_back(10);
		alignment_paths.back().mapqs.push_back(20);
		alignment_paths.back().scores.push_back(3);
		alignment_paths.back().scores.push_back(2);

		clustered_path_index.emplace(10, 2);
		clustered_path_index.emplace(50, 3);

		ReadPathProbs read_path_probs_3(4);
		read_path_probs_3.calcReadPathProbs(alignment_paths, clustered_path_index, fragment_length_dist);

		REQUIRE(doubleCompare(read_path_probs_3.noise_prob, 0.109));
		REQUIRE(read_path_probs_3.read_path_probs.size() == 4);
		REQUIRE(doubleCompare(read_path_probs_3.read_path_probs.front(), 0.2943204134206829));
		REQUIRE(doubleCompare(read_path_probs_3.read_path_probs.at(1), 0.2943204134206829));
		REQUIRE(doubleCompare(read_path_probs_3.read_path_probs.at(2), 0));
		REQUIRE(doubleCompare(read_path_probs_3.read_path_probs.back(), 0.4113591731586340));
	}
}



