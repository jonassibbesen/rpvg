
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
