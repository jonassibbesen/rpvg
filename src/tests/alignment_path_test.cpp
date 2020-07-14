
#include "catch.hpp"

#include "gbwt/gbwt.h"

#include "../alignment_path.hpp"
#include "../utils.hpp"


TEST_CASE("Multiple mapping qualities can be combined into single probability") {
    
	AlignmentSearchPath alignment_search_path;

	alignment_search_path.mapqs.push_back(10);
	REQUIRE(doubleCompare(alignment_search_path.mapqProb(), 0.1));

	alignment_search_path.mapqs.push_back(20);
	REQUIRE(doubleCompare(alignment_search_path.mapqProb(), 0.109));

    SECTION("Mapping quality of zero returns one") {

		alignment_search_path.mapqs.push_back(0);
		REQUIRE(doubleCompare(alignment_search_path.mapqProb(), 1));
	}

    SECTION("Mapping quality order not relevant") {

		auto alignment_path_mapq_rev = alignment_search_path;

    	swap(alignment_path_mapq_rev.mapqs.front(), alignment_path_mapq_rev.mapqs.back());
		REQUIRE(doubleCompare(alignment_path_mapq_rev.mapqProb(), alignment_search_path.mapqProb()));
	}
}

TEST_CASE("Empty AlignmentSearchPath is not complete") {
    
	AlignmentSearchPath alignment_search_path;
	REQUIRE(alignment_search_path.complete() == false);
}

TEST_CASE("AlignmentPath can be created from AlignmentSearchPath") {
    
	AlignmentSearchPath alignment_search_path;

	alignment_search_path.seq_length = 100;

	alignment_search_path.mapqs.push_back(10);
	alignment_search_path.mapqs.push_back(20);

	alignment_search_path.scores.push_back(50);
	alignment_search_path.scores.push_back(60);

	AlignmentPath alignment_path(alignment_search_path);
	
	REQUIRE(alignment_path.seq_length == 100);
	REQUIRE(alignment_path.mapq_comb == 10);
	REQUIRE(alignment_path.score_sum == 110);
	REQUIRE(alignment_path.search_state.empty());
}


