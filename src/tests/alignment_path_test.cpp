
#include "catch.hpp"

#include "gbwt/gbwt.h"

#include "../alignment_path.hpp"
#include "../utils.hpp"


TEST_CASE("Empty AlignmentSearchPath is not complete") {
    
	AlignmentSearchPath alignment_search_path;
	REQUIRE(alignment_search_path.complete() == false);
}

TEST_CASE("AlignmentPath can be created from AlignmentSearchPath") {
    
	AlignmentSearchPath alignment_search_path;

	alignment_search_path.seq_length = 100;
	alignment_search_path.min_mapq = 10;
	
	alignment_search_path.scores.emplace_back(50,60);
	alignment_search_path.scores.emplace_back(60,100);

	REQUIRE(doubleCompare(alignment_search_path.minRelativeScore(), 0.6));

	AlignmentPath alignment_path(alignment_search_path, false);
	
	REQUIRE(alignment_path.seq_length == 100);
	REQUIRE(alignment_path.min_mapq == 10);
	REQUIRE(alignment_path.score_sum == 110);
	REQUIRE(alignment_path.search_state.empty());
}


