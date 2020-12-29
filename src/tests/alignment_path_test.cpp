
#include "catch.hpp"

#include "gbwt/gbwt.h"

#include "../alignment_path.hpp"
#include "../utils.hpp"


TEST_CASE("Empty AlignmentSearchPath is not complete") {
    
	AlignmentSearchPath alignment_search_path;
	REQUIRE(alignment_search_path.isComplete() == false);
}

TEST_CASE("AlignmentPath can be created from AlignmentSearchPath") {
    
	AlignmentSearchPath alignment_search_path;

	alignment_search_path.seq_length = 100;
	alignment_search_path.min_mapq = 10;
	
	alignment_search_path.read_stats.emplace_back(ReadAlignmentStats());
	alignment_search_path.read_stats.back().score = 50;
	alignment_search_path.read_stats.back().length = 100;
	alignment_search_path.read_stats.back().left_softclip_length = 10;
	alignment_search_path.read_stats.back().right_softclip_length = 30;

	alignment_search_path.read_stats.emplace_back(ReadAlignmentStats());
	alignment_search_path.read_stats.back().score = 60;
	alignment_search_path.read_stats.back().length = 10;
	alignment_search_path.read_stats.back().left_softclip_length = 2;
	alignment_search_path.read_stats.back().right_softclip_length = 3;

	REQUIRE(doubleCompare(alignment_search_path.maxSoftclipFraction(), 0.5));

	AlignmentPath alignment_path(alignment_search_path, false);
	
	REQUIRE(alignment_path.seq_length == 100);
	REQUIRE(alignment_path.min_mapq == 10);
	REQUIRE(alignment_path.score_sum == 110);
	REQUIRE(alignment_path.search_state.empty());
}

