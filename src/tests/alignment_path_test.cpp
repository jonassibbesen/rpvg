
#include "catch.hpp"

#include "gbwt/gbwt.h"

#include "../alignment_path.hpp"
#include "../utils.hpp"


TEST_CASE("New AlignmentSearchPath is empty") {
    
	AlignmentSearchPath alignment_search_path;
	REQUIRE(alignment_search_path.isEmpty());
}

TEST_CASE("AlignmentPath can be created from AlignmentSearchPath") {
    
	AlignmentSearchPath alignment_search_path;
	alignment_search_path.insert_length = 100;
	
	alignment_search_path.read_stats.emplace_back(ReadAlignmentStats());
	alignment_search_path.read_stats.back().mapq = 10;
	alignment_search_path.read_stats.back().score = 50;
	alignment_search_path.read_stats.back().length = 100;
	alignment_search_path.read_stats.back().left_softclip_length = 10;
	alignment_search_path.read_stats.back().right_softclip_length = 30;

	alignment_search_path.read_stats.emplace_back(ReadAlignmentStats());
	alignment_search_path.read_stats.back().mapq = 20;
	alignment_search_path.read_stats.back().score = 4;
	alignment_search_path.read_stats.back().length = 10;
	alignment_search_path.read_stats.back().left_softclip_length = 2;
	alignment_search_path.read_stats.back().right_softclip_length = 3;

	REQUIRE(doubleCompare(alignment_search_path.minBestScoreFraction(), 0.4));
	REQUIRE(doubleCompare(alignment_search_path.maxSoftclipFraction(), 0.5));

	AlignmentPath alignment_path(alignment_search_path, false);
	
	REQUIRE(alignment_path.frag_length == 178);
	REQUIRE(alignment_path.min_mapq == 10);
	REQUIRE(alignment_path.score_sum == 54);
	REQUIRE(alignment_path.search_state.empty());

    SECTION("Insert length can be negative for overlapping paired-end alignments") {

    	alignment_search_path.insert_length = -8;
		AlignmentPath alignment_path_neg(alignment_search_path, false);
	
		REQUIRE(alignment_path_neg.frag_length == 70);
		REQUIRE(alignment_path_neg.min_mapq == alignment_path.min_mapq);
		REQUIRE(alignment_path_neg.score_sum == alignment_path.score_sum);
		REQUIRE(alignment_path_neg.search_state == alignment_path.search_state);
    }

    SECTION("Cleared AlignmentPath is empty") {

    	alignment_search_path.clear();
		REQUIRE(alignment_search_path.isEmpty());    	
    }    
}