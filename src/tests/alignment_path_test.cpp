
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
	alignment_search_path.read_stats.back().left_softclip_length = make_pair(10, true);
	alignment_search_path.read_stats.back().right_softclip_length = make_pair(30, true);

    alignment_search_path.read_stats.back().updateInternalStartOffset(20, true);
    alignment_search_path.read_stats.back().updateInternalEndOffset(20, false);

	REQUIRE(alignment_search_path.read_stats.back().clippedOffsetLeftBases() == 20);
	REQUIRE(alignment_search_path.read_stats.back().clippedOffsetRightBases() == 50);

	REQUIRE(alignment_search_path.read_stats.back().adjustedScore() == 20);
	REQUIRE(alignment_search_path.read_stats.back().clippedOffsetTotalBases() == 70);

	alignment_search_path.read_stats.emplace_back(ReadAlignmentStats());
	alignment_search_path.read_stats.back().mapq = 20;
	alignment_search_path.read_stats.back().score = 7;
	alignment_search_path.read_stats.back().length = 10;
	alignment_search_path.read_stats.back().left_softclip_length = make_pair(2, true);
	alignment_search_path.read_stats.back().right_softclip_length = make_pair(0, true);

	REQUIRE(alignment_search_path.read_stats.back().clippedOffsetLeftBases() == 2);
	REQUIRE(alignment_search_path.read_stats.back().clippedOffsetRightBases() == 0);

	REQUIRE(alignment_search_path.read_stats.back().adjustedScore() == 7);
	REQUIRE(alignment_search_path.read_stats.back().clippedOffsetTotalBases() == 2);

	REQUIRE(doubleCompare(alignment_search_path.fragmentLength(), 158));
	REQUIRE(doubleCompare(alignment_search_path.minMappingQuality(), 10));
	REQUIRE(doubleCompare(alignment_search_path.scoreSum(), 27));

	REQUIRE(doubleCompare(alignment_search_path.minBestScoreFraction(), 0.2));
	REQUIRE(doubleCompare(alignment_search_path.maxSoftclipFraction(), 0.4));

	AlignmentPath alignment_path(alignment_search_path, false);
	
	REQUIRE(alignment_path.frag_length == 158);
	REQUIRE(alignment_path.min_mapq == 10);
	REQUIRE(alignment_path.score_sum == 27);
	REQUIRE(alignment_path.gbwt_search.first.empty());

    SECTION("Insert length can be negative for overlapping paired-end alignments") {

    	alignment_search_path.insert_length = -8;
		AlignmentPath alignment_path_neg(alignment_search_path, false);
	
		REQUIRE(alignment_path_neg.frag_length == 50);
		REQUIRE(alignment_path_neg.min_mapq == alignment_path.min_mapq);
		REQUIRE(alignment_path_neg.score_sum == alignment_path.score_sum);
		REQUIRE(alignment_path_neg.gbwt_search == alignment_path.gbwt_search);
    }

    SECTION("Cleared AlignmentPath is empty") {

    	alignment_search_path.clear();
		REQUIRE(alignment_search_path.isEmpty());    	
    }    
}