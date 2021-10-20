
#include "catch.hpp"

#include "gbwt/gbwt.h"

#include "../alignment_path.hpp"
#include "../utils.hpp"


TEST_CASE("AlignmentPath can be created from AlignmentSearchPath") {
    
	AlignmentSearchPath alignment_search_path;
	alignment_search_path.insert_length = 100;
	
	alignment_search_path.read_align_stats.emplace_back(AlignmentStats());
	alignment_search_path.read_align_stats.back().score = 50;
	alignment_search_path.read_align_stats.back().length = 100;
	alignment_search_path.read_align_stats.back().left_softclip_length = 10;
	alignment_search_path.read_align_stats.back().right_softclip_length = 30;

	alignment_search_path.read_align_stats.back().internal_start.is_internal = true;
	alignment_search_path.read_align_stats.back().internal_start.penalty = 10;
	alignment_search_path.read_align_stats.back().internal_start.offset = 10;

	alignment_search_path.read_align_stats.back().internal_end.is_internal = true;
	alignment_search_path.read_align_stats.back().internal_end.penalty = 15;
	alignment_search_path.read_align_stats.back().internal_end.offset = 20;

	REQUIRE(alignment_search_path.read_align_stats.back().clippedOffsetLeftBases() == 20);
	REQUIRE(alignment_search_path.read_align_stats.back().clippedOffsetRightBases() == 50);

	REQUIRE(alignment_search_path.read_align_stats.back().adjustedScore() == 25);
	REQUIRE(alignment_search_path.read_align_stats.back().clippedOffsetTotalBases() == 70);

	alignment_search_path.read_align_stats.emplace_back(AlignmentStats());
	alignment_search_path.read_align_stats.back().score = 7;
	alignment_search_path.read_align_stats.back().length = 10;
	alignment_search_path.read_align_stats.back().left_softclip_length = 2;
	alignment_search_path.read_align_stats.back().right_softclip_length = 0;

	REQUIRE(alignment_search_path.read_align_stats.back().clippedOffsetLeftBases() == 2);
	REQUIRE(alignment_search_path.read_align_stats.back().clippedOffsetRightBases() == 0);

	REQUIRE(alignment_search_path.read_align_stats.back().adjustedScore() == 7);
	REQUIRE(alignment_search_path.read_align_stats.back().clippedOffsetTotalBases() == 2);

	REQUIRE(Utils::doubleCompare(alignment_search_path.fragmentLength(), 158));
	REQUIRE(Utils::doubleCompare(alignment_search_path.scoreSum(), 32));

	REQUIRE(Utils::doubleCompare(alignment_search_path.minOptimalScoreFraction(vector<int32_t>({100, 10})), 0.25));
	REQUIRE(Utils::doubleCompare(alignment_search_path.maxSoftclipFraction(), 0.4));

	AlignmentPath alignment_path(alignment_search_path, false, 10);
	
	REQUIRE(alignment_path.frag_length == 158);
	REQUIRE(alignment_path.min_mapq == 10);
	REQUIRE(alignment_path.score_sum == 32);
	REQUIRE(alignment_path.gbwt_search.first.empty());

    SECTION("Insert length can be negative for overlapping paired-end alignments") {

    	alignment_search_path.insert_length = -8;
		AlignmentPath alignment_path_neg(alignment_search_path, false, 10);
	
		REQUIRE(alignment_path_neg.frag_length == 50);
		REQUIRE(alignment_path_neg.min_mapq == alignment_path.min_mapq);
		REQUIRE(alignment_path_neg.score_sum == alignment_path.score_sum);
		REQUIRE(alignment_path_neg.gbwt_search == alignment_path.gbwt_search);
    }

    SECTION("Cleared AlignmentPath is empty") {

    	alignment_search_path.clear();

		REQUIRE(alignment_search_path.path.empty());    	
		REQUIRE(alignment_search_path.gbwt_search.first.empty());    	
    }    
}