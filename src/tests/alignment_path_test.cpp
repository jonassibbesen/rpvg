
#include "catch.hpp"

#include "../alignment_path.hpp"
#include "../utils.hpp"


TEST_CASE("Multiple mapping qualities can be combined into single probability") {
    
	AlignmentPath alignment_path;

	alignment_path.mapqs.push_back(10);
	REQUIRE(doubleCompare(alignment_path.mapqProb(), 0.1));

	alignment_path.mapqs.push_back(20);
	REQUIRE(doubleCompare(alignment_path.mapqProb(), 0.109));

    SECTION("Mapping quality of zero returns one") {

		alignment_path.mapqs.push_back(0);
		REQUIRE(doubleCompare(alignment_path.mapqProb(), 1));
	}

    SECTION("Mapping quality order not relevant") {

		auto alignment_path_mapq_rev = alignment_path;

    	swap(alignment_path_mapq_rev.mapqs.front(), alignment_path_mapq_rev.mapqs.back());
		REQUIRE(doubleCompare(alignment_path_mapq_rev.mapqProb(), alignment_path.mapqProb()));
	}
}

TEST_CASE("Empty alignmentPath is complete") {
    
	AlignmentPath alignment_path;
	REQUIRE(alignment_path.complete() == true);
}
