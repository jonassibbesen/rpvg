
#include "catch2/catch.hpp"

#include "../alignment_path.hpp"
#include "../utils.hpp"


TEST_CASE("Different AlignmentPaths can be equal") {
    
	AlignmentPath alignment_path_1;

	alignment_path_1.mapqs.push_back(10);
	alignment_path_1.mapqs.push_back(20);

	alignment_path_1.scores.push_back(1);
	alignment_path_1.scores.push_back(2);

	AlignmentPath alignment_path_2 = alignment_path_1;
    REQUIRE(alignment_path_1 == alignment_path_2);

    swap(alignment_path_1.mapqs.front(), alignment_path_1.mapqs.back());
	REQUIRE(alignment_path_1 == alignment_path_2);

	swap(alignment_path_1.scores.front(), alignment_path_1.scores.back());
    REQUIRE(alignment_path_1 == alignment_path_2);

	alignment_path_1.scores.push_back(3);
	alignment_path_1.scores.push_back(4);	    

	alignment_path_2.scores.push_back(4);
	alignment_path_2.scores.push_back(3);	   

    REQUIRE(alignment_path_1 == alignment_path_2);
}

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

    SECTION("Mapping quality order does not matter") {

		auto alignment_path_mapq_rev = alignment_path;

    	swap(alignment_path_mapq_rev.mapqs.front(), alignment_path_mapq_rev.mapqs.back());
		REQUIRE(doubleCompare(alignment_path_mapq_rev.mapqProb(), alignment_path.mapqProb()));
	}
}

