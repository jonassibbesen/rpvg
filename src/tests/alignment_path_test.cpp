
#include "catch2/catch.hpp"
#include "../alignment_path.hpp"


TEST_CASE("Alignment paths are equal", "[alignment_path]" ) {
    
	AlignmentPath alignment_path_1;
	AlignmentPath alignment_path_2;

    REQUIRE(alignment_path_1 == alignment_path_2);
}