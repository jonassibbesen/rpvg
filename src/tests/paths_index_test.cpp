
#include "catch.hpp"

#include "gbwt/dynamic_gbwt.h"

#include "../paths_index.hpp"
#include "../utils.hpp"


TEST_CASE("Path index can calculate path lengths") {

    const string graph_str = R"(
    	{
    		"node": [
    			{"id": 1, "sequence": "GGGG"},
    			{"id": 2, "sequence": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"},
    			{"id": 3, "sequence": "C"},
    			{"id": 4, "sequence": "TT"}
    		],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 3},
                {"from": 2, "to": 4},
                {"from": 3, "to": 4}   
            ]
    	}
    )";

	vg::Graph graph;
	json2pb(graph, graph_str);

	gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(4, true)));

    gbwt::vector_type gbwt_thread_1(3);
    gbwt::vector_type gbwt_thread_2(3);
   
    gbwt_thread_1[0] = gbwt::Node::encode(1, false);
    gbwt_thread_1[1] = gbwt::Node::encode(2, false);
    gbwt_thread_1[2] = gbwt::Node::encode(4, false);

    gbwt_thread_2[0] = gbwt::Node::encode(1, false);
    gbwt_thread_2[1] = gbwt::Node::encode(3, false);
    gbwt_thread_2[2] = gbwt::Node::encode(4, false);

    gbwt_builder.insert(gbwt_thread_1, false);
    gbwt_builder.insert(gbwt_thread_2, false);

    gbwt_builder.finish();

    std::stringstream gbwt_stream;
    gbwt_builder.index.serialize(gbwt_stream);

    gbwt::GBWT gbwt_index;
    gbwt_index.load(gbwt_stream);

    PathsIndex paths_index(gbwt_index, graph);
    REQUIRE(paths_index.index().bidirectional() == false);

    REQUIRE(paths_index.pathLength(0) == 38);
    REQUIRE(paths_index.pathLength(1) == 7);

	SECTION("Effective paths length are calculated using fragment length distribution") {

		FragmentLengthDist fragment_length_dist(5, 2);

    	REQUIRE(doubleCompare(paths_index.effectivePathLength(0, fragment_length_dist), 32.889504274642021));
    	REQUIRE(doubleCompare(paths_index.effectivePathLength(1, fragment_length_dist), 2.4592743581826583));

    	fragment_length_dist = FragmentLengthDist(20, 1);

    	REQUIRE(doubleCompare(paths_index.effectivePathLength(0, fragment_length_dist), 18));
    	REQUIRE(doubleCompare(paths_index.effectivePathLength(1, fragment_length_dist), 1));
	}
}

