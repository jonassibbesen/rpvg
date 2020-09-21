
#include "catch.hpp"

#include "gbwt/dynamic_gbwt.h"
#include "sparsepp/spp.h"

#include "../path_clusters.hpp"
#include "../utils.hpp"


TEST_CASE("GBWT paths can be clustered") {

	gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(7, true)));

    gbwt::vector_type gbwt_thread_1(3);
    gbwt::vector_type gbwt_thread_2(2);
    gbwt::vector_type gbwt_thread_3(1);
    gbwt::vector_type gbwt_thread_4(2);
   
    gbwt_thread_1[0] = gbwt::Node::encode(1, false);
    gbwt_thread_1[1] = gbwt::Node::encode(2, false);
    gbwt_thread_1[2] = gbwt::Node::encode(4, false);

    gbwt_thread_2[0] = gbwt::Node::encode(1, true);
    gbwt_thread_2[1] = gbwt::Node::encode(6, true);

    gbwt_thread_3[0] = gbwt::Node::encode(3, false);

    gbwt_thread_4[0] = gbwt::Node::encode(6, false);
    gbwt_thread_4[1] = gbwt::Node::encode(7, false);

    gbwt_builder.insert(gbwt_thread_1, false);
    gbwt_builder.insert(gbwt_thread_2, false);
    gbwt_builder.insert(gbwt_thread_3, false);
    gbwt_builder.insert(gbwt_thread_4, false);

    gbwt_builder.index.addMetadata();

    for (uint32_t i = 0; i < 4; ++i) {

    	gbwt_builder.index.metadata.addPath(gbwt::PathName());
    }    

    gbwt_builder.finish();

    std::stringstream gbwt_stream;
    gbwt_builder.index.serialize(gbwt_stream);

    gbwt::GBWT gbwt_index;
    gbwt_index.load(gbwt_stream);

    const string graph_str = R"(
    	{
    		"node": [
    			{"id": 1, "sequence": "A"},
    			{"id": 2, "sequence": "A"},
    			{"id": 3, "sequence": "A"},
    			{"id": 4, "sequence": "A"},
    			{"id": 5, "sequence": "A"},
    			{"id": 6, "sequence": "A"},
    			{"id": 7, "sequence": "A"}
    		],
    	}
    )";

	vg::Graph graph;
	json2pb(graph, graph_str);

    PathsIndex paths_index(gbwt_index, graph);
    REQUIRE(paths_index.index().metadata.paths() == 4);

	PathClusters path_clusters(paths_index, 1);

    REQUIRE(path_clusters.path_to_cluster_index.size() == 4);
    REQUIRE(path_clusters.path_to_cluster_index == vector<uint32_t>({0, 0, 1, 0}));
    REQUIRE(path_clusters.cluster_to_paths_index.size() == 2);
    REQUIRE(path_clusters.cluster_to_paths_index.at(0) == vector<uint32_t>({0, 1, 3}));
    REQUIRE(path_clusters.cluster_to_paths_index.at(1) == vector<uint32_t>({2}));

    REQUIRE(path_clusters.node_to_path_index.size() == 6);    
    REQUIRE(path_clusters.node_to_path_index.at(1) == 0);
    REQUIRE(path_clusters.node_to_path_index.at(2) == 0);
    REQUIRE(path_clusters.node_to_path_index.at(3) == 2);
    REQUIRE(path_clusters.node_to_path_index.at(4) == 0);
    REQUIRE(path_clusters.node_to_path_index.at(6) == 3);
    REQUIRE(path_clusters.node_to_path_index.at(7) == 3);

    SECTION("Bidirectionality does not affect clustering") {

    	gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
	    gbwt::GBWTBuilder gbwt_builder_bidi(gbwt::bit_length(gbwt::Node::encode(7, true)));

	    gbwt_builder_bidi.insert(gbwt_thread_1, true);
	    gbwt_builder_bidi.insert(gbwt_thread_2, true);
	    gbwt_builder_bidi.insert(gbwt_thread_3, true);
	    gbwt_builder_bidi.insert(gbwt_thread_4, true);

	    gbwt_builder_bidi.index.addMetadata();

	    for (uint32_t i = 0; i < 4; ++i) {

	    	gbwt_builder_bidi.index.metadata.addPath(gbwt::PathName());
	    }    

	    gbwt_builder_bidi.finish();

	    std::stringstream gbwt_stream_bidi;
	    gbwt_builder_bidi.index.serialize(gbwt_stream_bidi);

	    gbwt::GBWT gbwt_index_bidi;
	    gbwt_index_bidi.load(gbwt_stream_bidi);

	    PathsIndex paths_index_bidi(gbwt_index_bidi, graph);
	    REQUIRE(paths_index_bidi.index().metadata.paths() == 4);

		PathClusters path_clusters_bidi(paths_index, 1);

	    REQUIRE(path_clusters_bidi.path_to_cluster_index == path_clusters.path_to_cluster_index);
	    REQUIRE(path_clusters_bidi.path_to_cluster_index == path_clusters.path_to_cluster_index);
	    REQUIRE(path_clusters_bidi.node_to_path_index == path_clusters.node_to_path_index);
    }
}

