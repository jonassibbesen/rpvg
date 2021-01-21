
#include "catch.hpp"

#include "gbwt/dynamic_gbwt.h"
#include "gbwt/fast_locate.h"
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

    gbwt_thread_4[0] = gbwt::Node::encode(6, true);
    gbwt_thread_4[1] = gbwt::Node::encode(7, true);

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
	Utils::json2pb(graph, graph_str);

    gbwt::FastLocate r_index(gbwt_index);
    PathsIndex paths_index(gbwt_index, r_index, graph);

    REQUIRE(!paths_index.bidirectional());
    REQUIRE(paths_index.numberOfPaths() == 4);

    spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> align_paths_index;

    PathClusters path_clusters(1, paths_index, align_paths_index);
    path_clusters.addNodeClusters(paths_index);

    REQUIRE(path_clusters.path_to_cluster_index.size() == 4);
    REQUIRE(path_clusters.path_to_cluster_index == vector<uint32_t>({0, 1, 2, 1}));
    REQUIRE(path_clusters.cluster_to_paths_index.size() == 3);
    REQUIRE(path_clusters.cluster_to_paths_index.at(0) == vector<uint32_t>({0}));
    REQUIRE(path_clusters.cluster_to_paths_index.at(1) == vector<uint32_t>({1, 3}));
    REQUIRE(path_clusters.cluster_to_paths_index.at(2) == vector<uint32_t>({2}));

    SECTION("Bidirectionality affect clustering") {

    	gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
	    gbwt::GBWTBuilder gbwt_builder_bd(gbwt::bit_length(gbwt::Node::encode(7, true)));

	    gbwt_builder_bd.insert(gbwt_thread_1, true);
	    gbwt_builder_bd.insert(gbwt_thread_2, true);
	    gbwt_builder_bd.insert(gbwt_thread_3, true);
	    gbwt_builder_bd.insert(gbwt_thread_4, true);

	    gbwt_builder_bd.index.addMetadata();

	    for (uint32_t i = 0; i < 4; ++i) {

	    	gbwt_builder_bd.index.metadata.addPath(gbwt::PathName());
	    }    

	    gbwt_builder_bd.finish();

	    std::stringstream gbwt_stream_bd;
	    gbwt_builder_bd.index.serialize(gbwt_stream_bd);

	    gbwt::GBWT gbwt_index_bd;
	    gbwt_index_bd.load(gbwt_stream_bd);

        gbwt::FastLocate r_index_bd(gbwt_index_bd);
	    PathsIndex paths_index_bd(gbwt_index_bd, r_index_bd, graph);

        REQUIRE(paths_index_bd.bidirectional());
        REQUIRE(paths_index.numberOfPaths() == 4);

    	path_clusters.addNodeClusters(paths_index_bd);

	    REQUIRE(path_clusters.path_to_cluster_index.size() == 4);
	    REQUIRE(path_clusters.path_to_cluster_index == vector<uint32_t>({0, 0, 1, 0}));
	    REQUIRE(path_clusters.cluster_to_paths_index.size() == 2);
	    REQUIRE(path_clusters.cluster_to_paths_index.at(0) == vector<uint32_t>({0, 1, 3}));
	    REQUIRE(path_clusters.cluster_to_paths_index.at(1) == vector<uint32_t>({2}));
    }
}

