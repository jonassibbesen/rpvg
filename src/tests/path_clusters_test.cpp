
#include "catch.hpp"

#include "gbwt/dynamic_gbwt.h"
#include "sparsepp/spp.h"

#include "../path_clusters.hpp"
#include "../utils.hpp"


TEST_CASE("Connected paths can be clustered") {

    REQUIRE(false);

	// gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
 //    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(1, true)));

 //    gbwt_builder.index.addMetadata();

 //    for (uint32_t i = 0; i < 7; ++i) {

 //    	gbwt_builder.index.metadata.addPath(gbwt::PathName());
 //    }    

 //    gbwt_builder.finish();

 //    std::stringstream gbwt_stream;
 //    gbwt_builder.index.serialize(gbwt_stream);

 //    gbwt::GBWT gbwt_index;
 //    gbwt_index.load(gbwt_stream);

	// vg::Graph graph;

 //    PathsIndex paths_index(gbwt_index, graph);
 //    REQUIRE(paths_index.index().metadata.paths() == 7);

	// spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > connected_paths;
	
	// connected_paths[1].emplace(2);
	// connected_paths[1].emplace(5);
	// connected_paths[2].emplace(1);
	// connected_paths[5].emplace(1);

	// connected_paths[6].emplace(3);
	// connected_paths[3].emplace(6);

	// PathClusters path_clusters(1);
	// path_clusters.findPathClusters(&connected_paths, paths_index);

 //    REQUIRE(path_clusters.path_to_cluster_index.size() == 7);
 //    REQUIRE(path_clusters.path_to_cluster_index == vector<uint32_t>({0, 1, 1, 2, 3, 1, 2}));
 //    REQUIRE(path_clusters.cluster_to_paths_index.size() == 4);
 //    REQUIRE(path_clusters.cluster_to_paths_index.at(0) == vector<uint32_t>({0}));
 //    REQUIRE(path_clusters.cluster_to_paths_index.at(1) == vector<uint32_t>({1, 2, 5}));
 //    REQUIRE(path_clusters.cluster_to_paths_index.at(2) == vector<uint32_t>({3, 6}));
 //    REQUIRE(path_clusters.cluster_to_paths_index.at(3) == vector<uint32_t>({4}));
}

// TEST_CASE("GBWT paths can be clustered") {

// 	gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
//     gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(7, true)));

//     gbwt::vector_type gbwt_thread_1(3);
//     gbwt::vector_type gbwt_thread_2(2);
//     gbwt::vector_type gbwt_thread_3(1);
//     gbwt::vector_type gbwt_thread_4(2);
   
//     gbwt_thread_1[0] = gbwt::Node::encode(1, false);
//     gbwt_thread_1[1] = gbwt::Node::encode(2, false);
//     gbwt_thread_1[2] = gbwt::Node::encode(4, false);

//     gbwt_thread_2[0] = gbwt::Node::encode(1, true);
//     gbwt_thread_2[1] = gbwt::Node::encode(6, true);

//     gbwt_thread_3[0] = gbwt::Node::encode(3, false);

//     gbwt_thread_4[0] = gbwt::Node::encode(6, false);
//     gbwt_thread_4[1] = gbwt::Node::encode(7, false);

//     gbwt_builder.insert(gbwt_thread_1, false);
//     gbwt_builder.insert(gbwt_thread_2, false);
//     gbwt_builder.insert(gbwt_thread_3, false);
//     gbwt_builder.insert(gbwt_thread_4, false);

//     gbwt_builder.index.addMetadata();

//     for (uint32_t i = 0; i < 4; ++i) {

//     	gbwt_builder.index.metadata.addPath(gbwt::PathName());
//     }    

//     gbwt_builder.finish();

//     std::stringstream gbwt_stream;
//     gbwt_builder.index.serialize(gbwt_stream);

//     gbwt::GBWT gbwt_index;
//     gbwt_index.load(gbwt_stream);

//     const string graph_str = R"(
//     	{
//     		"node": [
//     			{"id": 1, "sequence": "A"},
//     			{"id": 2, "sequence": "A"},
//     			{"id": 3, "sequence": "A"},
//     			{"id": 4, "sequence": "A"},
//     			{"id": 5, "sequence": "A"},
//     			{"id": 6, "sequence": "A"},
//     			{"id": 7, "sequence": "A"}
//     		],
//     	}
//     )";

// 	vg::Graph graph;
// 	json2pb(graph, graph_str);

//     PathsIndex paths_index(gbwt_index, graph);
//     REQUIRE(paths_index.index().metadata.paths() == 4);

// 	PathClusters path_clusters(1);
// 	auto node_to_path_index = path_clusters.findPathNodeClusters(paths_index);

//     REQUIRE(path_clusters.path_to_cluster_index.size() == 4);
//     REQUIRE(path_clusters.path_to_cluster_index == vector<uint32_t>({0, 0, 1, 0}));
//     REQUIRE(path_clusters.cluster_to_paths_index.size() == 2);
//     REQUIRE(path_clusters.cluster_to_paths_index.at(0) == vector<uint32_t>({0, 1, 3}));
//     REQUIRE(path_clusters.cluster_to_paths_index.at(1) == vector<uint32_t>({2}));

//     REQUIRE(node_to_path_index.size() == 6);    
//     REQUIRE(node_to_path_index.at(1) == 0);
//     REQUIRE(node_to_path_index.at(2) == 0);
//     REQUIRE(node_to_path_index.at(3) == 2);
//     REQUIRE(node_to_path_index.at(4) == 0);
//     REQUIRE(node_to_path_index.at(6) == 3);
//     REQUIRE(node_to_path_index.at(7) == 3);

//     SECTION("Bidirectionality does not affect clustering") {

//     	gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
// 	    gbwt::GBWTBuilder gbwt_builder_bidi(gbwt::bit_length(gbwt::Node::encode(7, true)));

// 	    gbwt_builder_bidi.insert(gbwt_thread_1, true);
// 	    gbwt_builder_bidi.insert(gbwt_thread_2, true);
// 	    gbwt_builder_bidi.insert(gbwt_thread_3, true);
// 	    gbwt_builder_bidi.insert(gbwt_thread_4, true);

// 	    gbwt_builder_bidi.index.addMetadata();

// 	    for (uint32_t i = 0; i < 4; ++i) {

// 	    	gbwt_builder_bidi.index.metadata.addPath(gbwt::PathName());
// 	    }    

// 	    gbwt_builder_bidi.finish();

// 	    std::stringstream gbwt_stream_bidi;
// 	    gbwt_builder_bidi.index.serialize(gbwt_stream_bidi);

// 	    gbwt::GBWT gbwt_index_bidi;
// 	    gbwt_index_bidi.load(gbwt_stream_bidi);

// 	    PathsIndex paths_index_bidi(gbwt_index_bidi, graph);
// 	    REQUIRE(paths_index_bidi.index().metadata.paths() == 4);

// 		PathClusters path_clusters_bidi(1);
// 		auto node_to_path_index_bidi = path_clusters_bidi.findPathNodeClusters(paths_index_bidi);

// 	    REQUIRE(path_clusters_bidi.path_to_cluster_index == path_clusters.path_to_cluster_index);
// 	    REQUIRE(path_clusters_bidi.path_to_cluster_index == path_clusters.path_to_cluster_index);

//     	REQUIRE(node_to_path_index_bidi.size() == 6);    
//     	REQUIRE(node_to_path_index_bidi.at(1) == node_to_path_index.at(1));
//     	REQUIRE(node_to_path_index_bidi.at(2) == node_to_path_index.at(2));
//     	REQUIRE(node_to_path_index_bidi.at(3) == node_to_path_index.at(3));
//     	REQUIRE(node_to_path_index_bidi.at(4) == node_to_path_index.at(4));
//     	REQUIRE(node_to_path_index_bidi.at(6) == 1);
//     	REQUIRE(node_to_path_index_bidi.at(7) == node_to_path_index.at(7));
//     }
// }

