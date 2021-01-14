
#include "catch.hpp"

#include "gbwt/dynamic_gbwt.h"

#include "../alignment_path_finder.hpp"
#include "../utils.hpp"


TEST_CASE("Alignment path(s) can be found from a single-end alignment") {
    
    const string graph_str = R"(
    	{
    		"node": [
    			{"id": 1, "sequence": "AAAA"},
    			{"id": 2, "sequence": "A"},
    			{"id": 3, "sequence": "A"},
    			{"id": 4, "sequence": "AAAAAAAA"}
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

    vector<uint32_t> node_frag_lengths = {0, 4, 1, 1, 8};
    function<size_t(const uint32_t)> node_frag_length_func = [&](const uint32_t node_id) { return node_frag_lengths.at(node_id); };

	gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(4, true)));

    gbwt::vector_type gbwt_thread_1(3);
    gbwt::vector_type gbwt_thread_2(2);
   
    gbwt_thread_1[0] = gbwt::Node::encode(1, false);
    gbwt_thread_1[1] = gbwt::Node::encode(2, false);
    gbwt_thread_1[2] = gbwt::Node::encode(4, false);

    gbwt_thread_2[0] = gbwt::Node::encode(1, false);
    gbwt_thread_2[1] = gbwt::Node::encode(2, false);

    gbwt_builder.insert(gbwt_thread_1, true);
    gbwt_builder.insert(gbwt_thread_2, false);

    gbwt_builder.finish();

    std::stringstream gbwt_stream;
    gbwt_builder.index.serialize(gbwt_stream);

    gbwt::GBWT gbwt_index;
    gbwt_index.load(gbwt_stream);

    const string alignment_1_str = R"(
        {
            "path": {
            	"mapping": [
                	{
                		"position": {"node_id": 1, "offset": 2},
                    	"edit": [
                        	{"from_length": 2, "to_length": 2}
                    	]
                	},
                	{
                		"position": {"node_id": 2},
                    	"edit": [
                        	{"from_length": 1, "to_length": 1}
                    	]
                	},
                	{
                		"position": {"node_id": 4},
                    	"edit": [
                            {"from_length": 1, "to_length": 1},
                            {"from_length": 2, "to_length": 2, "sequence": "AA"},
                            {"from_length": 2, "to_length": 2}
                    	]
                	}
                ]
           	},
            "sequence": "AAAAAAAA",
           	"mapping_quality": 10,
           	"score": 4
        }
    )";

    vg::Alignment alignment_1;
    json2pb(alignment_1, alignment_1_str);

    PathsIndex paths_index(gbwt_index, graph);
    REQUIRE(!paths_index.index().bidirectional());

    AlignmentPathFinder<vg::Alignment> alignment_path_finder(paths_index, "unstranded", 1000, 0, 0, 0, 1);

    auto alignment_paths = alignment_path_finder.findAlignmentPaths(alignment_1);
    REQUIRE(alignment_paths.size() == 2);

    SECTION("Single-end read alignment finds alignment path(s)") {    

        REQUIRE(alignment_paths.front().frag_length == 8);
        REQUIRE(alignment_paths.front().min_mapq == 10);
        REQUIRE(alignment_paths.front().score_sum == 4);
        REQUIRE(paths_index.locatePathIds(alignment_paths.front().search_state) == vector<gbwt::size_type>({0}));

        REQUIRE(alignment_paths.back().frag_length == alignment_paths.front().frag_length);
        REQUIRE(alignment_paths.back().min_mapq == alignment_paths.front().min_mapq);
        REQUIRE(alignment_paths.back().score_sum == alignment_paths.front().score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths.back().search_state) == vector<gbwt::size_type>({1}));
    }

    SECTION("Reverse-complement single-end read alignment finds alignment path(s)") {

        auto alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_frag_length_func);
        alignment_1_rc.set_sequence("AAAAAAAA");
        
        auto alignment_paths_rc = alignment_path_finder.findAlignmentPaths(alignment_1_rc);
        REQUIRE(alignment_paths_rc.size() == 2);

        REQUIRE(alignment_paths_rc == alignment_paths);
    }

    SECTION("Soft-clipped single-end read alignment finds alignment path(s)") {

        alignment_1.mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_from_length(1);
        alignment_1.mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_to_length(1);

        auto new_edit = alignment_1.mutable_path()->mutable_mapping(0)->add_edit();
        new_edit->set_from_length(0);
        new_edit->set_to_length(1);
        new_edit->set_sequence("C");

        alignment_1.mutable_path()->mutable_mapping(2)->mutable_edit(2)->set_from_length(0);
        alignment_1.mutable_path()->mutable_mapping(2)->mutable_edit(2)->set_to_length(2);
        alignment_1.mutable_path()->mutable_mapping(2)->mutable_edit(2)->set_sequence("CC");

        auto alignment_paths_sc = alignment_path_finder.findAlignmentPaths(alignment_1);
        REQUIRE(alignment_paths_sc.size() == 2);

        REQUIRE(alignment_paths_sc == alignment_paths);
    }

    SECTION("Alternative single-end read alignment finds empty alignment path") {

        alignment_1.mutable_path()->mutable_mapping(1)->mutable_position()->set_node_id(3);
        
        auto alignment_paths_alt = alignment_path_finder.findAlignmentPaths(alignment_1);
        REQUIRE(alignment_paths_alt.empty());
    }

    SECTION("Single-end read alignment finds forward alignment path(s) in bidirectional index") {

        gbwt::GBWTBuilder gbwt_builder_bd(gbwt::bit_length(gbwt::Node::encode(4, true)));

        gbwt_builder_bd.insert(gbwt_thread_1, true);
        gbwt_builder_bd.insert(gbwt_thread_2, true);

        gbwt_builder_bd.finish();

        std::stringstream gbwt_stream_bd;
        gbwt_builder_bd.index.serialize(gbwt_stream_bd);

        gbwt::GBWT gbwt_index_bd;
        gbwt_index_bd.load(gbwt_stream_bd);

        PathsIndex paths_index_bd(gbwt_index_bd, graph);
        REQUIRE(paths_index_bd.index().bidirectional() == true);

        AlignmentPathFinder<vg::Alignment> alignment_path_finder_bd(paths_index_bd, "unstranded", 1000, 0, 0, 0, 1);
    
        auto alignment_paths_bd = alignment_path_finder_bd.findAlignmentPaths(alignment_1);
        REQUIRE(alignment_paths_bd.size() == 1);

        REQUIRE(alignment_paths_bd.front().frag_length == alignment_paths.front().frag_length);
        REQUIRE(alignment_paths_bd.front().min_mapq == alignment_paths.front().min_mapq);
        REQUIRE(alignment_paths_bd.front().score_sum == alignment_paths.front().score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths_bd.front().search_state) == vector<gbwt::size_type>({0}));
    }
}
    
TEST_CASE("Alignment path(s) can be found from a paired-end alignment") {
    
    const string graph_str = R"(
        {
            "node": [
                {"id": 1, "sequence": "AAAA"},
                {"id": 2, "sequence": "A"},
                {"id": 3, "sequence": "A"},
                {"id": 4, "sequence": "AAAAAAAA"},
                {"id": 5, "sequence": "AA"},
                {"id": 6, "sequence": "AAAAAAA"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 3},
                {"from": 2, "to": 4},
                {"from": 3, "to": 4},
                {"from": 4, "to": 5},
                {"from": 2, "to": 6},
                {"from": 4, "to": 6},
                {"from": 5, "to": 6}     
            ]
        }
    )";

    vg::Graph graph;
    json2pb(graph, graph_str);

    vector<uint32_t> node_frag_lengths = {0, 4, 1, 1, 8, 2, 7};
    function<size_t(const uint32_t)> node_frag_length_func = [&](const uint32_t node_id) { return node_frag_lengths.at(node_id); };

    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(6, true)));

    gbwt::vector_type gbwt_thread_1(5);
    gbwt::vector_type gbwt_thread_2(4);
    gbwt::vector_type gbwt_thread_3(3);
   
    gbwt_thread_1[0] = gbwt::Node::encode(1, false);
    gbwt_thread_1[1] = gbwt::Node::encode(2, false);
    gbwt_thread_1[2] = gbwt::Node::encode(4, false);
    gbwt_thread_1[3] = gbwt::Node::encode(5, false);
    gbwt_thread_1[4] = gbwt::Node::encode(6, false);
    
    gbwt_thread_2[0] = gbwt::Node::encode(6, true);
    gbwt_thread_2[1] = gbwt::Node::encode(4, true);
    gbwt_thread_2[2] = gbwt::Node::encode(2, true);
    gbwt_thread_2[3] = gbwt::Node::encode(1, true);

    gbwt_thread_3[0] = gbwt::Node::encode(1, false);
    gbwt_thread_3[1] = gbwt::Node::encode(2, false);
    gbwt_thread_3[2] = gbwt::Node::encode(6, false);

    gbwt_builder.insert(gbwt_thread_1, false);
    gbwt_builder.insert(gbwt_thread_2, true);    
    gbwt_builder.insert(gbwt_thread_3, false);

    gbwt_builder.finish();

    std::stringstream gbwt_stream;
    gbwt_builder.index.serialize(gbwt_stream);

    gbwt::GBWT gbwt_index;
    gbwt_index.load(gbwt_stream);

    const string alignment_1_str = R"(
        {
            "path": {
                "mapping": [
                    {
                        "position": {"node_id": 1, "offset": 2},
                        "edit": [
                            {"from_length": 2, "to_length": 2}
                        ]
                    },
                    {
                        "position": {"node_id": 2},
                        "edit": [
                            {"from_length": 1, "to_length": 1}
                        ]
                    },
                    {
                        "position": {"node_id": 4},
                        "edit": [
                            {"from_length": 5, "to_length": 5}
                        ]
                    }
                ]
            },
            "sequence": "AAAAAAAA",
            "mapping_quality": 10,
            "score": 8
        }
    )";

    vg::Alignment alignment_1;
    json2pb(alignment_1, alignment_1_str);

    const string alignment_2_str = R"(
        {
            "path": {
                "mapping": [
                    {
                        "position": {"node_id": 6, "offset": 1, "is_reverse": true},
                        "edit": [
                            {"from_length": 2, "to_length": 2},
                            {"from_length": 1, "to_length": 1, "sequence": "A"},
                            {"from_length": 1, "to_length": 1}
                        ]
                    }
                ]
            },
            "sequence": "AAAA",
            "mapping_quality": 20,
            "score": 2
        }
    )";

    vg::Alignment alignment_2;
    json2pb(alignment_2, alignment_2_str);

    PathsIndex paths_index(gbwt_index, graph);
    REQUIRE(!paths_index.index().bidirectional());

    AlignmentPathFinder<vg::Alignment> alignment_path_finder(paths_index, "unstranded", 1000, 0, 0, 0, 1);
    
    auto alignment_paths = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
    REQUIRE(alignment_paths.size() == 3);

    SECTION("Paired-end read alignment finds alignment path(s)") {

        REQUIRE(alignment_paths.front().frag_length == 19);
        REQUIRE(alignment_paths.front().min_mapq == 10);
        REQUIRE(alignment_paths.front().score_sum == 10);
        REQUIRE(paths_index.locatePathIds(alignment_paths.front().search_state) == vector<gbwt::size_type>({0}));

        REQUIRE(alignment_paths.at(1).frag_length == 17);
        REQUIRE(alignment_paths.at(1).min_mapq == alignment_paths.front().min_mapq);
        REQUIRE(alignment_paths.at(1).score_sum == alignment_paths.front().score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths.at(1).search_state) == vector<gbwt::size_type>({2}));

        REQUIRE(alignment_paths.back().frag_length == alignment_paths.at(1).frag_length);
        REQUIRE(alignment_paths.back().min_mapq == alignment_paths.at(1).min_mapq);
        REQUIRE(alignment_paths.back().score_sum == alignment_paths.at(1).score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths.back().search_state) == vector<gbwt::size_type>({1}));
    }

    SECTION("Incorrect oriented paired-end read alignment finds empty alignment path") {

        auto alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_frag_length_func);
        alignment_2_rc.set_sequence("AAAA");

        auto alignment_paths_rc = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2_rc);
        REQUIRE(alignment_paths_rc.empty());
    }

    SECTION("Extended paired-end read alignment finds alignment path(s)") {
  
        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_from_length(3);
        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_to_length(3);

        auto new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(5);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(2);
        new_edit->set_to_length(2);

        alignment_2.set_sequence(alignment_2.sequence() + "AAAA");        

        auto alignment_paths_ext = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ext.size() == 1);

        REQUIRE(alignment_paths_ext.front() == alignment_paths.front());

        new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(4);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_2.set_sequence(alignment_2.sequence() + "A");        

        alignment_paths_ext = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ext.size() == 1);

        REQUIRE(alignment_paths_ext.front() == alignment_paths.front());             
    }

    SECTION("Partial overlapping paired-end read alignment finds alignment path(s)") {

        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_from_length(3);
        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_to_length(3);

        auto new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(4);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(5);
        new_edit->set_to_length(5);

        alignment_2.set_sequence(alignment_2.sequence() + "AAAAAAA");

        auto alignment_paths_ov = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 2);

        REQUIRE(alignment_paths_ov.front() == alignment_paths.at(1));
        REQUIRE(alignment_paths_ov.back() == alignment_paths.back());

        new_edit->set_from_length(8);
        new_edit->set_to_length(8);

        alignment_2.set_sequence(alignment_2.sequence() + "AAA");

        new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(2);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_2.set_sequence(alignment_2.sequence() + "A");

        alignment_paths_ov = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 2);

        REQUIRE(alignment_paths_ov.front() == alignment_paths.at(1));
        REQUIRE(alignment_paths_ov.back() == alignment_paths.back());

        new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(1);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_2.set_sequence(alignment_2.sequence() + "A");

        alignment_paths_ov = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 2);

        REQUIRE(alignment_paths_ov.front() == alignment_paths.at(1));
        REQUIRE(alignment_paths_ov.back() == alignment_paths.back());
    }

    SECTION("Perfect overlapping paired-end read alignment finds alignment path(s)") {

        auto alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_frag_length_func);
        alignment_1_rc.set_sequence("AAAAAAAA");

        auto alignment_paths_ov_1 = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_1_rc);
        REQUIRE(alignment_paths_ov_1.size() == 2);

        REQUIRE(alignment_paths_ov_1.front().frag_length == 8);
        REQUIRE(alignment_paths_ov_1.front().min_mapq == 10);
        REQUIRE(alignment_paths_ov_1.front().score_sum == 16);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_1.front().search_state) == vector<gbwt::size_type>({0, 2}));

        REQUIRE(alignment_paths_ov_1.back().frag_length == alignment_paths_ov_1.front().frag_length);
        REQUIRE(alignment_paths_ov_1.back().min_mapq == alignment_paths_ov_1.front().min_mapq);
        REQUIRE(alignment_paths_ov_1.back().score_sum == alignment_paths_ov_1.front().score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_1.back().search_state) == vector<gbwt::size_type>({1}));

        auto alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_frag_length_func);
        alignment_2_rc.set_sequence("AAAA");

        auto alignment_paths_ov_2 = alignment_path_finder.findPairedAlignmentPaths(alignment_2, alignment_2_rc);
        REQUIRE(alignment_paths_ov_2.size() == 2);

        REQUIRE(alignment_paths_ov_2.front().frag_length == 4);
        REQUIRE(alignment_paths_ov_2.front().min_mapq == 20);
        REQUIRE(alignment_paths_ov_2.front().score_sum == 4);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_2.front().search_state) == vector<gbwt::size_type>({1}));

        REQUIRE(alignment_paths_ov_2.back().frag_length == alignment_paths_ov_2.front().frag_length);
        REQUIRE(alignment_paths_ov_2.back().min_mapq == alignment_paths_ov_2.front().min_mapq);
        REQUIRE(alignment_paths_ov_2.back().score_sum == alignment_paths_ov_2.front().score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_2.back().search_state) == vector<gbwt::size_type>({0, 2, 3}));
    }

    SECTION("Incorrect overlapping paired-end read alignment finds empty alignment path") {

        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_from_length(3);
        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_to_length(3);

        auto new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(2);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_2.set_sequence(alignment_2.sequence() + "AAA");

        auto alignment_paths_ov = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.empty());
    }

    SECTION("Paired-end read alignment finds forward alignment path(s) in bidirectional index") {

        gbwt::GBWTBuilder gbwt_builder_bd(gbwt::bit_length(gbwt::Node::encode(6, true)));

        gbwt_builder_bd.insert(gbwt_thread_1, true);
        gbwt_builder_bd.insert(gbwt_thread_2, true);    
        gbwt_builder_bd.insert(gbwt_thread_3, true);

        gbwt_builder_bd.finish();

        std::stringstream gbwt_stream_bd;
        gbwt_builder_bd.index.serialize(gbwt_stream_bd);

        gbwt::GBWT gbwt_index_bd;
        gbwt_index_bd.load(gbwt_stream_bd);

        PathsIndex paths_index_bd(gbwt_index_bd, graph);
        REQUIRE(paths_index_bd.index().bidirectional() == true);

        AlignmentPathFinder<vg::Alignment> alignment_path_finder_bd(paths_index_bd, "unstranded", 1000, 0, 0, 0, 1);
    
        auto alignment_paths_bd = alignment_path_finder_bd.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_bd.size() == 2);

        REQUIRE(alignment_paths_bd.front() == alignment_paths.front());
        REQUIRE(alignment_paths_bd.at(1) == alignment_paths.at(1));     
    }
}

TEST_CASE("Circular alignment path(s) can be found from a paired-end alignment") {
    
    const string graph_str = R"(
        {
            "node": [
                {"id": 1, "sequence": "AAAA"},
                {"id": 2, "sequence": "AAAA"},
                {"id": 3, "sequence": "AAAA"},
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 2, "to": 2},
                {"from": 2, "to": 3},
            ]
        }
    )";

    vg::Graph graph;
    json2pb(graph, graph_str);

    vector<uint32_t> node_frag_lengths = {0, 4, 4, 4, 4};
    function<size_t(const uint32_t)> node_frag_length_func = [&](const uint32_t node_id) { return node_frag_lengths.at(node_id); };

    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(6, true)));

    gbwt::vector_type gbwt_thread_1(3);
    gbwt::vector_type gbwt_thread_2(5);
   
    gbwt_thread_1[0] = gbwt::Node::encode(1, false);
    gbwt_thread_1[1] = gbwt::Node::encode(2, false);
    gbwt_thread_1[2] = gbwt::Node::encode(3, false);

    gbwt_thread_2[0] = gbwt::Node::encode(1, false);
    gbwt_thread_2[1] = gbwt::Node::encode(2, false);
    gbwt_thread_2[2] = gbwt::Node::encode(2, false);
    gbwt_thread_2[3] = gbwt::Node::encode(2, false);
    gbwt_thread_2[4] = gbwt::Node::encode(3, false);

    gbwt_builder.insert(gbwt_thread_1, false);
    gbwt_builder.insert(gbwt_thread_2, true);    

    gbwt_builder.finish();

    std::stringstream gbwt_stream;
    gbwt_builder.index.serialize(gbwt_stream);

    gbwt::GBWT gbwt_index;
    gbwt_index.load(gbwt_stream);

    const string alignment_1_str = R"(
        {
            "path": {
                "mapping": [
                    {
                        "position": {"node_id": 1, "offset": 2},
                        "edit": [
                            {"from_length": 2, "to_length": 2}
                        ]
                    }
                ]
            },
            "sequence": "AA",
            "mapping_quality": 10,
            "score": 2 
        }
    )";

    vg::Alignment alignment_1;
    json2pb(alignment_1, alignment_1_str);

    const string alignment_2_str = R"(
        {
            "path": {
                "mapping": [
                    {
                        "position": {"node_id": 3, "offset": 0, "is_reverse": true},
                        "edit": [
                            {"from_length": 2, "to_length": 2}
                        ]
                    }
                ]
            },
            "sequence": "AA",
            "mapping_quality": 20,
            "score": 2 
        }
    )";

    vg::Alignment alignment_2;
    json2pb(alignment_2, alignment_2_str);

    PathsIndex paths_index(gbwt_index, graph);
    REQUIRE(!paths_index.index().bidirectional());

    AlignmentPathFinder<vg::Alignment> alignment_path_finder(paths_index, "unstranded", 1000, 0, 0, 0, 1);

    auto alignment_paths = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
    REQUIRE(alignment_paths.size() == 3);

    SECTION("Paired-end read alignment finds circular alignment path(s)") {

        REQUIRE(alignment_paths.front().frag_length == 18);
        REQUIRE(alignment_paths.front().min_mapq == 10);
        REQUIRE(alignment_paths.front().score_sum == 4);
        REQUIRE(paths_index.locatePathIds(alignment_paths.front().search_state) == vector<gbwt::size_type>({1}));                

        REQUIRE(alignment_paths.at(1).frag_length == 10);
        REQUIRE(alignment_paths.at(1).min_mapq == alignment_paths.front().min_mapq);
        REQUIRE(alignment_paths.at(1).score_sum == alignment_paths.front().score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths.at(1).search_state) == vector<gbwt::size_type>({0}));        

        REQUIRE(alignment_paths.back().frag_length == alignment_paths.front().frag_length);
        REQUIRE(alignment_paths.back().min_mapq == alignment_paths.at(1).min_mapq);
        REQUIRE(alignment_paths.back().score_sum == alignment_paths.at(1).score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths.back().search_state) == vector<gbwt::size_type>({2}));   
    }

    SECTION("Non-circular paired-end read alignment finds non-circular alignment path(s)") {

        auto new_mapping = alignment_1.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(2);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(false);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(4);
        new_edit->set_to_length(4);

        alignment_1.set_sequence(alignment_1.sequence() + "AAAA");

        new_mapping = alignment_1.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(3);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(false);

        new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_1.set_sequence(alignment_1.sequence() + "A");

        auto alignment_paths_ncirc = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ncirc.size() == 1);

        REQUIRE(alignment_paths_ncirc.front() == alignment_paths.at(1));
    }

    SECTION("Circular paired-end read alignment finds circular alignment path(s)") {

        auto new_mapping = alignment_1.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(2);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(false);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(4);
        new_edit->set_to_length(4);

        alignment_1.set_sequence(alignment_1.sequence() + "AAAA");

        for (uint32_t i = 0; i < 2; i++) {

            new_mapping = alignment_1.mutable_path()->add_mapping();
            new_mapping->mutable_position()->set_node_id(2);
            new_mapping->mutable_position()->set_offset(0);
            new_mapping->mutable_position()->set_is_reverse(false);

            new_edit = new_mapping->add_edit();
            new_edit->set_from_length(4);
            new_edit->set_to_length(4);

            alignment_1.set_sequence(alignment_1.sequence() + "AAAA");

            auto alignment_paths_circ = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
            REQUIRE(alignment_paths_circ.size() == 2);

            REQUIRE(alignment_paths_circ.front() == alignment_paths.front());
            REQUIRE(alignment_paths_circ.back() == alignment_paths.back());
        }
    }

    SECTION("Partial overlapping non-circular paired-end read alignment finds non-circular alignment path(s)") {

        auto new_mapping = alignment_1.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(2);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(false);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(4);
        new_edit->set_to_length(4);

        alignment_1.set_sequence(alignment_1.sequence() + "AAAA");

        new_mapping = alignment_1.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(3);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(false);

        new_edit = new_mapping->add_edit();
        new_edit->set_from_length(4);
        new_edit->set_to_length(4);

        alignment_1.set_sequence(alignment_1.sequence() + "AAAA");

        auto alignment_paths_ncirc = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ncirc.size() == 1);

        REQUIRE(alignment_paths_ncirc.front() == alignment_paths.at(1));
    }

    SECTION("Partial overlapping circular paired-end read alignment finds circular alignment path(s)") {

        for (uint32_t i = 0; i < 2; i++) {

            auto new_mapping = alignment_1.mutable_path()->add_mapping();
            new_mapping->mutable_position()->set_node_id(2);
            new_mapping->mutable_position()->set_offset(0);
            new_mapping->mutable_position()->set_is_reverse(false);

            auto new_edit = new_mapping->add_edit();
            new_edit->set_from_length(4);
            new_edit->set_to_length(4);
    
            alignment_1.set_sequence(alignment_1.sequence() + "AAAA");
        }

        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_from_length(4);
        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_to_length(4);

        alignment_2.set_sequence(alignment_2.sequence() + "AA");

        for (uint32_t i = 0; i < 3; i++) {

            auto new_mapping = alignment_2.mutable_path()->add_mapping();
            new_mapping->mutable_position()->set_node_id(2);
            new_mapping->mutable_position()->set_offset(0);
            new_mapping->mutable_position()->set_is_reverse(true);

            auto new_edit = new_mapping->add_edit();
            new_edit->set_from_length(4);
            new_edit->set_to_length(4);

            alignment_2.set_sequence(alignment_2.sequence() + "AAAA");
        }

        auto alignment_paths_circ = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_circ.size() == 2);

        REQUIRE(alignment_paths_circ.front() == alignment_paths.front());
        REQUIRE(alignment_paths_circ.back() == alignment_paths.back());
    }

    SECTION("Circular paired-end read alignment finds forward alignment path(s) in bidirectional index") {

        gbwt::GBWTBuilder gbwt_builder_bd(gbwt::bit_length(gbwt::Node::encode(6, true)));

        gbwt_builder_bd.insert(gbwt_thread_1, true);
        gbwt_builder_bd.insert(gbwt_thread_2, true); 

        gbwt_builder_bd.finish();

        std::stringstream gbwt_stream_bd;
        gbwt_builder_bd.index.serialize(gbwt_stream_bd);

        gbwt::GBWT gbwt_index_bd;
        gbwt_index_bd.load(gbwt_stream_bd);

        PathsIndex paths_index_bd(gbwt_index_bd, graph);
        REQUIRE(paths_index_bd.index().bidirectional() == true);

        AlignmentPathFinder<vg::Alignment> alignment_path_finder_bd(paths_index_bd, "unstranded", 1000, 0, 0, 0, 1);
    
        auto alignment_paths_bd = alignment_path_finder_bd.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_bd.size() == 2);

        REQUIRE(alignment_paths_bd.front() == alignment_paths.front());

        REQUIRE(alignment_paths_bd.back().frag_length == alignment_paths.at(1).frag_length);
        REQUIRE(alignment_paths_bd.back().min_mapq == alignment_paths.at(1).min_mapq);
        REQUIRE(alignment_paths_bd.back().score_sum == alignment_paths.at(1).score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths_bd.back().search_state) == vector<gbwt::size_type>({0}));                
    }
}

TEST_CASE("Alignment path(s) can be found from a single-end multipath alignment") {

    const string graph_str = R"(
        {
            "node": [
                {"id": 1, "sequence": "A"},
                {"id": 2, "sequence": "A"},
                {"id": 3, "sequence": "AAA"},
                {"id": 4, "sequence": "AA"},
                {"id": 5, "sequence": "AAA"},
                {"id": 6, "sequence": "AAA"},
            ],
            "edge": [
                {"from": 1, "to": 3},
                {"from": 2, "to": 3},
                {"from": 3, "to": 4},
                {"from": 4, "to": 5},
                {"from": 5, "to": 6}
            ]
        }
    )";

    vg::Graph graph;
    json2pb(graph, graph_str);

    vector<uint32_t> node_frag_lengths = {0, 1, 1, 3, 2, 3, 3};
    function<size_t(const uint32_t)> node_frag_length_func = [&](const uint32_t node_id) { return node_frag_lengths.at(node_id); };

    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(6, true)));

    gbwt::vector_type gbwt_thread_1(4);
    gbwt::vector_type gbwt_thread_2(4);
   
    gbwt_thread_1[0] = gbwt::Node::encode(1, false);
    gbwt_thread_1[1] = gbwt::Node::encode(3, false);
    gbwt_thread_1[2] = gbwt::Node::encode(4, false);
    gbwt_thread_1[3] = gbwt::Node::encode(5, false);

    gbwt_thread_2[0] = gbwt::Node::encode(6, true);
    gbwt_thread_2[1] = gbwt::Node::encode(4, true);
    gbwt_thread_2[2] = gbwt::Node::encode(3, true);
    gbwt_thread_2[3] = gbwt::Node::encode(1, true);

    gbwt_builder.insert(gbwt_thread_1, false);
    gbwt_builder.insert(gbwt_thread_2, false);

    gbwt_builder.finish();

    std::stringstream gbwt_stream;
    gbwt_builder.index.serialize(gbwt_stream);

    gbwt::GBWT gbwt_index;
    gbwt_index.load(gbwt_stream);

    const string alignment_1_str = R"(
        {
            "start": [0,1],
            "subpath": [
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 1},
                                "edit": [
                                    {"from_length": 1, "to_length": 1}
                                ]
                            }
                        ]
                    },
                    "next": [2],
                    "score": 1
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 2},
                                "edit": [
                                    {"from_length": 1, "to_length": 1, "sequence": "A"}
                                ]
                            }
                        ]
                    },
                    "next": [2],
                    "score": -1
                },
                {                
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 3},
                                "edit": [
                                    {"from_length": 3, "to_length": 3}
                                ]
                            },
                            {
                                "position": {"node_id": 4},
                                "edit": [
                                    {"from_length": 2, "to_length": 2}
                                ]
                            }
                        ]
                    },
                    "next": [3,4],
                    "score": 5
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 5},
                                "edit": [
                                    {"from_length": 2, "to_length": 2}
                                ]
                            }
                        ]
                    },
                    "score": 2
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 6},
                                "edit": [
                                    {"from_length": 1, "to_length": 1, "sequence": "A"},
                                    {"from_length": 1, "to_length": 1}
                                ]
                            }
                        ]
                    },
                    "score": 0
                }
            ],
            "sequence": "AAAAAAAA",
            "mapping_quality": 10
        }
    )";

    vg::MultipathAlignment alignment_1;
    json2pb(alignment_1, alignment_1_str);

    PathsIndex paths_index(gbwt_index, graph);
    REQUIRE(!paths_index.index().bidirectional());

    AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder(paths_index, "unstranded", 1000, 0, 0, 0, 1);
    
    auto alignment_paths = alignment_path_finder.findAlignmentPaths(alignment_1);
    REQUIRE(alignment_paths.size() == 2);

    SECTION("Single-end multipath read alignment finds alignment path(s)") {

        REQUIRE(alignment_paths.front().frag_length == 8);
        REQUIRE(alignment_paths.front().min_mapq == 10);
        REQUIRE(alignment_paths.front().score_sum == 8);
        REQUIRE(paths_index.locatePathIds(alignment_paths.front().search_state) == vector<gbwt::size_type>({0}));

        REQUIRE(alignment_paths.back().frag_length == alignment_paths.front().frag_length);
        REQUIRE(alignment_paths.back().min_mapq == alignment_paths.front().min_mapq);
        REQUIRE(alignment_paths.back().score_sum == 6);
        REQUIRE(paths_index.locatePathIds(alignment_paths.back().search_state) == vector<gbwt::size_type>({1}));               
    }

    SECTION("Reverse-complement single-end multipath read alignment finds alignment path(s)") {

        auto alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_frag_length_func);
        alignment_1_rc.set_sequence("AAAAAAAA");
        
        auto alignment_paths_rc = alignment_path_finder.findAlignmentPaths(alignment_1_rc);
        REQUIRE(alignment_paths_rc.size() == 2);

        REQUIRE(alignment_paths_rc.front() == alignment_paths.front());
        REQUIRE(alignment_paths_rc.back() == alignment_paths.back());
    }

    SECTION("Soft-clipped single-end multipath read alignment finds alignment path(s)") {

        alignment_1.mutable_subpath(3)->mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_from_length(1);
        alignment_1.mutable_subpath(3)->mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_to_length(1);

        auto new_edit = alignment_1.mutable_subpath(3)->mutable_path()->mutable_mapping(0)->add_edit();
        new_edit->set_from_length(0);
        new_edit->set_to_length(1);
        new_edit->set_sequence("A");

        auto alignment_paths_sc = alignment_path_finder.findAlignmentPaths(alignment_1);
        REQUIRE(alignment_paths_sc.size() == 2);

        REQUIRE(alignment_paths_sc == alignment_paths);
    }

    SECTION("Single-end multipath read alignment finds forward alignment path(s) in bidirectional index") {

        gbwt::GBWTBuilder gbwt_builder_bd(gbwt::bit_length(gbwt::Node::encode(6, true)));

        gbwt_builder_bd.insert(gbwt_thread_1, true);
        gbwt_builder_bd.insert(gbwt_thread_2, true);

        gbwt_builder_bd.finish();

        std::stringstream gbwt_stream_bd;
        gbwt_builder_bd.index.serialize(gbwt_stream_bd);

        gbwt::GBWT gbwt_index_bd;
        gbwt_index_bd.load(gbwt_stream_bd);

        PathsIndex paths_index_bd(gbwt_index_bd, graph);
        REQUIRE(paths_index_bd.index().bidirectional() == true);

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_bd(paths_index_bd, "unstranded", 1000, 0, 0, 0, 1);
    
        auto alignment_paths_bd = alignment_path_finder_bd.findAlignmentPaths(alignment_1);
        REQUIRE(alignment_paths_bd.size() == 2);

        REQUIRE(alignment_paths_bd.front().frag_length == alignment_paths.back().frag_length);
        REQUIRE(alignment_paths_bd.front().min_mapq == alignment_paths.back().min_mapq);
        REQUIRE(alignment_paths_bd.front().score_sum == alignment_paths.back().score_sum);
        REQUIRE(paths_index_bd.locatePathIds(alignment_paths_bd.back().search_state) == vector<gbwt::size_type>({0}));
        REQUIRE(alignment_paths_bd.back() == alignment_paths.front());
    }
}

TEST_CASE("Alignment path(s) can be found from a paired-end multipath alignment") {

    const string graph_str = R"(
        {
            "node": [
                {"id": 1, "sequence": "A"},
                {"id": 2, "sequence": "AAAA"},
                {"id": 3, "sequence": "AA"},
                {"id": 4, "sequence": "AAAA"},
                {"id": 5, "sequence": "AA"},
                {"id": 6, "sequence": "A"},
                {"id": 7, "sequence": "AA"},
                {"id": 8, "sequence": "AAA"},
            ],
            "edge": [
                {"from": 1, "to": 3},
                {"from": 2, "to": 3},
                {"from": 3, "to": 4},
                {"from": 3, "to": 5},
                {"from": 4, "to": 5},
                {"from": 5, "to": 6},
                {"from": 5, "to": 7},
                {"from": 6, "to": 8},
                {"from": 7, "to": 8}
            ]
        }
    )";

    vg::Graph graph;
    json2pb(graph, graph_str);

    vector<uint32_t> node_frag_lengths = {0, 1, 4, 2, 4, 2, 1, 2, 3};
    function<size_t(const uint32_t)> node_frag_length_func = [&](const uint32_t node_id) { return node_frag_lengths.at(node_id); };

    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(8, true)));

    gbwt::vector_type gbwt_thread_1(5);
    gbwt::vector_type gbwt_thread_2(6);
   
    gbwt_thread_1[0] = gbwt::Node::encode(1, false);
    gbwt_thread_1[1] = gbwt::Node::encode(3, false);
    gbwt_thread_1[2] = gbwt::Node::encode(5, false);
    gbwt_thread_1[3] = gbwt::Node::encode(6, false);
    gbwt_thread_1[4] = gbwt::Node::encode(8, false);

    gbwt_thread_2[0] = gbwt::Node::encode(2, false);
    gbwt_thread_2[1] = gbwt::Node::encode(3, false);
    gbwt_thread_2[2] = gbwt::Node::encode(4, false);
    gbwt_thread_2[3] = gbwt::Node::encode(5, false);
    gbwt_thread_2[4] = gbwt::Node::encode(7, false);
    gbwt_thread_2[5] = gbwt::Node::encode(8, false);

    gbwt_builder.insert(gbwt_thread_1, false);
    gbwt_builder.insert(gbwt_thread_2, true);

    gbwt_builder.finish();

    std::stringstream gbwt_stream;
    gbwt_builder.index.serialize(gbwt_stream);

    gbwt::GBWT gbwt_index;
    gbwt_index.load(gbwt_stream);

    const string alignment_1_str = R"(
        {
            "start": [0,1],
            "subpath": [
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 1},
                                "edit": [
                                    {"to_length": 3, "sequence": "AAA"},
                                    {"from_length": 1, "to_length": 1}
                                ]
                            }
                        ]
                    },
                    "next": [2],
                    "score": 1
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 2},
                                "edit": [
                                    {"from_length": 4, "to_length": 4}
                                ]
                            }
                        ]
                    },
                    "next": [2],
                    "score": 4
                },
                {                
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 3},
                                "edit": [
                                    {"from_length": 2, "to_length": 2}
                                ]
                            },
                        ]
                    },
                    "score": 2
                }
            ],
            "sequence": "AAAAAA",
            "mapping_quality": 10
        }
    )";

    vg::MultipathAlignment alignment_1;
    json2pb(alignment_1, alignment_1_str);

    const string alignment_2_str = R"(
        {
            "start": [0],
            "subpath": [
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 8, "offset": 2, "is_reverse": true},
                                "edit": [
                                    {"from_length": 1, "to_length": 1}
                                ]
                            }
                        ]
                    },
                    "next": [1,4],
                    "score": 1
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 7, "is_reverse": true},
                                "edit": [
                                    {"from_length": 1, "to_length": 1}
                                ]
                            }
                        ]
                    },
                    "next": [2],
                    "score": 1
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 7, "offset": 1, "is_reverse": true},
                                "edit": [
                                    {"to_length": 1, "sequence": "A"}
                                ]
                            }
                        ]
                    },
                    "next": [3],
                    "score": -1
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 7, "offset": 1, "is_reverse": true},
                                "edit": [
                                    {"from_length": 1, "to_length": 1}
                                ]
                            }
                        ]
                    },
                    "next": [7],
                    "score": 1
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 6, "is_reverse": true},
                                "edit": [
                                    {"to_length": 2, "sequence": "AA"}
                                ]
                            }
                        ]
                    },
                    "next": [5],
                    "score": -2
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 6, "is_reverse": true},
                                "edit": [
                                    {"from_length": 1}
                                ]
                            }
                        ]
                    },
                    "next": [6],
                    "score": -1
                },
                                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 6, "offset": 1, "is_reverse": true},
                                "edit": [
                                    {"to_length": 1, "sequence": "A"}
                                ]
                            }
                        ]
                    },
                    "next": [7],
                    "score": -1
                },
                {                
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 5, "is_reverse": true},
                                "edit": [
                                    {"from_length": 1, "to_length": 1},
                                    {"to_length": 2, "sequence": "AA"}
                                ]
                            },
                        ]
                    },
                    "score": 1
                }
            ],
            "sequence": "AAAAAAA",
            "mapping_quality": 20
        }
    )";

    vg::MultipathAlignment alignment_2;
    json2pb(alignment_2, alignment_2_str);

    PathsIndex paths_index(gbwt_index, graph);
    REQUIRE(!paths_index.index().bidirectional());

    AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder(paths_index, "unstranded", 1000, 0, 0, 0, 1);

    auto alignment_paths = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
    REQUIRE(alignment_paths.size() == 3);

    SECTION("Paired-end multipath read alignment finds alignment path(s)") {

        REQUIRE(alignment_paths.front().frag_length == 16);
        REQUIRE(alignment_paths.front().min_mapq == 10);
        REQUIRE(alignment_paths.front().score_sum == 9);
        REQUIRE(paths_index.locatePathIds(alignment_paths.front().search_state) == vector<gbwt::size_type>({1}));                

        REQUIRE(alignment_paths.at(1).frag_length == 12);
        REQUIRE(alignment_paths.at(1).min_mapq == alignment_paths.front().min_mapq);
        REQUIRE(alignment_paths.at(1).score_sum == 1);
        REQUIRE(paths_index.locatePathIds(alignment_paths.at(1).search_state) == vector<gbwt::size_type>({0})); 

        REQUIRE(alignment_paths.back().frag_length == alignment_paths.front().frag_length);
        REQUIRE(alignment_paths.back().min_mapq == alignment_paths.front().min_mapq);
        REQUIRE(alignment_paths.back().score_sum == alignment_paths.front().score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths.back().search_state) == vector<gbwt::size_type>({2}));                              
    }

    SECTION("Incorrect oriented paired-end multipath read alignment finds empty alignment path") {

        auto alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_frag_length_func);
        alignment_2_rc.set_sequence("AAAAAAA");
        
        auto alignment_paths_rc = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2_rc);
        REQUIRE(alignment_paths_rc.empty());
    }

    SECTION("Extended paired-end multipath read alignment finds alignment path(s)") {

        alignment_1.mutable_subpath(2)->add_next(3);

        auto new_subpath = alignment_1.add_subpath();
        new_subpath->set_score(0);

        auto new_mapping = new_subpath->mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(4);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(false);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(2);
        new_edit->set_to_length(2);

        alignment_1.set_sequence(alignment_1.sequence() + "AA");

        auto alignment_paths_ext = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ext.size() == 2);

        REQUIRE(alignment_paths_ext.front() == alignment_paths.front());
        REQUIRE(alignment_paths_ext.back() == alignment_paths.back());
    }

    SECTION("Partial overlapping paired-end read alignment finds alignment path(s)") {

        alignment_1.mutable_subpath(2)->add_next(3);

        auto new_subpath = alignment_1.add_subpath();
        new_subpath->set_score(0);

        auto new_mapping = new_subpath->mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(5);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(false);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_1.set_sequence(alignment_1.sequence() + "A");

        auto alignment_paths_ov = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 1);

        REQUIRE(alignment_paths_ov.front() == alignment_paths.at(1));

        new_edit->set_from_length(2);
        new_edit->set_to_length(2);

        alignment_1.set_sequence(alignment_1.sequence() + "A");

        alignment_paths_ov = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 1);

        REQUIRE(alignment_paths_ov.front() == alignment_paths.at(1));

        alignment_1.mutable_subpath(3)->add_next(4);

        new_subpath = alignment_1.add_subpath();
        new_subpath->set_score(0);

        new_mapping = new_subpath->mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(6);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(false);

        new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_1.set_sequence(alignment_1.sequence() + "A");

        alignment_paths_ov = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 1);

        REQUIRE(alignment_paths_ov.front() == alignment_paths.at(1));

        new_edit->set_to_length(0);

        alignment_1.mutable_subpath(4)->add_next(5);

        new_subpath = alignment_1.add_subpath();
        new_subpath->set_score(0);

        new_mapping = new_subpath->mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(6);
        new_mapping->mutable_position()->set_offset(1);
        new_mapping->mutable_position()->set_is_reverse(false);

        new_edit = new_mapping->add_edit();
        new_edit->set_to_length(1);

        alignment_paths_ov = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 1);

        REQUIRE(alignment_paths_ov.front().frag_length == 11);
        REQUIRE(alignment_paths_ov.front().min_mapq == alignment_paths.at(1).min_mapq);
        REQUIRE(alignment_paths_ov.front().score_sum == alignment_paths.at(1).score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov.front().search_state) == vector<gbwt::size_type>({0}));   

        alignment_1.mutable_subpath(5)->add_next(6);

        new_subpath = alignment_1.add_subpath();
        new_subpath->set_score(0);

        new_mapping = new_subpath->mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(8);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(false);

        new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_1.set_sequence(alignment_1.sequence() + "A");

        alignment_paths_ov = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 1);

        REQUIRE(alignment_paths_ov.front() == alignment_paths.at(1));
    }

    SECTION("Perfect overlapping paired-end multipath read alignment finds alignment path(s)") {

        auto alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_frag_length_func);
        alignment_1_rc.set_sequence("AAAAAA");

        auto alignment_paths_ov_1 = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_1_rc);
        REQUIRE(alignment_paths_ov_1.size() == 3);

        REQUIRE(alignment_paths_ov_1.front().frag_length == 6);
        REQUIRE(alignment_paths_ov_1.front().min_mapq == 10);
        REQUIRE(alignment_paths_ov_1.front().score_sum == 12);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_1.front().search_state) == vector<gbwt::size_type>({1}));

        REQUIRE(alignment_paths_ov_1.at(1).frag_length == alignment_paths_ov_1.front().frag_length);
        REQUIRE(alignment_paths_ov_1.at(1).min_mapq == alignment_paths_ov_1.front().min_mapq);
        REQUIRE(alignment_paths_ov_1.at(1).score_sum == 6);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_1.at(1).search_state) == vector<gbwt::size_type>({0}));                

        REQUIRE(alignment_paths_ov_1.back().frag_length == alignment_paths_ov_1.at(1).frag_length);
        REQUIRE(alignment_paths_ov_1.back().min_mapq == alignment_paths_ov_1.at(1).min_mapq);
        REQUIRE(alignment_paths_ov_1.back().score_sum == alignment_paths_ov_1.front().score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_1.back().search_state) == vector<gbwt::size_type>({2}));                

        auto alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_frag_length_func);
        alignment_2_rc.set_sequence("AAAAAAA");

        auto alignment_paths_ov_2 = alignment_path_finder.findPairedAlignmentPaths(alignment_2, alignment_2_rc);
        REQUIRE(alignment_paths_ov_2.size() == 3);

        REQUIRE(alignment_paths_ov_2.front().frag_length == 8);
        REQUIRE(alignment_paths_ov_2.front().min_mapq == 20);
        REQUIRE(alignment_paths_ov_2.front().score_sum == 6);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_2.front().search_state) == vector<gbwt::size_type>({1}));                

        REQUIRE(alignment_paths_ov_2.at(1).frag_length == 9);
        REQUIRE(alignment_paths_ov_2.at(1).min_mapq == alignment_paths_ov_2.front().min_mapq);
        REQUIRE(alignment_paths_ov_2.at(1).score_sum == 0);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_2.at(1).search_state) == vector<gbwt::size_type>({0}));                

        REQUIRE(alignment_paths_ov_2.back().frag_length == alignment_paths_ov_2.front().frag_length);
        REQUIRE(alignment_paths_ov_2.back().min_mapq == alignment_paths_ov_2.at(1).min_mapq);
        REQUIRE(alignment_paths_ov_2.back().score_sum == alignment_paths_ov_2.front().score_sum);
        REQUIRE(paths_index.locatePathIds(alignment_paths_ov_2.back().search_state) == vector<gbwt::size_type>({2}));                
    }

    SECTION("Paired-end multipath read alignment finds forward alignment path(s) in bidirectional index") {

        gbwt::GBWTBuilder gbwt_builder_bd(gbwt::bit_length(gbwt::Node::encode(8, true)));

        gbwt_builder_bd.insert(gbwt_thread_1, true);
        gbwt_builder_bd.insert(gbwt_thread_2, true);

        gbwt_builder_bd.finish();

        std::stringstream gbwt_stream_bd;
        gbwt_builder_bd.index.serialize(gbwt_stream_bd);

        gbwt::GBWT gbwt_index_bd;
        gbwt_index_bd.load(gbwt_stream_bd);
        
        PathsIndex paths_index_bd(gbwt_index_bd, graph);
        REQUIRE(paths_index_bd.index().bidirectional() == true);

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_bd(paths_index_bd, "unstranded", 1000, 0, 0, 0, 1);
    
        auto alignment_paths_bd = alignment_path_finder_bd.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_bd.size() == 2);

        REQUIRE(alignment_paths_bd.front() == alignment_paths.front());
        REQUIRE(alignment_paths_bd.back() == alignment_paths.at(1));
    }

    SECTION("Strand-specific paired-end multipath read alignment finds unidirectional alignment path(s)") {

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_fr(paths_index, "fr", 1000, 0, 0, 0, 1);

        auto alignment_paths_fr = alignment_path_finder_fr.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_fr.size() == 2);

        REQUIRE(alignment_paths_fr.front() == alignment_paths.front());
        REQUIRE(alignment_paths_fr.back() == alignment_paths.at(1));

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_rf(paths_index, "rf", 1000, 0, 0, 0, 1);

        auto alignment_paths_rf = alignment_path_finder_rf.findPairedAlignmentPaths(alignment_1, alignment_2);
        REQUIRE(alignment_paths_rf.size() == 1);

        REQUIRE(alignment_paths_rf.front() == alignment_paths.back());
    }

    SECTION("Alignment pairs from a paired-end multipath alignment are filtered based on length") {

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_len16(paths_index, "unstranded", 16, 0, 0, 0, 1);

        auto alignment_paths_len16 = alignment_path_finder_len16.findPairedAlignmentPaths(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_len16.size() == 3);
        
        REQUIRE(alignment_paths_len16 == alignment_paths);

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_len12(paths_index, "unstranded", 12, 0, 0, 0, 1);

        auto alignment_paths_len12 = alignment_path_finder_len12.findPairedAlignmentPaths(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_len12.size() == 1);

        REQUIRE(alignment_paths_len12.front() == alignment_paths.at(1));
        
        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_len11(paths_index, "unstranded", 11, 0, 0, 0, 1);

        auto alignment_paths_len11 = alignment_path_finder_len11.findPairedAlignmentPaths(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_len11.empty());
    }

    SECTION("Alignment pairs from a paired-end multipath alignment are filtered based on mapping quality") {

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_mq10(paths_index, "unstranded", 1000, 0, 10, 0, 1);

        auto alignment_paths_mq10 = alignment_path_finder_mq10.findPairedAlignmentPaths(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_mq10.size() == 3);

        assert(alignment_paths_mq10 == alignment_paths);

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_mq11(paths_index, "unstranded", 1000, 0, 11, 0, 1);

        auto alignment_paths_mq11 = alignment_path_finder_mq11.findPairedAlignmentPaths(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_mq11.empty());
    }

    SECTION("Alignment pairs from a paired-end multipath alignment are filtered based on best score fraction") {

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_bs40(paths_index, "unstranded", 1000, 0, 0, 0.40, 1);

        auto alignment_paths_bs40 = alignment_path_finder_bs40.findPairedAlignmentPaths(alignment_1, alignment_2);    
        REQUIRE(alignment_paths_bs40.size() == 3);

        assert(alignment_paths_bs40 == alignment_paths);

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_bs45(paths_index, "unstranded", 1000, 0, 0, 0.45, 1);

        auto alignment_paths_bs45 = alignment_path_finder_bs45.findPairedAlignmentPaths(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_bs45.empty());
    }

    SECTION("Alignment pairs from a paired-end multipath alignment are filtered based on soft-clipping length") {

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_sc30(paths_index, "unstranded", 1000, 0, 0, 0, 0.30);

        auto alignment_paths_sc30 = alignment_path_finder_sc30.findPairedAlignmentPaths(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_sc30.size() == 3);

        assert(alignment_paths_sc30 == alignment_paths);

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder_sc25(paths_index, "unstranded", 1000, 0, 0, 0, 0.25);

        auto alignment_paths_sc25 = alignment_path_finder_sc25.findPairedAlignmentPaths(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_sc25.empty());
    }
}

