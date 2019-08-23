
#include "catch2/catch.hpp"

#include "gbwt/dynamic_gbwt.h"

#include "../alignment_path_finder.hpp"
#include "../utils.hpp"


TEST_CASE("Alignment paths can be found from single-end alignment") {
    
    const string graph_str = R"(
    	{
    		"node": [
    			{"id": 1, "sequence": "GGGG"},
    			{"id": 2, "sequence": "A"},
    			{"id": 3, "sequence": "C"},
    			{"id": 4, "sequence": "TTTTTTTT"}
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

    vector<int32_t> node_seq_lengths = {0, 4, 1, 1, 8};
    function<size_t(const int64_t)> node_seq_length_func = [&](const int64_t node_id) { return node_seq_lengths.at(node_id); };

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
                            {"from_length": 2, "to_length": 2, "sequence": "AG"},
                            {"from_length": 2, "to_length": 2}
                    	]
                	}
                ]
           	},
           	"mapping_quality": 10,
           	"score":1 
        }
    )";

    vg::Alignment alignment_1;
    json2pb(alignment_1, alignment_1_str);

    AlignmentPathFinder<vg::Alignment> alignment_path_finder(graph, gbwt_index, 1000);
    
    auto alignment_paths = alignment_path_finder.findAlignmentPathsIds(alignment_1);
    REQUIRE(alignment_paths.size() == 1);

    REQUIRE(alignment_paths.front().search.node == gbwt::Node::encode(4, false));
    REQUIRE(alignment_paths.front().search.size() == 1);   
    REQUIRE(alignment_paths.front().end_offset == 5);   
    REQUIRE(alignment_paths.front().ids == vector<gbwt::size_type>({0}));
    REQUIRE(alignment_paths.front().seq_length == 8);
    REQUIRE(alignment_paths.front().mapqs == vector<int32_t>({10}));
    REQUIRE(alignment_paths.front().scores == vector<int32_t>({1}));

    SECTION("Reverse-complement single-end read alignment finds alignment path") {

        alignment_1 = lazy_reverse_complement_alignment(alignment_1, node_seq_length_func);
        
        auto alignment_paths_rc = alignment_path_finder.findAlignmentPathsIds(alignment_1);
        REQUIRE(alignment_paths_rc.size() == 1);

        REQUIRE(alignment_paths_rc.front().search.node == gbwt::Node::encode(1, true));
        REQUIRE(alignment_paths_rc.front().search.size() == alignment_paths.front().search.size());
        REQUIRE(alignment_paths_rc.front().end_offset == 2);   
        REQUIRE(alignment_paths_rc.front().ids == vector<gbwt::size_type>({1}));
        REQUIRE(alignment_paths_rc.front().seq_length == alignment_paths.front().seq_length);
        REQUIRE(alignment_paths_rc.front().mapqs == alignment_paths.front().mapqs);
        REQUIRE(alignment_paths_rc.front().scores == alignment_paths.front().scores);
    }

    SECTION("Soft-clipped single-end read alignment finds alignment path") {

        alignment_1.mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_from_length(1);
        alignment_1.mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_to_length(1);

        auto new_edit = alignment_1.mutable_path()->mutable_mapping(0)->add_edit();
        new_edit->set_from_length(0);
        new_edit->set_to_length(1);
        new_edit->set_sequence("C");

        alignment_1.mutable_path()->mutable_mapping(2)->mutable_edit(2)->set_from_length(0);
        alignment_1.mutable_path()->mutable_mapping(2)->mutable_edit(2)->set_to_length(2);
        alignment_1.mutable_path()->mutable_mapping(2)->mutable_edit(2)->set_sequence("CC");

        auto alignment_paths_sc = alignment_path_finder.findAlignmentPathsIds(alignment_1);
        REQUIRE(alignment_paths_sc.size() == 1);

        REQUIRE(alignment_paths_sc.front().search.node == alignment_paths.front().search.node);
        REQUIRE(alignment_paths_sc.front().search.size() == alignment_paths.front().search.size());
        REQUIRE(alignment_paths_sc.front().end_offset == 3);   
        REQUIRE(alignment_paths_sc.front().ids == alignment_paths.front().ids);
        REQUIRE(alignment_paths_sc.front().seq_length == alignment_paths.front().seq_length);
        REQUIRE(alignment_paths_sc.front().mapqs == alignment_paths.front().mapqs);
        REQUIRE(alignment_paths_sc.front().scores == alignment_paths.front().scores);
    }

    SECTION("Alternative single-end read alignment does not find alignment path") {

        alignment_1.mutable_path()->mutable_mapping(1)->mutable_position()->set_node_id(3);
        
        auto alignment_paths_alt = alignment_path_finder.findAlignmentPathsIds(alignment_1);
        REQUIRE(alignment_paths_alt.empty());
    }
}

TEST_CASE("Alignment paths can be found from paired-end alignment") {
    
    const string graph_str = R"(
        {
            "node": [
                {"id": 1, "sequence": "GGGG"},
                {"id": 2, "sequence": "A"},
                {"id": 3, "sequence": "C"},
                {"id": 4, "sequence": "TTTTTTTT"},
                {"id": 5, "sequence": "CC"},
                {"id": 6, "sequence": "AAAAAAA"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 3},
                {"from": 2, "to": 4},
                {"from": 3, "to": 4},
                {"from": 4, "to": 5},
                {"from": 4, "to": 6},
                {"from": 5, "to": 6}     
            ]
        }
    )";

    vg::Graph graph;
    json2pb(graph, graph_str);

    vector<int32_t> node_seq_lengths = {0, 4, 1, 1, 8, 2, 7};
    function<size_t(const int64_t)> node_seq_length_func = [&](const int64_t node_id) { return node_seq_lengths.at(node_id); };

    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(6, true)));

    gbwt::vector_type gbwt_thread_1(5);
    gbwt::vector_type gbwt_thread_2(4);
   
    gbwt_thread_1[0] = gbwt::Node::encode(1, false);
    gbwt_thread_1[1] = gbwt::Node::encode(2, false);
    gbwt_thread_1[2] = gbwt::Node::encode(4, false);
    gbwt_thread_1[3] = gbwt::Node::encode(5, false);
    gbwt_thread_1[4] = gbwt::Node::encode(6, false);
    
    gbwt_thread_2[0] = gbwt::Node::encode(6, true);
    gbwt_thread_2[1] = gbwt::Node::encode(4, true);
    gbwt_thread_2[2] = gbwt::Node::encode(2, true);
    gbwt_thread_2[3] = gbwt::Node::encode(1, true);

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
            "mapping_quality": 10,
            "score":1 
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
                            {"from_length": 1, "to_length": 1, "sequence": "T"},
                            {"from_length": 1, "to_length": 1}
                        ]
                    }
                ]
            },
            "mapping_quality": 20,
            "score":2 
        }
    )";

    vg::Alignment alignment_2;
    json2pb(alignment_2, alignment_2_str);

    AlignmentPathFinder<vg::Alignment> alignment_path_finder(graph, gbwt_index, 1000);

    auto alignment_paths = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2);
    REQUIRE(alignment_paths.size() == 3);

    REQUIRE(alignment_paths.at(0).search.node == gbwt::Node::encode(6, false));
    REQUIRE(alignment_paths.at(0).search.size() == 1);   
    REQUIRE(alignment_paths.at(0).end_offset == 6);   
    REQUIRE(alignment_paths.at(0).ids == vector<gbwt::size_type>({2}));
    REQUIRE(alignment_paths.at(0).seq_length == 17);
    REQUIRE(alignment_paths.at(0).mapqs == vector<int32_t>({10, 20}));
    REQUIRE(alignment_paths.at(0).scores == vector<int32_t>({1, 2}));

    REQUIRE(alignment_paths.at(1).search.node == alignment_paths.at(0).search.node);
    REQUIRE(alignment_paths.at(1).search.size() == alignment_paths.at(0).search.size());    
    REQUIRE(alignment_paths.at(1).end_offset == alignment_paths.at(0).end_offset);   
    REQUIRE(alignment_paths.at(1).ids == vector<gbwt::size_type>({0}));
    REQUIRE(alignment_paths.at(1).seq_length == 19);
    REQUIRE(alignment_paths.at(1).mapqs == alignment_paths.at(0).mapqs);
    REQUIRE(alignment_paths.at(1).scores == alignment_paths.at(0).scores);

    REQUIRE(alignment_paths.at(2).search.node == gbwt::Node::encode(1, true));
    REQUIRE(alignment_paths.at(2).search.size() == alignment_paths.at(0).search.size());    
    REQUIRE(alignment_paths.at(2).end_offset == 2);   
    REQUIRE(alignment_paths.at(2).ids == vector<gbwt::size_type>({1}));
    REQUIRE(alignment_paths.at(2).seq_length == alignment_paths.at(0).seq_length);
    REQUIRE(alignment_paths.at(2).mapqs == vector<int32_t>({20, 10}));
    REQUIRE(alignment_paths.at(2).scores == vector<int32_t>({2, 1}));

    SECTION("Extended paired-end read alignment finds alignment path") {
  
        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_from_length(3);
        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_to_length(3);

        auto new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(5);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(2);
        new_edit->set_to_length(2);

        auto alignment_paths_ext = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ext.size() == 1);

        REQUIRE(alignment_paths_ext.front() == alignment_paths.at(1));

        new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(4);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_paths_ext = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ext.size() == 1);

        REQUIRE(alignment_paths_ext.front() == alignment_paths.at(1));             
    }

    SECTION("Incorrect oriented paired-end read alignment does not finds alignment path") {

        alignment_2 = lazy_reverse_complement_alignment(alignment_2, node_seq_length_func);
        
        auto alignment_paths_rc = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2);
        REQUIRE(alignment_paths_rc.empty());
    }

    SECTION("Partial overlapping paired-end read alignment finds alignment paths") {

        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_from_length(3);
        alignment_2.mutable_path()->mutable_mapping(0)->mutable_edit(2)->set_to_length(3);

        auto new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(4);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        auto new_edit = new_mapping->add_edit();
        new_edit->set_from_length(5);
        new_edit->set_to_length(5);

        auto alignment_paths_ov = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 2);

        REQUIRE(alignment_paths_ov.front() == alignment_paths.at(0));
        REQUIRE(alignment_paths_ov.back() == alignment_paths.at(2));

        alignment_2.mutable_path()->mutable_mapping(1)->mutable_edit(0)->set_from_length(8);
        alignment_2.mutable_path()->mutable_mapping(1)->mutable_edit(0)->set_to_length(8);

        new_mapping = alignment_2.mutable_path()->add_mapping();
        new_mapping->mutable_position()->set_node_id(2);
        new_mapping->mutable_position()->set_offset(0);
        new_mapping->mutable_position()->set_is_reverse(true);

        new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);

        alignment_paths_ov = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2);
        REQUIRE(alignment_paths_ov.size() == 2);

        REQUIRE(alignment_paths_ov.front() == alignment_paths.at(0));
        REQUIRE(alignment_paths_ov.back() == alignment_paths.at(2));
    }

    SECTION("Perfect overlapping paired-end read alignment finds alignment paths") {

        auto alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_seq_length_func);

        auto alignment_paths_ov_1 = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_1_rc);
        REQUIRE(alignment_paths_ov_1.size() == 2);

        REQUIRE(alignment_paths_ov_1.front().search.node == gbwt::Node::encode(4, false));
        REQUIRE(alignment_paths_ov_1.front().search.size() == 2);   
        REQUIRE(alignment_paths_ov_1.front().end_offset == 5);   
        REQUIRE(alignment_paths_ov_1.front().ids == vector<gbwt::size_type>({0, 2}));
        REQUIRE(alignment_paths_ov_1.front().seq_length == 8);
        REQUIRE(alignment_paths_ov_1.front().mapqs == vector<int32_t>({10, 10}));
        REQUIRE(alignment_paths_ov_1.front().scores == vector<int32_t>({1, 1}));

        REQUIRE(alignment_paths_ov_1.back().search.node == gbwt::Node::encode(1, true));
        REQUIRE(alignment_paths_ov_1.back().search.size() == 1);   
        REQUIRE(alignment_paths_ov_1.back().end_offset == 2);   
        REQUIRE(alignment_paths_ov_1.back().ids == vector<gbwt::size_type>({1}));
        REQUIRE(alignment_paths_ov_1.back().seq_length == alignment_paths_ov_1.front().seq_length);
        REQUIRE(alignment_paths_ov_1.back().mapqs == alignment_paths_ov_1.front().mapqs);
        REQUIRE(alignment_paths_ov_1.back().scores == alignment_paths_ov_1.front().scores);

        auto alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_seq_length_func);

        auto alignment_paths_ov_2 = alignment_path_finder.findPairedAlignmentPathsIds(alignment_2, alignment_2_rc);
        REQUIRE(alignment_paths_ov_2.size() == 2);

        REQUIRE(alignment_paths_ov_2.front().search.node == gbwt::Node::encode(6, true));
        REQUIRE(alignment_paths_ov_2.front().search.size() == 1);   
        REQUIRE(alignment_paths_ov_2.front().end_offset == 5);   
        REQUIRE(alignment_paths_ov_2.front().ids == vector<gbwt::size_type>({1}));
        REQUIRE(alignment_paths_ov_2.front().seq_length == 4);
        REQUIRE(alignment_paths_ov_2.front().mapqs == vector<int32_t>({20, 20}));
        REQUIRE(alignment_paths_ov_2.front().scores == vector<int32_t>({2, 2}));

        REQUIRE(alignment_paths_ov_2.back().search.node == gbwt::Node::encode(6, false));
        REQUIRE(alignment_paths_ov_2.back().search.size() == 2);   
        REQUIRE(alignment_paths_ov_2.back().end_offset == 6);   
        REQUIRE(alignment_paths_ov_2.back().ids == vector<gbwt::size_type>({0, 2}));
        REQUIRE(alignment_paths_ov_2.back().seq_length == alignment_paths_ov_2.front().seq_length);
        REQUIRE(alignment_paths_ov_2.back().mapqs == alignment_paths_ov_2.front().mapqs);
        REQUIRE(alignment_paths_ov_2.back().scores == alignment_paths_ov_2.front().scores);
    }

    SECTION("Maximum paired sequence length filters alignment pairs") {

        alignment_path_finder.setMaxPairSeqLength(19);

        auto alignment_paths_len = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_len.size() == 3);
        
        REQUIRE(alignment_paths_len == alignment_paths);

        alignment_path_finder.setMaxPairSeqLength(18);

        alignment_paths_len = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_len.size() == 2);
        
        REQUIRE(alignment_paths_len.front() == alignment_paths.at(0));
        REQUIRE(alignment_paths_len.back() == alignment_paths.at(2));

        alignment_path_finder.setMaxPairSeqLength(10);

        alignment_paths_len = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2);        
        REQUIRE(alignment_paths_len.empty());
    }
}

TEST_CASE("Alignment paths can be found from single-end multipath alignment") {

    const string graph_str = R"(
        {
            "node": [
                {"id": 1, "sequence": "A"},
                {"id": 2, "sequence": "C"},
                {"id": 3, "sequence": "TTT"},
                {"id": 4, "sequence": "TT"},
                {"id": 5, "sequence": "GGG"},
                {"id": 6, "sequence": "AGG"},
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

    vector<int32_t> node_seq_lengths = {0, 1, 1, 3, 2, 3, 3};
    function<size_t(const int64_t)> node_seq_length_func = [&](const int64_t node_id) { return node_seq_lengths.at(node_id); };

    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(6, true)));

    gbwt::vector_type gbwt_thread_1(4);
    gbwt::vector_type gbwt_thread_2(4);
    gbwt::vector_type gbwt_thread_3(4);
   
    gbwt_thread_1[0] = gbwt::Node::encode(1, false);
    gbwt_thread_1[1] = gbwt::Node::encode(3, false);
    gbwt_thread_1[2] = gbwt::Node::encode(4, false);
    gbwt_thread_1[3] = gbwt::Node::encode(5, false);

    gbwt_thread_2[0] = gbwt::Node::encode(1, false);
    gbwt_thread_2[1] = gbwt::Node::encode(3, false);
    gbwt_thread_2[2] = gbwt::Node::encode(4, false);
    gbwt_thread_2[3] = gbwt::Node::encode(6, false);

    gbwt_thread_3[0] = gbwt::Node::encode(5, true);
    gbwt_thread_3[1] = gbwt::Node::encode(4, true);
    gbwt_thread_3[2] = gbwt::Node::encode(3, true);
    gbwt_thread_3[3] = gbwt::Node::encode(2, true);

    gbwt_builder.insert(gbwt_thread_1, false);
    gbwt_builder.insert(gbwt_thread_2, false);
    gbwt_builder.insert(gbwt_thread_3, false);

    gbwt_builder.finish();

    std::stringstream gbwt_stream;
    gbwt_builder.index.serialize(gbwt_stream);

    gbwt::GBWT gbwt_index;
    gbwt_index.load(gbwt_stream);

    const string alignment_1_str = R"(
        {
            "start":[0,1],
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
                    "next":[2],
                    "score":4
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
                    "next":[2],
                    "score":1
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
                    "next":[3,4],
                    "score":6
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
                    "score":4
                },
                {
                    "path": {
                        "mapping": [
                            {
                                "position": {"node_id": 6},
                                "edit": [
                                    {"from_length": 1, "to_length": 1, "sequence": "G"},
                                    {"from_length": 1, "to_length": 1}
                                ]
                            }
                        ]
                    },
                    "score":2
                }
            ],
            "mapping_quality": 10
        }
    )";

    vg::MultipathAlignment alignment_1;
    json2pb(alignment_1, alignment_1_str);

    AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder(graph, gbwt_index, 1000);
    
    auto alignment_paths = alignment_path_finder.findAlignmentPathsIds(alignment_1);
    REQUIRE(alignment_paths.size() == 2);

    REQUIRE(alignment_paths.front().search.node == gbwt::Node::encode(5, false));
    REQUIRE(alignment_paths.front().search.size() == 1);   
    REQUIRE(alignment_paths.front().end_offset == 2);   
    REQUIRE(alignment_paths.front().ids == vector<gbwt::size_type>({0}));
    REQUIRE(alignment_paths.front().seq_length == 8);
    REQUIRE(alignment_paths.front().mapqs == vector<int32_t>({10}));
    REQUIRE(alignment_paths.front().scores == vector<int32_t>({14}));

    REQUIRE(alignment_paths.back().search.node == gbwt::Node::encode(6, false));
    REQUIRE(alignment_paths.back().search.size() == alignment_paths.front().search.size());   
    REQUIRE(alignment_paths.back().end_offset == alignment_paths.front().end_offset);   
    REQUIRE(alignment_paths.back().ids == vector<gbwt::size_type>({1}));
    REQUIRE(alignment_paths.back().seq_length == alignment_paths.front().seq_length);
    REQUIRE(alignment_paths.back().mapqs == alignment_paths.front().mapqs);
    REQUIRE(alignment_paths.back().scores == vector<int32_t>({12}));

    SECTION("Reverse-complement single-end multipath read alignment finds alignment path") {

        alignment_1 = lazy_reverse_complement_alignment(alignment_1, node_seq_length_func);
        
        auto alignment_paths_rc = alignment_path_finder.findAlignmentPathsIds(alignment_1);
        REQUIRE(alignment_paths_rc.size() == 1);

        REQUIRE(alignment_paths_rc.front().search.node == gbwt::Node::encode(2, true));
        REQUIRE(alignment_paths_rc.front().search.size() == alignment_paths.front().search.size());
        REQUIRE(alignment_paths_rc.front().end_offset == 1);   
        REQUIRE(alignment_paths_rc.front().ids == vector<gbwt::size_type>({2}));
        REQUIRE(alignment_paths_rc.front().seq_length == alignment_paths.front().seq_length);
        REQUIRE(alignment_paths_rc.front().mapqs == alignment_paths.front().mapqs);
        REQUIRE(alignment_paths_rc.front().scores == vector<int32_t>({11}));
    }

    SECTION("Soft-clipped single-end multipath read alignment finds alignment path") {

        alignment_1.mutable_subpath(3)->mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_from_length(1);
        alignment_1.mutable_subpath(3)->mutable_path()->mutable_mapping(0)->mutable_edit(0)->set_to_length(1);

        auto new_edit = alignment_1.mutable_subpath(3)->mutable_path()->mutable_mapping(0)->add_edit();
        new_edit->set_from_length(0);
        new_edit->set_to_length(1);
        new_edit->set_sequence("C");

        auto alignment_paths_sc = alignment_path_finder.findAlignmentPathsIds(alignment_1);
        REQUIRE(alignment_paths_sc.size() == 2);

        REQUIRE(alignment_paths_sc.front().search.node == alignment_paths.front().search.node);
        REQUIRE(alignment_paths_sc.front().search.size() == alignment_paths.front().search.size());
        REQUIRE(alignment_paths_sc.front().end_offset == 1);   
        REQUIRE(alignment_paths_sc.front().ids == alignment_paths.front().ids);
        REQUIRE(alignment_paths_sc.front().seq_length == alignment_paths.front().seq_length);
        REQUIRE(alignment_paths_sc.front().mapqs == alignment_paths.front().mapqs);
        REQUIRE(alignment_paths_sc.front().scores == alignment_paths.front().scores);
    
        REQUIRE(alignment_paths_sc.back() == alignment_paths.back());
    }

    SECTION("Alternative single-end multipath read alignment does not find alignment path") {

        alignment_1.clear_start();
        alignment_1.add_start(1);
        
        auto alignment_paths_alt = alignment_path_finder.findAlignmentPathsIds(alignment_1);
        REQUIRE(alignment_paths_alt.empty());
    }
}



