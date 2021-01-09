
#include "alignment_path_finder.hpp"

#include <assert.h>

#include "utils.hpp"

//#define debug

static const uint32_t max_internal_offset = 8;
static const uint32_t internal_start_end_penalty = 1;

static const uint32_t max_score_diff = scorePrecision(double_precision);


template<class AlignmentType>
AlignmentPathFinder<AlignmentType>::AlignmentPathFinder(const PathsIndex & paths_index_in, const string library_type_in, const uint32_t max_pair_frag_length_in, const uint32_t min_mapq_filter_in, const double min_best_score_filter_in, const double max_softclip_filter_in) : paths_index(paths_index_in), library_type(library_type_in), max_pair_frag_length(max_pair_frag_length_in), min_mapq_filter(min_mapq_filter_in), min_best_score_filter(min_best_score_filter_in), max_softclip_filter(max_softclip_filter_in) {}
        
template<class AlignmentType>
bool AlignmentPathFinder<AlignmentType>::alignmentHasPath(const vg::Alignment & alignment) const {

    return alignment.has_path();
}

template<class AlignmentType>
bool AlignmentPathFinder<AlignmentType>::alignmentHasPath(const vg::MultipathAlignment & alignment) const {

    return (alignment.subpath_size() > 0);
}

template<class AlignmentType>
bool AlignmentPathFinder<AlignmentType>::alignmentStartInGraph(const AlignmentType & alignment) const {

    auto alignment_start_nodes_index = getAlignmentStartNodesIndex(alignment, 0);

    for (auto & start_node: alignment_start_nodes_index) {

        if (!paths_index.hasNodeId(gbwt::Node::id(start_node.first))) {

            return false;
        } 
    } 

    return true;
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::findAlignmentPaths(const AlignmentType & alignment) const {

#ifdef debug

    cerr << endl;
    cerr << pb2json(alignment) << endl;

#endif

    if (!alignmentHasPath(alignment)) {

        return vector<AlignmentPath>();
    }

    if (!alignmentStartInGraph(alignment)) {

        return vector<AlignmentPath>();
    }

    vector<AlignmentSearchPath> align_search_paths;

    function<size_t(const uint32_t)> node_length_func = [&](const uint32_t node_id) { return paths_index.nodeLength(node_id); };

    if (library_type == "fr") {

        align_search_paths = extendAlignmentPath(AlignmentSearchPath(), alignment);

    } else if (library_type == "rf") {

        AlignmentType alignment_rc = lazy_reverse_complement_alignment(alignment, node_length_func);
        align_search_paths = extendAlignmentPath(AlignmentSearchPath(), alignment_rc);

    } else {

        assert(library_type == "unstranded");
        align_search_paths = extendAlignmentPath(AlignmentSearchPath(), alignment);

        if (!paths_index.index().bidirectional()) {

            AlignmentType alignment_rc = lazy_reverse_complement_alignment(alignment, node_length_func);
            auto align_search_paths_rc = extendAlignmentPath(AlignmentSearchPath(), alignment_rc);

            align_search_paths.reserve(align_search_paths.size() + align_search_paths_rc.size());
            align_search_paths.insert(align_search_paths.end(), align_search_paths_rc.begin(), align_search_paths_rc.end());
        }  
    }

    if (filterAlignmentSearchPaths(align_search_paths)) {

        return vector<AlignmentPath>();
    }

    auto align_paths = AlignmentPath::alignmentSearchPathsToAlignmentPaths(align_search_paths, max_score_diff, isAlignmentDisconnected(alignment));

#ifdef debug

    cerr << endl;
    cerr << align_search_paths << endl;
    cerr << align_paths << endl;
    cerr << endl;

#endif

    return align_paths;
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::Alignment & alignment) const {

    return extendAlignmentPath(align_search_path, alignment, 0, 0);
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::Alignment & alignment, const pair<uint32_t, uint32_t> start_node_path_idx) const {

    assert(alignment.mapping_quality() >= 0);

    vector<AlignmentSearchPath> extended_align_search_path(1, align_search_path);
    
    extended_align_search_path.front().read_stats.emplace_back(ReadAlignmentStats());
    extended_align_search_path.front().read_stats.back().mapq = alignment.mapping_quality();
    extended_align_search_path.front().read_stats.back().score = alignment.score();
    
    extendAlignmentPath(&extended_align_search_path, alignment.path(), start_node_path_idx.second);
    return extended_align_search_path;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPath(vector<AlignmentSearchPath> * align_search_paths, const vg::Path & path, const uint32_t path_start_idx) const {

    assert(align_search_paths->size() == 1);
    auto main_align_search_path = &(align_search_paths->front())

    assert(!main_align_search_path->read_stats.empty());
    auto main_read_stats = &(main_align_search_path->read_stats.back());

    assert(main_align_search_path->path_idx <= main_align_search_path->path.size());

    auto mapping_it = path.mapping().cbegin();
    assert(mapping_it != path.mapping().cend());

    const bool is_first_node = (main_read_stats->left_softclip_length == -1);

    if (!main_align_search_path->path.empty() && main_align_search_path->path_idx == 0) {

        assert(is_first_node);
 
        if (mapping_it->position().offset() < main_align_search_path->start_offset) {
 
            main_align_search_path->search_state = gbwt::SearchState();
            assert(main_align_search_path->search_state.empty());

            return;      
        } 
    }

    main_read_stats->updateSoftClippingLengths(path);

    assert(main_read_stats->left_softclip_length >= 0);
    assert(main_read_stats->right_softclip_length >= 0);

    assert(path_start_idx < path.mapping_size());

    for (uint32_t i = 0; i < path_start_idx; ++i) {

        main_align_search_path->path_idx++;
        main_align_search_path->path_offset = mapping_it->position().offset() + mapping_from_length(*mapping_it);



        ++mapping_it;
    }

    auto end_mapping_it = path.mapping().cend();
    --end_mapping_it;

    while (main_align_search_path->path_idx < main_align_search_path->path.size() && mapping_it != path.mapping().cend()) {

        auto cur_node = mapping_to_gbwt(*mapping_it);

        auto mapping_from_length = mapping_from_length(*mapping_it);
        auto mapping_to_length = mapping_to_length(*mapping_it);

        bool is_multi_visit = false;

        if (!is_first_node) {

            assert(main_align_search_path->path_idx > 0);

            if (main_align_search_path->path.at(main_align_search_path->path_idx - 1) == cur_node && main_align_search_path->path_offset == mapping_it->position().offset()) {

                is_multi_visit = true;
            }
        }

        if (!is_multi_visit) {

            ++main_align_search_path->path_idx;
        }     

        main_align_search_path->path_offset = mapping_it->position().offset() + mapping_from_length;
        main_read_stats->length += mapping_to_length;

        if (path_start_idx < mapping_it - path.mapping().cbegin()) {

            assert(main_align_search_path->path_idx < main_align_search_path->path.size());
            main_read_stats->updateInternalStartOffset(main_read_stats->length, true);

            if (main_read_stats->internal_start_offset > max_internal_offset) {

                main_align_search_path->search_state = gbwt::SearchState();
                assert(main_align_search_path->search_state.empty());

                return;
            }

        } else if (main_align_search_path->path.at(main_align_search_path->path_idx - 1) == cur_node || is_multi_visit) {

            if (main_align_search_path->path_idx == main_align_search_path->path.size()) {

                auto new_end_offset = mapping_it->position().offset() + mapping_from_length;

                if (new_end_offset < main_align_search_path->end_offset) {

                    main_align_search_path->search_state = gbwt::SearchState();
                    assert(main_align_search_path->search_state.empty());

                    return;
                }

                main_align_search_path->insert_length += (static_cast<int32_t>(mapping_it->position().offset()) - static_cast<int32_t>(main_align_search_path->end_offset));
                main_align_search_path->end_offset = new_end_offset;
            
            } else {

                main_align_search_path->insert_length -= mapping_from_length(*mapping_it);
            }

        } else {

            main_align_search_path->search_state = gbwt::SearchState();
            assert(main_align_search_path->search_state.empty());  
    
            return;  
        } 

        if (main_align_search_path->path_idx < main_align_search_path->path.size()) {

            for (auto & align_search_path: *align_search_paths) {

                if (align_search_path.read_stats.back().internal_end_offset > 0) {

                    align_search_path.read_stats.back().updateInternalEndOffset(mapping_to_length, mapping_it == end_mapping_it);

                    if (align_search_path.read_stats.back().internal_end_offset > max_internal_offset) {
       
                        align_search_path.search_state = gbwt::SearchState();
                        assert(align_search_path.search_state.empty());
                    }
                }

                align_search_path.read_stats.back().length += mapping_to_length;
            }

            auto internal_end_read_stats = *main_read_stats;
            internal_end_read_stats.updateInternalEndOffset(mapping_to_length, mapping_it == end_mapping_it);

            if (internal_end_read_stats.internal_end_offset <= max_internal_offset) {

                align_search_paths->emplace_back(align_search_paths->front());
                align_search_paths->back().read_stats.back() = internal_end_read_stats;
            }        
        }

        ++mapping_it;
    }

    while (mapping_it != path.mapping().cend()) {

        auto mapping_to_length = mapping_to_length(*mapping_it);

        if (main_read_stats->internal_end_offset > 0) {

            assert(align_search_paths->size() == 1);
            main_read_stats->updateInternalEndOffset(mapping_to_length, mapping_it == end_mapping_it);

            if (main_read_stats->internal_end_offset > max_internal_offset) {

                main_align_search_path->search_state = gbwt::SearchState();
                assert(main_align_search_path->search_state.empty());

                return;          
            }

        } else {

            auto internal_start_read_stats = *main_read_stats;
            internal_start_read_stats.internal_start_offset = std::numeric_limits<uint32_t>::max();

            if (!main_align_search_path->path.empty() && main_align_search_path->read_stats.size() == 1 && main_read_stats->internal_start_offset == 0) {

                internal_start_read_stats.internal_start_offset = 0;
                internal_start_read_stats.updateInternalStartOffset(main_read_stats->length, true);
            }

            for (auto & align_search_path: *align_search_paths) {

                if (align_search_path.read_stats.back().internal_end_offset == 0) {

                    extendAlignmentPath(&align_search_path, *mapping_it);

                } else {

                    align_search_path.read_stats.back().updateInternalEndOffset(mapping_to_length, mapping_it == end_mapping_it);

                    if (align_search_path.read_stats.back().internal_end_offset > max_internal_offset) {
       
                        align_search_path.search_state = gbwt::SearchState();
                        assert(align_search_path.search_state.empty());
                    }
                }

                align_search_path.read_stats.back().length += mapping_to_length;
            }

            auto internal_end_read_stats = *main_read_stats;
            internal_end_read_stats.updateInternalEndOffset(mapping_to_length, mapping_it == end_mapping_it);

            if (internal_end_read_stats.internal_end_offset <= max_internal_offset) {

                align_search_paths->emplace_back(align_search_paths->front());
                align_search_paths->back().read_stats.back() = internal_end_read_stats;
            }

            if (internal_start_read_stats.internal_start_offset <= max_internal_offset) {

                assert(internal_start_read_stats.internal_start_offset > 0);

                assert(main_align_search_path->insert_length == 0);
                assert(main_align_search_path->read_stats.size() == 1);

                AlignmentSearchPath internal_start_align_search_path;

                internal_start_align_search_path.start_offset = main_align_search_path->start_offset;
                internal_start_align_search_path.read_stats = internal_start_read_stats;

                extendAlignmentPath(&align_search_path, *mapping_it);
                assert(internal_start_align_search_path.search_state.node == main_align_search_path->search_state.node);

                if (!internal_start_align_search_path.search_state.empty() && internal_start_align_search_path.search_state.range != main_align_search_path->search_state.range) {

                    align_search_paths->emplace_back(move(internal_start_align_search_path));
                }
            }
        }

        ++mapping_it;
    }
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPath(AlignmentSearchPath * align_search_path, const vg::Mapping & mapping) const {

    assert(align_search_path->path_offset == align_search_path->end_offset);
    auto cur_node = mapping_to_gbwt(mapping);

    if (align_search_path->path.empty()) {

        assert(align_search_path->search_state.node == gbwt::ENDMARKER);
        assert(align_search_path->insert_length == 0);

        align_search_path->start_offset = mapping_it->position().offset();
        align_search_path->search_state = paths_index.index().find(cur_node);
  
        align_search_path->path.emplace_back(cur_node);
        ++align_search_path->path_idx;

    } else {

        bool is_cycle_visit = false;

        if (align_search_path->path.back() == cur_node && mapping_it->position().offset() != align_search_path->path_offset) {

            assert(mapping_it->position().offset() == 0);
            is_cycle_visit = true;      
        }

        if (align_search_path->path.back() != cur_node || is_cycle_visit) {

            align_search_path->search_state = paths_index.index().extend(align_search_path->search_state, cur_node);

            align_search_path->path.emplace_back(cur_node);
            ++align_search_path->path_idx;
        }
    }

    align_search_path->path_offset = mapping_it->position().offset() + mapping_from_length(mapping);
    align_search_path->end_offset = align_search_path->path_offset;
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::MultipathAlignment & alignment) const {

    vector<AlignmentSearchPath> extended_align_search_paths;

    for (auto & start_idx: alignment.start()) {

        auto cur_extended_align_search_paths = extendAlignmentPath(align_search_path, alignment, start_idx, 0);
        extended_align_search_paths.insert(extended_align_search_paths.end(), make_move_iterator(cur_extended_align_search_paths.begin()), make_move_iterator(cur_extended_align_search_paths.end()));
    }

    return extended_align_search_paths;
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::MultipathAlignment & alignment, const pair<uint32_t, uint32_t> start_node_path_idx) const {

    assert(alignment.mapping_quality() >= 0);

    vector<AlignmentSearchPath> extended_align_search_path(1, align_search_path);
    
    extended_align_search_path.front().read_stats.emplace_back(ReadAlignmentStats());
    extended_align_search_path.front().read_stats.back().mapq = alignment.mapping_quality();

    extendAlignmentPaths(&extended_align_search_path, alignment.subpath(), start_node_path_idx);
            
    return extended_align_search_path;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPaths(vector<AlignmentSearchPath> * align_search_paths, const google::protobuf::RepeatedPtrField<vg::Subpath> & subpaths, const pair<uint32_t, uint32_t> start_node_path_idx) const {

    std::queue<pair<AlignmentSearchPath, pair<uint32_t, uint32_t> > > align_search_paths_queue;

    for (auto & align_search_path: *align_search_paths) {

        align_search_paths_queue.push(make_pair(align_search_path, start_node_path_idx)));
    }

    align_search_paths->clear();

    // Perform depth-first alignment path extension.
    while (!align_search_paths_queue.empty()) {

        auto & cur_align_search_path = align_search_paths_queue.front();
        const vg::Subpath & subpath = subpaths.Get(cur_align_search_path.second.first);

        vector<AlignmentSearchPath> extended_align_search_path(1, cur_align_search_path.first);
        extended_align_search_path.front().read_stats.back().score += subpath.score();

        extendAlignmentPath(&extended_align_search_path, subpath.path(), cur_align_search_path.second.second);

        for (auto & align_search_path: extended_align_search_path) {

            if (align_search_path.path.empty() || !align_search_path.search_state.empty()) {

                if (subpath.next_size() > 0) {

                    for (auto & next_subpath_idx: subpath.next()) {

                        align_search_paths_queue.push(make_pair(move(align_search_path), make_pair(next_subpath_idx, 0)));
                    }

                } else {

                    align_search_paths->emplace_back(move(align_search_path));
                }
            }
        }

        align_search_paths_queue.pop();
    }
}

// Debug start

char quality_short_to_char(short i) {
    return static_cast<char>(i + 33);
}

string string_quality_short_to_char(const string& quality) {
    string buffer; buffer.resize(quality.size());
    for (int i = 0; i < quality.size(); ++i) {
        buffer[i] = quality_short_to_char(quality[i]);
    }
    return buffer;
}

// Debug end

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::findPairedAlignmentPaths(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const {

#ifdef debug

    cerr << endl;
    findAlignmentPaths(alignment_1);
    findAlignmentPaths(alignment_2);

#endif

    if (!alignmentHasPath(alignment_1) || !alignmentHasPath(alignment_2)) {

        return vector<AlignmentPath>();
    }

    if (!alignmentStartInGraph(alignment_1) || !alignmentStartInGraph(alignment_2)) {

        return vector<AlignmentPath>();
    }

    vector<AlignmentSearchPath> paired_align_search_paths;

    function<size_t(const uint32_t)> node_length_func = [&](const uint32_t node_id) { return paths_index.nodeLength(node_id); };
    AlignmentType alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_length_func);

    if (library_type == "fr") {

        pairAlignmentPaths(&paired_align_search_paths, alignment_1, alignment_2_rc);

    } else if (library_type == "rf") {

        AlignmentType alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_length_func);
        pairAlignmentPaths(&paired_align_search_paths, alignment_2, alignment_1_rc);

    } else {

        assert(library_type == "unstranded");

        AlignmentType alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_length_func);
        pairAlignmentPaths(&paired_align_search_paths, alignment_1, alignment_2_rc);

        if (!paths_index.index().bidirectional()) {

            AlignmentType alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_length_func);
            pairAlignmentPaths(&paired_align_search_paths, alignment_2, alignment_1_rc);
        }
    }

    if (filterAlignmentSearchPaths(paired_align_search_paths)) {

        return vector<AlignmentPath>();
    }

    auto paired_align_paths = AlignmentPath::alignmentSearchPathsToAlignmentPaths(paired_align_search_paths, max_score_diff, isAlignmentDisconnected(alignment_1) || isAlignmentDisconnected(alignment_2));

    // Debug start

    // string debug_paths = "";
    // int32_t debug_idx = -1;

    // string debug_paths2 = "";
    // int32_t debug_idx2 = -1;

    // for (size_t i = 0; i < paired_align_search_paths.size(); ++i) {

    //     if (paired_align_search_paths.at(i).complete()) {

    //         for (auto & path_id: paths_index.locatePathIds(paired_align_search_paths.at(i).search_state)) {

    //             auto path_name = paths_index.pathName(path_id);

    //             if (
    //                 path_name == "ENST00000580018.3_15" || 
    //                 path_name == "ENST00000374259.7" || 
    //                 path_name == "ENST00000325307.11" || 
    //                 path_name == "ENST00000216252.3_19" || 
    //                 path_name == "ENST00000271638.2"
    //             ) {   

    //                 debug_paths = path_name; 
    //                 debug_idx = i;         
                
    //             } else if (
    //                 path_name == "ENST00000580018.3_16" || 
    //                 path_name == "ENST00000374259.7_19" || 
    //                 path_name == "ENST00000325307.11_24" || 
    //                 path_name == "ENST00000216252.3_24" || 
    //                 path_name == "ENST00000271638.2_8"
    //                 ) {  

    //                 debug_paths2 = path_name; 
    //                 debug_idx2 = i;         
    //             }                
    //         }
    //     }
    // }

    // if (debug_idx != debug_idx2) {

    //     #pragma omp critical
    //     {
    //         cerr << "\n\n###" << endl;
    //         cerr << debug_paths << endl;
    //         cerr << debug_idx << endl;

    //         if (debug_idx >= 0) {

    //             cerr << paired_align_search_paths.at(debug_idx) << endl;
    //         }

    //         cerr << debug_paths2 << endl;
    //         cerr << debug_idx2 << endl;

    //         if (debug_idx2 >= 0) {

    //             cerr << paired_align_search_paths.at(debug_idx2) << endl;
    //         }
            
    //         cerr << endl;
    //         cerr << paired_align_search_paths << endl;
    //         cerr << endl;
    //         cerr << pb2json(alignment_1) << endl;
    //         cerr << string_quality_short_to_char(alignment_1.quality()) << endl;
    //         cerr << endl;
    //         cerr << pb2json(alignment_2) << endl;
    //         cerr << string_quality_short_to_char(alignment_2.quality()) << endl;
    //     }
    // }

    // Debug end

#ifdef debug

    cerr << endl;
    cerr << paired_align_search_paths << endl;
    cerr << paired_align_paths << endl;
    cerr << endl;

#endif

    return paired_align_paths;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::pairAlignmentPaths(vector<AlignmentSearchPath> * paired_align_search_paths, const AlignmentType & start_alignment, const AlignmentType & end_alignment) const {

    auto start_align_search_paths = extendAlignmentPath(AlignmentSearchPath(), start_alignment);
    auto end_alignment_start_nodes_index = getAlignmentStartNodesIndex(end_alignment, max_internal_offset);

    std::queue<AlignmentSearchPath> paired_align_search_path_queue;

    for (auto & align_search_path: start_align_search_paths) {

        assert(!align_search_path.search_state.empty());
        assert(!align_search_path.path.empty());

        auto node_length = paths_index.nodeLength(gbwt::Node::id(align_search_path.search_state.node));
        
        assert(align_search_path.path_offset <= node_length);
        assert(align_search_path.path_offset == align_search_path.end_offset);

        paired_align_search_path_queue.push(align_search_path);

        paired_align_search_path_queue.back().insert_length += (node_length - align_search_path.path_offset);
        paired_align_search_path_queue.back().path_offset = node_length;
        paired_align_search_path_queue.back().end_offset = node_length;

        for (auto & start_nodes: end_alignment_start_nodes_index) {

            auto path_it = find(align_search_path.path.begin(), align_search_path.path.end(), start_nodes.first); 
            
            auto align_search_path_end = align_search_path.path.end();
            --align_search_path_end;

            while (path_it != align_search_path.path.end()) {

                if (path_it == align_search_path_end) {

                    break;
                }

                align_search_path.path_idx = path_it - align_search_path.path.begin();
                auto complete_paired_align_search_paths = extendAlignmentPath(align_search_path, end_alignment, start_nodes.second);

                for (auto & complete_align_search_path: complete_paired_align_search_paths) {

                    if (!complete_align_search_path.search_state.empty() && complete_align_search_path.fragmentLength() <= max_pair_frag_length) {

                        paired_align_search_paths->emplace_back(complete_align_search_path);                         
                    }
                }

                ++path_it;
                path_it = find(path_it, align_search_path.path.end(), start_nodes.first); 
            }
        }
    }

    // Perform depth-first path extension.
    while (!paired_align_search_path_queue.empty()) {

        AlignmentSearchPath * cur_paired_align_search_path = &(paired_align_search_path_queue.front());
        assert(cur_paired_align_search_path->search_state.node != gbwt::ENDMARKER);

        auto end_alignment_start_nodes_index_itp = end_alignment_start_nodes_index.equal_range(cur_paired_align_search_path->search_state.node);

        if (end_alignment_start_nodes_index_itp.first != end_alignment_start_nodes_index_itp.second) {

            while (end_alignment_start_nodes_index_itp.first != end_alignment_start_nodes_index_itp.second) {

                AlignmentSearchPath cur_paired_align_search_path_end = *cur_paired_align_search_path;

                assert(cur_paired_align_search_path_end.path_idx == cur_paired_align_search_path_end.path.size());
                --cur_paired_align_search_path_end.path_idx;

                assert(cur_paired_align_search_path_end.path_offset == cur_paired_align_search_path_end.end_offset);
                cur_paired_align_search_path_end.insert_length -= cur_paired_align_search_path_end.path_offset;
                
                cur_paired_align_search_path_end.path_offset = 0;
                cur_paired_align_search_path_end.end_offset = 0;

                auto complete_paired_align_search_paths = extendAlignmentPath(cur_paired_align_search_path_end, end_alignment, end_alignment_start_nodes_index_itp.first->second);

                for (auto & complete_align_search_path: complete_paired_align_search_paths) {

                    if (!complete_align_search_path.search_state.empty() && complete_align_search_path.fragmentLength() <= max_pair_frag_length) {

                        paired_align_search_paths->emplace_back(complete_align_search_path);                         
                    }
                }

                ++end_alignment_start_nodes_index_itp.first;
            }
        }
           
        if (cur_paired_align_search_path->fragmentLength() + end_alignment.sequence().size() > max_pair_frag_length) {

            paired_align_search_path_queue.pop();
            continue;
        }

        auto out_edges = paths_index.index().edges(cur_paired_align_search_path->search_state.node);

        // End current extension if no outgoing edges exist.
        if (out_edges.empty()) {

            paired_align_search_path_queue.pop();
            continue;
        }

        auto out_edges_it = out_edges.begin(); 
        assert(out_edges_it != out_edges.end());
        
        ++out_edges_it;

        while (out_edges_it != out_edges.end()) {

            if (out_edges_it->first != gbwt::ENDMARKER) {

                auto extended_path = paths_index.index().extend(cur_paired_align_search_path->search_state, out_edges_it->first);

                // Add new extension to queue if not empty (path found).
                if (!extended_path.empty()) { 

                    paired_align_search_path_queue.push(*cur_paired_align_search_path);
                    paired_align_search_path_queue.back().path.emplace_back(extended_path.node);
                    ++paired_align_search_path_queue.back().path_idx;
                    paired_align_search_path_queue.back().path_offset = paths_index.nodeLength(gbwt::Node::id(extended_path.node));
                    paired_align_search_path_queue.back().search_state = extended_path;
                    paired_align_search_path_queue.back().end_offset = paired_align_search_path_queue.back().path_offset;
                    paired_align_search_path_queue.back().insert_length += paired_align_search_path_queue.back().path_offset;
                }
            }

            ++out_edges_it;
        }

        if (out_edges.begin()->first != gbwt::ENDMARKER) {
            
            cur_paired_align_search_path->search_state = paths_index.index().extend(cur_paired_align_search_path->search_state, out_edges.begin()->first);

            // End current extension if empty (no haplotypes found). 
            if (cur_paired_align_search_path->search_state.empty()) { 

                paired_align_search_path_queue.pop(); 

            } else {

                cur_paired_align_search_path->path.emplace_back(cur_paired_align_search_path->search_state.node);
                ++cur_paired_align_search_path->path_idx;
                cur_paired_align_search_path->path_offset = paths_index.nodeLength(gbwt::Node::id(cur_paired_align_search_path->search_state.node));
                cur_paired_align_search_path->end_offset = cur_paired_align_search_path->path_offset;
                cur_paired_align_search_path->insert_length += cur_paired_align_search_path->path_offset;
            }
    
        } else {

            paired_align_search_path_queue.pop();
        }
    }
}

template<class AlignmentType>
multimap<gbwt::node_type, pair<uint32_t, uint32_t> > AlignmentPathFinder<AlignmentType>::getAlignmentStartNodesIndex(const vg::Alignment & alignment, const uint32_t max_internal_offset) const {

    multimap<gbwt::node_type, uint32_t> alignment_start_nodes_index;

    assert(alignment.alignment.path().mapping_size() > 0);

    auto mapping_it = alignment.path().mapping().cbegin();
    assert(mapping_it != alignment.path().mapping().cend());

    uint32_t mapping_idx = 0;

    alignment_start_nodes_index.emplace(mapping_to_gbwt(*mapping_it), mapping_idx);
    ++mapping_idx;

    uint32_t internal_offset = 

    assert(mapping_it->edit_size() > 0);
    const vg::Edit & first_edit = mapping_it->edit(0);

    if (first_edit.from_length() == 0) {

         = first_edit.to_length();   
    } 

    while (mapping_it != path.mapping().cend()) {

        ++mapping_it;

        alignment_start_nodes_index.emplace(mapping_to_gbwt(*mapping_it), mapping_idx);        
        ++mapping_idx;
    }

    return alignment_start_nodes_index;
}

template<class AlignmentType>
multimap<gbwt::node_type, pair<uint32_t, uint32_t> > AlignmentPathFinder<AlignmentType>::getAlignmentStartNodesIndex(const vg::MultipathAlignment & alignment, const uint32_t max_internal_offset) const {

    multimap<gbwt::node_type, uint32_t> alignment_start_nodes_index;

    for (auto & start_idx: alignment.start()) {

        assert(alignment.subpath(start_idx).path().mapping_size() > 0);
        alignment_start_nodes_index.emplace(mapping_to_gbwt(alignment.subpath(start_idx).path().mapping(0)), start_idx);
    }

    return alignment_start_nodes_index;
}

template<class AlignmentType>
bool AlignmentPathFinder<AlignmentType>::isAlignmentDisconnected(const vg::Alignment & alignment) const {

    return false;
}

template<class AlignmentType>
bool AlignmentPathFinder<AlignmentType>::isAlignmentDisconnected(const vg::MultipathAlignment & alignment) const {

    bool is_connected = false;

    if (alignment.has_annotation()) {

        auto annotation_it = alignment.annotation().fields().find("disconnected");

        if (annotation_it != alignment.annotation().fields().end()) {

            assert(annotation_it->second.bool_value());
            is_connected = true;
        }
    }

    return is_connected;
}

template<class AlignmentType>
bool AlignmentPathFinder<AlignmentType>::filterAlignmentSearchPaths(const vector<AlignmentSearchPath> & align_search_paths) const {

    double max_min_best_score_frac = 0;
    double min_max_softclip_frac = 1;

    for (auto & align_search_path: align_search_paths) {

        if (align_search_path.isComplete()) {

            if (align_search_path.minMappingQuality() < min_mapq_filter) {

                return true;
            }

            max_min_best_score_frac = max(max_min_best_score_frac, align_search_path.minBestScoreFraction());
            min_max_softclip_frac = min(min_max_softclip_frac, align_search_path.maxSoftclipFraction());
        }
    }

    if (max_min_best_score_frac < min_best_score_filter || min_max_softclip_frac > max_softclip_filter) {

        return true;
    
    } else {

        return false;
    }
}

template class AlignmentPathFinder<vg::Alignment>;
template class AlignmentPathFinder<vg::MultipathAlignment>;

