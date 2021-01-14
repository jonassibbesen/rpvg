
#include "alignment_path_finder.hpp"

#include <assert.h>

#include "utils.hpp"

//#define debug

static const uint32_t max_score_diff = scorePrecision(double_precision);


template<class AlignmentType>
AlignmentPathFinder<AlignmentType>::AlignmentPathFinder(const PathsIndex & paths_index_in, const string library_type_in, const uint32_t max_pair_frag_length_in, const uint32_t max_internal_offset_in, const uint32_t min_mapq_filter_in, const double min_best_score_filter_in, const double max_softclip_filter_in) : paths_index(paths_index_in), library_type(library_type_in), max_pair_frag_length(max_pair_frag_length_in), max_internal_offset(max_internal_offset_in), min_mapq_filter(min_mapq_filter_in), min_best_score_filter(min_best_score_filter_in), max_softclip_filter(max_softclip_filter_in) {}
        
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

    auto alignment_start_nodes = getAlignmentStartNodes(alignment);

    for (auto & start_node: alignment_start_nodes) {

        if (!paths_index.hasNodeId(gbwt::Node::id(start_node))) {

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

    const uint32_t alignment_length = alignment.sequence().size();
    assert(alignment_length > 0);

    if (library_type == "fr") {

        align_search_paths = extendAlignmentPath(AlignmentSearchPath(), alignment, alignment_length);

    } else if (library_type == "rf") {

        AlignmentType alignment_rc = lazy_reverse_complement_alignment(alignment, node_length_func);
        align_search_paths = extendAlignmentPath(AlignmentSearchPath(), alignment_rc, alignment_length);

    } else {

        assert(library_type == "unstranded");
        align_search_paths = extendAlignmentPath(AlignmentSearchPath(), alignment, alignment_length);

        if (!paths_index.bidirectional()) {

            AlignmentType alignment_rc = lazy_reverse_complement_alignment(alignment, node_length_func);
            auto align_search_paths_rc = extendAlignmentPath(AlignmentSearchPath(), alignment_rc, alignment_length);

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
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::Alignment & alignment, const uint32_t alignment_length) const {

    return extendAlignmentPath(align_search_path, alignment, alignment_length, 0);
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::Alignment & alignment, const uint32_t alignment_length, const uint32_t subpath_idx) const {

    assert(alignment.mapping_quality() >= 0);

    vector<AlignmentSearchPath> extended_align_search_path(1, align_search_path);
    
    extended_align_search_path.front().read_stats.emplace_back(ReadAlignmentStats());
    extended_align_search_path.front().read_stats.back().mapq = alignment.mapping_quality();
    extended_align_search_path.front().read_stats.back().score = alignment.score();

    const uint32_t max_right_softclip_length = getMaxAlignmentEndSoftClip(alignment);
    assert(max_right_softclip_length <= alignment_length);

    extended_align_search_path.front().read_stats.back().internal_end_offset.first = alignment_length - max_right_softclip_length;

    extendAlignmentPath(&extended_align_search_path, alignment.path(), true, true);
    return extended_align_search_path;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPath(vector<AlignmentSearchPath> * align_search_paths, const vg::Path & path, const bool is_first_path, const bool is_last_path) const {

    assert(align_search_paths->size() == 1);
    assert(!align_search_paths->front().read_stats.empty());

    auto mapping_it = path.mapping().cbegin();
    assert(mapping_it != path.mapping().cend());

    if (is_first_path) {

        align_search_paths->front().read_stats.back().updateLeftSoftClipLength(path);
        assert(align_search_paths->front().read_stats.back().left_softclip_length.second);
    }

    if (is_last_path) {

        align_search_paths->front().read_stats.back().updateRightSoftClipLength(path);
        assert(align_search_paths->front().read_stats.back().right_softclip_length.second);
    }

    auto end_mapping_it = path.mapping().cend();
    --end_mapping_it;

    while (mapping_it != path.mapping().cend()) {

        auto mapping_read_length = mapping_to_length(*mapping_it);
        const bool is_last_mapping = (is_last_path && mapping_it == end_mapping_it);

        if (align_search_paths->front().read_stats.back().internal_end_offset.second) {

            assert(align_search_paths->size() == 1);
            align_search_paths->front().read_stats.back().updateInternalEndOffset(mapping_read_length, is_last_mapping);

            if (align_search_paths->front().read_stats.back().internal_end_offset.first > max_internal_offset) {

                align_search_paths->front().clear();
                return;          
            }

        } else {

            AlignmentSearchPath main_align_search_paths;

            if (!main_align_search_paths.isEmpty() && max_internal_offset > 0 && !main_align_search_paths.read_stats.back().internal_start_offset.second) {

                main_align_search_paths = align_search_paths->front();
            }

            for (auto & align_search_path: *align_search_paths) {

                if (align_search_path.read_stats.back().internal_end_offset.second) {

                    align_search_path.read_stats.back().updateInternalEndOffset(mapping_read_length, is_last_mapping);

                    if (align_search_path.read_stats.back().internal_end_offset.first > max_internal_offset) {
                        
                        align_search_path.clear();
                    }
                
                } else {

                    extendAlignmentPath(&align_search_path, *mapping_it);
                }
            }

            if (!main_align_search_paths.isEmpty()) {

                assert(!main_align_search_paths.read_stats.back().internal_end_offset.second);
                assert(main_align_search_paths.search_state.size() >= align_search_paths->front().search_state.size());

                if (main_align_search_paths.search_state.size() > align_search_paths->front().search_state.size()) {

                    const uint32_t align_length_left = max(int32_t(0), static_cast<int32_t>(main_align_search_paths.read_stats.back().internal_end_offset.first) - static_cast<int32_t>(main_align_search_paths.read_stats.back().length));

                    if (align_length_left <= max_internal_offset) {

                        auto internal_end_read_stats = main_align_search_paths.read_stats.back();
                        internal_end_read_stats.updateInternalEndOffset(mapping_read_length, is_last_mapping);

                        if (internal_end_read_stats.internal_end_offset.first <= max_internal_offset) {

                            align_search_paths->emplace_back(main_align_search_paths);
                            align_search_paths->back().read_stats.back() = internal_end_read_stats;               
                        }
                    }
                }
            }

            if (align_search_paths->front().path.size() > 1 && max_internal_offset > 0 && !align_search_paths->front().read_stats.back().internal_start_offset.second) {

                auto internal_start_read_stats = align_search_paths->front().read_stats.back();
                internal_start_read_stats.updateInternalStartOffset(align_search_paths->front().read_stats.back().length, true);

                if (internal_start_read_stats.internal_start_offset.first <= max_internal_offset) {

                    align_search_paths->emplace_back(AlignmentSearchPath());
                    align_search_paths->back().read_stats = vector<ReadAlignmentStats>(1, internal_start_read_stats);

                    extendAlignmentPath(&(align_search_paths->back()), *mapping_it);

                    if (align_search_paths->back().isEmpty()) {

                        align_search_paths->pop_back();

                    } else {

                        assert(align_search_paths->back().search_state.size() >= align_search_paths->front().search_state.size());

                        if (align_search_paths->back().search_state.size() == align_search_paths->front().search_state.size()) {

                            align_search_paths->pop_back();
                        }
                    }
                }
            }
        }

        for (auto & align_search_path: *align_search_paths) {

            align_search_path.read_stats.back().length += mapping_read_length;
        }

        ++mapping_it;
    }
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPath(AlignmentSearchPath * align_search_path, const vg::Mapping & mapping) const {

    auto cur_node = mapping_to_gbwt(mapping);
    auto mapping_read_length = mapping_to_length(mapping);

    if (align_search_path->path.empty()) {

        if (mapping_read_length > 0) {

            assert(align_search_path->search_state.node == gbwt::ENDMARKER);

            align_search_path->path.emplace_back(cur_node);
            align_search_path->search_state = paths_index.find(cur_node);
      
            align_search_path->start_offset = mapping.position().offset();
            align_search_path->end_offset = mapping.position().offset() + mapping_from_length(mapping);
        }

    } else {

        if (align_search_path->path.back() == cur_node || mapping_read_length > 0) {

            bool is_cycle_visit = false;

            if (align_search_path->path.back() == cur_node && mapping.position().offset() != align_search_path->end_offset) {

                assert(mapping.position().offset() == 0);
                is_cycle_visit = true;      
            }

            if (align_search_path->path.back() != cur_node || is_cycle_visit) {

                align_search_path->path.emplace_back(cur_node);
                align_search_path->search_state = paths_index.extend(align_search_path->search_state, cur_node);
            } 

            align_search_path->end_offset = mapping.position().offset() + mapping_from_length(mapping);
        }
    }
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::MultipathAlignment & alignment, const uint32_t alignment_length) const {

    vector<AlignmentSearchPath> extended_align_search_paths;

    for (auto & subpath_start_idx: alignment.start()) {

        auto cur_extended_align_search_paths = extendAlignmentPath(align_search_path, alignment, alignment_length, subpath_start_idx);
        extended_align_search_paths.insert(extended_align_search_paths.end(), make_move_iterator(cur_extended_align_search_paths.begin()), make_move_iterator(cur_extended_align_search_paths.end()));
    }

    return extended_align_search_paths;
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::MultipathAlignment & alignment, const uint32_t alignment_length, const uint32_t subpath_idx) const {

    assert(alignment.mapping_quality() >= 0);

    vector<AlignmentSearchPath> extended_align_search_path(1, align_search_path);
    
    extended_align_search_path.front().read_stats.emplace_back(ReadAlignmentStats());
    extended_align_search_path.front().read_stats.back().mapq = alignment.mapping_quality();

    const uint32_t max_right_softclip_length = getMaxAlignmentEndSoftClip(alignment);
    assert(max_right_softclip_length <= alignment_length);

    extended_align_search_path.front().read_stats.back().internal_end_offset.first = alignment_length - max_right_softclip_length;

    extendAlignmentPaths(&extended_align_search_path, alignment.subpath(), subpath_idx);
            
    return extended_align_search_path;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPaths(vector<AlignmentSearchPath> * align_search_paths, const google::protobuf::RepeatedPtrField<vg::Subpath> & subpaths, const uint32_t subpath_idx) const {

    std::queue<pair<AlignmentSearchPath, uint32_t> > align_search_paths_queue;

    for (auto & align_search_path: *align_search_paths) {

        align_search_paths_queue.push(make_pair(align_search_path, subpath_idx));
    }

    align_search_paths->clear();

    // Perform depth-first alignment path extension.
    while (!align_search_paths_queue.empty()) {

        auto & cur_align_search_path = align_search_paths_queue.front();
        const vg::Subpath & subpath = subpaths.Get(cur_align_search_path.second);

        vector<AlignmentSearchPath> extended_align_search_path(1, cur_align_search_path.first);
        extended_align_search_path.front().read_stats.back().score += subpath.score();

        extendAlignmentPath(&extended_align_search_path, subpath.path(), cur_align_search_path.second == subpath_idx, subpath.next_size() == 0);

        for (auto & align_search_path: extended_align_search_path) {

            assert(align_search_path.read_stats.back().left_softclip_length.second);
            assert(align_search_path.read_stats.back().left_softclip_length.first <= align_search_path.read_stats.back().length);

            if (!align_search_path.isEmpty() || (align_search_path.read_stats.back().length - align_search_path.read_stats.back().left_softclip_length.first) <= max_internal_offset) {

                if (subpath.next_size() > 0) {

                    for (auto & next_subpath_idx: subpath.next()) {

                        align_search_paths_queue.push(make_pair(align_search_path, next_subpath_idx));
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

    const uint32_t alignment_1_length = alignment_1.sequence().size();
    assert(alignment_1_length > 0);

    const uint32_t alignment_2_length = alignment_2.sequence().size();
    assert(alignment_2_length > 0);

    if (library_type == "fr") {

        pairAlignmentPaths(&paired_align_search_paths, alignment_1, alignment_1_length, alignment_2_rc, alignment_2_length);

    } else if (library_type == "rf") {

        AlignmentType alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_length_func);
        pairAlignmentPaths(&paired_align_search_paths, alignment_2, alignment_2_length, alignment_1_rc, alignment_1_length);

    } else {

        assert(library_type == "unstranded");

        AlignmentType alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_length_func);
        pairAlignmentPaths(&paired_align_search_paths, alignment_1, alignment_1_length, alignment_2_rc, alignment_2_length);

        if (!paths_index.bidirectional()) {

            AlignmentType alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_length_func);
            pairAlignmentPaths(&paired_align_search_paths, alignment_2, alignment_2_length, alignment_1_rc, alignment_1_length);
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
void AlignmentPathFinder<AlignmentType>::mergeAlignmentPaths(AlignmentSearchPath * main_align_search_path, uint32_t main_path_start_idx, const AlignmentSearchPath & second_align_search_path) const {

    assert(!main_align_search_path->isEmpty());
    assert(!second_align_search_path.isEmpty());

    assert(main_path_start_idx < main_align_search_path->path.size());

    assert(main_align_search_path->read_stats.size() == 1);
    assert(second_align_search_path.read_stats.size() == 1);

    assert(main_align_search_path->read_stats.back().internal_start_offset.first <= max_internal_offset);
    assert(second_align_search_path.read_stats.back().internal_start_offset.first <= max_internal_offset);

    assert(!main_align_search_path->read_stats.back().internal_start_offset.second || main_align_search_path->read_stats.back().internal_start_offset.first <= max_internal_offset);
    assert(!second_align_search_path.read_stats.back().internal_start_offset.second || second_align_search_path.read_stats.back().internal_start_offset.first <= max_internal_offset);

    assert(!main_align_search_path->read_stats.back().internal_end_offset.second || main_align_search_path->read_stats.back().internal_end_offset.first <= max_internal_offset);
    assert(!second_align_search_path.read_stats.back().internal_end_offset.second || second_align_search_path.read_stats.back().internal_end_offset.first <= max_internal_offset);

    if (second_align_search_path.path.size() < main_align_search_path->path.size() - main_path_start_idx) {

        main_align_search_path->clear();
        return;  
    }

    if (main_path_start_idx == 0) {

        const int32_t main_read_left_offset = static_cast<int32_t>(main_align_search_path->start_offset) - static_cast<int32_t>(main_align_search_path->read_stats.back().clippedOffsetLeftBases());
        const int32_t second_read_left_offset = static_cast<int32_t>(second_align_search_path.start_offset) - static_cast<int32_t>(second_align_search_path.read_stats.back().clippedOffsetLeftBases());

        if (second_read_left_offset < main_read_left_offset) {
                
            main_align_search_path->clear();
            return;    
        } 
    }

    uint32_t second_path_start_idx = 0;

    while (main_path_start_idx < main_align_search_path->path.size()) {

        assert(second_path_start_idx < second_align_search_path.path.size());

        if (main_align_search_path->path.at(main_path_start_idx) != second_align_search_path.path.at(second_path_start_idx)) {

            main_align_search_path->clear();
            return; 
        }   

        if (main_path_start_idx + 1 == main_align_search_path->path.size()) {
            
            if (second_path_start_idx + 1 == second_align_search_path.path.size()) {

                const uint32_t main_read_right_offset = main_align_search_path->end_offset + main_align_search_path->read_stats.back().clippedOffsetRightBases();
                const uint32_t second_read_right_offset = second_align_search_path.end_offset + second_align_search_path.read_stats.back().clippedOffsetRightBases();

                if (second_read_right_offset < main_read_right_offset) {
                        
                    main_align_search_path->clear();
                    return;   
                }

                if (main_path_start_idx == 0) {

                    assert(second_path_start_idx == 0);

                    main_align_search_path->insert_length += (static_cast<int32_t>(max(main_align_search_path->start_offset, second_align_search_path.start_offset)) - static_cast<int32_t>(min(main_align_search_path->end_offset, second_align_search_path.end_offset)));             
                
                } else if (second_path_start_idx == 0) {

                    main_align_search_path->insert_length += (static_cast<int32_t>(second_align_search_path.start_offset) - static_cast<int32_t>(min(main_align_search_path->end_offset, second_align_search_path.end_offset)));             

                } else {

                    main_align_search_path->insert_length -= min(main_align_search_path->end_offset, second_align_search_path.end_offset);
                }

            } else if (second_path_start_idx == 0) {

                main_align_search_path->insert_length += (static_cast<int32_t>(second_align_search_path.start_offset) - static_cast<int32_t>(main_align_search_path->end_offset));             
            
            } else {

                main_align_search_path->insert_length -= main_align_search_path->end_offset;
            } 

        } else if (second_path_start_idx == 0) {

            assert(main_align_search_path->path.size() > 1);
            assert(second_align_search_path.path.size() > 1);

            const uint32_t node_length = paths_index.nodeLength(gbwt::Node::id(main_align_search_path->path.at(main_path_start_idx)));
            assert(second_align_search_path.start_offset <= node_length);

            if (main_path_start_idx == 0) {

                assert(main_align_search_path->start_offset <= node_length);
                main_align_search_path->insert_length -= (node_length - max(main_align_search_path->start_offset, second_align_search_path.start_offset));

            } else {

                main_align_search_path->insert_length -= (node_length - second_align_search_path.start_offset);             
            }
   
        } else {

            main_align_search_path->insert_length -= paths_index.nodeLength(gbwt::Node::id(main_align_search_path->path.at(main_path_start_idx)));
        } 

        ++main_path_start_idx;
        ++second_path_start_idx;
    }

    main_align_search_path->end_offset = second_align_search_path.end_offset;
    main_align_search_path->read_stats.emplace_back(second_align_search_path.read_stats.front());
    
    assert(main_path_start_idx == main_align_search_path->path.size());
    assert(second_path_start_idx <= second_align_search_path.path.size());

    while (second_path_start_idx < second_align_search_path.path.size()) {

        main_align_search_path->path.emplace_back(second_align_search_path.path.at(second_path_start_idx));
        main_align_search_path->search_state = paths_index.extend(main_align_search_path->search_state, main_align_search_path->path.back());

        if (main_align_search_path->isEmpty()) {

            break;            
        }

        ++second_path_start_idx;
    }
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::pairAlignmentPaths(vector<AlignmentSearchPath> * paired_align_search_paths, const AlignmentType & start_alignment, const uint32_t start_alignment_length, const AlignmentType & end_alignment, const uint32_t end_alignment_length) const {

    auto start_align_search_paths = extendAlignmentPath(AlignmentSearchPath(), start_alignment, start_alignment_length);
    auto end_align_search_paths = extendAlignmentPath(AlignmentSearchPath(), end_alignment, end_alignment_length);

    if (start_align_search_paths.empty() || end_align_search_paths.empty()) {

        return;
    }

    sort(start_align_search_paths.begin(), start_align_search_paths.end());
    sort(end_align_search_paths.begin(), end_align_search_paths.end());

    spp::sparse_hash_map<gbwt::node_type, vector<uint32_t> > end_alignment_start_nodes;

    for (size_t i = 0; i < end_align_search_paths.size(); ++i) {

        if (!end_align_search_paths.at(i).isEmpty()) {

            if (i > 0) {

                if (end_align_search_paths.at(i) == end_align_search_paths.at(i - 1)) {

                    continue;
                }
            }

            const AlignmentSearchPath & end_align_search_path = end_align_search_paths.at(i);

            assert(end_align_search_path.read_stats.size() == 1);
            assert(end_align_search_path.read_stats.back().length == end_alignment_length);

            auto end_alignment_start_nodes_it = end_alignment_start_nodes.emplace(end_align_search_path.path.front(), vector<uint32_t>());
            end_alignment_start_nodes_it.first->second.emplace_back(i);
        }
    }

    std::queue<pair<AlignmentSearchPath, bool> > paired_align_search_path_queue;

    for (size_t i = 0; i < start_align_search_paths.size(); ++i) {

        if (start_align_search_paths.at(i).isEmpty()) {

            continue;
        }

        if (i > 0) {

            if (start_align_search_paths.at(i) == start_align_search_paths.at(i - 1)) {

                continue;
            }
        }

        const AlignmentSearchPath & start_align_search_path = start_align_search_paths.at(i);

        assert(start_align_search_path.read_stats.size() == 1);
        assert(start_align_search_path.read_stats.back().length == start_alignment_length);

        auto node_length = paths_index.nodeLength(gbwt::Node::id(start_align_search_path.search_state.node));
        assert(start_align_search_path.end_offset <= node_length);

        paired_align_search_path_queue.push(make_pair(start_align_search_path, false));

        paired_align_search_path_queue.back().first.insert_length += (node_length - start_align_search_path.end_offset);
        paired_align_search_path_queue.back().first.end_offset = node_length;

        for (auto & end_alignment_start_nodes: end_alignment_start_nodes) {

            auto path_it = find(start_align_search_path.path.begin(), start_align_search_path.path.end(), end_alignment_start_nodes.first); 

            while (path_it != start_align_search_path.path.end()) {

                auto main_path_start_idx = path_it - start_align_search_path.path.begin();

                for (auto end_alignment_idx: end_alignment_start_nodes.second) {

                    AlignmentSearchPath complete_paired_align_search_path = start_align_search_path;
                    mergeAlignmentPaths(&complete_paired_align_search_path, main_path_start_idx, end_align_search_paths.at(end_alignment_idx));

                    if (!complete_paired_align_search_path.isEmpty() && complete_paired_align_search_path.fragmentLength() <= max_pair_frag_length) {

                        paired_align_search_paths->emplace_back(complete_paired_align_search_path);                         
                    }
                }

                ++path_it;
                path_it = find(path_it, start_align_search_path.path.end(), end_alignment_start_nodes.first); 
            }
        }
    }

    auto max_left_softclip_length = getMaxAlignmentStartSoftClip(end_alignment);
    assert(max_left_softclip_length <= end_alignment_length);

    // Perform depth-first path extension.
    while (!paired_align_search_path_queue.empty()) {

        AlignmentSearchPath * cur_paired_align_search_path = &(paired_align_search_path_queue.front().first);
        
        assert(cur_paired_align_search_path->search_state.node != gbwt::ENDMARKER);
        assert(!cur_paired_align_search_path->isEmpty());

        if (paired_align_search_path_queue.front().second) {

            auto end_alignment_start_nodes_it = end_alignment_start_nodes.find(cur_paired_align_search_path->search_state.node);

            if (end_alignment_start_nodes_it != end_alignment_start_nodes.end()) {

                for (auto end_alignment_idx: end_alignment_start_nodes_it->second) {

                    AlignmentSearchPath complete_paired_align_search_path = *cur_paired_align_search_path;
                    complete_paired_align_search_path.insert_length -= complete_paired_align_search_path.end_offset;

                    complete_paired_align_search_path.end_offset = end_align_search_paths.at(end_alignment_idx).start_offset;
                    complete_paired_align_search_path.insert_length += complete_paired_align_search_path.end_offset;

                    mergeAlignmentPaths(&complete_paired_align_search_path, cur_paired_align_search_path->path.size() - 1, end_align_search_paths.at(end_alignment_idx));

                    if (!complete_paired_align_search_path.isEmpty() && complete_paired_align_search_path.fragmentLength() <= max_pair_frag_length) {

                        paired_align_search_paths->emplace_back(complete_paired_align_search_path);                         
                    }
                }
            }
        }

        paired_align_search_path_queue.front().second = true;
           
        if (cur_paired_align_search_path->fragmentLength() + end_alignment_length - max_left_softclip_length > max_pair_frag_length) {

            paired_align_search_path_queue.pop();
            continue;
        }

        auto out_edges = paths_index.edges(cur_paired_align_search_path->search_state.node);

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

                auto extended_path = paths_index.extend(cur_paired_align_search_path->search_state, out_edges_it->first);

                // Add new extension to queue if not empty (path found).
                if (!extended_path.empty()) { 

                    paired_align_search_path_queue.push(make_pair(*cur_paired_align_search_path, true));

                    paired_align_search_path_queue.back().first.path.emplace_back(extended_path.node);
                    paired_align_search_path_queue.back().first.search_state = extended_path;
                    paired_align_search_path_queue.back().first.end_offset = paths_index.nodeLength(gbwt::Node::id(extended_path.node));
                    paired_align_search_path_queue.back().first.insert_length += paired_align_search_path_queue.back().first.end_offset;
                }
            }

            ++out_edges_it;
        }

        if (out_edges.begin()->first != gbwt::ENDMARKER) {
            
            cur_paired_align_search_path->search_state = paths_index.extend(cur_paired_align_search_path->search_state, out_edges.begin()->first);

            // End current extension if empty (no haplotypes found). 
            if (cur_paired_align_search_path->search_state.empty()) { 

                paired_align_search_path_queue.pop(); 

            } else {

                cur_paired_align_search_path->path.emplace_back(cur_paired_align_search_path->search_state.node);
                cur_paired_align_search_path->end_offset = paths_index.nodeLength(gbwt::Node::id(cur_paired_align_search_path->search_state.node));
                cur_paired_align_search_path->insert_length += cur_paired_align_search_path->end_offset;
            }
    
        } else {

            paired_align_search_path_queue.pop();
        }
    }
}

template<class AlignmentType>
vector<gbwt::node_type> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodes(const vg::Alignment & alignment) const {

    vector<gbwt::node_type> alignment_start_nodes;

    assert(alignment.path().mapping_size() > 0);
    alignment_start_nodes.emplace_back(mapping_to_gbwt(alignment.path().mapping(0)));

    return alignment_start_nodes;
}

template<class AlignmentType>
vector<gbwt::node_type> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodes(const vg::MultipathAlignment & alignment) const {

    vector<gbwt::node_type> alignment_start_nodes;

    for (auto & start_idx: alignment.start()) {

        assert(alignment.subpath(start_idx).path().mapping_size() > 0);
        alignment_start_nodes.emplace_back(mapping_to_gbwt(alignment.subpath(start_idx).path().mapping(0)));
    }

    return alignment_start_nodes;
}

template<class AlignmentType>
uint32_t AlignmentPathFinder<AlignmentType>::getMaxAlignmentStartSoftClip(const vg::Alignment & alignment) const {

    ReadAlignmentStats read_stats;
    read_stats.updateLeftSoftClipLength(alignment.path());

    assert(read_stats.left_softclip_length.second);
    return read_stats.left_softclip_length.first;
}

template<class AlignmentType>
uint32_t AlignmentPathFinder<AlignmentType>::getMaxAlignmentStartSoftClip(const vg::MultipathAlignment & alignment) const {

    uint32_t max_left_softclip_length = 0;
    ReadAlignmentStats read_stats;

    for (auto & start_idx: alignment.start()) {

        read_stats.left_softclip_length.second = false;
        read_stats.updateLeftSoftClipLength(alignment.subpath(start_idx).path());

        assert(read_stats.left_softclip_length.second);
        max_left_softclip_length = max(max_left_softclip_length, read_stats.left_softclip_length.first);
    }

    return max_left_softclip_length;
}

template<class AlignmentType>
uint32_t AlignmentPathFinder<AlignmentType>::getMaxAlignmentEndSoftClip(const vg::Alignment & alignment) const {

    ReadAlignmentStats read_stats;
    read_stats.updateRightSoftClipLength(alignment.path());

    assert(read_stats.right_softclip_length.second);
    return read_stats.right_softclip_length.first;
}

template<class AlignmentType>
uint32_t AlignmentPathFinder<AlignmentType>::getMaxAlignmentEndSoftClip(const vg::MultipathAlignment & alignment) const {

    uint32_t max_right_softclip_length = 0;
    ReadAlignmentStats read_stats;

    for (auto & subpath: alignment.subpath()) {

        if (subpath.next_size() == 0) {

            read_stats.right_softclip_length.second = false;
            read_stats.updateRightSoftClipLength(subpath.path());

            assert(read_stats.right_softclip_length.second);
            max_right_softclip_length = max(max_right_softclip_length, read_stats.right_softclip_length.first);
        }
    }

    return max_right_softclip_length;
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

        if (!align_search_path.isEmpty()) {

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

