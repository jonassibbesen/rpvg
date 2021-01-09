
#include "alignment_path_finder.hpp"

#include <assert.h>

#include "utils.hpp"

//#define debug

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

    auto alignment_start_nodes_index = getAlignmentStartNodesIndex(alignment);

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

    return extendAlignmentPath(align_search_path, alignment, 0);
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::Alignment & alignment, const uint32_t subpath_start_idx) const {

    assert(alignment.mapping_quality() >= 0);

    vector<AlignmentSearchPath> extended_align_search_path(1, align_search_path);
    
    extended_align_search_path.front().read_stats.emplace_back(ReadAlignmentStats());
    extended_align_search_path.front().read_stats.back().mapq = alignment.mapping_quality();
    extended_align_search_path.front().read_stats.back().score = alignment.score();
    
    extendAlignmentPath(&extended_align_search_path.front(), alignment.path());

    if (extended_align_search_path.front().search_state.empty()) {

        return vector<AlignmentSearchPath>();
    
    } else {

        return extended_align_search_path;
    }
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPath(AlignmentSearchPath * align_search_path, const vg::Path & path) const {

    assert(align_search_path->path_idx <= align_search_path->path.size());
    
    auto mapping_it = path.mapping().cbegin();
    assert(mapping_it != path.mapping().cend());

    const bool is_first_node = (align_search_path->read_stats.back().left_softclip_length == -1);

    if (!align_search_path->path.empty() && align_search_path->path_idx == 0) {

        assert(is_first_node);
 
        if (mapping_it->position().offset() < align_search_path->start_offset) {
 
            align_search_path->search_state = gbwt::SearchState();
            assert(align_search_path->search_state.empty());

            return;      
        } 
    }

    if (is_first_node) {

        align_search_path->read_stats.back().left_softclip_length = 0;

        assert(mapping_it->edit_size() > 0);
        const vg::Edit & first_edit = mapping_it->edit(0);

        if (first_edit.from_length() == 0) {

            align_search_path->read_stats.back().left_softclip_length = first_edit.to_length();   
        } 
    }

    auto mapping_rit = path.mapping().rbegin();
    assert(mapping_rit != path.mapping().rend());   

    align_search_path->read_stats.back().right_softclip_length = 0;

    assert(mapping_rit->edit_size() > 0);
    const vg::Edit & last_edit = mapping_rit->edit(mapping_rit->edit_size() - 1);

    if (last_edit.from_length() == 0) {

        align_search_path->read_stats.back().right_softclip_length = last_edit.to_length();   
    }

    while (align_search_path->path_idx < align_search_path->path.size() && mapping_it != path.mapping().cend()) {

        auto cur_node = mapping_to_gbwt(*mapping_it);

        bool is_multi_visit = false;

        if (!is_first_node) {

            assert(align_search_path->path_idx > 0);

            if (align_search_path->path.at(align_search_path->path_idx - 1) == cur_node && align_search_path->path_offset == mapping_it->position().offset()) {

                is_multi_visit = true;
            }
        }

        if (align_search_path->path.at(align_search_path->path_idx) == cur_node || is_multi_visit) {

            if (!is_multi_visit) {

                ++align_search_path->path_idx;
            }

            if (align_search_path->path_idx == align_search_path->path.size()) {

                auto new_end_offset = mapping_it->position().offset() + mapping_from_length(*mapping_it);

                if (new_end_offset < align_search_path->end_offset) {

                    align_search_path->search_state = gbwt::SearchState();
                    assert(align_search_path->search_state.empty());

                    return;           
                }

                align_search_path->insert_length += (static_cast<int32_t>(mapping_it->position().offset()) - static_cast<int32_t>(align_search_path->end_offset));
                align_search_path->end_offset = new_end_offset;
            
            } else if (!is_multi_visit) {

                align_search_path->insert_length -= mapping_from_length(*mapping_it);
            }
    
            align_search_path->path_offset = mapping_it->position().offset() + mapping_from_length(*mapping_it);
            align_search_path->read_stats.back().length += mapping_to_length(*mapping_it);

        } else {

            align_search_path->search_state = gbwt::SearchState();
            assert(align_search_path->search_state.empty());  
    
            return;  
        } 

        ++mapping_it;
    }

    while (mapping_it != path.mapping().cend()) {

        assert(align_search_path->path_offset == align_search_path->end_offset);

        auto cur_node = mapping_to_gbwt(*mapping_it);

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

        align_search_path->path_offset = mapping_it->position().offset() + mapping_from_length(*mapping_it);
        align_search_path->end_offset = align_search_path->path_offset;

        align_search_path->read_stats.back().length += mapping_to_length(*mapping_it);

        if (align_search_path->search_state.empty()) {

            break;
        }

        ++mapping_it;
    }
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::MultipathAlignment & alignment) const {

    vector<AlignmentSearchPath> extended_align_search_paths;

    for (auto & start_idx: alignment.start()) {

        auto cur_extended_align_search_paths = extendAlignmentPath(align_search_path, alignment, start_idx);
        extended_align_search_paths.insert(extended_align_search_paths.end(), cur_extended_align_search_paths.begin(), cur_extended_align_search_paths.end());
    }

    return extended_align_search_paths;
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentSearchPath & align_search_path, const vg::MultipathAlignment & alignment, const uint32_t subpath_start_idx) const {

    assert(alignment.mapping_quality() >= 0);

    vector<AlignmentSearchPath> extended_align_search_path(1, align_search_path);
    
    extended_align_search_path.front().read_stats.emplace_back(ReadAlignmentStats());
    extended_align_search_path.front().read_stats.back().mapq = alignment.mapping_quality();

    extendAlignmentPaths(&extended_align_search_path, alignment.subpath(), subpath_start_idx);
            
    return extended_align_search_path;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPaths(vector<AlignmentSearchPath> * align_search_paths, const google::protobuf::RepeatedPtrField<vg::Subpath> & subpaths, const uint32_t subpath_start_idx) const {

    std::queue<pair<AlignmentSearchPath, uint32_t> > align_search_paths_queue;

    for (auto & align_search_path: *align_search_paths) {

        align_search_paths_queue.push(make_pair(align_search_path, subpath_start_idx));
    }

    align_search_paths->clear();

    // Perform depth-first alignment path extension.
    while (!align_search_paths_queue.empty()) {

        auto & cur_align_search_path = align_search_paths_queue.front();

        const vg::Subpath & subpath = subpaths.Get(cur_align_search_path.second);

        cur_align_search_path.first.read_stats.back().score += subpath.score();
        extendAlignmentPath(&cur_align_search_path.first, subpath.path());

        if (cur_align_search_path.first.path.empty() || !cur_align_search_path.first.search_state.empty()) {

            if (subpath.next_size() > 0) {

                for (auto & next_subpath_idx: subpath.next()) {

                    align_search_paths_queue.push(make_pair(cur_align_search_path.first, next_subpath_idx));
                }

            } else {

                align_search_paths->emplace_back(cur_align_search_path.first);
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

    string debug_paths = "";
    int32_t debug_idx = -1;

    string debug_paths2 = "";
    int32_t debug_idx2 = -1;

    for (size_t i = 0; i < paired_align_search_paths.size(); ++i) {

        if (paired_align_search_paths.at(i).isComplete()) {

            for (auto & path_id: paths_index.locatePathIds(paired_align_search_paths.at(i).search_state)) {

                auto path_name = paths_index.pathName(path_id);

                if (
                    path_name == "ENST00000346234.6_28" || 
                     path_name == "ENST00000461096.6_28" || 
                     path_name == "ENST00000317897.4_26" || 
                     path_name == "ENST00000594159.1_9" || 
                     path_name == "ENST00000396062.3_21" || 
                     path_name == "ENST00000296677.4_79" || 
                     path_name == "ENST00000568280.1_20" || 
                     path_name == "ENST00000370206.8_3" || 
                     path_name == "ENST00000379375.5_66" || 
                     path_name == "ENST00000275766.1_44" || 
                     path_name == "ENST00000317811.5" || 
                     path_name == "ENST00000378045.4" || 
                     path_name == "ENST00000329235.6_9" || 
                     path_name == "ENST00000596580.2_72" || 
                     path_name == "ENST00000368847.4_38" ||
                     path_name == "ENST00000580018.3_15" || 
                     path_name == "ENST00000374259.7" || 
                     path_name == "ENST00000325307.11" || 
                     path_name == "ENST00000216252.3_19" || 
                     path_name == "ENST00000271638.2" ||
                     path_name == "ENST00000253788.11_9" ||
                     path_name == "ENST00000323699.8_55" ||
                     path_name == "ENST00000378119.8_158" ||
                     path_name == "ENST00000228506.7_157" ||
                     path_name == "ENST00000340913.10" ||
                     path_name == "ENST00000592588.6" ||
                     path_name == "ENST00000584577.5" ||
                     path_name == "ENST00000223641.4" ||
                     path_name == "ENST00000591776.5" ||
                     path_name == "ENST00000221975.6" 
                ) {   

                    debug_paths = path_name; 
                    debug_idx = i;         
                
                } else if (
                    path_name == "ENST00000346234.6_38" || 
                     path_name == "ENST00000461096.6" || 
                     path_name == "ENST00000317897.4_29" || 
                     path_name == "ENST00000594159.1_55" || 
                     path_name == "ENST00000396062.3" || 
                     path_name == "ENST00000296677.4_88" || 
                     path_name == "ENST00000568280.1_45" || 
                     path_name == "ENST00000370206.8_55" || 
                     path_name == "ENST00000379375.5_43" || 
                     path_name == "ENST00000275766.1" || 
                     path_name == "ENST00000317811.5_46" || 
                     path_name == "ENST00000378045.4_11" || 
                     path_name == "ENST00000329235.6_14" || 
                     path_name == "ENST00000596580.2_175" || 
                     path_name == "ENST00000368847.4_110" ||
                     path_name == "ENST00000580018.3_16" || 
                     path_name == "ENST00000374259.7_19" || 
                     path_name == "ENST00000325307.11_24" || 
                     path_name == "ENST00000216252.3_24" || 
                     path_name == "ENST00000271638.2_8" || 
                     path_name == "ENST00000253788.11_6" ||
                     path_name == "ENST00000323699.8_56" ||
                     path_name == "ENST00000378119.8_161" ||
                     path_name == "ENST00000228506.7_171" ||
                     path_name == "ENST00000340913.10_40" ||
                     path_name == "ENST00000592588.6_3" ||
                     path_name == "ENST00000584577.5_6" ||
                     path_name == "ENST00000223641.4_1" ||
                     path_name == "ENST00000591776.5_4" ||
                     path_name == "ENST00000221975.6_4" 
                    ) {  

                    debug_paths2 = path_name; 
                    debug_idx2 = i;         
                }                
            }
        }
    }

    if (debug_idx != debug_idx2) {

        #pragma omp critical
        {
            cerr << "\n\n###" << endl;
            cerr << debug_paths << endl;
            cerr << debug_idx << endl;

            if (debug_idx >= 0) {

                cerr << paired_align_search_paths.at(debug_idx) << endl;
            }

            cerr << debug_paths2 << endl;
            cerr << debug_idx2 << endl;

            if (debug_idx2 >= 0) {

                cerr << paired_align_search_paths.at(debug_idx2) << endl;
            }
            
            cerr << endl;
            cerr << paired_align_search_paths << endl;
            cerr << endl;
            cerr << pb2json(alignment_1) << endl;
            cerr << string_quality_short_to_char(alignment_1.quality()) << endl;
            cerr << endl;
            cerr << pb2json(alignment_2) << endl;
            cerr << string_quality_short_to_char(alignment_2.quality()) << endl;
        }
    }

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
    auto end_alignment_start_nodes_index = getAlignmentStartNodesIndex(end_alignment);

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
multimap<gbwt::node_type, uint32_t> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodesIndex(const vg::Alignment & alignment) const {

    multimap<gbwt::node_type, uint32_t> alignment_start_nodes_index;

    assert(alignment.path().mapping_size() > 0);
    alignment_start_nodes_index.emplace(mapping_to_gbwt(alignment.path().mapping(0)), 0);

    return alignment_start_nodes_index;
}

template<class AlignmentType>
multimap<gbwt::node_type, uint32_t> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodesIndex(const vg::MultipathAlignment & alignment) const {

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

