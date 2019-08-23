
#include "alignment_path_finder.hpp"

#include <assert.h>

#include "utils.hpp"

// #define debug


template<class AlignmentType>
AlignmentPathFinder<AlignmentType>::AlignmentPathFinder(const vg::Graph & graph, const gbwt::GBWT & paths_index_in, const int32_t max_pair_seq_length_in) : paths_index(paths_index_in), max_pair_seq_length(max_pair_seq_length_in) {

    node_seq_lengths = vector<int32_t>(graph.node_size() + 1, 0);

    for (auto & node: graph.node()) {

        if (node.id() >= node_seq_lengths.size()) {

            node_seq_lengths.resize(node.id() + 1, 0);
        }

        assert(node_seq_lengths.at(node.id()) == 0);
        node_seq_lengths.at(node.id()) = node.sequence().size();
    }
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::setMaxPairSeqLength(const int32_t max_pair_seq_length_in) {

    max_pair_seq_length = max_pair_seq_length_in;
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::findAlignmentPaths(const AlignmentType & alignment) const {
            
    AlignmentPath align_path;
    return extendAlignmentPath(align_path, alignment);
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::findAlignmentPathsIds(const AlignmentType & alignment) const {

    auto align_paths = findAlignmentPaths(alignment);

    for (auto & align_path: align_paths) {

        align_path.ids = paths_index.locate(align_path.search);
    }

    return align_paths;
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment, const vector<gbwt::node_type> & stop_nodes) const {

    return extendAlignmentPath(align_path, alignment, 0, stop_nodes);
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment, const int32_t sub_path_start_idx, const vector<gbwt::node_type> & stop_nodes) const {

    vector<AlignmentPath> extended_align_path(1, align_path);
    extended_align_path.front().mapqs.emplace_back(alignment.mapping_quality());
    extended_align_path.front().scores.emplace_back(alignment.score());
    
    extendAlignmentPath(&extended_align_path.front(), alignment.path(), stop_nodes);

    if (extended_align_path.front().search.empty()) {

        return vector<AlignmentPath>();
    
    } else {

        return extended_align_path;
    }
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPath(AlignmentPath * align_path, const vg::Path & path, const vector<gbwt::node_type> & stop_nodes) const {
    
    auto mapping_it = path.mapping().cbegin();
    assert(mapping_it != path.mapping().cend());

    if (align_path->search.node == mapping_to_gbwt(*mapping_it)) {
        
        align_path->end_offset = mapping_it->position().offset() + mapping_from_length(*mapping_it);
        align_path->seq_length += mapping_it->position().offset() + mapping_to_length(*mapping_it);

        if (find(stop_nodes.begin(), stop_nodes.end(), mapping_to_gbwt(*mapping_it)) != stop_nodes.end()) {

            mapping_it = path.mapping().cend();
        
        } else {

            ++mapping_it;
        }
    } 

    while (mapping_it != path.mapping().cend()) {

        if (align_path->seq_length == 0) {

            assert(align_path->search.node == gbwt::ENDMARKER);
            align_path->search = paths_index.find(mapping_to_gbwt(*mapping_it));
        
        } else {

            align_path->search = paths_index.extend(align_path->search, mapping_to_gbwt(*mapping_it));                
        }

        align_path->end_offset = mapping_it->position().offset() + mapping_from_length(*mapping_it);
        align_path->seq_length += mapping_to_length(*mapping_it);

        if (find(stop_nodes.begin(), stop_nodes.end(), mapping_to_gbwt(*mapping_it)) != stop_nodes.end()) {

            break;
        }

        ++mapping_it;
    }
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment, const vector<gbwt::node_type> & stop_nodes) const {

    vector<AlignmentPath> extended_align_paths;

    for (auto & sub_path_start_idx: alignment.start()) {

        auto cur_extended_align_paths = extendAlignmentPath(align_path, alignment, sub_path_start_idx, stop_nodes);
        extended_align_paths.insert(extended_align_paths.end(), cur_extended_align_paths.begin(), cur_extended_align_paths.end());
    }

    return extended_align_paths;
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment, const int32_t sub_path_start_idx, const vector<gbwt::node_type> & stop_nodes) const {

    vector<AlignmentPath> extended_align_path(1, align_path);
    extended_align_path.front().mapqs.emplace_back(alignment.mapping_quality());
    extended_align_path.front().scores.emplace_back(0);

    extendAlignmentPaths(&extended_align_path, alignment.subpath(), sub_path_start_idx, stop_nodes);
            
    return extended_align_path;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPaths(vector<AlignmentPath> * align_paths, const google::protobuf::RepeatedPtrField<vg::Subpath> & subpaths, const int32_t sub_path_start_idx, const vector<gbwt::node_type> & stop_nodes) const {

    std::queue<pair<AlignmentPath, int32_t> > align_paths_queue;

    for (auto & align_path: *align_paths) {

        align_paths_queue.push(make_pair(align_path, sub_path_start_idx));
    }

    align_paths->clear();

    // Perform depth-first alignment path extension.
    while (!align_paths_queue.empty()) {

        auto & cur_align_path = align_paths_queue.front();

        const vg::Subpath & subpath = subpaths[cur_align_path.second];

        cur_align_path.first.scores.back() += subpath.score();
        extendAlignmentPath(&cur_align_path.first, subpath.path(), stop_nodes);

        if (find(stop_nodes.begin(), stop_nodes.end(), cur_align_path.first.search.node) == stop_nodes.end() && subpath.next_size() > 0) {

            for (auto & next_subpath_idx: subpath.next()) {

                align_paths_queue.push(make_pair(cur_align_path.first, next_subpath_idx));
            }

        } else if (!cur_align_path.first.search.empty()) {

            align_paths->emplace_back(cur_align_path.first);
        }

        align_paths_queue.pop();
    }
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::findPairedAlignmentPaths(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const {

    function<size_t(const int64_t)> node_seq_length_func = [&](const int64_t node_id) { return node_seq_lengths.at(node_id); };

    vector<AlignmentPath> paired_align_paths;

    AlignmentType alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_seq_length_func);
    
    auto alignment_2_rc_start_nodes = getAlignmentStartNodes(alignment_2_rc);
    assert(!alignment_2_rc_start_nodes.empty());
    
    for (auto & align_path: extendAlignmentPath(AlignmentPath(), alignment_1, alignment_2_rc_start_nodes)) {

        pairAlignmentPaths(&paired_align_paths, align_path, alignment_2_rc);
    }

    AlignmentType alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_seq_length_func);

    auto alignment_1_rc_start_nodes = getAlignmentStartNodes(alignment_1_rc);
    assert(!alignment_1_rc_start_nodes.empty());

    for (auto & align_path: extendAlignmentPath(AlignmentPath(), alignment_2, alignment_1_rc_start_nodes)) {

        pairAlignmentPaths(&paired_align_paths, align_path, alignment_1_rc);
    }

    return paired_align_paths;
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::findPairedAlignmentPathsIds(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const {

#ifdef debug
        
    printDebug(alignment_1, alignment_2);

#endif

    auto paired_align_paths = findPairedAlignmentPaths(alignment_1, alignment_2);

    for (auto & align_path: paired_align_paths) {

        align_path.ids = paths_index.locate(align_path.search);
    }

    return paired_align_paths;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::pairAlignmentPaths(vector<AlignmentPath> * paired_align_paths, const AlignmentPath & start_align_path, const AlignmentType & end_alignment) const {

    assert(!start_align_path.search.empty());

    std::queue<AlignmentPath> paired_align_path_queue;
    paired_align_path_queue.push(start_align_path);

    paired_align_path_queue.front().seq_length += (node_seq_lengths.at(gbwt::Node::id(paired_align_path_queue.front().search.node)) - paired_align_path_queue.front().end_offset);
    paired_align_path_queue.front().end_offset = node_seq_lengths.at(gbwt::Node::id(paired_align_path_queue.front().search.node));

    auto end_alignment_start_nodes_index = getAlignmentStartNodesIndex(end_alignment);

    // Perform depth-first path extension.
    while (!paired_align_path_queue.empty()) {

        AlignmentPath * cur_paired_align_path = &(paired_align_path_queue.front());
        assert(cur_paired_align_path->search.node != gbwt::ENDMARKER);

        auto end_alignment_start_nodes_index_it = end_alignment_start_nodes_index.find(cur_paired_align_path->search.node);

        if (end_alignment_start_nodes_index_it != end_alignment_start_nodes_index.end()) {

            cur_paired_align_path->seq_length -= cur_paired_align_path->end_offset;

            auto complete_paired_align_paths = extendAlignmentPath(*cur_paired_align_path, end_alignment, end_alignment_start_nodes_index_it->second);

            for (auto & complete_align_path: complete_paired_align_paths) {

                if (!complete_align_path.search.empty() && complete_align_path.seq_length <= max_pair_seq_length) {

                    paired_align_paths->emplace_back(complete_align_path);                         
                }
            }

            paired_align_path_queue.pop();
            continue;
        }
       
        if (cur_paired_align_path->seq_length > max_pair_seq_length) {

            paired_align_path_queue.pop();
            continue;
        }

        auto out_edges = paths_index.edges(cur_paired_align_path->search.node);

        // End current extension if no outgoing edges exist.
        if (out_edges.empty()) {

            paired_align_path_queue.pop();
            continue;
        }

        auto out_edges_it = out_edges.begin(); 

        while (out_edges_it != out_edges.end()) {

            if (out_edges_it->first != gbwt::ENDMARKER) {

                auto extended_path = paths_index.extend(cur_paired_align_path->search, out_edges_it->first);

                // Add new extension to queue if not empty (path found).
                if (!extended_path.empty()) { 

                    paired_align_path_queue.push(*cur_paired_align_path);
                    paired_align_path_queue.back().search = extended_path;
                    paired_align_path_queue.back().end_offset = node_seq_lengths.at(gbwt::Node::id(extended_path.node));
                    paired_align_path_queue.back().seq_length += paired_align_path_queue.back().end_offset;
                }
            }

            ++out_edges_it;
        }

        paired_align_path_queue.pop();
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

    for (size_t i = 0; i < alignment.subpath_size(); ++i) {

        assert(alignment.subpath(i).path().mapping_size() > 0);
        alignment_start_nodes.emplace_back(mapping_to_gbwt(alignment.subpath(i).path().mapping(0)));
    }

    return alignment_start_nodes;
}

template<class AlignmentType>
map<gbwt::node_type, int32_t> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodesIndex(const vg::Alignment & alignment) const {

    map<gbwt::node_type, int32_t> alignment_start_nodes_index;

    assert(alignment.path().mapping_size() > 0);
    assert(alignment_start_nodes_index.emplace(mapping_to_gbwt(alignment.path().mapping(0)), 0).second);

    return alignment_start_nodes_index;
}

template<class AlignmentType>
map<gbwt::node_type, int32_t> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodesIndex(const vg::MultipathAlignment & alignment) const {

    map<gbwt::node_type, int32_t> alignment_start_nodes_index;

    for (size_t i = 0; i < alignment.subpath_size(); ++i) {

        assert(alignment.subpath(i).path().mapping_size() > 0);
        assert(alignment_start_nodes_index.emplace(mapping_to_gbwt(alignment.subpath(i).path().mapping(0)), i).second);
    }

    return alignment_start_nodes_index;
}


template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::printDebug(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const {

    function<size_t(const int64_t)> node_seq_length_func = [&](const int64_t node_id) { return node_seq_lengths.at(node_id); };

    auto alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_seq_length_func);
    auto alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_seq_length_func);

    cout << pb2json(alignment_1) << endl;
    cout << pb2json(alignment_2_rc) << endl;
    cout << pb2json(alignment_2) << endl;
    cout << pb2json(alignment_1_rc) << endl;
    cout << endl;

    cout << findAlignmentPathsIds(alignment_1) << endl;
    cout << findAlignmentPathsIds(alignment_2_rc) << endl;
    cout << findAlignmentPathsIds(alignment_2) << endl;
    cout << findAlignmentPathsIds(alignment_1_rc) << endl;

    auto paired_align_paths = findPairedAlignmentPaths(alignment_1, alignment_2);

    for (auto & align_path: paired_align_paths) {

        align_path.ids = paths_index.locate(align_path.search);
    }

    cout << paired_align_paths << endl;
    cout << endl;
}

template class AlignmentPathFinder<vg::Alignment>;
template class AlignmentPathFinder<vg::MultipathAlignment>;

