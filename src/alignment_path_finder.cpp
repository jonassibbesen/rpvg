
#include <assert.h>

#include "alignment_path_finder.hpp"
#include "utils.hpp"

//#define debug


template<class AlignmentType>
AlignmentPathFinder<AlignmentType>::AlignmentPathFinder(const vg::Graph & graph, const gbwt::GBWT & paths_index_in, const uint32_t & max_pair_distance_in) : paths_index(paths_index_in), max_pair_distance(max_pair_distance_in) {

    node_seq_lengths = vector<uint32_t>(graph.node_size() + 1, 0);

    for (auto & node: graph.node()) {

        if (node.id() >= node_seq_lengths.size()) {

            node_seq_lengths.resize(node.id() + 1, 0);
        }

        assert(node_seq_lengths.at(node.id()) == 0);
        node_seq_lengths.at(node.id()) = node.sequence().size();
    }
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

        align_path.path_ids = paths_index.locate(align_path.path);
    }

    return align_paths;
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment) const {

    return extendAlignmentPath(align_path, alignment, make_pair(0, 0));
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment, const pair<uint32_t, uint32_t> offset) const {

    vector<AlignmentPath> extended_align_path(1, align_path);
    extended_align_path.front().mapqs.emplace_back(alignment.mapping_quality());
    extended_align_path.front().scores.emplace_back(alignment.score());

    assert(offset.first == 0);
    extendAlignmentPath(&extended_align_path.front(), alignment.path(), offset.second);

    if (extended_align_path.front().path.empty()) {

        return vector<AlignmentPath>();
    
    } else {

        return extended_align_path;
    }
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPath(AlignmentPath * align_path, const vg::Path & path, const uint32_t offset) const {
    
    if (offset < path.mapping().size()) {
    
        auto mapping_it = path.mapping().cbegin() + offset;
        assert(mapping_it != path.mapping().cend());

        if (align_path->node_length == 0) {

            assert(align_path->seq_length == 0);

            align_path->path = paths_index.find(mapping_to_gbwt(*mapping_it));
            align_path->node_length++;
            align_path->seq_length = mapping_to_length(*mapping_it);
            ++mapping_it;
        
        } else if (align_path->path.node == mapping_to_gbwt(*mapping_it) && offset == 0) {

            align_path->seq_length += mapping_to_length(*mapping_it);
            ++mapping_it;            
        } 

        while (mapping_it != path.mapping().cend()) {

            align_path->path = paths_index.extend(align_path->path, mapping_to_gbwt(*mapping_it));
            align_path->node_length++;
            align_path->seq_length += mapping_to_length(*mapping_it);
            ++mapping_it;
        }
    }
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment) const {

    vector<AlignmentPath> extended_align_paths;

    for (auto & start_idx: alignment.start()) {

        auto cur_extended_align_paths = extendAlignmentPath(align_path, alignment, make_pair(start_idx, 0));
        extended_align_paths.insert(extended_align_paths.end(), cur_extended_align_paths.begin(), cur_extended_align_paths.end());
    }

    return extended_align_paths;
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment, const pair<uint32_t, uint32_t> offset) const {

    vector<AlignmentPath> extended_align_path(1, align_path);
    extended_align_path.front().mapqs.emplace_back(alignment.mapping_quality());
    extended_align_path.front().scores.emplace_back(0);

    assert(offset.first < alignment.subpath_size());
    extendAlignmentPaths(&extended_align_path, alignment.subpath(), offset);
            
    return extended_align_path;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPaths(vector<AlignmentPath> * align_paths, const google::protobuf::RepeatedPtrField<vg::Subpath> & sub_path, const pair<uint32_t, uint32_t> offset) const {

    std::queue<pair<AlignmentPath, pair<int32_t, int32_t> > > align_paths_queue;

    for (auto & align_path: *align_paths) {

        align_paths_queue.push(make_pair(align_path, offset));
    }

    align_paths->clear();

    // Perform depth-first alignment path extension.
    while (!align_paths_queue.empty()) {

        auto & cur_align_path = align_paths_queue.front();

        const vg::Subpath & subpath = sub_path[cur_align_path.second.first];

        cur_align_path.first.scores.back() += subpath.score();
        extendAlignmentPath(&cur_align_path.first, subpath.path(), cur_align_path.second.second);

        if (subpath.next_size() > 0) {

            for (auto & next_subpath_idx: subpath.next()) {

                align_paths_queue.push(make_pair(cur_align_path.first, make_pair(next_subpath_idx, 0)));
            }

        } else if (!cur_align_path.first.path.empty()) {

            align_paths->emplace_back(cur_align_path.first);
        }

        align_paths_queue.pop();
    }
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::findPairedAlignmentPaths(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const {

    function<size_t(const int64_t)> node_seq_length_func = [&](const int64_t node_id) { return node_seq_lengths.at(node_id); };

    vector<AlignmentPath> paired_align_paths;

    auto align_paths_1 = findAlignmentPaths(alignment_1);

    if (!align_paths_1.empty()) {

        AlignmentType alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_seq_length_func);
        
        for (auto & align_path: align_paths_1) {

            pairAlignmentPaths(&paired_align_paths, align_path, alignment_2_rc);
        }
    }

    auto align_paths_2 = findAlignmentPaths(alignment_2);

    if (!align_paths_2.empty()) {

        AlignmentType alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_seq_length_func);

        for (auto & align_path: align_paths_2) {

            pairAlignmentPaths(&paired_align_paths, align_path, alignment_1_rc);
        }
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

        align_path.path_ids = paths_index.locate(align_path.path);
    }

    return paired_align_paths;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::pairAlignmentPaths(vector<AlignmentPath> * paired_align_paths, const AlignmentPath & start_align_path, const AlignmentType & end_alignment) const {

    assert(!start_align_path.path.empty());

    std::queue<AlignmentPath> paired_align_path_queue;
    paired_align_path_queue.push(start_align_path);

    auto end_alignment_node_index = getAlignmentNodeIndex(end_alignment);

    // Perform depth-first path extension.
    while (!paired_align_path_queue.empty()) {

        auto & cur_paired_align_path = paired_align_path_queue.front();

        if (cur_paired_align_path.seq_length > max_pair_distance) {

            paired_align_path_queue.pop();
            continue;
        }

        assert(cur_paired_align_path.path.node != gbwt::ENDMARKER);
        auto end_alignment_node_index_it = end_alignment_node_index.equal_range(cur_paired_align_path.path.node);

        // Stop current extension if end node is reached.
        if (end_alignment_node_index_it.first != end_alignment_node_index_it.second) {

            while (end_alignment_node_index_it.first != end_alignment_node_index_it.second) {

                auto node_length = cur_paired_align_path.node_length;

                auto cur_offset = end_alignment_node_index_it.first->second;
                cur_offset.second++;

                auto ended_paired_align_paths = extendAlignmentPath(cur_paired_align_path, end_alignment, cur_offset);

                for (auto & ended_align_path: ended_paired_align_paths) {

                    if (!ended_align_path.path.empty() && ended_align_path.seq_length <= max_pair_distance) {

                        // Was not extended;
                        if (ended_align_path.node_length == node_length) {

                            auto end_mapping = getMapping(end_alignment, end_alignment_node_index_it.first->second);
                            assert(mapping_to_gbwt(end_mapping) == ended_align_path.path.node);
                            
                            ended_align_path.seq_length -= (node_seq_lengths.at(gbwt::Node::id(ended_align_path.path.node)) - mapping_to_length(end_mapping));
                        }
                      
                        paired_align_paths->emplace_back(ended_align_path);                         
                    }
                }

                end_alignment_node_index_it.first++;
            }

            paired_align_path_queue.pop();
            continue;
        }
        
        auto out_edges = paths_index.edges(cur_paired_align_path.path.node);

        // End current extension if no outgoing edges exist.
        if (out_edges.empty()) {

            paired_align_path_queue.pop();
            continue;
        }

        auto out_edges_it = out_edges.begin(); 

        while (out_edges_it != out_edges.end()) {

            if (out_edges_it->first != gbwt::ENDMARKER) {

                auto extended_path = paths_index.extend(cur_paired_align_path.path, out_edges_it->first);

                // Add new extension to queue if not empty (path found).
                if (!extended_path.empty()) { 

                    AlignmentPath new_paired_align_path = cur_paired_align_path;
                    new_paired_align_path.path = extended_path;
                    new_paired_align_path.node_length++;
                    new_paired_align_path.seq_length += node_seq_lengths.at(gbwt::Node::id(out_edges_it->first));

                    paired_align_path_queue.push(new_paired_align_path);
                }
            }

            ++out_edges_it;
        }

        paired_align_path_queue.pop();
    }
}

template<class AlignmentType>
multimap<gbwt::node_type, pair<int32_t, int32_t> > AlignmentPathFinder<AlignmentType>::getAlignmentNodeIndex(const vg::Alignment & alignment) const {

    multimap<gbwt::node_type, pair<int32_t, int32_t> > alignment_node_index;

    for (size_t i = 0; i < alignment.path().mapping_size(); ++i) {

        alignment_node_index.emplace(mapping_to_gbwt(alignment.path().mapping(i)), make_pair(0, i));
    }

    return alignment_node_index;
}

template<class AlignmentType>
multimap<gbwt::node_type, pair<int32_t, int32_t> > AlignmentPathFinder<AlignmentType>::getAlignmentNodeIndex(const vg::MultipathAlignment & alignment) const {

    multimap<gbwt::node_type, pair<int32_t, int32_t> > alignment_node_index;

    for (size_t i = 0; i < alignment.subpath_size(); ++i) {

        for (size_t j = 0; j < alignment.subpath(i).path().mapping_size(); ++j) {

            alignment_node_index.emplace(mapping_to_gbwt(alignment.subpath(i).path().mapping(j)), make_pair(i, j));
        }
    }

    return alignment_node_index;
}

template<class AlignmentType>
vg::Mapping AlignmentPathFinder<AlignmentType>::getMapping(const vg::Alignment & alignment, const pair<uint32_t, uint32_t> offset) const {

    assert(offset.first == 0);
    assert(offset.second < alignment.path().mapping_size());

    return alignment.path().mapping(offset.second);
}

template<class AlignmentType>
vg::Mapping AlignmentPathFinder<AlignmentType>::getMapping(const vg::MultipathAlignment & alignment, const pair<uint32_t, uint32_t> offset) const {

    assert(offset.first < alignment.subpath_size());
    assert(offset.second < alignment.subpath(offset.first).path().mapping_size());

    return alignment.subpath(offset.first).path().mapping(offset.second);
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

        align_path.path_ids = paths_index.locate(align_path.path);
    }

    cout << paired_align_paths << endl;
    cout << endl;
}

template class AlignmentPathFinder<vg::Alignment>;
template class AlignmentPathFinder<vg::MultipathAlignment>;

