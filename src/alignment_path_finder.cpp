
#include "alignment_path_finder.hpp"

#include <assert.h>

#include "utils.hpp"

//#define debug


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

#ifdef debug

    cerr << endl;
    cerr << pb2json(alignment) << endl;

#endif

    auto align_paths = extendAlignmentPath(AlignmentPath(), alignment);

    if (!paths_index.bidirectional()) {

        function<size_t(const int64_t)> node_seq_length_func = [&](const int64_t node_id) { return node_seq_lengths.at(node_id); };
        AlignmentType alignment_rc = lazy_reverse_complement_alignment(alignment, node_seq_length_func);

        auto align_paths_rc = extendAlignmentPath(AlignmentPath(), alignment_rc);

        align_paths.reserve(align_paths.size() + align_paths_rc.size());
        align_paths.insert(align_paths.end(), align_paths_rc.begin(), align_paths_rc.end());
    }       

    for (auto & align_path: align_paths) {

        align_path.ids = paths_index.locate(align_path.search);

        if (paths_index.bidirectional()) {

            for (auto & id: align_path.ids) {

                id = gbwt::Path::id(id);
            }
        }
    }

#ifdef debug

    cerr << endl;
    cerr << align_paths << endl;
    cerr << endl;

#endif

    return align_paths;
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment) const {

    return extendAlignmentPath(align_path, alignment, 0);
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment, const int32_t subpath_start_idx) const {

    vector<AlignmentPath> extended_align_path(1, align_path);
    extended_align_path.front().mapqs.emplace_back(alignment.mapping_quality());
    extended_align_path.front().scores.emplace_back(alignment.score());
    
    extendAlignmentPath(&extended_align_path.front(), alignment.path());

    if (extended_align_path.front().search.empty()) {

        return vector<AlignmentPath>();
    
    } else {

        return extended_align_path;
    }
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPath(AlignmentPath * align_path, const vg::Path & path) const {

    assert(align_path->path_end_pos <= align_path->path.size());
    
    auto mapping_it = path.mapping().cbegin();
    assert(mapping_it != path.mapping().cend());

    while (align_path->path_end_pos < align_path->path.size() && mapping_it != path.mapping().cend()) {

        ++align_path->path_end_pos;

        if (align_path->path.at(align_path->path_end_pos - 1) != mapping_to_gbwt(*mapping_it)) {

            align_path->search = gbwt::SearchState();
            assert(align_path->search.empty());  
        
        } else {

            align_path->seq_length -= align_path->seq_end_offset;
            align_path->seq_end_offset = mapping_it->position().offset() + mapping_from_length(*mapping_it);
            align_path->seq_length += mapping_it->position().offset() + mapping_to_length(*mapping_it);
        } 

        ++mapping_it;
    }

    while (mapping_it != path.mapping().cend()) {

        align_path->path.emplace_back(mapping_to_gbwt(*mapping_it));
        ++align_path->path_end_pos;
        align_path->seq_end_offset = mapping_it->position().offset() + mapping_from_length(*mapping_it);

        if (align_path->path.size() == 1) {

            assert(align_path->search.node == gbwt::ENDMARKER);
            assert(align_path->seq_length == 0);

            align_path->search = paths_index.find(align_path->path.back());
        
        } else {

            align_path->search = paths_index.extend(align_path->search, align_path->path.back());                
        }

        align_path->seq_length += mapping_to_length(*mapping_it);
        ++mapping_it;
    }
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment) const {

    vector<AlignmentPath> extended_align_paths;

    for (auto & start_idx: alignment.start()) {

        auto cur_extended_align_paths = extendAlignmentPath(align_path, alignment, start_idx);
        extended_align_paths.insert(extended_align_paths.end(), cur_extended_align_paths.begin(), cur_extended_align_paths.end());
    }

    return extended_align_paths;
}

template<class AlignmentType>
vector<AlignmentPath> AlignmentPathFinder<AlignmentType>::extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment, const int32_t subpath_start_idx) const {

    vector<AlignmentPath> extended_align_path(1, align_path);
    extended_align_path.front().mapqs.emplace_back(alignment.mapping_quality());
    extended_align_path.front().scores.emplace_back(0);

    extendAlignmentPaths(&extended_align_path, alignment.subpath(), subpath_start_idx);
            
    return extended_align_path;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentPaths(vector<AlignmentPath> * align_paths, const google::protobuf::RepeatedPtrField<vg::Subpath> & subpaths, const int32_t subpath_start_idx) const {

    std::queue<pair<AlignmentPath, int32_t> > align_paths_queue;

    for (auto & align_path: *align_paths) {

        align_paths_queue.push(make_pair(align_path, subpath_start_idx));
    }

    align_paths->clear();

    // Perform depth-first alignment path extension.
    while (!align_paths_queue.empty()) {

        auto & cur_align_path = align_paths_queue.front();

        const vg::Subpath & subpath = subpaths[cur_align_path.second];

        cur_align_path.first.scores.back() += subpath.score();
        extendAlignmentPath(&cur_align_path.first, subpath.path());

        if (subpath.next_size() > 0) {

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

#ifdef debug

    cerr << endl;
    findAlignmentPaths(alignment_1);
    findAlignmentPaths(alignment_2);

#endif

    vector<AlignmentPath> paired_align_paths;

    function<size_t(const int64_t)> node_seq_length_func = [&](const int64_t node_id) { return node_seq_lengths.at(node_id); };
    AlignmentType alignment_2_rc = lazy_reverse_complement_alignment(alignment_2, node_seq_length_func);

    pairAlignmentPaths(&paired_align_paths, alignment_1, alignment_2_rc);

    if (!paths_index.bidirectional()) {

        AlignmentType alignment_1_rc = lazy_reverse_complement_alignment(alignment_1, node_seq_length_func);
        pairAlignmentPaths(&paired_align_paths, alignment_2, alignment_1_rc);
    }

    for (auto & align_path: paired_align_paths) {

        align_path.ids = paths_index.locate(align_path.search);

        if (paths_index.bidirectional()) {

            for (auto & id: align_path.ids) {

                id = gbwt::Path::id(id);
            }
        }
    }

#ifdef debug

    cerr << endl;
    cerr << paired_align_paths << endl;
    cerr << endl;

#endif

    return paired_align_paths;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::pairAlignmentPaths(vector<AlignmentPath> * paired_align_paths, const AlignmentType & start_alignment, const AlignmentType & end_alignment) const {

    auto start_align_paths = extendAlignmentPath(AlignmentPath(), start_alignment);
    auto end_alignment_start_nodes_index = getAlignmentStartNodesIndex(end_alignment);

    std::queue<AlignmentPath> paired_align_path_queue;

    for (auto & align_path: start_align_paths) {

        assert(!align_path.search.empty());
        assert(!align_path.path.empty());

        align_path.seq_length += (node_seq_lengths.at(gbwt::Node::id(align_path.search.node)) - align_path.seq_end_offset);
        align_path.seq_end_offset = node_seq_lengths.at(gbwt::Node::id(align_path.search.node));

        bool found_overlap = false;

        for (auto & start_nodes: end_alignment_start_nodes_index) {

            auto path_it = find(align_path.path.begin(), align_path.path.end(), start_nodes.first); 

            while (path_it != align_path.path.end()) {

                align_path.path_end_pos = path_it - align_path.path.begin();
                auto complete_paired_align_paths = extendAlignmentPath(align_path, end_alignment, start_nodes.second);

                for (auto & complete_align_path: complete_paired_align_paths) {

                    if (!complete_align_path.search.empty() && complete_align_path.seq_length <= max_pair_seq_length) {

                        paired_align_paths->emplace_back(complete_align_path);                         
                    }
                }

                found_overlap = true;

                ++path_it;
                path_it = find(path_it, align_path.path.end(), start_nodes.first); 
            }
        }

        if (!found_overlap) {

            paired_align_path_queue.push(align_path);
        }
    }

    // Perform depth-first path extension.
    while (!paired_align_path_queue.empty()) {

        AlignmentPath * cur_paired_align_path = &(paired_align_path_queue.front());
        assert(cur_paired_align_path->search.node != gbwt::ENDMARKER);

        auto end_alignment_start_nodes_index_itp = end_alignment_start_nodes_index.equal_range(cur_paired_align_path->search.node);

        if (end_alignment_start_nodes_index_itp.first != end_alignment_start_nodes_index_itp.second) {

            while (end_alignment_start_nodes_index_itp.first != end_alignment_start_nodes_index_itp.second) {

                AlignmentPath cur_paired_align_path_end = *cur_paired_align_path;

                while (true) {

                    assert(cur_paired_align_path_end.path_end_pos == cur_paired_align_path_end.path.size());
                    --cur_paired_align_path_end.path_end_pos;

                    auto complete_paired_align_paths = extendAlignmentPath(cur_paired_align_path_end, end_alignment, end_alignment_start_nodes_index_itp.first->second);

                    for (auto & complete_align_path: complete_paired_align_paths) {

                        if (!complete_align_path.search.empty() && complete_align_path.seq_length <= max_pair_seq_length) {

                            paired_align_paths->emplace_back(complete_align_path);                         
                        }
                    }

                    cur_paired_align_path_end.search = paths_index.extend(cur_paired_align_path_end.search, cur_paired_align_path_end.search.node);

                    if (!cur_paired_align_path_end.search.empty()) { 

                        cur_paired_align_path_end.path.emplace_back(cur_paired_align_path_end.search.node);
                        cur_paired_align_path_end.path_end_pos = cur_paired_align_path_end.path.size();
                        cur_paired_align_path_end.seq_length += cur_paired_align_path_end.seq_end_offset;
                    
                    } else {

                        break;
                    }
                }

                ++end_alignment_start_nodes_index_itp.first;
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
                    paired_align_path_queue.back().path.emplace_back(extended_path.node);
                    ++paired_align_path_queue.back().path_end_pos;
                    paired_align_path_queue.back().seq_end_offset = node_seq_lengths.at(gbwt::Node::id(extended_path.node));
                    paired_align_path_queue.back().search = extended_path;
                    paired_align_path_queue.back().seq_length += paired_align_path_queue.back().seq_end_offset;
                }
            }

            ++out_edges_it;
        }

        paired_align_path_queue.pop();
    }
}

template<class AlignmentType>
multimap<gbwt::node_type, int32_t> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodesIndex(const vg::Alignment & alignment) const {

    multimap<gbwt::node_type, int32_t> alignment_start_nodes_index;

    assert(alignment.path().mapping_size() > 0);
    alignment_start_nodes_index.emplace(mapping_to_gbwt(alignment.path().mapping(0)), 0);

    return alignment_start_nodes_index;
}

template<class AlignmentType>
multimap<gbwt::node_type, int32_t> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodesIndex(const vg::MultipathAlignment & alignment) const {

    multimap<gbwt::node_type, int32_t> alignment_start_nodes_index;

    for (auto & start_idx: alignment.start()) {

        assert(alignment.subpath(start_idx).path().mapping_size() > 0);
        alignment_start_nodes_index.emplace(mapping_to_gbwt(alignment.subpath(start_idx).path().mapping(0)), start_idx);
    }

    return alignment_start_nodes_index;
}

template class AlignmentPathFinder<vg::Alignment>;
template class AlignmentPathFinder<vg::MultipathAlignment>;

