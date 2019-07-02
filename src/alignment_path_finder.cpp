
#include <assert.h>

#include "alignment_path_finder.hpp"
#include "utils.hpp"

// #define debug


AlignmentPathFinder::AlignmentPathFinder(const vg::Graph & graph, const gbwt::GBWT & paths_index_in) : paths_index(paths_index_in) {

    node_seq_lengths = vector<uint32_t>(graph.node_size() + 1, 0);

    for (auto & node: graph.node()) {

        if (node.id() >= node_seq_lengths.size()) {

            node_seq_lengths.resize(node.id() + 1, 0);
        }

        assert(node_seq_lengths.at(node.id()) == 0);
        node_seq_lengths.at(node.id()) = node.sequence().size();
    }
}

AlignmentPath AlignmentPathFinder::findAlignmentPath(const vg::Alignment & alignment) const {
            
    AlignmentPath align_path;
    align_path.scores.first = alignment.score();
    align_path.mapqs.first = alignment.mapping_quality();

    extendAlignmentPath(&align_path, alignment.path(), 0);

    return align_path;
}

AlignmentPath AlignmentPathFinder::findAlignmentPathIds(const vg::Alignment & alignment) const {

    auto align_path = findAlignmentPath(alignment);
    align_path.path_ids = paths_index.locate(align_path.path);

    return align_path;
}

void AlignmentPathFinder::extendAlignmentPath(AlignmentPath * align_path, const vg::Path & extend_path, const uint32_t & node_offset) const {
    
    assert(node_offset < extend_path.mapping().size());
    auto mapping_it = extend_path.mapping().cbegin() + node_offset;

    if (align_path->node_length == 0) {

        assert(align_path->seq_length == 0);

        align_path->path = paths_index.find(mapping_to_gbwt(*mapping_it));
        align_path->node_length++;
        align_path->seq_length = mapping_to_length(*mapping_it);
        ++mapping_it;
    } 

    while (mapping_it != extend_path.mapping().cend()) {

        align_path->path = paths_index.extend(align_path->path, mapping_to_gbwt(*mapping_it));
        align_path->node_length++;
        align_path->seq_length += mapping_to_length(*mapping_it);
        ++mapping_it;
    }
}

vector<AlignmentPath> AlignmentPathFinder::findPairedAlignmentPaths(const vg::Alignment & alignment_1, const vg::Alignment & alignment_2, const int32_t max_pair_distance) const {

    function<size_t(const int64_t)> node_seq_length_func = [&](const int64_t node_id) { return node_seq_lengths.at(node_id); };

    vector<AlignmentPath> paired_align_paths;

    auto align_path_1 = findAlignmentPath(alignment_1);

    if (!align_path_1.path.empty()) {

        vg::Alignment alignment_2_rc = alignment_2; 
        *alignment_2_rc.mutable_path() = lazy_reverse_complement_path(alignment_2.path(), node_seq_length_func);
        pairAlignmentPaths(&paired_align_paths, align_path_1, alignment_2_rc, max_pair_distance);
    }

    auto align_path_2 = findAlignmentPath(alignment_2);

    if (!align_path_2.path.empty()) {

        vg::Alignment alignment_1_rc = alignment_1; 
        *alignment_1_rc.mutable_path() = lazy_reverse_complement_path(alignment_1.path(), node_seq_length_func);
        pairAlignmentPaths(&paired_align_paths, align_path_2, alignment_1_rc, max_pair_distance);
    }

#ifdef debug
        cout << endl;
        cout << pb2json(alignment_1) << endl;
        cout << pb2json(alignment_2) << endl;

        vg::Alignment alignment_1_rc = alignment_1; 
        *alignment_1_rc.mutable_path() = lazy_reverse_complement_path(alignment_1.path(), node_seq_length_func);

        vg::Alignment alignment_2_rc = alignment_2; 
        *alignment_2_rc.mutable_path() = lazy_reverse_complement_path(alignment_2.path(), node_seq_length_func);

        cout << pb2json(alignment_1_rc) << endl;
        cout << pb2json(alignment_2_rc) << endl;

        cout << findAlignmentPathIds(alignment_1) << endl;
        cout << findAlignmentPathIds(alignment_1_rc) << endl;
        cout << findAlignmentPathIds(alignment_2) << endl;
        cout << findAlignmentPathIds(alignment_2_rc) << endl;       

        cout << paired_align_paths << endl;
#endif 

    return paired_align_paths;
}

vector<AlignmentPath> AlignmentPathFinder::findPairedAlignmentPathsIds(const vg::Alignment & alignment_1, const vg::Alignment & alignment_2, const int32_t max_pair_distance) const {

    auto paired_align_paths = findPairedAlignmentPaths(alignment_1, alignment_2, max_pair_distance);

    for (auto & align_path: paired_align_paths) {

        align_path.path_ids = paths_index.locate(align_path.path);
    }

    return paired_align_paths;
}

void AlignmentPathFinder::pairAlignmentPaths(vector<AlignmentPath> * paired_align_paths, const AlignmentPath & start_align_path, const vg::Alignment & end_alignment, const int32_t max_pair_distance) const {

    assert(!start_align_path.path.empty());

    std::queue<AlignmentPath> paired_align_path_queue;
    paired_align_path_queue.push(start_align_path);

    assert(end_alignment.path().mapping_size() > 0);

    auto end_alignment_node_index = getAlignmentNodeIndex(end_alignment);

    // Perform depth-first path extension.
    while (!paired_align_path_queue.empty()) {

        auto & cur_paired_align_path = paired_align_path_queue.front();

        if (cur_paired_align_path.seq_length > max_pair_distance) {

            paired_align_path_queue.pop();
            continue;                
        }

        auto end_alignment_node_index_it = end_alignment_node_index.equal_range(cur_paired_align_path.path.node);

        // Stop current extension if end node is reached.
        if (end_alignment_node_index_it.first != end_alignment_node_index_it.second) {

            while (end_alignment_node_index_it.first != end_alignment_node_index_it.second) {

                auto start_mapping_idx = end_alignment_node_index_it.first->second.second + 1;

                for (size_t i = start_mapping_idx; i < end_alignment.path().mapping_size(); ++i) {

                    auto gbwt_node_id = mapping_to_gbwt(end_alignment.path().mapping(i));

                    cur_paired_align_path.path = paths_index.extend(cur_paired_align_path.path, gbwt_node_id);
                    cur_paired_align_path.node_length++;
                    cur_paired_align_path.seq_length += node_seq_lengths.at(gbwt::Node::id(gbwt_node_id));
                }

                cur_paired_align_path.seq_length -= (node_seq_lengths.at(end_alignment.path().mapping().rbegin()->position().node_id()) - mapping_to_length(*end_alignment.path().mapping().rbegin()));

                if (!cur_paired_align_path.path.empty() && cur_paired_align_path.seq_length <= max_pair_distance) {

                    paired_align_paths->emplace_back(cur_paired_align_path);

                    paired_align_paths->back().path_ids = paths_index.locate(paired_align_paths->back().path);                                
                    paired_align_paths->back().scores.second = end_alignment.score();
                    paired_align_paths->back().mapqs.second = end_alignment.mapping_quality();           
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

multimap<gbwt::node_type, pair<int32_t, int32_t> > AlignmentPathFinder::getAlignmentNodeIndex(const vg::Alignment & alignment) const {

    multimap<gbwt::node_type, pair<int32_t, int32_t> > alignment_node_index;

    for (size_t i = 0; i < alignment.path().mapping_size(); ++i) {

        alignment_node_index.emplace(mapping_to_gbwt(alignment.path().mapping(i)), make_pair(0, i));
    }

    return alignment_node_index;
}


// vector<AlignmentPath> find_align_paths(const vg::MultipathAlignment & mp_alignment, const gbwt::GBWT & paths_index) {

//     vector<AlignmentPath> align_paths;

//     std::queue<pair<AlignmentPath, int32_t> > mp_align_path_queue;
//     assert(mp_alignment.start_size() > 0); 

//     for (auto & start_subpath_idx: mp_alignment.start()) {

//         mp_align_path_queue.push(make_pair(AlignmentPath(), start_subpath_idx));
//     }

//     // Perform depth-first alignment path extension.
//     while (!mp_align_path_queue.empty()) {

//         auto & cur_mp_align_path = mp_align_path_queue.front();

//         const vg::Subpath & subpath = mp_alignment.subpath(cur_mp_align_path.second);

//         cur_mp_align_path.first.scores.first += subpath.score();
//         cur_mp_align_path.first.extendPath(subpath.path(), 0, paths_index);

//         if (subpath.next_size() > 0) {

//             for (auto & next_subpath_idx: subpath.next()) {

//                 mp_align_path_queue.push(make_pair(cur_mp_align_path.first, next_subpath_idx));
//             }

//         } else {

//             align_paths.emplace_back(cur_mp_align_path.first);
//         }

//         mp_align_path_queue.pop();
//     }

//     return align_paths;
// }








