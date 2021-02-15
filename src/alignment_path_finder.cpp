
#include "alignment_path_finder.hpp"

#include <assert.h>
#include <stack>

#include "utils.hpp"

//#define debug

static const int32_t max_score_diff = (Utils::default_match + Utils::default_mismatch) * 4;
static const int32_t max_noise_score_diff = (Utils::default_match + Utils::default_mismatch) * 2;


template<class AlignmentType>
AlignmentPathFinder<AlignmentType>::AlignmentPathFinder(const PathsIndex & paths_index_in, const string library_type_in, const uint32_t max_pair_frag_length_in, const uint32_t max_partial_offset_in, const bool est_missing_noise_prob_in, const double min_best_score_filter_in) : paths_index(paths_index_in), library_type(library_type_in), max_pair_frag_length(max_pair_frag_length_in), max_partial_offset(max_partial_offset_in), est_missing_noise_prob(est_missing_noise_prob_in), min_best_score_filter(min_best_score_filter_in) {}
        
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
int32_t AlignmentPathFinder<AlignmentType>::alignmentScore(const char & quality) const {

    return Utils::qual_score_matrix.at(25 * quality);
}

template<class AlignmentType>
int32_t AlignmentPathFinder<AlignmentType>::alignmentScore(const string & quality, const uint32_t & start_offset, const uint32_t & length) const {

    if (quality.empty()) {

        return length;
    }

    assert(start_offset + length <= quality.size());
    int32_t optimal_score = 0;

    for (size_t i = start_offset; i < start_offset + length; ++i) {

        optimal_score += alignmentScore(quality.at(i));
    } 

    return optimal_score;
}

template<class AlignmentType>
int32_t AlignmentPathFinder<AlignmentType>::optimalAlignmentScore(const string & quality, const uint32_t seq_length) const {

    if (quality.empty()) {

        return seq_length * Utils::default_match + 2 * Utils::default_full_length_bonus;
    } 

    assert(quality.size() == seq_length);

    int32_t optimal_score = alignmentScore(quality, 0, seq_length);
    optimal_score += Utils::qual_full_length_bonuses.at(quality.front()) + Utils::qual_full_length_bonuses.at(quality.back());

    return optimal_score;
}

template<class AlignmentType>
int32_t AlignmentPathFinder<AlignmentType>::optimalAlignmentScore(const vg::Alignment & alignment) const {

    return optimalAlignmentScore(alignment.quality(), alignment.sequence().size());
}

template<class AlignmentType>
int32_t AlignmentPathFinder<AlignmentType>::optimalAlignmentScore(const vg::MultipathAlignment & alignment) const {

    return optimalAlignmentScore(alignment.quality(), alignment.sequence().size());
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

        findAlignmentSearchPaths(&align_search_paths, alignment);

    } else if (library_type == "rf") {

        AlignmentType alignment_rc = Utils::lazy_reverse_complement_alignment(alignment, node_length_func);
        findAlignmentSearchPaths(&align_search_paths, alignment_rc);

    } else {

        assert(library_type == "unstranded");
        findAlignmentSearchPaths(&align_search_paths, alignment);

        if (!paths_index.bidirectional()) {

            AlignmentType alignment_rc = Utils::lazy_reverse_complement_alignment(alignment, node_length_func);
            findAlignmentSearchPaths(&align_search_paths, alignment_rc);
        }  
    }

    auto align_paths = AlignmentPath::alignmentSearchPathsToAlignmentPaths(align_search_paths, isAlignmentDisconnected(alignment));

#ifdef debug

    cerr << endl;
    cerr << align_search_paths << endl;
    cerr << align_paths << endl;
    cerr << endl;

#endif

    return align_paths;
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentSearchPath(const AlignmentSearchPath & align_search_path, const vg::Alignment & alignment) const {

    assert(alignment.mapping_quality() >= 0);
    auto optimal_score = optimalAlignmentScore(alignment);

    vector<AlignmentSearchPath> extended_align_search_paths(1, align_search_path);
    
    extended_align_search_paths.back().read_align_stats.emplace_back(AlignmentStats());
    AlignmentStats * read_align_stats = &(extended_align_search_paths.back().read_align_stats.back());

    read_align_stats->mapq = alignment.mapping_quality();
    read_align_stats->score = alignment.score();

    read_align_stats->internal_start.max_offset = min(read_align_stats->left_softclip_length + max_partial_offset, static_cast<uint32_t>(alignment.sequence().size()));
    read_align_stats->internal_end.max_offset = min(read_align_stats->right_softclip_length + max_partial_offset, static_cast<uint32_t>(alignment.sequence().size()));

    extendAlignmentSearchPath(&extended_align_search_paths, alignment.path(), true, true, alignment.quality(), alignment.sequence().size(), true);

    int32_t max_align_path_score = 0;

    for (auto & align_search_path: extended_align_search_paths) {

        assert(align_search_path.read_align_stats.back().length <= alignment.sequence().size());
        assert(!align_search_path.read_align_stats.back().complete);

        if ((align_search_path.isInternal() || !est_missing_noise_prob) && align_search_path.gbwt_search.first.empty()) {

            continue;
        }

        if (align_search_path.read_align_stats.back().length == alignment.sequence().size()) {

            align_search_path.read_align_stats.back().complete = true;
            max_align_path_score = max(max_align_path_score, align_search_path.scoreSum());
        }
    }

    assert(max_align_path_score <= optimal_score);

    for (auto & align_search_path: extended_align_search_paths) {

        if (align_search_path.read_align_stats.back().complete) {

            const int32_t align_path_score = align_search_path.scoreSum();
            assert(align_path_score <= max_align_path_score);

            if (max_align_path_score - align_path_score > max_score_diff) {

                align_search_path.read_align_stats.back().complete = false;
            }
        }
    }

    if (filterAlignmentSearchPaths(extended_align_search_paths, vector<int32_t>({optimal_score}))) {

        extended_align_search_paths.emplace_back(AlignmentSearchPath());
        extended_align_search_paths.back().path.emplace_back(gbwt::ENDMARKER);

        extended_align_search_paths.back().read_align_stats.emplace_back(AlignmentStats());
        AlignmentStats * error_read_align_stats = &(extended_align_search_paths.back().read_align_stats.back());

        error_read_align_stats->mapq = alignment.mapping_quality();
        error_read_align_stats->score = numeric_limits<int32_t>::max();

        error_read_align_stats->length = alignment.sequence().size();
        error_read_align_stats->complete = true;
    }

    return extended_align_search_paths;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentSearchPath(vector<AlignmentSearchPath> * align_search_paths, const vg::Path & path, const bool is_first_path, const bool is_last_path, const string & quality, const uint32_t seq_length, const bool add_internal_start) const {

    assert(align_search_paths->size() == 1);
    assert(!align_search_paths->front().read_align_stats.empty());

    if (is_first_path) {

        align_search_paths->front().read_align_stats.back().updateLeftSoftclipLength(path);
    }

    if (is_last_path) {

        align_search_paths->front().read_align_stats.back().updateRightSoftclipLength(path);
    }

    uint32_t last_internal_start_idx = 0;

    auto mapping_it = path.mapping().cbegin();
    assert(mapping_it != path.mapping().cend());

    auto end_mapping_it = path.mapping().cend();
    --end_mapping_it;

    while (mapping_it != path.mapping().cend()) {

        auto cur_node = Utils::mapping_to_gbwt(*mapping_it);
        auto mapping_read_length = Utils::mapping_to_length(*mapping_it);

        const bool is_last_mapping = (is_last_path && mapping_it == end_mapping_it);

        AlignmentSearchPath main_align_search_path;

        if (max_partial_offset > 0 && !align_search_paths->front().gbwt_search.first.empty() && !align_search_paths->front().read_align_stats.back().internal_end.is_internal) {

            assert(align_search_paths->front().read_align_stats.back().internal_end.offset == 0);
            assert(align_search_paths->front().read_align_stats.back().length <= seq_length);

            if (seq_length - align_search_paths->front().read_align_stats.back().length <= align_search_paths->front().read_align_stats.back().internal_end.max_offset) {

                main_align_search_path = align_search_paths->front();
            }
        }

        for (auto & align_search_path: *align_search_paths) {

            if (align_search_path.read_align_stats.back().internal_end.is_internal) {

                assert(max_partial_offset > 0);

                AlignmentStats * internal_end_read_align_stats = &(align_search_path.read_align_stats.back());
                auto internal_end_new_offset = mapping_read_length;

                if (is_last_mapping) {

                    assert(internal_end_read_align_stats->right_softclip_length <= internal_end_new_offset);
                    internal_end_new_offset -= internal_end_read_align_stats->right_softclip_length;
                } 

                internal_end_read_align_stats->internal_end.offset += internal_end_new_offset;

                if (internal_end_read_align_stats->internal_end.offset <= max_partial_offset) {
                    
                    internal_end_read_align_stats->internal_end.penalty += alignmentScore(quality, internal_end_read_align_stats->length, internal_end_new_offset);
                
                } else {

                    align_search_path.clear();
                }
            
            } else {

                extendAlignmentSearchPath(&align_search_path, *mapping_it);
            }
        }

        if (max_partial_offset > 0 && !main_align_search_path.gbwt_search.first.empty()) {

            assert(main_align_search_path.gbwt_search.first.size() >= align_search_paths->front().gbwt_search.first.size());

            if (main_align_search_path.gbwt_search.first.size() > align_search_paths->front().gbwt_search.first.size()) {

                main_align_search_path.read_align_stats.back().internal_end.is_internal = true;
                main_align_search_path.read_align_stats.back().internal_end.offset = mapping_read_length;

                if (is_last_mapping) {

                    assert(main_align_search_path.read_align_stats.back().right_softclip_length <= main_align_search_path.read_align_stats.back().internal_end.offset);
                    main_align_search_path.read_align_stats.back().internal_end.offset -= main_align_search_path.read_align_stats.back().right_softclip_length;
                } 

                if (main_align_search_path.read_align_stats.back().internal_end.offset <= max_partial_offset) {

                    main_align_search_path.read_align_stats.back().internal_end_next_node = cur_node;
                    main_align_search_path.read_align_stats.back().internal_end.penalty = alignmentScore(quality, main_align_search_path.read_align_stats.back().length, main_align_search_path.read_align_stats.back().internal_end.offset);

                    align_search_paths->emplace_back(main_align_search_path);
                }
            }
        }

        if (max_partial_offset > 0 && add_internal_start && align_search_paths->at(last_internal_start_idx).path.size() > 1 && !align_search_paths->at(last_internal_start_idx).read_align_stats.back().internal_end.is_internal) {

            if (align_search_paths->at(last_internal_start_idx).read_align_stats.back().length <= align_search_paths->at(last_internal_start_idx).read_align_stats.back().internal_start.max_offset) {

                auto internal_start_read_align_stats = align_search_paths->at(last_internal_start_idx).read_align_stats.back();
                assert(internal_start_read_align_stats.left_softclip_length <= internal_start_read_align_stats.length);

                internal_start_read_align_stats.internal_start.is_internal = true;
                internal_start_read_align_stats.internal_start.offset = internal_start_read_align_stats.length - internal_start_read_align_stats.left_softclip_length;

                if (internal_start_read_align_stats.internal_start.offset <= max_partial_offset) {

                    AlignmentSearchPath new_start_align_search_path;
                    extendAlignmentSearchPath(&new_start_align_search_path, *mapping_it);

                    if (!new_start_align_search_path.gbwt_search.first.empty()) {

                        assert(new_start_align_search_path.gbwt_search.first.size() >= align_search_paths->at(last_internal_start_idx).gbwt_search.first.size());

                        if (new_start_align_search_path.gbwt_search.first.size() > align_search_paths->at(last_internal_start_idx).gbwt_search.first.size()) {

                            align_search_paths->emplace_back(new_start_align_search_path);
                            last_internal_start_idx = align_search_paths->size() - 1;

                            internal_start_read_align_stats.internal_start.penalty = alignmentScore(quality, internal_start_read_align_stats.left_softclip_length, internal_start_read_align_stats.internal_start.offset);
                            align_search_paths->back().read_align_stats = vector<AlignmentStats>(1, internal_start_read_align_stats);
                        }
                    }
                }
            }
        }

        for (auto & align_search_path: *align_search_paths) {

            align_search_path.read_align_stats.back().length += mapping_read_length;
        }

        ++mapping_it;
    }
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentSearchPath(AlignmentSearchPath * align_search_path, const vg::Mapping & mapping) const {

    auto cur_node = Utils::mapping_to_gbwt(mapping);

    if (align_search_path->path.empty()) {

        assert(align_search_path->gbwt_search.first.node == gbwt::ENDMARKER);

        align_search_path->path.emplace_back(cur_node);
        paths_index.find(&(align_search_path->gbwt_search), cur_node);
  
        align_search_path->start_offset = mapping.position().offset();

    } else {

        bool is_cycle_visit = false;

        if (align_search_path->path.back() == cur_node && mapping.position().offset() != align_search_path->end_offset) {

            assert(mapping.position().offset() == 0);
            is_cycle_visit = true;      
        }

        if (align_search_path->path.back() != cur_node || is_cycle_visit) {

            align_search_path->path.emplace_back(cur_node);

            if (!align_search_path->gbwt_search.first.empty()) {

                paths_index.extend(&(align_search_path->gbwt_search), cur_node);
            }
        } 
    }

    align_search_path->end_offset = mapping.position().offset() + Utils::mapping_from_length(mapping);
}

template<class AlignmentType>
vector<AlignmentSearchPath> AlignmentPathFinder<AlignmentType>::extendAlignmentSearchPath(const AlignmentSearchPath & align_search_path, const vg::MultipathAlignment & alignment) const {

    assert(alignment.mapping_quality() >= 0);
    auto optimal_score = optimalAlignmentScore(alignment);

    vector<AlignmentSearchPath> extended_align_search_paths;

    auto left_softclip_lengths = getAlignmentStartSoftclipLengths(alignment);
    auto right_softclip_lengths = getAlignmentEndSoftclipLengths(alignment);

    auto min_right_softclip_length = *min_element(right_softclip_lengths.begin(), right_softclip_lengths.end());
    assert(min_right_softclip_length <= alignment.sequence().size());

    auto max_right_softclip_length = *max_element(right_softclip_lengths.begin(), right_softclip_lengths.end());
    assert(max_right_softclip_length <= alignment.sequence().size());

    vector<pair<int32_t, uint32_t> > start_score_indexes;

    for (auto & start_subpath_idx: alignment.start()) {

        start_score_indexes.emplace_back(alignment.subpath().Get(start_subpath_idx).score(), start_subpath_idx);
    }

    sort(start_score_indexes.rbegin(), start_score_indexes.rend());

    spp::sparse_hash_map<pair<uint32_t, uint32_t>, int32_t> internal_node_subpaths;
    int32_t best_align_score = floor(optimal_score * min_best_score_filter);

    for (auto & start_score_idx: start_score_indexes) {

        AlignmentSearchPath init_align_search_path;

        init_align_search_path.read_align_stats.emplace_back(AlignmentStats());
        AlignmentStats * init_read_align_stats = &(init_align_search_path.read_align_stats.back());

        init_read_align_stats->mapq = alignment.mapping_quality();

        AlignmentStats tpm_read_align_stats;
        tpm_read_align_stats.updateLeftSoftclipLength(alignment.subpath(start_score_idx.second).path());
        assert(tpm_read_align_stats.left_softclip_length <= alignment.sequence().size());

        init_read_align_stats->internal_start.max_offset = min(tpm_read_align_stats.left_softclip_length + max_partial_offset, static_cast<uint32_t>(alignment.sequence().size()));
        init_read_align_stats->internal_end.max_offset = min(max_right_softclip_length + max_partial_offset, static_cast<uint32_t>(alignment.sequence().size()));

        extendAlignmentSearchPaths(&extended_align_search_paths, init_align_search_path, alignment.subpath(), start_score_idx.second, alignment.quality(), alignment.sequence().size(), &internal_node_subpaths, &best_align_score, min_right_softclip_length == 0);
    }

    assert(best_align_score <= optimal_score);

    for (auto & align_search_path: extended_align_search_paths) {

        assert(align_search_path.read_align_stats.back().complete);

        const int32_t align_path_score = align_search_path.scoreSum();
        assert(align_path_score <= best_align_score);

        if (best_align_score - align_path_score > max_score_diff) {

            align_search_path.read_align_stats.back().complete = false;
        }
    }

    if (filterAlignmentSearchPaths(extended_align_search_paths, vector<int32_t>({optimal_score}))) {

        extended_align_search_paths.emplace_back(AlignmentSearchPath());
        extended_align_search_paths.back().path.emplace_back(gbwt::ENDMARKER);

        extended_align_search_paths.back().read_align_stats.emplace_back(AlignmentStats());
        AlignmentStats * error_read_align_stats = &(extended_align_search_paths.back().read_align_stats.back());

        error_read_align_stats->mapq = alignment.mapping_quality();
        error_read_align_stats->score = numeric_limits<int32_t>::max();

        error_read_align_stats->length = alignment.sequence().size();
        error_read_align_stats->complete = true;
    }

    return extended_align_search_paths;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::extendAlignmentSearchPaths(vector<AlignmentSearchPath> * align_search_paths, const AlignmentSearchPath & init_align_search_path, const google::protobuf::RepeatedPtrField<vg::Subpath> & subpaths, const uint32_t start_subpath_idx, const string & quality, const uint32_t seq_length, spp::sparse_hash_map<pair<uint32_t, uint32_t>, int32_t> * internal_node_subpaths, int32_t * best_align_score, const bool has_right_bonus) const {

    stack<pair<AlignmentSearchPath, uint32_t> > align_search_paths_stack;
    align_search_paths_stack.emplace(init_align_search_path, start_subpath_idx);

    // Perform depth-first alignment path extension.
    while (!align_search_paths_stack.empty()) {

        vector<AlignmentSearchPath> extended_align_search_paths(1, align_search_paths_stack.top().first);
        const uint32_t subpath_idx = align_search_paths_stack.top().second;

        align_search_paths_stack.pop();

        const vg::Subpath & subpath = subpaths.Get(subpath_idx);
        AlignmentSearchPath * extended_align_search_path = &(extended_align_search_paths.front());
        extended_align_search_path->read_align_stats.back().score += subpath.score();

        uint32_t subpath_length = 0;

        for (auto & mapping: subpath.path().mapping()) {

            subpath_length += Utils::mapping_to_length(mapping);
        }

        assert(extended_align_search_path->read_align_stats.back().length + subpath_length <= seq_length);
        const int32_t seq_length_left = seq_length - (extended_align_search_path->read_align_stats.back().length + subpath_length);

        int32_t max_score = extended_align_search_path->read_align_stats.back().score + seq_length_left;

        if (has_right_bonus && subpath.next_size() > 0) {

            max_score += Utils::default_full_length_bonus;
        }

        if (*best_align_score - max_score > max_score_diff) {

            continue;
        } 

        bool add_internal_start = false;

        if (max_partial_offset > 0 && extended_align_search_path->read_align_stats.back().length <= extended_align_search_path->read_align_stats.back().internal_start.max_offset) {

            add_internal_start = true;

            assert(extended_align_search_path->read_align_stats.back().left_softclip_length <= extended_align_search_path->read_align_stats.back().length);
            auto internal_node_subpaths_it = internal_node_subpaths->emplace(make_pair(subpath_idx, extended_align_search_path->read_align_stats.back().length - extended_align_search_path->read_align_stats.back().left_softclip_length), extended_align_search_path->read_align_stats.back().score);

            if (!internal_node_subpaths_it.second) {

                if (extended_align_search_path->read_align_stats.back().score <= internal_node_subpaths_it.first->second) {

                    add_internal_start = false;

                } else {

                    internal_node_subpaths_it.first->second = extended_align_search_path->read_align_stats.back().score;
                }
            }
        
        } else if (extended_align_search_path->gbwt_search.first.empty()) {

            if (*best_align_score - max_score > max_noise_score_diff) {

                continue;
            } 
        }

        extendAlignmentSearchPath(&extended_align_search_paths, subpath.path(), subpath_idx == start_subpath_idx, subpath.next_size() == 0, quality, seq_length, add_internal_start);        

        for (auto & align_search_path: extended_align_search_paths) {

            if (align_search_path.gbwt_search.first.empty()) {

                if (align_search_path.isInternal()) {

                    continue;
                }

                if (!est_missing_noise_prob && max_partial_offset == 0) {

                    continue;
                }

                if (!est_missing_noise_prob && align_search_path.read_align_stats.back().length > align_search_path.read_align_stats.back().internal_start.max_offset) {

                    continue;
                }
            }

            assert(!align_search_path.path.empty());

            if (subpath.next_size() > 0) {

                vector<pair<int32_t, uint32_t> > next_score_indexes;

                for (auto & next_subpath_idx: subpath.next()) {

                    next_score_indexes.emplace_back(subpaths.Get(next_subpath_idx).score(), next_subpath_idx);
                }

                sort(next_score_indexes.begin(), next_score_indexes.end());

                for (auto & next_score_idx: next_score_indexes) {

                    align_search_paths_stack.emplace(align_search_path, next_score_idx.second);
                }

            } else if (subpath.connection_size() == 0) {

                *best_align_score = max(*best_align_score, align_search_path.scoreSum());

                assert(align_search_path.read_align_stats.back().length == seq_length);
                assert(!align_search_path.read_align_stats.back().complete);

                align_search_path.read_align_stats.back().complete = true;
                align_search_paths->emplace_back(move(align_search_path));
            }
        }
    }
}

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
    AlignmentType alignment_2_rc = Utils::lazy_reverse_complement_alignment(alignment_2, node_length_func);

    if (library_type == "fr") {

        findPairedAlignmentSearchPaths(&paired_align_search_paths, alignment_1, alignment_2_rc);

    } else if (library_type == "rf") {

        AlignmentType alignment_1_rc = Utils::lazy_reverse_complement_alignment(alignment_1, node_length_func);
        findPairedAlignmentSearchPaths(&paired_align_search_paths, alignment_2, alignment_1_rc);

    } else {

        assert(library_type == "unstranded");

        AlignmentType alignment_2_rc = Utils::lazy_reverse_complement_alignment(alignment_2, node_length_func);
        findPairedAlignmentSearchPaths(&paired_align_search_paths, alignment_1, alignment_2_rc);

        if (!paths_index.bidirectional()) {

            AlignmentType alignment_1_rc = Utils::lazy_reverse_complement_alignment(alignment_1, node_length_func);
            findPairedAlignmentSearchPaths(&paired_align_search_paths, alignment_2, alignment_1_rc);
        }
    }

    auto paired_align_paths = AlignmentPath::alignmentSearchPathsToAlignmentPaths(paired_align_search_paths, isAlignmentDisconnected(alignment_1) || isAlignmentDisconnected(alignment_2));

#ifdef debug

    cerr << endl;
    cerr << paired_align_search_paths << endl;
    cerr << paired_align_paths << endl;
    cerr << endl;

#endif

    return paired_align_paths;
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::findAlignmentSearchPaths(vector<AlignmentSearchPath> * align_search_paths, const AlignmentType & alignment) const {

    auto single_align_search_paths = extendAlignmentSearchPath(AlignmentSearchPath(), alignment);

    if (single_align_search_paths.empty()) {

        return;
    }

    sort(single_align_search_paths.rbegin(), single_align_search_paths.rend());

    double joint_single_align_score = numeric_limits<int32_t>::lowest();
    double joint_empty_single_align_score = numeric_limits<int32_t>::lowest();

    for (size_t i = 0; i < single_align_search_paths.size(); ++i) {

        const AlignmentSearchPath & single_search_path = single_align_search_paths.at(i);
        assert(single_search_path.read_align_stats.size() == 1);

        if (!single_search_path.isComplete()) {

            continue;
        }

        assert(!single_search_path.path.empty());
        assert(single_search_path.read_align_stats.back().length == alignment.sequence().size());

        if (i > 0) {

            if (single_search_path.path == single_align_search_paths.at(i - 1).path) {

                assert(single_search_path.gbwt_search == single_align_search_paths.at(i - 1).gbwt_search);
                assert(single_search_path.scoreSum() <= single_align_search_paths.at(i - 1).scoreSum());

                continue;
            }
        }

        const int32_t score_sum = single_search_path.scoreSum();

        if (single_search_path.gbwt_search.first.empty()) {

            assert(!single_search_path.isInternal());
            joint_empty_single_align_score  = Utils::add_log(joint_empty_single_align_score, score_sum * Utils::score_log_base);
            
            continue;
        }

        if (!single_search_path.isInternal()) {

            joint_single_align_score = Utils::add_log(joint_single_align_score, score_sum * Utils::score_log_base);
        }

        align_search_paths->emplace_back(move(single_search_path));
    }

    align_search_paths->emplace_back(AlignmentSearchPath());

    align_search_paths->back().read_align_stats.emplace_back(AlignmentStats());
    align_search_paths->back().read_align_stats.back().score = Utils::doubleToInt((joint_single_align_score - joint_empty_single_align_score) / Utils::noise_score_log_base);
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::findPairedAlignmentSearchPaths(vector<AlignmentSearchPath> * paired_align_search_paths, const AlignmentType & start_alignment, const AlignmentType & end_alignment) const {

    auto start_align_search_paths = extendAlignmentSearchPath(AlignmentSearchPath(), start_alignment);
    auto end_align_search_paths = extendAlignmentSearchPath(AlignmentSearchPath(), end_alignment);

    if (start_align_search_paths.empty() || end_align_search_paths.empty()) {

        return;
    }

    sort(start_align_search_paths.rbegin(), start_align_search_paths.rend());
    sort(end_align_search_paths.rbegin(), end_align_search_paths.rend());

    uint32_t num_unique_end_search_paths = 0;
    uint32_t end_max_left_softclip_length = 0;

    spp::sparse_hash_map<gbwt::node_type, uint32_t> end_search_paths_nodes;
    spp::sparse_hash_map<gbwt::node_type, vector<uint32_t> > end_search_paths_start_nodes_index;

    double joint_end_align_score = numeric_limits<int32_t>::lowest();
    double joint_empty_end_align_score = numeric_limits<int32_t>::lowest();

    for (size_t i = 0; i < end_align_search_paths.size(); ++i) {

        const AlignmentSearchPath & end_search_path = end_align_search_paths.at(i);
        assert(end_search_path.read_align_stats.size() == 1);

        if (!end_search_path.isComplete()) {

            continue;
        }

        assert(!end_search_path.path.empty());
        assert(end_search_path.read_align_stats.back().length == end_alignment.sequence().size());

        if (i > 0) {

            if (end_search_path.path == end_align_search_paths.at(i - 1).path) {

                assert(end_search_path.gbwt_search.first == end_align_search_paths.at(i - 1).gbwt_search.first);
                assert(end_search_path.scoreSum() <= end_align_search_paths.at(i - 1).scoreSum());

                continue;
            }
        }

        const int32_t score_sum = end_search_path.scoreSum();

        if (end_search_path.gbwt_search.first.empty()) {

            assert(!end_search_path.isInternal());
            joint_empty_end_align_score  = Utils::add_log(joint_empty_end_align_score, score_sum * Utils::score_log_base);
            
            continue;
        }

        if (!end_search_path.isInternal()) {

            joint_end_align_score = Utils::add_log(joint_end_align_score, score_sum * Utils::score_log_base);
        }

        num_unique_end_search_paths++;
        end_max_left_softclip_length = max(end_max_left_softclip_length, end_search_path.read_align_stats.back().left_softclip_length);

        for (auto & path_id: end_search_path.path) {

            auto end_search_paths_nodes_it = end_search_paths_nodes.emplace(path_id, 0);
            end_search_paths_nodes_it.first->second++;
        }

        auto end_search_paths_start_nodes_index_it = end_search_paths_start_nodes_index.emplace(end_search_path.path.front(), vector<uint32_t>());
        end_search_paths_start_nodes_index_it.first->second.emplace_back(i);
    }

    assert(end_max_left_softclip_length <= end_alignment.sequence().size());

    bool end_alignment_in_cycle = false;

    for (auto end_search_paths_start_node: end_search_paths_start_nodes_index) {

        pair<gbwt::SearchState, gbwt::size_type> start_node_gbwt_search;
        paths_index.find(&start_node_gbwt_search, end_search_paths_start_node.first);

        const uint32_t num_start_node_paths = paths_index.locatePathIds(start_node_gbwt_search).size();
        assert(num_start_node_paths <= start_node_gbwt_search.first.size());

        if (num_start_node_paths < start_node_gbwt_search.first.size()) {

            end_alignment_in_cycle = true;
            break;
        }
    }

    stack<pair<AlignmentSearchPath, bool> > paired_align_search_path_stack;

    double joint_start_align_score = numeric_limits<int32_t>::lowest();
    double joint_empty_start_align_score = numeric_limits<int32_t>::lowest();

    for (size_t i = 0; i < start_align_search_paths.size(); ++i) {

        const AlignmentSearchPath & start_search_path = start_align_search_paths.at(i);
        assert(start_search_path.read_align_stats.size() == 1);

        if (!start_search_path.isComplete()) {

            continue;
        }

        assert(!start_search_path.path.empty());
        assert(start_search_path.read_align_stats.back().length == start_alignment.sequence().size());

        if (i > 0) {

            if (start_search_path.path == start_align_search_paths.at(i - 1).path) {

                assert(start_search_path.gbwt_search == start_align_search_paths.at(i - 1).gbwt_search);
                assert(start_search_path.scoreSum() <= start_align_search_paths.at(i - 1).scoreSum());

                continue;
            }
        }

        const int32_t score_sum = start_search_path.scoreSum();

        if (start_search_path.gbwt_search.first.empty()) {

            assert(!start_search_path.isInternal());
            joint_empty_start_align_score = Utils::add_log(joint_empty_start_align_score, score_sum * Utils::score_log_base);

            continue;
        } 
        
        if (!start_search_path.isInternal()) {

            joint_start_align_score = Utils::add_log(joint_start_align_score, score_sum * Utils::score_log_base);
        }

        auto node_length = paths_index.nodeLength(gbwt::Node::id(start_search_path.gbwt_search.first.node));
        assert(start_search_path.end_offset <= node_length);

        for (auto & end_search_paths_start_node: end_search_paths_start_nodes_index) {

            auto path_it = find(start_search_path.path.begin(), start_search_path.path.end(), end_search_paths_start_node.first); 

            while (path_it != start_search_path.path.end()) {

                auto main_path_start_idx = path_it - start_search_path.path.begin();

                for (auto end_alignment_idx: end_search_paths_start_node.second) {

                    AlignmentSearchPath complete_paired_align_search_path = start_search_path;
                    mergeAlignmentSearchPath(&complete_paired_align_search_path, main_path_start_idx, end_align_search_paths.at(end_alignment_idx));

                    if (!complete_paired_align_search_path.gbwt_search.first.empty() && complete_paired_align_search_path.fragmentLength() <= max_pair_frag_length) {

                        paired_align_search_paths->emplace_back(complete_paired_align_search_path);                         
                    }
                }

                ++path_it;
                path_it = find(path_it, start_search_path.path.end(), end_search_paths_start_node.first); 
            }
        }

        paired_align_search_path_stack.emplace(start_search_path, false);

        paired_align_search_path_stack.top().first.insert_length += (node_length - start_search_path.end_offset);
        paired_align_search_path_stack.top().first.end_offset = node_length;
    }

    // Perform depth-first path extension.
    while (!paired_align_search_path_stack.empty()) {

        const pair<AlignmentSearchPath, bool> cur_paired_align_search_path = paired_align_search_path_stack.top();
        paired_align_search_path_stack.pop();
   
        assert(!cur_paired_align_search_path.first.gbwt_search.first.empty());
        assert(cur_paired_align_search_path.first.path.back() == cur_paired_align_search_path.first.gbwt_search.first.node);

        if (cur_paired_align_search_path.second) {

            auto end_search_paths_start_nodes_index_it = end_search_paths_start_nodes_index.find(cur_paired_align_search_path.first.path.back());

            if (end_search_paths_start_nodes_index_it != end_search_paths_start_nodes_index.end()) {

                for (auto end_alignment_idx: end_search_paths_start_nodes_index_it->second) {

                    AlignmentSearchPath complete_paired_align_search_path = cur_paired_align_search_path.first;
                    complete_paired_align_search_path.insert_length -= complete_paired_align_search_path.end_offset;

                    complete_paired_align_search_path.end_offset = end_align_search_paths.at(end_alignment_idx).start_offset;
                    complete_paired_align_search_path.insert_length += complete_paired_align_search_path.end_offset;

                    mergeAlignmentSearchPath(&complete_paired_align_search_path, cur_paired_align_search_path.first.path.size() - 1, end_align_search_paths.at(end_alignment_idx));

                    if (!complete_paired_align_search_path.gbwt_search.first.empty() && complete_paired_align_search_path.fragmentLength() <= max_pair_frag_length) {

                        paired_align_search_paths->emplace_back(complete_paired_align_search_path);                         
                    }
                }
            }
        }

        if (!end_alignment_in_cycle) {

            auto end_search_paths_nodes_it = end_search_paths_nodes.find(cur_paired_align_search_path.first.path.back());

            if (end_search_paths_nodes_it != end_search_paths_nodes.end()) {

                if (end_search_paths_nodes_it->second == num_unique_end_search_paths) {

                    continue;  
                }
            }
        }
           
        if (cur_paired_align_search_path.first.fragmentLength() + end_alignment.sequence().size() - end_max_left_softclip_length > max_pair_frag_length) {

            continue;
        }

        auto out_edges = paths_index.edges(cur_paired_align_search_path.first.gbwt_search.first.node);

        // End current extension if no outgoing edges exist.
        if (out_edges.empty()) {

            continue;
        }

        auto out_edges_it = out_edges.begin(); 
        assert(out_edges_it != out_edges.end());
        
        while (out_edges_it != out_edges.end()) {

            if (out_edges_it->first != gbwt::ENDMARKER && out_edges_it->first != cur_paired_align_search_path.first.read_align_stats.back().internal_end_next_node) {

                auto extended_gbwt_search = cur_paired_align_search_path.first.gbwt_search;
                paths_index.extend(&extended_gbwt_search, out_edges_it->first);

                // Add new extension to queue if not empty (path found).
                if (!extended_gbwt_search.first.empty()) { 

                    paired_align_search_path_stack.emplace(cur_paired_align_search_path.first, true);

                    paired_align_search_path_stack.top().first.path.emplace_back(extended_gbwt_search.first.node);
                    paired_align_search_path_stack.top().first.gbwt_search = extended_gbwt_search;
                    paired_align_search_path_stack.top().first.end_offset = paths_index.nodeLength(gbwt::Node::id(paired_align_search_path_stack.top().first.path.back()));
                    paired_align_search_path_stack.top().first.insert_length += paired_align_search_path_stack.top().first.end_offset;
                    paired_align_search_path_stack.top().first.read_align_stats.back().internal_end_next_node = gbwt::ENDMARKER;
                }
            }

            ++out_edges_it;
        }
    }

    paired_align_search_paths->emplace_back(AlignmentSearchPath());

    paired_align_search_paths->back().read_align_stats.emplace_back(AlignmentStats());
    paired_align_search_paths->back().read_align_stats.back().score = Utils::doubleToInt((joint_start_align_score - joint_empty_start_align_score) / Utils::noise_score_log_base);

    paired_align_search_paths->back().read_align_stats.emplace_back(AlignmentStats());
    paired_align_search_paths->back().read_align_stats.back().score = Utils::doubleToInt((joint_end_align_score - joint_empty_end_align_score) / Utils::noise_score_log_base);
}

template<class AlignmentType>
void AlignmentPathFinder<AlignmentType>::mergeAlignmentSearchPath(AlignmentSearchPath * main_align_search_path, uint32_t main_path_start_idx, const AlignmentSearchPath & second_align_search_path) const {

    assert(!main_align_search_path->gbwt_search.first.empty());
    assert(!second_align_search_path.gbwt_search.first.empty());

    assert(main_align_search_path->isComplete());
    assert(second_align_search_path.isComplete());

    assert(main_path_start_idx < main_align_search_path->path.size());

    assert(main_align_search_path->read_align_stats.size() == 1);
    assert(second_align_search_path.read_align_stats.size() == 1);

    assert(main_align_search_path->read_align_stats.back().maxInternalOffset() <= max_partial_offset);
    assert(second_align_search_path.read_align_stats.back().maxInternalOffset() <= max_partial_offset);

    if (second_align_search_path.path.size() < main_align_search_path->path.size() - main_path_start_idx) {

        main_align_search_path->clear();
        return;  
    }

    if (main_path_start_idx == 0) {

        const int32_t main_read_left_offset = static_cast<int32_t>(main_align_search_path->start_offset) - static_cast<int32_t>(main_align_search_path->read_align_stats.back().clippedOffsetLeftBases());
        const int32_t second_read_left_offset = static_cast<int32_t>(second_align_search_path.start_offset) - static_cast<int32_t>(second_align_search_path.read_align_stats.back().clippedOffsetLeftBases());

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

                const uint32_t main_read_right_offset = main_align_search_path->end_offset + main_align_search_path->read_align_stats.back().clippedOffsetRightBases();
                const uint32_t second_read_right_offset = second_align_search_path.end_offset + second_align_search_path.read_align_stats.back().clippedOffsetRightBases();

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
    main_align_search_path->read_align_stats.emplace_back(second_align_search_path.read_align_stats.front());
    
    assert(main_path_start_idx == main_align_search_path->path.size());
    assert(second_path_start_idx <= second_align_search_path.path.size());

    while (second_path_start_idx < second_align_search_path.path.size()) {

        main_align_search_path->path.emplace_back(second_align_search_path.path.at(second_path_start_idx));
        paths_index.extend(&(main_align_search_path->gbwt_search), main_align_search_path->path.back());

        if (main_align_search_path->gbwt_search.first.empty()) {

            break;            
        }

        ++second_path_start_idx;
    }
}

template<class AlignmentType>
vector<gbwt::node_type> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodes(const vg::Alignment & alignment) const {

    vector<gbwt::node_type> alignment_start_nodes;

    assert(alignment.path().mapping_size() > 0);
    alignment_start_nodes.emplace_back(Utils::mapping_to_gbwt(alignment.path().mapping(0)));

    return alignment_start_nodes;
}

template<class AlignmentType>
vector<gbwt::node_type> AlignmentPathFinder<AlignmentType>::getAlignmentStartNodes(const vg::MultipathAlignment & alignment) const {

    vector<gbwt::node_type> alignment_start_nodes;

    for (auto & start_idx: alignment.start()) {

        assert(alignment.subpath(start_idx).path().mapping_size() > 0);
        alignment_start_nodes.emplace_back(Utils::mapping_to_gbwt(alignment.subpath(start_idx).path().mapping(0)));
    }

    return alignment_start_nodes;
}

template<class AlignmentType>
vector<uint32_t> AlignmentPathFinder<AlignmentType>::getAlignmentStartSoftclipLengths(const vg::MultipathAlignment & alignment) const {

    vector<uint32_t> start_softclip_lengths;
    AlignmentStats read_align_stats;

    for (auto & start_idx: alignment.start()) {

        read_align_stats.updateLeftSoftclipLength(alignment.subpath(start_idx).path());
        start_softclip_lengths.emplace_back(read_align_stats.left_softclip_length);
    }

    assert(!start_softclip_lengths.empty());
    return start_softclip_lengths;
}

template<class AlignmentType>
vector<uint32_t> AlignmentPathFinder<AlignmentType>::getAlignmentEndSoftclipLengths(const vg::MultipathAlignment & alignment) const {

    vector<uint32_t> end_softclip_lengths;
    AlignmentStats read_align_stats;

    for (auto & subpath: alignment.subpath()) {

        if (subpath.next_size() == 0) {

            read_align_stats.updateRightSoftclipLength(subpath.path());
            end_softclip_lengths.emplace_back(read_align_stats.right_softclip_length);
        }
    }

    assert(!end_softclip_lengths.empty());
    return end_softclip_lengths;
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
bool AlignmentPathFinder<AlignmentType>::filterAlignmentSearchPaths(const vector<AlignmentSearchPath> & align_search_paths, const vector<int32_t> & optimal_align_scores) const {

    double max_min_optim_score_frac = 0;

    for (auto & align_search_path: align_search_paths) {

        if (align_search_path.isComplete()) {

            max_min_optim_score_frac = max(max_min_optim_score_frac, align_search_path.minOptimalScoreFraction(optimal_align_scores));
        }
    }

    if (max_min_optim_score_frac < min_best_score_filter) {

        return true;
    
    } else {

        return false;
    }
}

template class AlignmentPathFinder<vg::Alignment>;
template class AlignmentPathFinder<vg::MultipathAlignment>;

