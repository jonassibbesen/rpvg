
#include "alignment_path.hpp"

#include <algorithm>
#include <numeric>

AlignmentPath::AlignmentPath(const uint32_t frag_length_in, const uint32_t min_mapq_in, const uint32_t score_sum_in, const bool is_multimap_in, const gbwt::SearchState & search_state_in) : frag_length(frag_length_in), min_mapq(min_mapq_in), is_multimap(is_multimap_in), score_sum(score_sum_in), search_state(search_state_in) {}

AlignmentPath::AlignmentPath(const AlignmentSearchPath & align_path_in, const bool is_multimap_in) : frag_length(align_path_in.fragmentLength()), min_mapq(align_path_in.minMappingQuality()), score_sum(align_path_in.scoreSum()), search_state(align_path_in.search_state), is_multimap(is_multimap_in) {}

vector<AlignmentPath> AlignmentPath::alignmentSearchPathsToAlignmentPaths(const vector<AlignmentSearchPath> & align_search_paths, const uint32_t max_score_diff, const bool is_multimap) {

    uint32_t max_score = 0;

    for (auto & align_search_path: align_search_paths) {

        if (align_search_path.isComplete()) {

            max_score = max(max_score, align_search_path.scoreSum());
        }
    }    

    vector<AlignmentPath> align_paths;
    align_paths.reserve(align_search_paths.size());

    for (auto & align_search_path: align_search_paths) {

        if (align_search_path.isComplete()) {

            const uint32_t score_sum = align_search_path.scoreSum();
            assert(score_sum <= max_score);

            if (max_score - score_sum <= max_score_diff) {

                align_paths.emplace_back(align_search_path, is_multimap);
            }
        }
    }

    align_paths.shrink_to_fit();

    return align_paths;
}

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.frag_length == rhs.frag_length && lhs.min_mapq == rhs.min_mapq && lhs.score_sum == rhs.score_sum && lhs.is_multimap == rhs.is_multimap && lhs.search_state == rhs.search_state);
}

bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    if (lhs.frag_length != rhs.frag_length) {

        return (lhs.frag_length < rhs.frag_length);    
    } 

    if (lhs.min_mapq != rhs.min_mapq) {

        return (lhs.min_mapq < rhs.min_mapq);    
    } 

    if (lhs.score_sum != rhs.score_sum) {

        return (lhs.score_sum < rhs.score_sum);    
    } 

    if (lhs.is_multimap != rhs.is_multimap) {

        return (lhs.is_multimap < rhs.is_multimap);    
    } 

    if (lhs.search_state.node != rhs.search_state.node) {

        return (lhs.search_state.node < rhs.search_state.node);    
    } 

    if (lhs.search_state.range != rhs.search_state.range) {

        return (lhs.search_state.range < rhs.search_state.range);    
    } 

    return false;
}

ostream & operator<<(ostream & os, const AlignmentPath & align_path) {

    os << align_path.frag_length;
    os << " | " << align_path.min_mapq;
    os << " | " << align_path.score_sum;
    os << " | " << align_path.is_multimap;
    os << " | " << gbwt::Node::id(align_path.search_state.node);
    os << " | " << align_path.search_state.range.first << " " << align_path.search_state.range.second;

    return os;
}

ostream & operator<<(ostream & os, const vector<AlignmentPath> & align_path) {

    os << "# " << align_path.size() << endl;

    for (auto & align_path: align_path) {

        os << align_path << endl;
    }

    return os;
}

ostream & operator<<(ostream & os, const ReadAlignmentStats & read_stats) {

    os << read_stats.mapq;
    os << "," << read_stats.score;
    os << "," << read_stats.length;
    os << "," << read_stats.left_softclip_length;
    os << "," << read_stats.right_softclip_length;
    os << "," << read_stats.internal_start_offset;
    os << "," << read_stats.internal_end_offset;

    return os;
}


ReadAlignmentStats::ReadAlignmentStats() {

    mapq = 0;
    score = 0;
    length = 0;

    left_softclip_length = -1;
    right_softclip_length = -1; 

    internal_start_offset = 0;
    internal_end_offset = 0;
}

void ReadAlignmentStats::updateSoftClippingLengths(const vg::Path & path) {

    auto mapping_it = path.mapping().cbegin();
    assert(mapping_it != path.mapping().cend());

    if (left_softclip_length == -1) {

        left_softclip_length = 0;

        assert(mapping_it->edit_size() > 0);
        const vg::Edit & first_edit = mapping_it->edit(0);

        if (first_edit.from_length() == 0) {

            left_softclip_length = first_edit.to_length();   
        } 
    }

    auto mapping_rit = path.mapping().rbegin();
    assert(mapping_rit != path.mapping().rend());   

    right_softclip_length = 0;

    assert(mapping_rit->edit_size() > 0);
    const vg::Edit & last_edit = mapping_rit->edit(mapping_rit->edit_size() - 1);

    if (last_edit.from_length() == 0) {

        right_softclip_length = last_edit.to_length();   
    }
}

void ReadAlignmentStats::updateInternalStartOffset(const uint32_t offset, const bool is_first) {

    internal_start_offset += offset;

    if (is_first) {

        assert(left_softclip_length <= offset);
        internal_start_offset -= left_softclip_length;
    }
}

void ReadAlignmentStats::updateInternalEndOffset(const uint32_t offset, const bool is_last) {

    internal_end_offset += offset;

    if (is_last) {

        assert(right_softclip_length <= offset);
        internal_end_offset -= right_softclip_length;
    }
}

uint32_t ReadAlignmentStats::clippedOffsetLeftBases() {

    return (left_softclip_length + internal_start_offset);
}

uint32_t ReadAlignmentStats::clippedOffsetRightBases() {

    return (right_softclip_length + internal_end_offset);
}

AlignmentSearchPath::AlignmentSearchPath() {

    path_idx = 0;
    path_offset = 0;

    start_offset = 0;
    end_offset = 0;

    insert_length = 0;
}

uint32_t AlignmentSearchPath::fragmentLength() const {

    assert(!read_stats.empty());
    assert(read_stats.size() <= 2);

    if (read_stats.size() == 1) {

        assert(insert_length >= 0);
        return (read_stats.front().length + insert_length);

    } else {

        assert(read_stats.size() == 2);

        int32_t frag_length = read_stats.front().length + read_stats.back().length + insert_length;
        assert(frag_length >= 0);

        assert(read_stats.front().right_softclip_length >= 0);
        assert(read_stats.back().left_softclip_length >= 0);

        assert(read_stats.front().clippedOffsetRightBases() + read_stats.back().clippedOffsetLeftBases() <= frag_length);
        return (frag_length - read_stats.front().clippedOffsetRightBases() - read_stats.back().clippedOffsetLeftBases());
    }
}

uint32_t AlignmentSearchPath::minMappingQuality() const {

    uint32_t min_mapq = std::numeric_limits<uint32_t>::max();
    assert(!read_stats.empty());

    for (auto & stats: read_stats) {

        min_mapq = min(min_mapq, stats.mapq);
    }

    return min_mapq;
}

uint32_t AlignmentSearchPath::scoreSum() const {

    int32_t score_sum = 0;
    assert(!read_stats.empty());

    for (auto & stats: read_stats) {

        score_sum += (stats.score - stats.internal_start_offset - stats.internal_end_offset);
    }

    return max(0, score_sum);
}

double AlignmentSearchPath::minBestScoreFraction() const {

    double min_best_score_frac = 1;
    assert(!read_stats.empty());

    for (auto & stats: read_stats) {

        assert(stats.score <= static_cast<int32_t>(stats.length));
        min_best_score_frac = min(min_best_score_frac, max(0, stats.score) / static_cast<double>(stats.length));
    }

    return min_best_score_frac;
}

double AlignmentSearchPath::maxSoftclipFraction() const {

    double max_softclip_frac = 0;
    assert(!read_stats.empty());

    for (auto & stats: read_stats) {

        assert(stats.left_softclip_length >= 0);
        assert(stats.right_softclip_length >= 0);

        assert(stats.left_softclip_length + stats.right_softclip_length <= stats.length);
        max_softclip_frac = max(max_softclip_frac, (stats.left_softclip_length + stats.right_softclip_length) / static_cast<double>(stats.length));
    }

    return max_softclip_frac;
}

bool AlignmentSearchPath::isComplete() const {

    if (path.empty() || path_idx != path.size()) {

        return false;
    }

    assert(search_state.node == path.back());

    return true;
}

ostream & operator<<(ostream & os, const AlignmentSearchPath & align_search_path) {

    os << "(" << align_search_path.path << ")";
    os << " | " << align_search_path.path_idx;
    os << " | " << align_search_path.path_offset;
    os << " | " << gbwt::Node::id(align_search_path.search_state.node);
    os << " | " << align_search_path.search_state.size();
    os << " | " << align_search_path.start_offset;
    os << " | " << align_search_path.end_offset;
    os << " | " << align_search_path.insert_length;
    os << " | (" << align_search_path.read_stats << ")";

    return os;
}

ostream & operator<<(ostream & os, const vector<AlignmentSearchPath> & align_search_path) {

    os << "# " << align_search_path.size() << endl;

    for (auto & align_search_path: align_search_path) {

        os << align_search_path << endl;
    }

    return os;
}

