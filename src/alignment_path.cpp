
#include "alignment_path.hpp"

#include <algorithm>
#include <numeric>

AlignmentPath::AlignmentPath(const uint32_t seq_length_in, const uint32_t min_mapq_in, const uint32_t score_sum_in, const bool is_multimap_in, const gbwt::SearchState & search_state_in) : seq_length(seq_length_in), min_mapq(min_mapq_in), is_multimap(is_multimap_in), score_sum(score_sum_in), search_state(search_state_in) {}

AlignmentPath::AlignmentPath(const AlignmentSearchPath & align_path_in, const bool is_multimap_in) : seq_length(align_path_in.seq_length), min_mapq(align_path_in.min_mapq), score_sum(align_path_in.scoreSum()), search_state(align_path_in.search_state), is_multimap(is_multimap_in) {}

vector<AlignmentPath> AlignmentPath::alignmentSearchPathsToAlignmentPaths(const vector<AlignmentSearchPath> & align_search_paths, const bool is_multimap) {

    vector<AlignmentPath> align_paths;
    align_paths.reserve(align_search_paths.size());

    for (auto & align_search_path: align_search_paths) {

        if (align_search_path.isComplete()) {

            align_paths.emplace_back(align_search_path, is_multimap);
        }
    }

    return align_paths;
}

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.seq_length == rhs.seq_length && lhs.min_mapq == rhs.min_mapq && lhs.score_sum == rhs.score_sum && lhs.is_multimap == rhs.is_multimap && lhs.search_state == rhs.search_state);
}

bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    if (lhs.seq_length != rhs.seq_length) {

        return (lhs.seq_length < rhs.seq_length);    
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

    os << align_path.seq_length;
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

    os << read_stats.score;
    os << "," << read_stats.length;
    os << "," << read_stats.left_softclip_length;
    os << "," << read_stats.right_softclip_length;

    return os;
}


AlignmentSearchPath::AlignmentSearchPath() {

    path_end_idx = 0;
    seq_start_offset = 0;
    seq_end_offset = 0;
    seq_length = 0;
    min_mapq = std::numeric_limits<uint32_t>::max();
}

void AlignmentSearchPath::softclipAdjustSeqLength() {

    assert(!read_stats.empty());
    assert(read_stats.size() <= 2);

    if (read_stats.size() == 2) {

        assert(read_stats.front().right_softclip_length >= 0);
        assert(read_stats.back().left_softclip_length >= 0);

        assert(read_stats.front().right_softclip_length + read_stats.back().left_softclip_length <= seq_length);
        seq_length -= (read_stats.front().right_softclip_length + read_stats.back().left_softclip_length);
    }
}

uint32_t AlignmentSearchPath::scoreSum() const {

    int32_t score_sum = 0;
    assert(!read_stats.empty());

    for (auto & stats: read_stats) {

        score_sum += stats.score;
    }

    return max(0, score_sum);
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

    if (path.empty() || path_end_idx != path.size()) {

        return false;
    }

    assert(search_state.node == path.back());

    return true;
}

ostream & operator<<(ostream & os, const AlignmentSearchPath & align_search_path) {

    os << "(" << align_search_path.path << ")";
    os << " | " << align_search_path.path_end_idx;
    os << " | " << align_search_path.seq_start_offset;
    os << " | " << align_search_path.seq_end_offset;
    os << " | " << gbwt::Node::id(align_search_path.search_state.node);
    os << " | " << align_search_path.search_state.size();
    os << " | " << align_search_path.seq_length;
    os << " | " << align_search_path.min_mapq;
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

