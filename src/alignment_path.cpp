
#include "alignment_path.hpp"

#include <algorithm>
#include <numeric>

AlignmentPath::AlignmentPath(const pair<gbwt::SearchState, gbwt::size_type> & gbwt_search_in, const bool is_multimap_in, const uint32_t frag_length_in, const uint32_t min_mapq_in, const uint32_t score_sum_in) : gbwt_search(gbwt_search_in), is_multimap(is_multimap_in), frag_length(frag_length_in), min_mapq(min_mapq_in), score_sum(score_sum_in) {}

AlignmentPath::AlignmentPath(const AlignmentSearchPath & align_path_in, const bool is_multimap_in) : gbwt_search(align_path_in.gbwt_search), is_multimap(is_multimap_in), frag_length(align_path_in.fragmentLength()), min_mapq(align_path_in.minMappingQuality()), score_sum(align_path_in.scoreSum()) {}

vector<AlignmentPath> AlignmentPath::alignmentSearchPathsToAlignmentPaths(const vector<AlignmentSearchPath> & align_search_paths, const uint32_t max_score_diff, const bool is_multimap) {

    uint32_t max_score = 0;

    for (auto & align_search_path: align_search_paths) {

        if (!align_search_path.isEmpty()) {

            max_score = max(max_score, align_search_path.scoreSum());
        }
    }    

    vector<AlignmentPath> align_paths;
    align_paths.reserve(align_search_paths.size());

    for (auto & align_search_path: align_search_paths) {

        if (!align_search_path.isEmpty()) {

            const uint32_t score_sum = align_search_path.scoreSum();
            assert(score_sum <= max_score);

            if (max_score - score_sum <= max_score_diff) {

                align_paths.emplace_back(align_search_path, is_multimap);
            }
        }
    }

    align_paths.shrink_to_fit();
    sort(align_paths.rbegin(), align_paths.rend());

    return align_paths;
}

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.gbwt_search == rhs.gbwt_search && lhs.is_multimap == rhs.is_multimap && lhs.frag_length == rhs.frag_length && lhs.min_mapq == rhs.min_mapq && lhs.score_sum == rhs.score_sum);
}

bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    if (lhs.gbwt_search.first.node != rhs.gbwt_search.first.node) {

        return (lhs.gbwt_search.first.node < rhs.gbwt_search.first.node);    
    } 

    if (lhs.gbwt_search.first.range != rhs.gbwt_search.first.range) {

        return (lhs.gbwt_search.first.range < rhs.gbwt_search.first.range);    
    } 

    if (lhs.gbwt_search.second != rhs.gbwt_search.second) {

        return (lhs.gbwt_search.second < rhs.gbwt_search.second);    
    } 

    if (lhs.is_multimap != rhs.is_multimap) {

        return (lhs.is_multimap < rhs.is_multimap);    
    } 

    if (lhs.frag_length != rhs.frag_length) {

        return (lhs.frag_length < rhs.frag_length);    
    } 

    if (lhs.min_mapq != rhs.min_mapq) {

        return (lhs.min_mapq < rhs.min_mapq);    
    } 

    if (lhs.score_sum != rhs.score_sum) {

        return (lhs.score_sum < rhs.score_sum);    
    } 

    return false;
}

ostream & operator<<(ostream & os, const AlignmentPath & align_path) {

    os << gbwt::Node::id(align_path.gbwt_search.first.node);
    os << " | " << align_path.gbwt_search.first.range;
    os << " | " << align_path.gbwt_search.second;
    os << " | " << align_path.is_multimap;
    os << " | " << align_path.frag_length;
    os << " | " << align_path.min_mapq;
    os << " | " << align_path.score_sum;

    return os;
}

ostream & operator<<(ostream & os, const vector<AlignmentPath> & align_path) {

    os << "# " << align_path.size() << endl;

    for (auto & align_path: align_path) {

        os << align_path << endl;
    }

    return os;
}


ReadAlignmentStats::ReadAlignmentStats() {

    mapq = 0;
    score = 0;
    length = 0;

    left_softclip_length = make_pair(0, false);
    right_softclip_length = make_pair(0, false);

    internal_start_offset = make_pair(0, false);
    internal_end_offset = make_pair(0, false);
    internal_end_next_node = gbwt::ENDMARKER;
}

void ReadAlignmentStats::updateLeftSoftClipLength(const vg::Path & path) {

    assert(!left_softclip_length.second);
    left_softclip_length.second = true;

    auto mapping_it = path.mapping().cbegin();
    assert(mapping_it != path.mapping().cend());

    assert(mapping_it->edit_size() > 0);
    const vg::Edit & first_edit = mapping_it->edit(0);

    if (first_edit.from_length() == 0) {

        left_softclip_length.first = first_edit.to_length();   
    } 
}

void ReadAlignmentStats::updateRightSoftClipLength(const vg::Path & path) {

    assert(!right_softclip_length.second);
    right_softclip_length.second = true;

    auto mapping_rit = path.mapping().rbegin();
    assert(mapping_rit != path.mapping().rend());   

    assert(mapping_rit->edit_size() > 0);
    const vg::Edit & last_edit = mapping_rit->edit(mapping_rit->edit_size() - 1);

    if (last_edit.from_length() == 0) {

        right_softclip_length.first = last_edit.to_length();   
    }
}

void ReadAlignmentStats::updateInternalStartOffset(const uint32_t offset, const bool is_first) {

    if (!internal_start_offset.second) {

        internal_start_offset.first = 0;
    }

    internal_start_offset.second = true;
    internal_start_offset.first += offset;

    if (is_first) {

        assert(left_softclip_length.second);
        assert(left_softclip_length.first <= offset);

        internal_start_offset.first -= left_softclip_length.first;
    }
}

void ReadAlignmentStats::updateInternalEndOffset(const uint32_t offset, const bool is_last) {

    if (!internal_end_offset.second) {

        internal_end_offset.first = 0;
    }

    internal_end_offset.second = true;
    internal_end_offset.first += offset;

    if (is_last) {

        assert(right_softclip_length.second);
        assert(right_softclip_length.first <= offset);

        internal_end_offset.first -= right_softclip_length.first;
    }
}

int32_t ReadAlignmentStats::adjustedScore() const {

    int32_t adjust_score = score;

    if (internal_start_offset.second) {

        adjust_score -= internal_start_offset.first;
    }

    if (internal_end_offset.second) {

        adjust_score -= internal_end_offset.first;
    }

    return adjust_score;
}

uint32_t ReadAlignmentStats::clippedOffsetLeftBases() const {

    assert(left_softclip_length.second);

    if (internal_start_offset.second) {

        return (left_softclip_length.first + internal_start_offset.first);

    } else {

        return left_softclip_length.first;
    }
}

uint32_t ReadAlignmentStats::clippedOffsetRightBases() const {

    assert(right_softclip_length.second);

    if (internal_end_offset.second) {

        return (right_softclip_length.first + internal_end_offset.first);

    } else {

        return right_softclip_length.first;
    }
}

uint32_t ReadAlignmentStats::clippedOffsetTotalBases() const {

    return (clippedOffsetLeftBases() + clippedOffsetRightBases());
}

bool operator==(const ReadAlignmentStats & lhs, const ReadAlignmentStats & rhs) { 

    return (lhs.mapq == rhs.mapq && lhs.score == rhs.score && lhs.length == rhs.length && lhs.left_softclip_length == rhs.left_softclip_length && lhs.right_softclip_length == rhs.right_softclip_length && lhs.internal_start_offset == rhs.internal_start_offset && lhs.internal_end_offset == rhs.internal_end_offset);
}

bool operator!=(const ReadAlignmentStats & lhs, const ReadAlignmentStats & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const ReadAlignmentStats & lhs, const ReadAlignmentStats & rhs) { 

    if (lhs.mapq != rhs.mapq) {

        return (lhs.mapq < rhs.mapq);    
    } 

    if (lhs.score != rhs.score) {

        return (lhs.score < rhs.score);    
    } 

    if (lhs.length != rhs.length) {

        return (lhs.length < rhs.length);    
    } 

    if (lhs.left_softclip_length != rhs.left_softclip_length) {

        return (lhs.left_softclip_length < rhs.left_softclip_length);    
    } 

   if (lhs.right_softclip_length != rhs.right_softclip_length) {

        return (lhs.right_softclip_length < rhs.right_softclip_length);    
    } 

    if (lhs.internal_start_offset != rhs.internal_start_offset) {

        return (lhs.internal_start_offset < rhs.internal_start_offset);    
    } 

    if (lhs.internal_end_offset != rhs.internal_end_offset) {

        return (lhs.internal_end_offset < rhs.internal_end_offset);    
    } 

    return false;
}

ostream & operator<<(ostream & os, const ReadAlignmentStats & read_stats) {

    os << read_stats.mapq;
    os << "," << read_stats.score;
    os << "," << read_stats.length;
    os << "," << read_stats.left_softclip_length;
    os << "," << read_stats.right_softclip_length;
    os << "," << read_stats.internal_start_offset;
    os << "," << read_stats.internal_end_offset;
    os << "," << read_stats.internal_end_next_node;

    return os;
}


AlignmentSearchPath::AlignmentSearchPath() {
    
    start_offset = 0;
    end_offset = 0;

    insert_length = 0;
}

uint32_t AlignmentSearchPath::fragmentLength() const {

    assert(!read_stats.empty());
    assert(read_stats.size() <= 2);

    if (read_stats.size() == 1) {

        assert(insert_length >= 0);

        if (insert_length == 0) {

            return read_stats.front().length;
        
        } else {

            int32_t frag_length = read_stats.front().length + insert_length;
            assert(frag_length >= 0);

            assert(read_stats.front().clippedOffsetRightBases() <= frag_length);
            return frag_length - read_stats.front().clippedOffsetRightBases();
        }

    } else {

        assert(read_stats.size() == 2);

        int32_t frag_length = read_stats.front().length + read_stats.back().length + insert_length;
        assert(frag_length >= 0);

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

        score_sum += stats.adjustedScore();
    }

    return max(0, score_sum);
}

double AlignmentSearchPath::minBestScoreFraction() const {

    double min_best_score_frac = 1;
    assert(!read_stats.empty());

    for (auto & stats: read_stats) {

        assert(stats.adjustedScore() <= static_cast<int32_t>(stats.length));
        min_best_score_frac = min(min_best_score_frac, stats.adjustedScore() / static_cast<double>(stats.length));
    }

    return max(0.0, min_best_score_frac);
}

double AlignmentSearchPath::maxSoftclipFraction() const {

    double max_softclip_frac = 0;
    assert(!read_stats.empty());

    for (auto & stats: read_stats) {

        assert(stats.left_softclip_length.second);
        assert(stats.right_softclip_length.second);

        assert(stats.left_softclip_length.first + stats.right_softclip_length.first <= stats.length);
        max_softclip_frac = max(max_softclip_frac, (stats.left_softclip_length.first + stats.right_softclip_length.first) / static_cast<double>(stats.length));
    }

    return max_softclip_frac;
}

bool AlignmentSearchPath::isEmpty() const {

    if (path.empty() || gbwt_search.first.empty()) {

        return true;
    }

    assert(gbwt_search.first.node == path.back());
    return false;
}

void AlignmentSearchPath::clear() {

    path.clear();

    gbwt_search.first = gbwt::SearchState();
    assert(gbwt_search.first.empty());
}

bool operator==(const AlignmentSearchPath & lhs, const AlignmentSearchPath & rhs) { 

    return (lhs.path == rhs.path && lhs.gbwt_search == rhs.gbwt_search && lhs.start_offset == rhs.start_offset && lhs.end_offset == rhs.end_offset && lhs.insert_length == rhs.insert_length && lhs.read_stats == rhs.read_stats);
}

bool operator!=(const AlignmentSearchPath & lhs, const AlignmentSearchPath & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const AlignmentSearchPath & lhs, const AlignmentSearchPath & rhs) { 

    if (lhs.gbwt_search.first.node != rhs.gbwt_search.first.node) {

        return (lhs.gbwt_search.first.node < rhs.gbwt_search.first.node);    
    } 

    if (lhs.gbwt_search.first.range != rhs.gbwt_search.first.range) {

        return (lhs.gbwt_search.first.range < rhs.gbwt_search.first.range);    
    } 

    if (lhs.gbwt_search.second != rhs.gbwt_search.second) {

        return (lhs.gbwt_search.second < rhs.gbwt_search.second);    
    } 

    if (lhs.start_offset != rhs.start_offset) {

        return (lhs.start_offset < rhs.start_offset);    
    } 

    if (lhs.end_offset != rhs.end_offset) {

        return (lhs.end_offset < rhs.end_offset);    
    } 

    if (lhs.insert_length != rhs.insert_length) {

        return (lhs.insert_length < rhs.insert_length);    
    } 

    if (lhs.read_stats != rhs.read_stats) {

        return (lhs.read_stats < rhs.read_stats);    
    } 

    if (lhs.path.size() != rhs.path.size()) {

        return (lhs.path.size() < rhs.path.size());    
    }

    for (size_t i = 0; i < lhs.path.size(); ++i) {

        if (lhs.path.at(i) != rhs.path.at(i)) {

            return (lhs.path.at(i) < rhs.path.at(i));    
        } 
    }

    return false;
}

ostream & operator<<(ostream & os, const AlignmentSearchPath & align_search_path) {

    os << "(" << align_search_path.path << ")";
    os << " | " << gbwt::Node::id(align_search_path.gbwt_search.first.node);
    os << " | " << align_search_path.gbwt_search.first.size();
    os << " | " << align_search_path.gbwt_search.second;
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

