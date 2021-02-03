
#include "gbwt/fast_locate.h"

#include "alignment_path.hpp"

#include <algorithm>
#include <numeric>

AlignmentPath::AlignmentPath(const pair<gbwt::SearchState, gbwt::size_type> & gbwt_search_in, const bool is_multimap_in, const uint32_t frag_length_in, const uint32_t min_mapq_in, const int32_t score_sum_in) : gbwt_search(gbwt_search_in), is_multimap(is_multimap_in), frag_length(frag_length_in), min_mapq(min_mapq_in), score_sum(score_sum_in) {}

AlignmentPath::AlignmentPath(const AlignmentSearchPath & align_path_in, const bool is_multimap_in) : gbwt_search(align_path_in.gbwt_search), is_multimap(is_multimap_in), frag_length(align_path_in.fragmentLength()), min_mapq(align_path_in.minMappingQuality()), score_sum(align_path_in.scoreSum()) {}

vector<AlignmentPath> AlignmentPath::alignmentSearchPathsToAlignmentPaths(const vector<AlignmentSearchPath> & align_search_paths, const bool is_multimap) {

    vector<AlignmentPath> align_paths;
    align_paths.reserve(align_search_paths.size());

    if (!align_search_paths.empty()) {

        for (size_t i = 0; i < align_search_paths.size() - 1; ++i) {

            if (align_search_paths.at(i).isComplete()) {
                
                assert(!align_search_paths.at(i).gbwt_search.first.empty());

                align_paths.emplace_back(align_search_paths.at(i), is_multimap);
                assert(align_paths.front().min_mapq == align_paths.back().min_mapq);
            }
        }
    }

    sort(align_paths.rbegin(), align_paths.rend());

    if (!align_paths.empty()) {

        assert(align_search_paths.back().gbwt_search.first.empty());
        assert(align_search_paths.back().insert_length == 0);

        assert(align_search_paths.back().read_align_stats.size() == 2);
        assert(!align_search_paths.back().isInternal());

        const double read_1_error_prob = 1 / (1 + exp(align_search_paths.back().read_align_stats.front().score * Utils::noise_log_base));
        const double read_2_error_prob = 1 / (1 + exp(align_search_paths.back().read_align_stats.back().score * Utils::noise_log_base));

        const double noise_prob = 1 - ((1 - read_1_error_prob) * (1 - read_2_error_prob));

        if (Utils::doubleCompare(noise_prob, 0)) {

            align_paths.emplace_back(make_pair(gbwt::SearchState(), gbwt::FastLocate::NO_POSITION), is_multimap, 0, align_paths.front().min_mapq, numeric_limits<int32_t>::lowest());

        } else {

            assert(noise_prob > 0);
            assert(noise_prob < 1 || Utils::doubleCompare(noise_prob, 1));

            align_paths.emplace_back(make_pair(gbwt::SearchState(), gbwt::FastLocate::NO_POSITION), is_multimap, 0, align_paths.front().min_mapq, Utils::doubleToInt(log(noise_prob) / Utils::noise_log_base));
        }

        assert(align_paths.back().score_sum <= 0);
    }

    align_paths.shrink_to_fit();

    return align_paths;
}

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.gbwt_search == rhs.gbwt_search && 
            lhs.is_multimap == rhs.is_multimap && 
            lhs.frag_length == rhs.frag_length && 
            lhs.min_mapq == rhs.min_mapq && 
            lhs.score_sum == rhs.score_sum);
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


InternalAlignment::InternalAlignment() {

    is_internal = false;
    penalty = 0;
    offset = 0;

    max_offset = 0;    
}

bool operator==(const InternalAlignment & lhs, const InternalAlignment & rhs) { 

    return (lhs.is_internal == rhs.is_internal && 
            lhs.penalty == rhs.penalty && 
            lhs.offset == rhs.offset && 
            lhs.max_offset == rhs.max_offset);
}

bool operator!=(const InternalAlignment & lhs, const InternalAlignment & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const InternalAlignment & lhs, const InternalAlignment & rhs) { 

    if (lhs.is_internal != rhs.is_internal) {

        return (lhs.is_internal < rhs.is_internal);    
    } 

    if (lhs.penalty != rhs.penalty) {

        return (lhs.penalty < rhs.penalty);    
    } 

    if (lhs.offset != rhs.offset) {

        return (lhs.offset < rhs.offset);    
    } 

    if (lhs.max_offset != rhs.max_offset) {

        return (lhs.max_offset < rhs.max_offset);    
    } 

    return false;
}

ostream & operator<<(ostream & os, const InternalAlignment & internal_align) {

    os << internal_align.is_internal;
    os << "," << internal_align.penalty;
    os << "," << internal_align.offset;
    os << "," << internal_align.max_offset;

    return os;
}


AlignmentStats::AlignmentStats() {

    mapq = 0;
    score = 0;

    length = 0;
    complete = false;

    left_softclip_length = 0;
    right_softclip_length = 0;

    internal_end_next_node = gbwt::ENDMARKER;
}

void AlignmentStats::updateLeftSoftclipLength(const vg::Path & path) {

    auto mapping_it = path.mapping().cbegin();
    assert(mapping_it != path.mapping().cend());

    assert(mapping_it->edit_size() > 0);
    const vg::Edit & first_edit = mapping_it->edit(0);

    if (first_edit.from_length() == 0) {

        left_softclip_length = first_edit.to_length();   
    
    } else {

        left_softclip_length = 0;
    } 
}

void AlignmentStats::updateRightSoftclipLength(const vg::Path & path) {

    auto mapping_rit = path.mapping().rbegin();
    assert(mapping_rit != path.mapping().rend());   

    assert(mapping_rit->edit_size() > 0);
    const vg::Edit & last_edit = mapping_rit->edit(mapping_rit->edit_size() - 1);

    if (last_edit.from_length() == 0) {

        right_softclip_length = last_edit.to_length();   

    } else {

        right_softclip_length = 0;
    }
}

bool AlignmentStats::isInternal() const {

    return (internal_start.is_internal || internal_end.is_internal);
}

uint32_t AlignmentStats::internalPenalty() const {

    return (internal_start.penalty + internal_end.penalty);
}

uint32_t AlignmentStats::maxInternalOffset() const {

    return max(internal_start.offset, internal_end.offset);
}

int32_t AlignmentStats::adjustedScore() const {

    return (score - static_cast<int32_t>(internalPenalty()));
}

uint32_t AlignmentStats::clippedOffsetLeftBases() const {

    return (left_softclip_length + internal_start.offset);
}

uint32_t AlignmentStats::clippedOffsetRightBases() const {

    return (right_softclip_length + internal_end.offset);
}

uint32_t AlignmentStats::clippedOffsetTotalBases() const {

    return (clippedOffsetLeftBases() + clippedOffsetRightBases());
}

bool operator==(const AlignmentStats & lhs, const AlignmentStats & rhs) { 

    return (lhs.mapq == rhs.mapq && 
            lhs.score == rhs.score && 
            lhs.length == rhs.length && 
            lhs.complete == rhs.complete && 
            lhs.left_softclip_length == rhs.left_softclip_length && 
            lhs.right_softclip_length == rhs.right_softclip_length && 
            lhs.internal_start == rhs.internal_start && 
            lhs.internal_end == rhs.internal_end && 
            lhs.internal_end_next_node == rhs.internal_end_next_node);
}

bool operator!=(const AlignmentStats & lhs, const AlignmentStats & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const AlignmentStats & lhs, const AlignmentStats & rhs) { 

    if (lhs.mapq != rhs.mapq) {

        return (lhs.mapq < rhs.mapq);    
    } 

    if (lhs.score != rhs.score) {

        return (lhs.score < rhs.score);    
    } 

    if (lhs.length != rhs.length) {

        return (lhs.length < rhs.length);    
    } 

    if (lhs.complete != rhs.complete) {

        return (lhs.complete < rhs.complete);    
    } 

    if (lhs.left_softclip_length != rhs.left_softclip_length) {

        return (lhs.left_softclip_length < rhs.left_softclip_length);    
    } 

   if (lhs.right_softclip_length != rhs.right_softclip_length) {

        return (lhs.right_softclip_length < rhs.right_softclip_length);    
    } 

    if (lhs.internal_start != rhs.internal_start) {

        return (lhs.internal_start < rhs.internal_start);    
    } 

    if (lhs.internal_end != rhs.internal_end) {

        return (lhs.internal_end < rhs.internal_end);    
    } 

    if (lhs.internal_end_next_node != rhs.internal_end_next_node) {

        return (lhs.internal_end_next_node < rhs.internal_end_next_node);    
    } 

    return false;
}

ostream & operator<<(ostream & os, const AlignmentStats & read_align_stats) {

    os << read_align_stats.mapq;
    os << "," << read_align_stats.score;
    os << "," << read_align_stats.length;
    os << "," << read_align_stats.complete;
    os << "," << read_align_stats.left_softclip_length;
    os << "," << read_align_stats.right_softclip_length;
    os << ",(" << read_align_stats.internal_start << ")";
    os << ",(" << read_align_stats.internal_end << ")";
    os << "," << read_align_stats.internal_end_next_node;

    return os;
}


AlignmentSearchPath::AlignmentSearchPath() {

    assert(gbwt_search.first.empty());
    gbwt_search.second = gbwt::FastLocate::NO_POSITION;
    
    start_offset = 0;
    end_offset = 0;

    insert_length = 0;
}

uint32_t AlignmentSearchPath::fragmentLength() const {

    assert(!read_align_stats.empty());
    assert(read_align_stats.size() <= 2);

    if (read_align_stats.size() == 1) {

        assert(insert_length >= 0);

        if (insert_length == 0) {

            return read_align_stats.front().length;
        
        } else {

            int32_t frag_length = read_align_stats.front().length + insert_length;
            assert(frag_length >= 0);

            assert(read_align_stats.front().clippedOffsetRightBases() <= frag_length);
            return frag_length - read_align_stats.front().clippedOffsetRightBases();
        }

    } else {

        assert(read_align_stats.size() == 2);

        int32_t frag_length = read_align_stats.front().length + read_align_stats.back().length + insert_length;
        assert(frag_length >= 0);

        assert(read_align_stats.front().clippedOffsetRightBases() + read_align_stats.back().clippedOffsetLeftBases() <= frag_length);
        return (frag_length - read_align_stats.front().clippedOffsetRightBases() - read_align_stats.back().clippedOffsetLeftBases());
    }
}

uint32_t AlignmentSearchPath::minMappingQuality() const {

    uint32_t min_mapq = std::numeric_limits<uint32_t>::max();
    assert(!read_align_stats.empty());

    for (auto & stats: read_align_stats) {

        min_mapq = min(min_mapq, stats.mapq);
    }

    return min_mapq;
}

int32_t AlignmentSearchPath::scoreSum() const {

    int32_t score_sum = 0;
    assert(!read_align_stats.empty());

    for (auto & stats: read_align_stats) {

        score_sum += stats.adjustedScore();
    }

    return score_sum;
}

double AlignmentSearchPath::minOptimalScoreFraction(const vector<int32_t> & optimal_align_scores) const {

    double min_optim_score_frac = 1;

    assert(!read_align_stats.empty());
    assert(optimal_align_scores.size() == read_align_stats.size());

    for (size_t i = 0; i < read_align_stats.size(); ++i) {

        assert(read_align_stats.at(i).adjustedScore() <= optimal_align_scores.at(i));
        min_optim_score_frac = min(min_optim_score_frac, read_align_stats.at(i).adjustedScore() / static_cast<double>(optimal_align_scores.at(i)));
    }

    return max(0.0, min_optim_score_frac);
}

double AlignmentSearchPath::maxSoftclipFraction() const {

    double max_softclip_frac = 0;
    assert(!read_align_stats.empty());

    for (auto & stats: read_align_stats) {

        assert(stats.left_softclip_length + stats.right_softclip_length <= stats.length);
        max_softclip_frac = max(max_softclip_frac, (stats.left_softclip_length + stats.right_softclip_length) / static_cast<double>(stats.length));
    }

    return max_softclip_frac;
}

bool AlignmentSearchPath::isComplete() const {

    for (auto & stats: read_align_stats) {

        if (!stats.complete) {

            return false;
        }
    }

    return true;
}

bool AlignmentSearchPath::isInternal() const {

    for (auto & stats: read_align_stats) {

        if (stats.isInternal()) {

            return true;
        }
    }

    return false;
}


void AlignmentSearchPath::clear() {

    path.clear();

    gbwt_search.first = gbwt::SearchState();
    assert(gbwt_search.first.empty());

    gbwt_search.second = gbwt::FastLocate::NO_POSITION;
}

bool operator==(const AlignmentSearchPath & lhs, const AlignmentSearchPath & rhs) { 

    return (lhs.path == rhs.path && 
            lhs.gbwt_search == rhs.gbwt_search && 
            lhs.start_offset == rhs.start_offset && 
            lhs.end_offset == rhs.end_offset && 
            lhs.insert_length == rhs.insert_length && 
            lhs.read_align_stats == rhs.read_align_stats);
}

bool operator!=(const AlignmentSearchPath & lhs, const AlignmentSearchPath & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const AlignmentSearchPath & lhs, const AlignmentSearchPath & rhs) { 

    if (lhs.path.size() != rhs.path.size()) {

        return (lhs.path.size() < rhs.path.size());    
    }

    for (size_t i = 0; i < lhs.path.size(); ++i) {

        if (lhs.path.at(i) != rhs.path.at(i)) {

            return (lhs.path.at(i) < rhs.path.at(i));    
        } 
    }

    if (lhs.gbwt_search.first.node != rhs.gbwt_search.first.node) {

        return (lhs.gbwt_search.first.node < rhs.gbwt_search.first.node);    
    } 

    if (lhs.gbwt_search.first.range != rhs.gbwt_search.first.range) {

        return (lhs.gbwt_search.first.range < rhs.gbwt_search.first.range);    
    } 

    if (lhs.gbwt_search.second != rhs.gbwt_search.second) {

        return (lhs.gbwt_search.second < rhs.gbwt_search.second);    
    } 

    if (lhs.insert_length != rhs.insert_length) {

        return (lhs.insert_length < rhs.insert_length);    
    } 

    if (lhs.scoreSum() != rhs.scoreSum()) {

        return (lhs.scoreSum() < rhs.scoreSum());    
    } 

    if (lhs.read_align_stats != rhs.read_align_stats) {

        return (lhs.read_align_stats < rhs.read_align_stats);    
    } 

    if (lhs.start_offset != rhs.start_offset) {

        return (lhs.start_offset < rhs.start_offset);    
    } 

    if (lhs.end_offset != rhs.end_offset) {

        return (lhs.end_offset < rhs.end_offset);    
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
    os << " | (" << align_search_path.read_align_stats << ")";

    return os;
}

ostream & operator<<(ostream & os, const vector<AlignmentSearchPath> & align_search_path) {

    os << "# " << align_search_path.size() << endl;

    for (auto & align_search_path: align_search_path) {

        os << align_search_path << endl;
    }

    return os;
}

