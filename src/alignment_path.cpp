
#include "alignment_path.hpp"

#include <algorithm>
#include <numeric>

AlignmentPath::AlignmentPath(const AlignmentSearchPath & align_search_path) : seq_length(align_search_path.seq_length), mapq_comb(align_search_path.mapqComb()), score_sum(align_search_path.scoreSum()), search(align_search_path.search) {}

vector<AlignmentPath> AlignmentPath::alignmentSearchPathsToAlignmentPaths(const vector<AlignmentSearchPath> & align_search_paths) {

    vector<AlignmentPath> align_paths;
    align_paths.reserve(align_search_paths.size());

    for (auto & align_search_path: align_search_paths) {

        if (align_search_path.complete()) {

            align_paths.emplace_back(align_search_path);
        }
    }

    align_paths.shrink_to_fit();

    return align_paths;
}

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.seq_length == rhs.seq_length && lhs.mapq_comb == rhs.mapq_comb && lhs.score_sum == rhs.score_sum && lhs.search == rhs.search);
}

bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    if (lhs.seq_length != rhs.seq_length) {

        return (lhs.seq_length < rhs.seq_length);    
    } 

    if (lhs.mapq_comb != rhs.mapq_comb) {

        return (lhs.mapq_comb < rhs.mapq_comb);    
    } 

    if (lhs.score_sum != rhs.score_sum) {

        return (lhs.score_sum < rhs.score_sum);    
    } 

    if (lhs.search.node != rhs.search.node) {

        return (lhs.search.node < rhs.search.node);    
    } 

    if (lhs.search.range.first != rhs.search.range.first) {

        return (lhs.search.range.first < rhs.search.range.first);    
    } 

    if (lhs.search.range.second != rhs.search.range.second) {

        return (lhs.search.range.second < rhs.search.range.second);    
    }

    return false;
}

ostream & operator<<(ostream & os, const AlignmentPath & align_path) {

    os << align_path.seq_length;
    os << " | " << align_path.mapq_comb;
    os << " | " << align_path.score_sum;
    os << " | (" << align_path.search << ")";

    return os;
}

ostream & operator<<(ostream & os, const vector<AlignmentPath> & align_path) {

    os << "# " << align_path.size() << endl;

    for (auto & align_path: align_path) {

        os << align_path << endl;
    }

    return os;
}

AlignmentSearchPath::AlignmentSearchPath() {

    path_end_pos = 0;
    seq_start_offset = 0;
    seq_end_offset = 0;
    seq_length = 0;
}

double AlignmentSearchPath::mapqProb() const {

    double prob = 1;

    for (auto & mapq: mapqs) {

        if (mapq > 0) {

            prob *= (1 - phred_to_prob(mapq));

        } else {

            return 1;        
        }
    }

    return (1 - prob);
}

uint32_t AlignmentSearchPath::mapqComb() const {

    return round(prob_to_phred(mapqProb()));
}

uint32_t AlignmentSearchPath::scoreSum() const {

    return accumulate(scores.begin(), scores.end(), 0);
}

bool AlignmentSearchPath::complete() const {

    if (path.empty() || path_end_pos != path.size()) {

        return false;
    }

    assert(search.node == path.back());

    return true;
}

ostream & operator<<(ostream & os, const AlignmentSearchPath & align_search_path) {

    os << "(" << align_search_path.path << ")";
    os << " | " << align_search_path.path_end_pos;
    os << " | " << align_search_path.seq_start_offset;
    os << " | " << align_search_path.seq_end_offset;
    os << " | " << gbwt::Node::id(align_search_path.search.node);
    os << " | " << align_search_path.search.size();
    os << " | " << align_search_path.seq_length;
    os << " | (" << align_search_path.mapqs << ")";
    os << " | (" << align_search_path.scores << ")";

    return os;
}

ostream & operator<<(ostream & os, const vector<AlignmentSearchPath> & align_search_path) {

    os << "# " << align_search_path.size() << endl;

    for (auto & align_search_path: align_search_path) {

        os << align_search_path << endl;
    }

    return os;
}

