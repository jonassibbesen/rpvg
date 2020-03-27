
#include "alignment_path.hpp"

#include <algorithm>
#include <numeric>

AlignmentPath::AlignmentPath(const uint32_t seq_length_in, const uint32_t mapq_comb_in, const uint32_t score_sum_in, const vector<gbwt::size_type> & ids_in) : seq_length(seq_length_in), mapq_comb(mapq_comb_in), score_sum(score_sum_in), ids(ids_in) {}

AlignmentPath::AlignmentPath(const AlignmentSearchPath & align_path_in, const vector<gbwt::size_type> & ids_in) : seq_length(align_path_in.seq_length), mapq_comb(align_path_in.mapqComb()), score_sum(align_path_in.scoreSum()), ids(ids_in) {}

vector<AlignmentPath> AlignmentPath::alignmentSearchPathsToAlignmentPaths(const vector<AlignmentSearchPath> & align_search_paths, const PathsIndex & paths_index) {

    vector<AlignmentPath> align_paths;
    align_paths.reserve(align_search_paths.size());

    for (auto & align_search_path: align_search_paths) {

        if (align_search_path.complete()) {

            auto align_search_path_ids = paths_index.index().locate(align_search_path.search);

            if (paths_index.index().bidirectional()) {

                for (auto & id: align_search_path_ids) {

                    id = gbwt::Path::id(id);
                }
            }

            auto align_paths_it = align_paths.begin();

            while (align_paths_it != align_paths.end()) {

                if (align_paths_it->seq_length == align_search_path.seq_length && align_paths_it->score_sum == align_search_path.scoreSum()) {

                    assert(align_paths_it->mapq_comb == align_search_path.mapqComb());
                    break;
                }

                ++align_paths_it;
            }

            if (align_paths_it == align_paths.end()) {

                align_paths.emplace_back(align_search_path, align_search_path_ids);

            } else {

                align_paths_it->ids.insert(align_paths_it->ids.end(), align_search_path_ids.begin(), align_search_path_ids.end());
            }
        }
    }

    align_paths.shrink_to_fit();

    for (auto & align_path: align_paths) {

        sort(align_path.ids.begin(), align_path.ids.end());
    }

    return align_paths;
}

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.seq_length == rhs.seq_length && lhs.mapq_comb == rhs.mapq_comb && lhs.score_sum == rhs.score_sum && lhs.ids == rhs.ids);
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

    if (lhs.ids.size() != rhs.ids.size()) {

        return (lhs.ids.size() < rhs.ids.size());    
    } 

    for (size_t i = 0; i < lhs.ids.size(); ++i) {

        if (lhs.ids.at(i) != rhs.ids.at(i)) {

            return (lhs.ids.at(i) < rhs.ids.at(i));    
        }         
    }   

    return false;
}

ostream & operator<<(ostream & os, const AlignmentPath & align_path) {

    os << align_path.seq_length;
    os << " | " << align_path.mapq_comb;
    os << " | " << align_path.score_sum;
    os << " | (" << align_path.ids << ")";

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

