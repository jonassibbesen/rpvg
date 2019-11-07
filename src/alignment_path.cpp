
#include "alignment_path.hpp"

#include <algorithm>
#include <numeric>


AlignmentPath::AlignmentPath() {

    path_end_pos = 0;
    seq_end_offset = 0;
    seq_length = 0;
}

int32_t AlignmentPath::mapqMin() const {

    return *min_element(mapqs.begin(), mapqs.end());
}

double AlignmentPath::mapqProb() const {

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

int32_t AlignmentPath::scoreSum() const {

    return accumulate(scores.begin(), scores.end(), 0);
}

bool AlignmentPath::complete() const {

    if (path_end_pos != path.size()) {

        return false;
    }

    if (search.size() != ids.size()) {

        return false;
    }

    if (path.size() == 0) {

        return true;
    }

    if (search.node != path.back()) {

        return false;
    }

    return true;
}

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.path == rhs.path && lhs.path_end_pos == rhs.path_end_pos && lhs.seq_end_offset == rhs.seq_end_offset && lhs.search.node == rhs.search.node && lhs.search.size() == rhs.search.size() && lhs.ids == rhs.ids && lhs.seq_length == rhs.seq_length && lhs.mapqs == rhs.mapqs && lhs.scores == rhs.scores);
}

bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return !(lhs == rhs);
}

ostream & operator<<(ostream & os, const AlignmentPath & align_path) {

    os << "(" << align_path.path << ")";
    os << " | " << align_path.path_end_pos;
    os << " | " << align_path.seq_end_offset;
    os << " | " << gbwt::Node::id(align_path.search.node);
    os << " | " << align_path.search.size();
    os << " | (" << align_path.ids << ")";
    os << " | " << align_path.seq_length;
    os << " | (" << align_path.mapqs << ")";
    os << " | (" << align_path.scores << ")";

    return os;
}

ostream & operator<<(ostream & os, const vector<AlignmentPath> & align_paths) {

    os << "# " << align_paths.size() << endl;

    for (auto & align_path: align_paths) {

        os << align_path << endl;
    }

    return os;
}

