
#include "alignment_path.hpp"

#include <algorithm>
#include <numeric>


AlignmentPath::AlignmentPath() {

    end_offset = 0;
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

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.search.node == rhs.search.node && lhs.search.size() == rhs.search.size() && lhs.end_offset == rhs.end_offset && lhs.ids == rhs.ids && lhs.seq_length == rhs.seq_length && lhs.mapqProb() == rhs.mapqProb() && lhs.scoreSum() == rhs.scoreSum());
}

bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return !(lhs == rhs);
}

ostream & operator<<(ostream & os, const AlignmentPath & align_path) {

    os << gbwt::Node::id(align_path.search.node);
    os << " | " << align_path.search.size();
    os << " | " << align_path.end_offset;
    os << " | (" << align_path.ids << ")";
    os << " | " << align_path.seq_length;
    os << " | (" << align_path.scores << ")";
    os << " | (" << align_path.mapqs << ")";

    return os;
}

ostream & operator<<(ostream & os, const vector<AlignmentPath> & align_paths) {

    os << "# " << align_paths.size() << endl;

    for (auto & align_path: align_paths) {

        os << align_path << endl;
    }

    return os;
}

