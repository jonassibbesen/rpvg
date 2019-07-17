
#include "alignment_path.hpp"

#include <algorithm>
#include <numeric>

#include "utils.hpp"


AlignmentPath::AlignmentPath() {
    
    node_length = 0;
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

    return (lhs.path_ids == rhs.path_ids && lhs.node_length == rhs.node_length && lhs.seq_length == rhs.seq_length && lhs.mapqProb() == rhs.mapqProb() && lhs.scoreSum() == rhs.scoreSum());
}

bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    if (lhs.path_ids.size() != rhs.path_ids.size()) {

        return (lhs.path_ids.size() < rhs.path_ids.size());    
    } 

    for (size_t i = 0; i < lhs.path_ids.size(); ++i) {

        if (lhs.path_ids.at(i) != rhs.path_ids.at(i)) {

            return (lhs.path_ids.at(i) < rhs.path_ids.at(i));    
        }         
    }

    if (lhs.node_length != rhs.node_length) {

        return (lhs.node_length < rhs.node_length);
    }

    if (lhs.seq_length != rhs.seq_length) {

        return (lhs.seq_length < rhs.seq_length);
    }

    if (lhs.mapqProb() != rhs.mapqProb()) {

        return (lhs.mapqProb() < rhs.mapqProb());
    }

    if (lhs.scoreSum() != rhs.scoreSum()) {

        return (lhs.scoreSum() < rhs.scoreSum());
    }

    return false;
}

ostream& operator<<(ostream& os, const vector<int32_t> & values) {

    auto values_it = values.cbegin();

    if (values_it == values.cend()) {

        return os;
    }

    os << *values_it;
    ++values_it;

    while (values_it != values.cend()) {

        os << " " << *values_it;
        ++values_it;
    }

    return os;
}

ostream& operator<<(ostream& os, const AlignmentPath & align_path) {

    for (auto & id: align_path.path_ids) {

        os << id << " ";
    }

    os << "| " << align_path.node_length;
    os << " | " << align_path.seq_length;
    os << " | (" << align_path.scores << ")";
    os << " | (" << align_path.mapqs << ")";

    return os;
}

ostream& operator<<(ostream& os, const vector<AlignmentPath> & align_paths) {

    os << "# " << align_paths.size() << endl;

    for (auto & align_path: align_paths) {

        os << align_path << endl;
    }

    return os;
}

