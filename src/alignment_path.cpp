
#include "alignment_path.hpp"
#include "utils.hpp"


AlignmentPath::AlignmentPath() {
    
    node_length = 0;
    seq_length = 0;

    scores = make_pair(0,0);
    mapqs = make_pair(0,0);
}

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.path_ids == rhs.path_ids && lhs.node_length == rhs.node_length && lhs.seq_length == rhs.seq_length && mapqsToProb(lhs.mapqs) == mapqsToProb(rhs.mapqs) && lhs.scores.first + lhs.scores.second == rhs.scores.first + rhs.scores.second);
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

    if (mapqsToProb(lhs.mapqs) != mapqsToProb(rhs.mapqs)) {

        return (mapqsToProb(lhs.mapqs) < mapqsToProb(rhs.mapqs));
    }

    if (lhs.scores.first + lhs.scores.second != rhs.scores.first + rhs.scores.second) {

        return (lhs.scores.first + lhs.scores.second < rhs.scores.first + rhs.scores.second);
    }

    return false;
}

ostream& operator<<(ostream& os, const AlignmentPath & align_path) {

    for (auto & id: align_path.path_ids) {

        os << id << " ";
    }

    os << "| " << align_path.node_length;
    os << " | " << align_path.seq_length;
    os << " | (" << align_path.scores.first << ", " << align_path.scores.second << ")";
    os << " | (" << align_path.mapqs.first << ", " << align_path.mapqs.second << ")";
    os << endl;

    return os;
}

ostream& operator<<(ostream& os, const vector<AlignmentPath> & align_paths) {

    os << "# " << align_paths.size() << endl;

    for (auto & align_path: align_paths) {

        for (auto & id: align_path.path_ids) {

            os << id << " ";
        }

        os << "| " << align_path.node_length;
        os << " | " << align_path.seq_length;
        os << " | (" << align_path.scores.first << ", " << align_path.scores.second << ")";
        os << " | (" << align_path.mapqs.first << ", " << align_path.mapqs.second << ")";
        os << endl;
    }

    return os;
}

