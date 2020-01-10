
#include "paths_index.hpp"

#include <sstream>
#include <math.h>

#include "utils.hpp"


PathsIndex::PathsIndex(const gbwt::GBWT & gbwt_index, const vg::Graph & graph) : index_(gbwt_index) {

    node_lengths = vector<int32_t>(graph.node_size() + 1, 0);

    for (auto & node: graph.node()) {

        if (node.id() >= node_lengths.size()) {

            node_lengths.resize(node.id() + 1, 0);
        }

        assert(node_lengths.at(node.id()) == 0);
        node_lengths.at(node.id()) = node.sequence().size();
    }
}

PathsIndex::PathsIndex(const gbwt::GBWT & gbwt_index, const handlegraph::PathPositionHandleGraph & graph) : index_(gbwt_index) {

    node_lengths = vector<int32_t>(graph.get_node_count() + 1, 0);

    assert(graph.for_each_handle([&](const handlegraph::handle_t & handle) {

        auto id = graph.get_id(handle);

        while (id >= node_lengths.size()) {

            node_lengths.resize(node_lengths.size() * 2, 0);
        }

        assert(node_lengths.at(id) == 0);
        node_lengths.at(id) =  graph.get_length(handle);
    }));
} 

const gbwt::GBWT & PathsIndex::index() const {

    return index_;
}
        
int32_t PathsIndex::nodeLength(const int32_t node_id) const {

    return node_lengths.at(node_id);
}

string PathsIndex::pathName(const int32_t path_id) const {

    stringstream sstream;

    if (!index_.hasMetadata() || !index_.metadata.hasPathNames() || index_.metadata.paths() <= path_id || !index_.metadata.hasSampleNames()) {
        
        sstream << path_id + 1;
    
    } else {

        const gbwt::PathName& path_name = index_.metadata.path(path_id);

        sstream << index_.metadata.sample(path_name.sample);

        if (index_.metadata.hasContigNames()) {

            sstream << "_" << index_.metadata.contig(path_name.contig);
            sstream << "_" << path_name.phase;
            sstream << "_" << path_name.count;
        }
    }

    return sstream.str();
}

int32_t PathsIndex::pathLength(const int32_t path_id) const {

    int32_t path_length = 0;
    
    for (auto & node: index_.extract(path_id)) {

        path_length += nodeLength(gbwt::Node::id(node));
    }

    return path_length;
}

double PathsIndex::effectivePathLength(const int32_t path_id, const FragmentLengthDist & fragment_length_dist) const {

    const int32_t path_length = pathLength(path_id); 

    if (path_length == 0) {

        return 0;
    }

    // https://en.wikipedia.org/wiki/Truncated_normal_distribution
    const double alpha = (1 - fragment_length_dist.mean()) / fragment_length_dist.sd();
    const double beta = (path_length - fragment_length_dist.mean()) / fragment_length_dist.sd();

    const double trunc_fragment_length_mean = fragment_length_dist.mean() + fragment_length_dist.sd() * (calculateLowerPhi(alpha) - calculateLowerPhi(beta)) / (calculateUpperPhi(beta) - calculateUpperPhi(alpha));

    if (!isfinite(trunc_fragment_length_mean)) {

        return 1;
    }

    return max(1.0, path_length - trunc_fragment_length_mean);
}

double PathsIndex::calculateLowerPhi(const double value) const {

    return (exp(-0.5 * pow(value, 2)) / sqrt(2 * acos(-1)));
}

double PathsIndex::calculateUpperPhi(const double value) const {

    return (0.5 * (1 + erf(value / sqrt(2))));
}


