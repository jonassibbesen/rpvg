
#include "paths_index.hpp"

#include <sstream>
#include <math.h>

#include "utils.hpp"


PathsIndex::PathsIndex(const gbwt::GBWT & gbwt_index, const vg::Graph & graph) : index_(gbwt_index) {

    node_lengths = vector<int32_t>(graph.node_size() + 1, -1);
    uint32_t max_node_id = 0;

    for (auto & node: graph.node()) {

        auto id = node.id();
        max_node_id = max(max_node_id, static_cast<uint32_t>(id));

        while (id >= node_lengths.size()) {

            node_lengths.resize(node_lengths.size() * 2, -1);
        }

        assert(node_lengths.at(id) == -1);
        node_lengths.at(node.id()) = node.sequence().size();
    }

    assert(node_lengths.size() > max_node_id);
    node_lengths.resize(max_node_id + 1);
}

PathsIndex::PathsIndex(const gbwt::GBWT & gbwt_index, const handlegraph::HandleGraph & graph) : index_(gbwt_index) {

    node_lengths = vector<int32_t>(graph.get_node_count() + 1, -1);
    uint32_t max_node_id = 0;

    assert(graph.for_each_handle([&](const handlegraph::handle_t & handle) {

        auto id = graph.get_id(handle);
        max_node_id = max(max_node_id, static_cast<uint32_t>(id));

        while (id >= node_lengths.size()) {

            node_lengths.resize(node_lengths.size() * 2, -1);
        }

        assert(node_lengths.at(id) == -1);
        node_lengths.at(id) = graph.get_length(handle);
    }));

    assert(node_lengths.size() > max_node_id);
    node_lengths.resize(max_node_id + 1);
} 

const gbwt::GBWT & PathsIndex::index() const {

    return index_;
}

bool PathsIndex::hasNodeId(const uint32_t node_id) const {

    if (node_id >= node_lengths.size()) {

        return false;
    }

    return (node_lengths.at(node_id) != -1);
}
        
uint32_t PathsIndex::nodeLength(const uint32_t node_id) const {

    assert(hasNodeId(node_id));
    return node_lengths.at(node_id);
}

vector<gbwt::size_type> PathsIndex::locatePathIds(const gbwt::SearchState & search) const {

    auto path_ids = index_.locate(search);

    if (index_.bidirectional()) {

        for (auto & id: path_ids) {

            id = gbwt::Path::id(id);
        }
    }

    return path_ids;
}

string PathsIndex::pathName(const uint32_t path_id) const {

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

uint32_t PathsIndex::pathLength(const uint32_t path_id) const {

    uint32_t path_length = 0;
    
    for (auto & node: index_.extract(path_id)) {

        path_length += nodeLength(gbwt::Node::id(node));
    }

    return path_length;
}

double PathsIndex::effectivePathLength(const uint32_t path_id, const FragmentLengthDist & fragment_length_dist) const {

    const uint32_t path_length = pathLength(path_id); 

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


