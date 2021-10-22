
#include "paths_index.hpp"

#include <sstream>
#include <math.h>

#include "utils.hpp"


PathsIndex::PathsIndex(const gbwt::GBWT & gbwt_index_in, const gbwt::FastLocate & r_index_in, const vg::Graph & graph) : gbwt_index(gbwt_index_in), r_index(r_index_in) {

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

PathsIndex::PathsIndex(const gbwt::GBWT & gbwt_index_in, const gbwt::FastLocate & r_index_in,  const handlegraph::HandleGraph & graph) : gbwt_index(gbwt_index_in), r_index(r_index_in) {

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

uint32_t PathsIndex::numberOfNodes() const {

    return node_lengths.size();
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

vector<gbwt::edge_type> PathsIndex::edges(const gbwt::node_type gbwt_node) const {

    return gbwt_index.edges(gbwt_node);
}

bool PathsIndex::bidirectional() const {

    return gbwt_index.bidirectional();
}

uint32_t PathsIndex::numberOfPaths() const {

    if (bidirectional()) {

        assert(gbwt_index.sequences() % 2 == 0);
        return (gbwt_index.sequences() / 2);

    } else {

        return gbwt_index.sequences();
    }
}

void PathsIndex::find(pair<gbwt::SearchState, gbwt::size_type> * gbwt_search, const gbwt::node_type gbwt_node) const {

    if (r_index.empty()) {

        gbwt_search->first = gbwt_index.find(gbwt_node);

    } else {

        gbwt_search->first = r_index.find(gbwt_node, gbwt_search->second);
    }
}

void PathsIndex::extend(pair<gbwt::SearchState, gbwt::size_type> * gbwt_search, const gbwt::node_type gbwt_node) const {

    if (r_index.empty()) {

        gbwt_search->first = gbwt_index.extend(gbwt_search->first, gbwt_node);

    } else {

        gbwt_search->first = r_index.extend(gbwt_search->first, gbwt_node, gbwt_search->second);        
    }
}

vector<gbwt::size_type> PathsIndex::locatePathIds(const pair<gbwt::SearchState, gbwt::size_type> & gbwt_search) const {

    vector<gbwt::size_type> path_ids;

    if (r_index.empty()) {

        path_ids = gbwt_index.locate(gbwt_search.first);

    } else {

        path_ids = r_index.locate(gbwt_search.first, gbwt_search.second);        
    }

    if (bidirectional()) {

        for (auto & id: path_ids) {

            id = gbwt::Path::id(id);
        }
    }

    return path_ids;
}

string PathsIndex::pathName(const uint32_t path_id) const {

    stringstream sstream;

    if (!gbwt_index.hasMetadata() || !gbwt_index.metadata.hasPathNames() || gbwt_index.metadata.paths() <= path_id || !gbwt_index.metadata.hasSampleNames()) {
        
        sstream << path_id + 1;
    
    } else {

        const gbwt::PathName& path_name = gbwt_index.metadata.path(path_id);

        sstream << gbwt_index.metadata.sample(path_name.sample);

        if (gbwt_index.metadata.hasContigNames()) {

            sstream << "_" << gbwt_index.metadata.contig(path_name.contig);
            sstream << "_" << path_name.phase;
            sstream << "_" << path_name.count;
        }
    }

    return sstream.str();
}

uint32_t PathsIndex::pathLength(uint32_t path_id) const {

    if (bidirectional()) {

        path_id = gbwt::Path::encode(path_id, false);
    }

    uint32_t path_length = 0;
    
    for (auto & node: gbwt_index.extract(path_id)) {

        path_length += nodeLength(gbwt::Node::id(node));
    }

    return path_length;
}

double PathsIndex::effectivePathLength(const uint32_t path_id, const FragmentLengthDist & fragment_length_dist) const {

    const uint32_t path_length = pathLength(path_id); 

    if (path_length == 0) {

        return 0;
    }
    
    double trunc_fragment_length_mean = 0.0;
    if (Utils::doubleCompare(fragment_length_dist.shape(), 0.0)) {
        // https://en.wikipedia.org/wiki/Truncated_normal_distribution
        const double alpha = (1.0 - fragment_length_dist.loc()) / fragment_length_dist.scale();
        const double beta = (path_length - fragment_length_dist.loc()) / fragment_length_dist.scale();
        
        trunc_fragment_length_mean = fragment_length_dist.loc() + fragment_length_dist.scale() * (calculateLowerPhi(alpha) - calculateLowerPhi(beta)) / (calculateUpperPhi(beta) - calculateUpperPhi(alpha));
    }
    else {
        trunc_fragment_length_mean = Utils::truncated_skew_normal_expected_value<double>(fragment_length_dist.loc(),
                                                                                         fragment_length_dist.scale(),
                                                                                         fragment_length_dist.shape(),
                                                                                         1.0, path_length);
    }
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

