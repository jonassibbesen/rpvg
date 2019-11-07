
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
    
    for (auto & node: index_.extract(gbwt::Path::encode(path_id, false))) {

        path_length += nodeLength(gbwt::Node::id(node));
    }

    return path_length;
}

int32_t PathsIndex::effectivePathLength(const int32_t path_id, const FragmentLengthDist & fragment_length_dist) const {

    const int32_t path_length = pathLength(path_id); 

    if (path_length == 0) {

        return 0;
    }

    const double pi = acos(-1);

    // https://en.wikipedia.org/wiki/Truncated_normal_distribution
    const double beta = (path_length - fragment_length_dist.mean()) / fragment_length_dist.sd();
    const double lower_phi = exp(-0.5 * pow(beta, 2)) / sqrt(2 * pi);
    const double upper_phi = 0.5 * (1 + erf(beta / sqrt(2)));

    const double trunc_frag_length_mean = fragment_length_dist.mean() - fragment_length_dist.sd() * lower_phi / upper_phi;

    return max(0.0, path_length - trunc_frag_length_mean);
}

