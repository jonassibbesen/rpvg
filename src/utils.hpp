
#ifndef VGPROB_UTILS_HPP
#define VGPROB_UTILS_HPP

#include <string>
#include <vector>
#include <limits>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <assert.h>

#include <google/protobuf/util/json_util.h>
#include <vg/io/basic_stream.hpp>

using namespace std;


//------------------------------------------------------------------------------

/*
The following code was copied and modified from https://github.com/vgteam/vg
*/

// Convert integer Phred quality score to probability of wrongness.
inline double phred_to_prob(int phred) {
    return pow(10, -((double)phred) / 10);
}

// normal pdf, from http://stackoverflow.com/a/10848293/238609
template <typename T>
inline T normal_pdf(T x, T m, T s)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return inv_sqrt_2pi / s * exp(-T(0.5) * a * a);
}

/// Convert vg::Mapping to gbwt::node_type.
inline gbwt::node_type mapping_to_gbwt(const vg::Mapping & mapping) {
    return gbwt::Node::encode(mapping.position().node_id(), mapping.position().is_reverse());
}

inline int mapping_to_length(const vg::Mapping & m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const vg::Edit& e = m.edit(i);
        l += e.to_length();
    }
    return l;
}

inline int mapping_from_length(const vg::Mapping & m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const vg::Edit& e = m.edit(i);
        l += e.from_length();
    }
    return l;

}

// Note that edit sequences are not reverse complemented.
// Original function in vg repo: reverse_complement_mapping().
inline vg::Mapping lazy_reverse_complement_mapping(const vg::Mapping& mapping,
                                   const function<int64_t(id_t)> & node_length) {
    // Make a new reversed mapping
    vg::Mapping mapping_rc = mapping;

    // switching around to the reverse strand requires us to change offsets
    // that are nonzero to count the unused bases on the other side of the block
    // of used bases.
    if(mapping.has_position() && mapping.position().node_id() != 0) {
        vg::Position* p = mapping_rc.mutable_position();
        
        // How many node bases are used by the mapping?
        size_t used_bases = mapping_from_length(mapping);
        // How many are taken up by the offset on the other strand?
        size_t unused_bases_after = p->offset();
        // The remainder ought to be taken up by the offset on this strand.
        size_t unused_bases_before = node_length(p->node_id()) - used_bases - unused_bases_after;
            
        // Adopt the new offset
        p->set_offset(unused_bases_before);
        // Toggle the reversed-ness flag
        p->set_is_reverse(!p->is_reverse());
    }

    // Clear out all the edits. TODO: we wasted time copying them
    mapping_rc.clear_edit();

    for (int64_t i = mapping.edit_size() - 1; i >= 0; i--) {
        // For each edit in reverse order, put it in reverse complemented
        *mapping_rc.add_edit() = mapping.edit(i);
    }

    return mapping_rc;
}

// Reverse complements path. Note that edit sequences are not reverse complemented. 
// Original function in vg repo:: reverse_complement_path().
inline vg::Path lazy_reverse_complement_path(const vg::Path& path,
                             const function<int64_t(id_t)> & node_length) {

    vg::Path path_rc;

    for(int64_t i = path.mapping_size() - 1; i >= 0; i--) {
        // For each mapping in reverse order, put it in reverse complemented and
        // measured from the other end of the node.
        *path_rc.add_mapping() = lazy_reverse_complement_mapping(path.mapping(i), node_length);
    }
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        path_rc.mutable_mapping(i)->set_rank(i+1);
    }

    return path_rc;
}

// Reverse complements subpats. Note that edit sequences are not reverse complemented.
// Original name in vg repo:: rev_comp_multipath_alignment().
inline google::protobuf::RepeatedPtrField<vg::Subpath> lazy_reverse_complement_subpaths(const google::protobuf::RepeatedPtrField<vg::Subpath> & subpaths, google::protobuf::RepeatedField<google::protobuf::uint32> * subpath_starts, const function<int64_t(id_t)> & node_length) {

    google::protobuf::RepeatedPtrField<vg::Subpath> subpaths_rc;
    
    vector<vector<size_t> > reverse_edge_lists(subpaths.size());
    vector<size_t> reverse_starts;
    
    // add subpaths in reverse order to maintain topological ordering
    for (int64_t i = subpaths.size() - 1; i >= 0; i--) {
        const vg::Subpath & subpath = subpaths[i];
        vg::Subpath * rc_subpath = subpaths_rc.Add();
        
        *rc_subpath->mutable_path() = lazy_reverse_complement_path(subpath.path(), node_length);
        rc_subpath->set_score(subpath.score());

        if (subpath.next_size() > 0) {
            // collect edges by their target (for reversing)
            for (size_t j = 0; j < subpath.next_size(); j++) {
                reverse_edge_lists[subpath.next(j)].push_back(i);
            }
        }
        else {
            // sink subpaths become sources in reverse
            reverse_starts.push_back(i);
        }
    }
    
    // add reversed edges
    for (size_t i = 0; i < subpaths.size(); i++) {
        vg::Subpath * rc_subpath = &(subpaths_rc[i]);
        vector<size_t>& reverse_edge_list = reverse_edge_lists[subpaths.size() - i - 1];
        for (size_t j = 0; j < reverse_edge_list.size(); j++) {
            rc_subpath->add_next(subpaths.size() - reverse_edge_list[j] - 1);
        }
    }
    
    assert(subpath_starts->empty());
    
    // assume that if the original multipath alignment had its starts labeled they want them
    // labeled in the reverse complement too
    if (subpaths.size() > 0) {
        for (size_t i = 0; i < reverse_starts.size(); i++) {
            google::protobuf::uint32 * subpath_start = subpath_starts->Add();
            *subpath_start = subpaths.size() - reverse_starts[i] - 1;
        }
    }

    return subpaths_rc;
}

inline string pb2json(const google::protobuf::Message &msg) {
    // Set options to preserve field names and not camel case them
    google::protobuf::util::JsonPrintOptions opts;
    opts.preserve_proto_field_names = true;

	string buffer;
    auto status = google::protobuf::util::MessageToJsonString(msg, &buffer, opts);
    
    if (!status.ok()) {
        throw runtime_error("Could not serialize " + msg.GetTypeName() + ": " + status.ToString());
    }
    
    return buffer;
}

//------------------------------------------------------------------------------


// Convert paired mapping quality to probability.
inline double mapqsToProb(const pair<int32_t, int32_t> & mapqs) {

    if (mapqs.first > 0 && mapqs.second > 0) {

        return (1 - (1 - phred_to_prob(mapqs.first)) * (1 - phred_to_prob(mapqs.second)));
    
    } else {

        return 1;        
    }
}

// Precision used when comparing double variables.
static const double double_precision = numeric_limits<double>::epsilon() * 100;

// Compare double variables using above precision.
inline bool doubleCompare(const double a, const double b) {

    assert(isfinite(a));
    assert(isfinite(b));

    return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision));
}

// Get path name from GBWT index using path id.
inline string getPathName(const gbwt::GBWT & paths_index, size_t path_id) {

    if (paths_index.bidirectional()) {

        path_id = gbwt::Path::id(path_id);
    }

    stringstream sstream;

    if (!paths_index.hasMetadata() || !paths_index.metadata.hasPathNames() || paths_index.metadata.paths() <= path_id || !paths_index.metadata.hasSampleNames()) {
        
        sstream << path_id + 1;
    
    } else {

        sstream << paths_index.metadata.sample(paths_index.metadata.path(path_id).sample);
    }

    return sstream.str();
}

#endif






