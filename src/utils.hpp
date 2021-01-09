
#ifndef RPVG_SRC_UTILS_HPP
#define RPVG_SRC_UTILS_HPP

#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <limits>
#include <sstream>
#include <algorithm>
#include <cmath>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "google/protobuf/util/json_util.h"
#include "vg/io/basic_stream.hpp"
#include "gbwt/gbwt.h"
#include "handlegraph/handle_graph.hpp"

using namespace std;


namespace Eigen {

    typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, 1, Eigen::ColMajor> ColVectorXui;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> ColVectorXd;

    typedef Eigen::Matrix<uint32_t, 1, Eigen::Dynamic, Eigen::RowMajor> RowVectorXui;
    typedef Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> RowVectorXd;
    
    typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMatrixXb;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMatrixXd;

    typedef Eigen::SparseMatrix<bool, Eigen::ColMajor> ColSparseMatrixXb;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor> ColSparseMatrixXd;
}

//------------------------------------------------------------------------------

/*
The following code have been copied and modified from https://github.com/vgteam/vg
*/

// Convert integer Phred quality score to probability of wrongness.
inline double phred_to_prob(uint32_t phred) {
    return pow(10, -((double)phred) / 10);
}

// Convert probability of wrongness to integer Phred quality score.
inline double prob_to_phred(double prob) {
    return -10.0 * log10(prob);
}

// log normal pdf, from http://stackoverflow.com/a/10848293/238609
template <typename T>
inline T log_normal_pdf(T x, T m, T s)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return log(inv_sqrt_2pi) - log(s) - T(0.5) * a * a;
}

/*
 * Return the log of the sum of two log-transformed values without taking them
 * out of log space.
 */
inline double add_log(double log_x, double log_y) {
    return log_x > log_y ? log_x + log1p(exp(log_y - log_x)) : log_y + log1p(exp(log_x - log_y));
}

/// Convert vg::Mapping to gbwt::node_type.
inline gbwt::node_type mapping_to_gbwt(const vg::Mapping & mapping) {
    return gbwt::Node::encode(mapping.position().node_id(), mapping.position().is_reverse());
}

inline uint32_t mapping_to_length(const vg::Mapping & m) {
    uint32_t l = 0;
    for (uint32_t i = 0; i < m.edit_size(); ++i) {
        const vg::Edit& e = m.edit(i);
        l += e.to_length();
    }
    return l;
}

inline uint32_t mapping_from_length(const vg::Mapping & m) {
    uint32_t l = 0;
    for (uint32_t i = 0; i < m.edit_size(); ++i) {
        const vg::Edit& e = m.edit(i);
        l += e.from_length();
    }
    return l;
}

// Note that edit sequences are not reverse complemented.
// Original function in vg repo: reverse_complement_mapping().
inline vg::Mapping lazy_reverse_complement_mapping(const vg::Mapping& mapping,
                                   const function<int64_t(int64_t)> & node_length) {
    // Make a new reversed mapping
    vg::Mapping mapping_rc;
    *mapping_rc.mutable_position() = mapping.position();

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

    for (int64_t i = mapping.edit_size() - 1; i >= 0; i--) {
        // For each edit in reverse order, put it in reverse complemented
        *mapping_rc.add_edit() = mapping.edit(i);
    }

    return mapping_rc;
}

// Reverse complements path. Note that edit sequences are not reverse complemented. 
// Original function in vg repo: reverse_complement_path().
inline vg::Path lazy_reverse_complement_path(const vg::Path& path,
                             const function<int64_t(int64_t)> & node_length) {

    vg::Path path_rc;

    for(int64_t i = path.mapping_size() - 1; i >= 0; i--) {
        // For each mapping in reverse order, put it in reverse complemented and
        // measured from the other end of the node.
        *path_rc.add_mapping() = lazy_reverse_complement_mapping(path.mapping(i), node_length);
    }

    return path_rc;
}

// Reverse complements alignment. Note that sequences, paths and edit sequences 
// are not reverse complemented. Original function in vg repo: reverse_complement_alignment().
inline vg::Alignment lazy_reverse_complement_alignment(const vg::Alignment& aln,
                                       const function<int64_t(int64_t)>& node_length) {
    // We're going to reverse the alignment and all its mappings.
    // TODO: should we/can we do this in place?
    
    vg::Alignment aln_rc;

    aln_rc.set_score(aln.score());
    aln_rc.set_mapping_quality(aln.mapping_quality());

    *aln_rc.mutable_path() = lazy_reverse_complement_path(aln.path(), node_length);
    
    return aln_rc;
}


// Reverse complements multipath alignment. Note that sequences, paths and edit sequences 
// are not reverse complemented. Original name in vg repo: rev_comp_multipath_alignment().
inline vg::MultipathAlignment lazy_reverse_complement_alignment(const vg::MultipathAlignment& multipath_aln, const function<int64_t(int64_t)>& node_length) {
    
    vg::MultipathAlignment multipath_aln_rc;

    multipath_aln_rc.set_mapping_quality(multipath_aln.mapping_quality());

    vector<vector<size_t> > reverse_edge_lists(multipath_aln.subpath_size());
    vector<vector<pair<size_t, int32_t> > > reverse_connection_lists(multipath_aln.subpath_size());

    vector<size_t> reverse_starts;
    
    // remove subpaths to avoid duplicating
    // multipath_aln_rc.clear_subpath();
    
    // add subpaths in reverse order to maintain topological ordering
    for (int64_t i = multipath_aln.subpath_size() - 1; i >= 0; i--) {
        const vg::Subpath& subpath = multipath_aln.subpath(i);
        vg::Subpath* rc_subpath = multipath_aln_rc.add_subpath();
        
        *(rc_subpath->mutable_path()) = lazy_reverse_complement_path(subpath.path(), node_length);
        rc_subpath->set_score(subpath.score());

        if (subpath.next_size() > 0 || subpath.connection_size() > 0) {
            // collect edges by their target (for reversing)
            for (size_t j = 0; j < subpath.next_size(); j++) {
                reverse_edge_lists[subpath.next(j)].push_back(i);
            }
            for (auto & connection: subpath.connection()) {
                reverse_connection_lists[connection.next()].emplace_back(i, connection.score());
            }
        }
        else {
            // sink subpaths become sources in reverse
            reverse_starts.push_back(i);
        }
    }
    
    // add reversed edges
    for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
        vg::Subpath* rc_subpath = multipath_aln_rc.mutable_subpath(i);
        vector<size_t>& reverse_edge_list = reverse_edge_lists[multipath_aln.subpath_size() - i - 1];
        for (size_t j = 0; j < reverse_edge_list.size(); j++) {
            rc_subpath->add_next(multipath_aln.subpath_size() - reverse_edge_list[j] - 1);
        }
        vector<pair<size_t, int32_t> >& reverse_connection_list = reverse_connection_lists[multipath_aln.subpath_size() - i - 1];
        for (size_t j = 0; j < reverse_connection_list.size(); ++j) {
            vg::Connection* connection = rc_subpath->add_connection();
            connection->set_next(multipath_aln.subpath_size() - reverse_connection_list[j].first - 1);
            connection->set_score(reverse_connection_list[j].second);
        }
    }
    
    // remove start nodes that are invalid in reverse
    // multipath_aln_rc.clear_start();
    
    // assume that if the original multipath alignment had its starts labeled they want them
    // labeled in the reverse complement too
    if (multipath_aln.start_size() > 0) {
        for (size_t i = 0; i < reverse_starts.size(); i++) {
            multipath_aln_rc.add_start(multipath_aln.subpath_size() - reverse_starts[i] - 1);
        }
    }

    return multipath_aln_rc;
}

inline string pb2json(const google::protobuf::Message &msg) {

	string buffer;
    auto status = google::protobuf::util::MessageToJsonString(msg, &buffer);
    
    if (!status.ok()) {
        throw runtime_error("Could not serialize " + msg.GetTypeName() + ": " + status.ToString());
    }
    
    return buffer;
}

inline void json2pb(google::protobuf::Message &msg, const string& buf) {
    auto status = google::protobuf::util::JsonStringToMessage(buf, &msg);
    
    if (!status.ok()) {
        // This generally will happen if someone feeds in the wrong type of JSON.
        // TODO: It would be nice to be able to find the neme of the offending non-existent field.
        throw runtime_error("Could not deserialize " + msg.GetTypeName() + ": " + status.ToString());
    }
}


//------------------------------------------------------------------------------


inline vector<string> splitString(const string & str, const char delim) {

    stringstream ss(str);
    vector<string> elems;

    for (string item; getline(ss, item, delim);) {

        elems.push_back(item);
    }

    return elems;
}

static const double score_log_base = 1.383325268738;

inline uint32_t scorePrecision(const double prob_precision) {

    return ceil(-1 * log(prob_precision) / score_log_base);
}

// Precision used when comparing double variables.
static const double double_precision = numeric_limits<double>::epsilon() * 100;

// Compare double variables using above precision.
inline bool doubleCompare(const double a, const double b) {

    assert(isfinite(a));
    assert(isfinite(b));

    return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision));
}

inline uint32_t numPermutations(vector<uint32_t> values) {

    assert(!values.empty());

    if (values.size() == 1) {

        return 1;
    }

    sort(values.begin(), values.end());

    uint32_t num_unique_values = 1;

    for (size_t i = 1; i < values.size(); ++i) {

        if (values.at(i - 1) != values.at(i)) {

            num_unique_values++;
        }
    }

    return (tgamma(values.size() + 1) / tgamma(values.size() - num_unique_values + 2));
}

template<class T>
inline ostream & operator<<(ostream & os, const pair<T,T> & values) {

    os << "(" << values.first << "," << values.second << ")";
    return os;
}

template<class T>
inline ostream & operator<<(ostream & os, const vector<T> & values) {

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


#endif
