
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


template<class T, class T2>
inline ostream & operator<<(ostream & os, const pair<T,T2> & values) {

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

namespace Utils {

    typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, 1, Eigen::ColMajor> ColVectorXui;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> ColVectorXd;

    typedef Eigen::Matrix<uint32_t, 1, Eigen::Dynamic, Eigen::RowMajor> RowVectorXui;
    typedef Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> RowVectorXd;
    
    typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMatrixXb;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMatrixXd;

    typedef Eigen::SparseMatrix<bool, Eigen::ColMajor> ColSparseMatrixXb;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor> ColSparseMatrixXd;

    inline vector<string> splitString(const string & str, const char delim) {

        stringstream ss(str);
        vector<string> elems;

        for (string item; getline(ss, item, delim);) {

            elems.push_back(item);
        }

        return elems;
    }

    static const double score_log_base = 1.383325268738;

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

    //------------------------------------------------------------------------------

    /*
    All the following code have been copied and modified from https://github.com/vgteam/vg
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

    inline char quality_short_to_char(short i) {
        return static_cast<char>(i + 33);
    }

    inline string string_quality_short_to_char(const string& quality) {
        string buffer; buffer.resize(quality.size());
        for (int i = 0; i < quality.size(); ++i) {
            buffer[i] = quality_short_to_char(quality[i]);
        }
        return buffer;
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

        aln_rc.set_sequence(string(aln.sequence().rbegin(), aln.sequence().rend()));
        aln_rc.set_quality(string(aln.quality().rbegin(), aln.quality().rend()));

        aln_rc.set_score(aln.score());
        aln_rc.set_mapping_quality(aln.mapping_quality());

        *aln_rc.mutable_path() = lazy_reverse_complement_path(aln.path(), node_length);
        
        return aln_rc;
    }

    // Reverse complements multipath alignment. Note that sequences, paths and edit sequences 
    // are not reverse complemented. Original name in vg repo: rev_comp_multipath_alignment().
    inline vg::MultipathAlignment lazy_reverse_complement_alignment(const vg::MultipathAlignment& multipath_aln, const function<int64_t(int64_t)>& node_length) {
        
        vg::MultipathAlignment multipath_aln_rc;

        multipath_aln_rc.set_sequence(string(multipath_aln.sequence().rbegin(), multipath_aln.sequence().rend()));
        multipath_aln_rc.set_quality(string(multipath_aln.quality().rbegin(), multipath_aln.quality().rend()));

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

    static const int8_t default_match = 1;
    static const int8_t default_mismatch = 4;
    static const int8_t default_full_length_bonus = 5;

    static const int8_t score_matrix[16] = {
         default_match,    -default_mismatch, -default_mismatch, -default_mismatch,
        -default_mismatch,  default_match,    -default_mismatch, -default_mismatch,
        -default_mismatch, -default_mismatch,  default_match,    -default_mismatch,
        -default_mismatch, -default_mismatch, -default_mismatch,  default_match
    };

    inline vector<int8_t> qual_adjusted_matrix(double gc_content = 0.5, uint32_t max_qual = 255) {
        
        // TODO: duplicative with GSSWAligner()
        double* nt_freqs = (double*) malloc(sizeof(double) * 4);
        nt_freqs[0] = 0.5 * (1 - gc_content);
        nt_freqs[1] = 0.5 * gc_content;
        nt_freqs[2] = 0.5 * gc_content;
        nt_freqs[3] = 0.5 * (1 - gc_content);
        
        // recover the emission probabilities of the align state of the HMM
        double* align_prob = (double*) malloc(sizeof(double) * 16);
        
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                align_prob[i * 4 + j] = (exp(score_log_base * score_matrix[i * 4 + j])
                                         * nt_freqs[i] * nt_freqs[j]);
            }
        }
        
        // compute the sum of the emission probabilities under a base error
        double* align_complement_prob = (double*) malloc(sizeof(double) * 16);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                align_complement_prob[i * 4 + j] = 0.0;
                for (int k = 0; k < 4; k++) {
                    if (k != j) {
                        align_complement_prob[i * 4 + j] += align_prob[i * 4 + k];
                    }
                }
            }
        }
        
        // quality score of random guessing
        int lowest_meaningful_qual = ceil(-10.0 * log10(0.75));
        
        // compute the adjusted alignment scores for each quality level
        vector<int8_t> qual_adj_mat(25 * (max_qual + 1));
        for (int q = 0; q <= max_qual; q++) {
            double err = pow(10.0, -q / 10.0);
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) {
                    int8_t score;
                    if (i == 4 || j == 4 || q < lowest_meaningful_qual) {
                        score = 0;
                    }
                    else {
                        score = round(log(((1.0 - err) * align_prob[i * 4 + j] + (err / 3.0) * align_complement_prob[i * 4 + j])
                                          / (nt_freqs[i] * ((1.0 - err) * nt_freqs[j] + (err / 3.0) * (1.0 - nt_freqs[j])))) / score_log_base);
                    }
                    qual_adj_mat.at(q * 25 + i * 5 + j) = round(score);
                }
            }
        }
        
        free(align_complement_prob);
        free(align_prob);
        free(nt_freqs);
        
        return qual_adj_mat;
    }

    inline vector<int8_t> qual_adjusted_bonuses(uint32_t max_qual = 255) {
        
        
        double p_full_len = exp(score_log_base * default_full_length_bonus) / (1.0 + exp(score_log_base * default_full_length_bonus));
        
        vector<int8_t> qual_adj_bonuses(max_qual + 1);
        
        int lowest_meaningful_qual = ceil(-10.0 * log10(0.75));
        // hack because i want the minimum qual value from illumina (2) to have zero score, but phred
        // values are spaced out in a way to approximate this singularity well
        ++lowest_meaningful_qual;
        
        for (int q = lowest_meaningful_qual; q <= max_qual; ++q) {
            double err = pow(10.0, -q / 10.0);
            double score = log(((1.0 - err * 4.0 / 3.0) * p_full_len + (err * 4.0 / 3.0) * (1.0 - p_full_len)) / (1.0 - p_full_len)) / score_log_base;
            qual_adj_bonuses[q] = round(score);
        }
        
        return qual_adj_bonuses;
    }

    static const vector<int8_t> qual_score_matrix = qual_adjusted_matrix();
    static const vector<int8_t> qual_full_length_bonuses = qual_adjusted_bonuses();

    //------------------------------------------------------------------------------
}


#endif
