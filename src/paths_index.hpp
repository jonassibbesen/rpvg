
#ifndef RPVG_SRC_PATHSINDEX_HPP
#define RPVG_SRC_PATHSINDEX_HPP

#include <vector>
#include <string>

#include "gbwt/gbwt.h"
#include "gbwt/fast_locate.h"
#include "handlegraph/handle_graph.hpp"
#include "vg/io/basic_stream.hpp"
#include "fragment_length_dist.hpp"

using namespace std;


class PathsIndex {

    public: 
    	
        PathsIndex(const gbwt::GBWT & gbwt_index_in, const gbwt::FastLocate & r_index_in, const vg::Graph & graph);
        PathsIndex(const gbwt::GBWT & gbwt_index_in, const gbwt::FastLocate & r_index_in, const handlegraph::HandleGraph & graph);

        uint32_t numberOfNodes() const;
        bool hasNodeId(const uint32_t node_id) const;
        uint32_t nodeLength(const uint32_t node_id) const;

        vector<gbwt::edge_type> edges(const gbwt::node_type gbwt_node) const;

        bool bidirectional() const;
        uint32_t numberOfPaths() const;

        void find(pair<gbwt::SearchState, gbwt::size_type> * gbwt_search, const gbwt::node_type gbwt_node) const;
        void extend(pair<gbwt::SearchState, gbwt::size_type> * gbwt_search, const gbwt::node_type gbwt_node) const;
        vector<gbwt::size_type> locatePathIds(const pair<gbwt::SearchState, gbwt::size_type> & gbwt_search) const;

        string pathName(const uint32_t path_id) const;
        uint32_t pathLength(uint32_t path_id) const;
        double effectivePathLength(const uint32_t path_id, const FragmentLengthDist & fragment_length_dist) const;

    private:

        const gbwt::GBWT & gbwt_index;
        const gbwt::FastLocate & r_index;

        vector<int32_t> node_lengths;

        double calculateLowerPhi(const double value) const;
        double calculateUpperPhi(const double value) const;
};


#endif
