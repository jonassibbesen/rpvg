
#ifndef RPVG_SRC_PATHSINDEX_HPP
#define RPVG_SRC_PATHSINDEX_HPP

#include <vector>
#include <string>

#include "gbwt/gbwt.h"
#include "handlegraph/handle_graph.hpp"
#include "vg/io/basic_stream.hpp"
#include "fragment_length_dist.hpp"

using namespace std;


class PathsIndex {

    public: 
    	
        PathsIndex(const gbwt::GBWT & gbwt_index, const vg::Graph & graph);
        PathsIndex(const gbwt::GBWT & gbwt_index, const handlegraph::HandleGraph & graph);

        const gbwt::GBWT & index() const;

        bool hasNodeId(const uint32_t node_id) const;
        uint32_t nodeLength(const uint32_t node_id) const;

        string pathName(const uint32_t path_id) const;
        uint32_t pathLength(uint32_t path_id) const;
        double effectivePathLength(const uint32_t path_id, const FragmentLengthDist & fragment_length_dist) const;

    private:

        const gbwt::GBWT & index_;
        vector<int32_t> node_lengths;

        double calculateLowerPhi(const double value) const;
        double calculateUpperPhi(const double value) const;
};


#endif
