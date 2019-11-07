
#ifndef RPVG_SRC_PATHSINDEX_HPP
#define RPVG_SRC_PATHSINDEX_HPP

#include <vector>
#include <string>

#include "gbwt/gbwt.h"
#include "vg/io/basic_stream.hpp"
#include "fragment_length_dist.hpp"

using namespace std;


class PathsIndex {

    public: 
    	
        PathsIndex(const gbwt::GBWT & gbwt_index, const vg::Graph & graph);

        const gbwt::GBWT & index() const;
        int32_t nodeLength(const int32_t node_id) const;

        string pathName(const int32_t path_id) const;
        int32_t pathLength(const int32_t path_id) const;
        int32_t effectivePathLength(const int32_t path_id, const FragmentLengthDist & fragment_length_dist) const;

    private:

        const gbwt::GBWT & index_;
        vector<int32_t> node_lengths;
;
};


#endif
