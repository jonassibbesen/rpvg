
#ifndef RPVG_SRC_PATHINFO_HPP
#define RPVG_SRC_PATHINFO_HPP

#include <vector>

#include "utils.hpp"

using namespace std;


struct PathInfo {
        
    string name;
    string origin;
    uint32_t length;
    double effective_length;
    
    PathInfo() {

    	name = "";
    	origin = "";
        length = 0;
        effective_length = 0;
    }
};


#endif
