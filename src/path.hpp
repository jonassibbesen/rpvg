
#ifndef FERSKEN_SRC_PATHINFO_HPP
#define FERSKEN_SRC_PATHINFO_HPP

#include <vector>

#include "utils.hpp"

using namespace std;


struct Path {
        
    string name;
    string origin;
    uint32_t length;
    double effective_length;
    
    Path() {

    	name = "";
    	origin = "";
        length = 0;
        effective_length = 0;
    }
};


#endif
