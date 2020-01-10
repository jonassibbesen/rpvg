/**
 * \file register_libvg_io.hpp
 * Includes calls to register all libvg types with libvgio.
 */

// Keep these includes in alphabetical order.

#include "register_loader_saver_gbwt.hpp"
#include "register_loader_saver_xg.hpp"

#include "register_libvg_io.hpp"


namespace vg {

namespace io {

using namespace std;

bool register_libvg_io() {
    register_loader_saver_gbwt();
    register_loader_saver_xg();
    return true;
}
    
}

}
