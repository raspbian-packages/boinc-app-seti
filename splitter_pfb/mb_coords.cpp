#include "setilib.h"
#include "mb_splitter.h"

void get_coord_history() {
// insert telescope coordinates into the coordinate history.
// this should be converted to a more accurate routine.

    std::vector<dr2_compact_block_t>::iterator i=tapebuffer.begin();

    for (;i!=tapebuffer.end();i++) {
        coord_history[i->header.coord_time].ra   = i->header.ra;
        coord_history[i->header.coord_time].dec  = i->header.dec;
        coord_history[i->header.coord_time].time = i->header.coord_time.jd().uval();
    }
}
