#ifndef NUMERICS_H_
#define NUMERICS_H_

#include "precision.h"

namespace Numerics {
    extern std::string predef_case;
    
    extern real xwidth, ywidth, zwidth;
    extern int nx, ny, nz;
    extern int lx, ly, lz;
    extern int istep0, nsteps;

    extern real dt;
    extern real dx, dy, dz;

    extern int noutput2;
    extern int noutput3;
    extern int ndump;

    extern int restart;
    extern int nrestart;
    extern int perturb;

    extern real visc;

    extern real ambient_fx, ambient_fy, ambient_fz;
    
    void init_numerics(void );

}
 
#endif
