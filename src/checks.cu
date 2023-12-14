#include <mpi.h>
#include <cassert>

#include "precision.h"
#include "numerics.h"
#include "parallel.h"

//===========================================================================
// Performs necessary parameter checks to avoid meaningless simulation
// cases. For example, this function should check whether nx is
// divisible by dims[0], to be sure that each subdomain has the same size.
//===========================================================================
void perform_checks(void)
{
    using namespace Numerics;
    
    assert(nx%dims[0]==0);
    assert(ny%dims[1]==0);
    assert(nz%dims[2]==0);
    
    return;
}
