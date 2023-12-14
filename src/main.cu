#include <iostream>
#include <iomanip>
#include <mpi.h>

#include "precision.h"
#include "parallel.h"
#include "cuda.h"
#include "numerics.h"
#include "checks.h"
#include "dugks.h"
#include "output.h"

int main(int argc, char* argv[])
{
    // Step 1: Perform initializations, readin parameters.
    Numerics::init_numerics();
    init_parallel(argc, argv);
    init_cuda(myrank);
    perform_checks();
    init_dugks();

    // Step 2: Start the simulation
    timestart1=MPI_Wtime();
    
    // Main timeloop
    for (int istep=Numerics::istep0; istep<=Numerics::istep0+Numerics::nsteps-1; istep++) {
        if (root()) {
            std::cout << "istep = " << std::setw(8) << istep << std::endl;
        }

        dugks_step(istep);

        // output, dump, postprocess
        if (istep==Numerics::istep0 || istep%Numerics::noutput2==0) {
            upload_macro();
            copy_vel();
            output_2d(istep, h_rho, h_ux, h_uy, h_uz);
        }

        if (istep%Numerics::ndump==0) {
            copy_ddf();
            dump(istep, h_f);
        }
    }

    // Step 3: Output runtime details and free memory and MPI
    timeend1=MPI_Wtime();

    if (root()) {
        std::cout << "Main loop time: " << timeend1-timestart1 << " seconds." << std::endl;
    }

    free_parallel();

    return 0;
}
