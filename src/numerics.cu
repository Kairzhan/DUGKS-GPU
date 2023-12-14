#include <iostream>
#include <string>

#include "precision.h"
#include "lbm_d3q19.h"
#include "config.h"

//***************************************************************************
// Numerical parameters which control the simulation are stored/
// computed here
//***************************************************************************
namespace Numerics {
    // String name for predefined cases
    std::string predef_case;
    
    // Simulation domain size
    real xwidth;
    real ywidth;
    real zwidth;

    // Mesh resolution
    int nx, ny, nz;

    // Cell sizes
    real dx, dy, dz;

    // CFL number and timestep
    real CFL, dt;

    // Subdomain mesh resolution
    int lx, ly, lz;

    // Timestepping
    int istep0, nsteps;

    // Monitoring point (thread 0)
    int imon, jmon, kmon;

    // Output and save frequency
    // noutput3 (3D) and poutput (particles)
    // will be used in future extensions
    int noutput2;
    int noutput3;
    int poutput;
    int ndump;
    int ndiag, nstat;

    // Restart flag (0 - do not restart, 1 - do restart)
    int restart;
    // Load restart data from this timestep:
    // input filename would be: ddf00[nrestart].dat
    int nrestart;

    // Set to 1 to activate force perturbations.
    // Used to bootstrap turbulence for channel
    // flow simulation.
    int perturb;
    
    real visc;
    real ambient_fx, ambient_fy, ambient_fz;
    
    void init_numerics()
    {
	Config config;

        imon=1;
        jmon=1;
        kmon=1;

        predef_case=config.get("case");

        std::cout << "Predefined case: #" << predef_case << "#\n";
        
	CFL=std::stof(config.get("CFL"));
	xwidth=std::stof(config.get("xwidth"));
	ywidth=std::stof(config.get("ywidth"));
	zwidth=std::stof(config.get("zwidth"));

	nx=std::stoi(config.get("nx"));
	ny=std::stoi(config.get("ny"));
	nz=std::stoi(config.get("nz"));

	istep0=std::stoi(config.get("istep0"));
	nsteps=std::stoi(config.get("nsteps"));

	noutput2=std::stoi(config.get("noutput2"));
	noutput3=std::stoi(config.get("noutput3"));
	ndump=std::stoi(config.get("ndump"));

	restart=std::stoi(config.get("restart"));
	nrestart=std::stoi(config.get("nrestart"));

	perturb=std::stoi(config.get("perturb"));

	visc=std::stof(config.get("visc"));

        ambient_fx=std::stof(config.get("ambient_fx"));
        ambient_fy=std::stof(config.get("ambient_fy"));
        ambient_fz=std::stof(config.get("ambient_fz"));

	dx=xwidth/nx;
	dy=ywidth/ny;
	dz=zwidth/nz;
    
	dt=CFL*dx/sqrt(6*DVM::RT);
    }
}
