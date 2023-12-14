#include <iostream>
#include <cstring>
#include <mpi.h>

#include "precision.h"
#include "constants.h"
#include "numerics.h"
#include "lbm_d3q19.h"
#include "parallel.h"
#include "output.h"
#include "cuda.h"

// Multidimensional arrays are emulated throught
// C macros with Fortran array layout
#define F(ip,i,j,k)       f[(k)*(ly+2)*(lx+2)*NPOP + (j)*(lx+2)*NPOP + (i)*NPOP + ip]
#define FB(ip,i,j,k)     fb[(k)*(ly+2)*(lx+2)*NPOP + (j)*(lx+2)*NPOP + (i)*NPOP + ip]
#define H_F(ip,i,j,k)   h_f[(k)*(ly+2)*(lx+2)*NPOP + (j)*(lx+2)*NPOP + (i)*NPOP + ip]
#define H_FB(ip,i,j,k) h_fb[(k)*(ly+2)*(lx+2)*NPOP + (j)*(lx+2)*NPOP + (i)*NPOP + ip]

#define FORCEX(i,j,k) forcex[(k)*(ly+2)*(lx+2) + (j)*(lx+2) + (i)]
#define FORCEY(i,j,k) forcey[(k)*(ly+2)*(lx+2) + (j)*(lx+2) + (i)]
#define FORCEZ(i,j,k) forcez[(k)*(ly+2)*(lx+2) + (j)*(lx+2) + (i)]

#define FBs(ip, i, j, k) fbs[(k)*(blockDim.y+2)*(blockDim.x+2)*NPOP + (j)*(blockDim.x+2)*NPOP + (i)*NPOP+ip]

// Host variables
real *h_rho, *h_ux, *h_uy, *h_uz;
real *h_f, *h_fb;
real *h_tp;

// Device variables
real *rho, *ux, *uy, *uz;
real *f, *fb;
real *forcex, *forcey, *forcez;

struct __align__(8) DugksParameters {
    real visc;
    real RT;
    real hdt;
    real bw0, bw1, bw2;
    real fw0, fw1, fw2;
    real qw1, qw2;
    real gw1, gw2;
    real *tp;
    int *cixyzo;
};

struct __align__(8) HydroParameters {
    real ambient_fx, ambient_fy, ambient_fz;
};

struct __align__(8) NumericalParameters {
    int npop;
    int lx, ly, lz;
    int nx, ny, nz;
    real dx, dy, dz;
    real dt;
    real xwidth, ywidth, zwidth;
    real zoffset;
    int periods[3];
};

DugksParameters *params;
NumericalParameters *numparams;
HydroParameters *hydroparams;

// predeclare functions
real feq(int ip, real *tp, real r, real u, real v, real w);

__global__ void test_params(NumericalParameters *numparams, DugksParameters *parameters);

__global__ void calc_macro(real r[], real u[], real v[], real w[], const real f[], const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp);
void upload_macro();

__global__ void calc_fb(real fb[], real f[], real r[], real u[], real v[], real w[], const NumericalParameters *n, const DugksParameters *p, int ip);
__global__ void calc_fb1(real fb[], real f[], real forcex[], real forcey[], real forcez[], const NumericalParameters *n, const DugksParameters *p);

__global__ void calc_fx(real f[], const real fb[], const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp);
__global__ void calc_fx_tri(real f[], const real fb[], const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp);
__device__ real trilinear(real a, real b, real c, real d, real e, real f, real g, real h, real mx, real my, real mz);

__global__ void calc_perturb_force(int step, real forcex[], real forcey[], real forcez[], const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp);
__global__ void calc_force(int step, real forcex[], real forcey[], real forcez[], const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp);

__device__  real feq_s(real r, real u, real v, real w, real tp, int cix, int ciy, int ciz);

__global__ void extend_fb_x(real *fb, NumericalParameters *n, DugksParameters *p);
__global__ void extend_fb_y(real *fb, NumericalParameters *n, DugksParameters *p);
__global__ void extend_fb_z(real *fb, NumericalParameters *n, DugksParameters *p);
__global__ void extend_fb_x_wall(real *fb, NumericalParameters *n, DugksParameters *p);

__device__ void f2ruvw(real *r, real *u, real *v, real *w, const real *f, const DugksParameters *p);
__device__ void f2ruvwh(real *r, real *u, real *v, real *w, const real f[], real force1[], const DugksParameters *p, const HydroParameters *hp);

void init_dugks()
{
    using DVM::RT;
    using Numerics::dt;
    using DVM::npop;

    const real visc=Numerics::visc;
    const real tau=visc/RT;
    real hdt=dt/2.;

    // Setup DUGKS parameters
    cudaMallocManaged(&params, sizeof(DugksParameters));
    std::cout << "alloc. success\n";

    params->visc=visc;
    params->RT=RT;
    params->hdt=hdt;

    params->bw0=(3.*tau*hdt)/(2.*tau+dt);
    params->bw1=(2.*tau-hdt)/(2.*tau+dt);
    params->bw2=3.*hdt/(2.*tau+dt);

    params->gw1=2.*tau/(2.*tau+dt);
    params->gw2=dt/(2.*tau+dt);

    params->qw1=(2.*tau+dt)/(2.*tau);
    params->qw2=-dt/(2.*tau);

    params->fw0=tau*hdt/(2.*tau+hdt);
    params->fw1=2.*tau/(2.*tau+hdt);
    params->fw2=hdt/(2.*tau+hdt);

    // set up discrete veloicty vectors
    cudaMallocManaged(&(params->cixyzo), 4*npop*sizeof(int));
    for (int ip=0; ip<npop; ip++) {
        params->cixyzo[ip+0*npop]=DVM::opp[ip];
        params->cixyzo[ip+1*npop]=DVM::cix[ip];
        params->cixyzo[ip+2*npop]=DVM::ciy[ip];
        params->cixyzo[ip+3*npop]=DVM::ciz[ip];
    }

    // set up directional weights
    cudaMallocManaged(&(params->tp), npop*sizeof(real));
    params->tp[0]=DVM::ww0;
    for (int ip=1; ip<=6; ip++)
        params->tp[ip]=DVM::ww1;
    for (int ip=7; ip<npop; ip++)
        params->tp[ip]=DVM::ww2;

    h_tp=(real *)malloc(npop*sizeof(real));
    h_tp[0]=DVM::ww0;
    for (int ip=1; ip<=6; ip++)
        h_tp[ip]=DVM::ww1;
    for (int ip=7; ip<npop; ip++)
        h_tp[ip]=DVM::ww2;

    // Setup numerical parameters
    using Numerics::lx;
    using Numerics::ly;
    using Numerics::lz;

    using Numerics::dx;
    using Numerics::dy;
    using Numerics::dz;

    cudaMallocManaged(&numparams, sizeof(NumericalParameters));
    numparams->npop=npop;
    numparams->lx=lx;
    numparams->ly=ly;
    numparams->lz=lz;
    numparams->dx=dx;
    numparams->dy=dy;
    numparams->dz=dz;
    numparams->dt=dt;

    numparams->nx=Numerics::nx;
    numparams->ny=Numerics::ny;
    numparams->nz=Numerics::nz;
    numparams->xwidth=Numerics::xwidth;
    numparams->ywidth=Numerics::ywidth;
    numparams->zwidth=Numerics::zwidth;

    numparams->periods[0]=periods[0];
    numparams->periods[1]=periods[1];
    numparams->periods[2]=periods[2];
    
    int lx2ly2lz2=(lx+2)*(ly+2)*(lz+2);
    cudaMalloc(&rho,lx2ly2lz2*sizeof(real));
    cudaMalloc(&ux, lx2ly2lz2*sizeof(real));
    cudaMalloc(&uy, lx2ly2lz2*sizeof(real));
    cudaMalloc(&uz, lx2ly2lz2*sizeof(real));
    cudaMalloc( &f, npop*lx2ly2lz2*sizeof(real));
    cudaMalloc(&fb, npop*lx2ly2lz2*sizeof(real));

    cudaMalloc(&forcex, lx2ly2lz2*sizeof(real));
    cudaMalloc(&forcey, lx2ly2lz2*sizeof(real));
    cudaMalloc(&forcez, lx2ly2lz2*sizeof(real));

    cudaMemset(rho, 0, lx2ly2lz2*sizeof(real));
    cudaMemset(ux, 0, lx2ly2lz2*sizeof(real));
    cudaMemset(uy, 0, lx2ly2lz2*sizeof(real));
    cudaMemset(uz, 0, lx2ly2lz2*sizeof(real));
    cudaMemset( f, 0, npop*lx2ly2lz2*sizeof(real));
    cudaMemset(fb, 0, npop*lx2ly2lz2*sizeof(real));

    cudaMemset(forcex, 0, lx2ly2lz2*sizeof(real));
    cudaMemset(forcey, 0, lx2ly2lz2*sizeof(real));
    cudaMemset(forcez, 0, lx2ly2lz2*sizeof(real));

    h_rho=(real *) malloc(lx2ly2lz2*sizeof(real));
    h_ux=(real *) malloc(lx2ly2lz2*sizeof(real));
    h_uy=(real *) malloc(lx2ly2lz2*sizeof(real));
    h_uz=(real *) malloc(lx2ly2lz2*sizeof(real));
    h_f =(real *) malloc(npop*lx2ly2lz2*sizeof(real));
    h_fb=(real *) malloc(npop*lx2ly2lz2*sizeof(real));

    // Setup hydro parameters
    cudaMallocManaged(&hydroparams, sizeof(HydroParameters));
    hydroparams->ambient_fx=Numerics::ambient_fx;
    hydroparams->ambient_fy=Numerics::ambient_fy;
    hydroparams->ambient_fz=Numerics::ambient_fz;
    
    using Numerics::xwidth;
    using Numerics::ywidth;
    using Numerics::zwidth;

    // Initial conditions for velocity and density
    if (Numerics::predef_case=="tg_eq") {
        for (int k=0; k<lz+2; k++) {
            real z=(k-0.5)*dz+zoffset;

            for (int j=0; j<ly+2; j++) {
                real y=(j-0.5)*dy+yoffset;

                for (int i=0; i<lx+2; i++) {
                    real x=(i-0.5)*dx+xoffset;

                    int ijk=k*(ly+2)*(lx+2)+j*(lx+2)+i;
                    real vel_scale=0.1;
                    h_rho[ijk]=1.0+1.0*vel_scale*vel_scale/16*(cos(2*pi2/xwidth*x)+cos(2*pi2/ywidth*y))*(cos(2*pi2/zwidth*z)+2)/RT;
                    h_ux[ijk]= vel_scale*sin(pi2/xwidth*x)*cos(pi2/ywidth*y)*cos(pi2/zwidth*z);
                    h_uy[ijk]=-vel_scale*cos(pi2/xwidth*x)*sin(pi2/ywidth*y)*cos(pi2/zwidth*z);
                    h_uz[ijk]=0;
                }
            }
        }
    } else if (Numerics::predef_case=="channel"){
        for (int k=0; k<lz+2; k++) {
            real z=(k-0.5)*dz+zoffset;

            for (int j=0; j<ly+2; j++) {
                real y=(j-0.5)*dy+yoffset;

                for (int i=0; i<lx+2; i++) {
                    real x=(i-0.5)*dx+xoffset;

                    int ijk=k*(ly+2)*(lx+2)+j*(lx+2)+i;
                    real vel_scale=0.1;
                    h_rho[ijk]=1.0;
                    h_ux[ijk]=0.0;
                    h_uy[ijk]=0.0;
                    h_uz[ijk]=vel_scale;
                }
            }
        }
    } else if (Numerics::predef_case=="general") {
        ;
    }
    else {
        std::cout << "ERROR: Unknown case.\n";
    }

    if (Numerics::restart) {
        read_dump(Numerics::nrestart, h_f);
    }
    else {
        for (int k=0; k<lz+2; k++) {
            for (int j=0; j<ly+2; j++) {
                for (int i=0; i<lx+2; i++) {
                    for (int ip=0; ip<npop; ip++) {
                        int ijk=k*(ly+2)*(lx+2)+j*(lx+2)+i;
                        H_F(ip,i,j,k)=feq(ip, h_tp, h_rho[ijk], h_ux[ijk], h_uy[ijk], h_uz[ijk]);
                    }
                }
            }
        }
    }

    cudaMemcpy(f, h_f, npop*(lx+2)*(ly+2)*(lz+2)*sizeof(real), cudaMemcpyHostToDevice);

    cudaMemcpy(rho, h_rho, lx2ly2lz2*sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(ux, h_ux, lx2ly2lz2*sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(uy, h_uy, lx2ly2lz2*sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(uz, h_uz, lx2ly2lz2*sizeof(real), cudaMemcpyHostToDevice);

    // cleanup host macroscopic (density-velocity) arrays
    std::memset(h_rho, 0, lx2ly2lz2*sizeof(real));
    std::memset(h_ux, 0, lx2ly2lz2*sizeof(real));
    std::memset(h_uy, 0, lx2ly2lz2*sizeof(real));
    std::memset(h_uz, 0, lx2ly2lz2*sizeof(real));

    // Perform test output from single thread
    dim3 threadsPerBlock(16, 16, 1);
    dim3 numBlocks(lx/threadsPerBlock.x, ly/threadsPerBlock.y, lz/threadsPerBlock.z);
    test_params<<<numBlocks, threadsPerBlock>>>(numparams, params);

    return ;
}

real feq(int ip, real tp[], real r, real u, real v, real w)
{
    real eu=(DVM::cix[ip]*u+DVM::ciy[ip]*v+DVM::ciz[ip]*w)/DVM::RT;
    real uv=(u*u+v*v+w*w)/DVM::RT;
    return tp[ip]*r*(1.0+eu+0.5*(eu*eu-uv));
}

void upload_macro()
{
    using Numerics::lx;
    using Numerics::ly;
    using Numerics::lz;

    const dim3 threadsPerBlock(32, 8, 2);
    const dim3 numBlocks(lx/threadsPerBlock.x, ly/threadsPerBlock.y, lz/threadsPerBlock.z);
    calc_macro<<<numBlocks, threadsPerBlock>>>(rho, ux, uy, uz, f, numparams, params, hydroparams);
}

void dugks_step(int step)
{
    using DVM::npop;

    using Numerics::lx;
    using Numerics::ly;
    using Numerics::lz;

    //-------------------------------------------
    // Stage 1.0: Compute f^bar
    //-------------------------------------------
    const dim3 threadsPerBlock(16, 8, 1);
    const dim3 numBlocks(lx/threadsPerBlock.x, ly/threadsPerBlock.y, lz/threadsPerBlock.z);

    if (Numerics::perturb) {
        calc_perturb_force<<<numBlocks, threadsPerBlock>>>(step, forcex, forcey, forcez, numparams, params, hydroparams);
    }
    else {
        calc_force<<<numBlocks, threadsPerBlock>>>(step, forcex, forcey, forcez, numparams, params, hydroparams);
    }
    
    // calc fb
    calc_fb1<<<numBlocks, threadsPerBlock>>>(fb, f, forcex, forcey, forcez, numparams, params);

    cudaDeviceSynchronize();

    //-------------------------------------------
    // Stage 1.1: Apply BC's to f^bar
    //-------------------------------------------
    dim3 threadsPerBlockX(1, 16, 16);
    dim3 blocksX(1, (ly+2+16)/16, (lz+2+16)/16);
    if (xnext==MPI_PROC_NULL)
        extend_fb_x_wall<<<blocksX, threadsPerBlockX>>>(fb, numparams, params);
    else
        extend_fb_x<<<blocksX, threadsPerBlockX>>>(fb, numparams, params);

    dim3 threadsPerBlockY(16, 1, 16);
    dim3 blocksY((lx+2+16)/16, 1, (lz+2+16)/16);
    if (ynext==MPI_PROC_NULL) {
        std::cout << "ERROR: wall BC is not supported in y-direction" << std::endl;
        std::exit(2);
    }
    else {
        extend_fb_y<<<blocksY, threadsPerBlockY>>>(fb, numparams, params);
    }
        
    // dim3 threadsPerBlockZ(16, 16, 1);
    // dim3 blocksZ((lx+2+16)/16, (ly+2+16)/16, 1);
    // extend_fb_z<<<blocksZ, threadsPerBlockZ>>>(fb, numparams, params);
    if (znext==MPI_PROC_NULL) {
        std::cout << "ERROR: wall BC is not supported in z-direction" << std::endl;
        std::exit(2);
    }
    else {
        exchange(fb);
    }
    
    // Shared memory based version
    const dim3 threadsPerBlockXX(8, 4, 2);
    const dim3 numBlocksXX(lx/threadsPerBlockXX.x, ly/threadsPerBlockXX.y, lz/threadsPerBlockXX.z);
    const int mysize=NPOP*(threadsPerBlockXX.x+2)*(threadsPerBlockXX.y+2)*(threadsPerBlockXX.z+2)*sizeof(real);
    if (mysize>SHMEM_LIMIT) {
        std::cout << "Shared memory size: " << mysize << "\n\n";
        std::cout << "ERROR: shmem. error, SIZE: " << mysize << std::endl;
        exit(1);
    }
    calc_fx_tri<<<numBlocksXX, threadsPerBlockXX, mysize>>>(f, fb, numparams, params, hydroparams);

    cudaDeviceSynchronize();
}

void copy_vel()
{
    using Numerics::lx;
    using Numerics::ly;
    using Numerics::lz;

    int lx2ly2lz2=(lx+2)*(ly+2)*(lz+2);

    cudaMemcpy(h_rho, rho, lx2ly2lz2*sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_ux, ux, lx2ly2lz2*sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_uy, uy, lx2ly2lz2*sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_uz, uz, lx2ly2lz2*sizeof(real), cudaMemcpyDeviceToHost);
}

void copy_ddf()
{
    using Numerics::lx;
    using Numerics::ly;
    using Numerics::lz;

    using DVM::npop;

    cudaMemcpy(h_f, f, npop*(lx+2)*(ly+2)*(lz+2)*sizeof(real), cudaMemcpyDeviceToHost);
}

//===========================================================================
// Computes fluxes and evolves distribution function.
//===========================================================================
__global__ void calc_fx(real f[], const real fb[], const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp)
{
    int i=threadIdx.x+blockIdx.x*blockDim.x+1;
    int j=threadIdx.y+blockIdx.y*blockDim.y+1;
    int k=threadIdx.z+blockIdx.z*blockDim.z+1;

    int i0=threadIdx.x+1;
    int j0=threadIdx.y+1;
    int k0=threadIdx.z+1;

    const real& hdt=p->hdt;
    const int npop=NPOP;
    //         0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
    const int cix[]= {0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
    const int ciy[]= {0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1};
    const int ciz[]= {0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1};

    const int& lx=n->lx;
    const int& ly=n->ly;
    const int& lz=n->lz;

    const real dx=n->dx;
    const real dy=n->dy;
    const real dz=n->dz;

    extern __shared__ real fbs[];

#include "block_copy.h"    
    __syncthreads();

    real f_temp[6][19];
#define FTEMP(ip, s) f_temp[s][ip]

    // compute only for internal cells
    real temp, aa, bb, gradx, grady, gradz;
    aa=bb=0;
    gradx=grady=gradz=temp=0;

    // ++ X-direction
    for (int ip=0; ip<npop; ip++) {
        temp=0.5*(FBs(ip,i0+1, j0, k0)+FBs(ip,i0, j0, k0));

        gradx=(FBs(ip,i0+1, j0, k0)-FBs(ip,i0, j0, k0))/(dx);

        aa=0.5*(FBs(ip,i0,j0+1,k0)+FBs(ip,i0+1,j0+1,k0));
        bb=0.5*(FBs(ip,i0,j0-1,k0)+FBs(ip,i0+1,j0-1,k0));
        grady=(aa-bb)/(2*dy);

        aa=0.5*(FBs(ip,i0,j0,k0+1)+FBs(ip,i0+1,j0,k0+1));
        bb=0.5*(FBs(ip,i0,j0,k0-1)+FBs(ip,i0+1,j0,k0-1));
        gradz=(aa-bb)/(2*dz);

        FTEMP(ip, 1)=temp-(hdt)*(cix[ip]*gradx+ciy[ip]*grady+ciz[ip]*gradz);
    }
    gradx=grady=gradz=temp=0;

    // -- X-direction
    for (int ip=0; ip<npop; ip++) {
        temp=0.5*(FBs(ip,i0, j0, k0)+FBs(ip,i0-1, j0, k0));

        gradx=(FBs(ip,i0, j0, k0)-FBs(ip,i0-1, j0, k0))/(dx);

        aa=0.5*(FBs(ip,i0,j0+1,k0)+FBs(ip,i0-1,j0+1,k0));
        bb=0.5*(FBs(ip,i0,j0-1,k0)+FBs(ip,i0-1,j0-1,k0));
        grady=(aa-bb)/(2*dy);

        aa=0.5*(FBs(ip,i0,j0,k0+1)+FBs(ip,i0-1,j0,k0+1));
        bb=0.5*(FBs(ip,i0,j0,k0-1)+FBs(ip,i0-1,j0,k0-1));
        gradz=(aa-bb)/(2*dz);

        FTEMP(ip, 0)=temp-(hdt)*(cix[ip]*gradx+ciy[ip]*grady+ciz[ip]*gradz);
    }
    gradx=grady=gradz=temp=0;

    // Wall boundary conditions
    if (n->periods[0]==0) {
        if (i==1) {
            FTEMP(1, 0)=FTEMP(2, 0);
            FTEMP(7, 0)=FTEMP(10,0);
            FTEMP(9, 0)=FTEMP(8, 0);
            FTEMP(11, 0)=FTEMP(14, 0);
            FTEMP(13, 0)=FTEMP(12, 0);
        }
        if (i==lx) {
            FTEMP(2, 1)=FTEMP(1, 1);
            FTEMP(10,1)=FTEMP(7,1);
            FTEMP(8, 1)=FTEMP(9, 1);
            FTEMP(14, 1)=FTEMP(11, 1);
            FTEMP(12, 1)=FTEMP(13, 1);
        }
    }
    
    // ++ Y-direction
    for (int ip=0; ip<npop; ip++) {
        temp=0.5*(FBs(ip,i0, j0+1, k0)+FBs(ip,i0, j0, k0));

        aa=0.5*(FBs(ip,i0+1,j0,k0)+FBs(ip,i0+1,j0+1,k0));
        bb=0.5*(FBs(ip,i0-1,j0,k0)+FBs(ip,i0-1,j0+1,k0));
        gradx=(aa-bb)/(2*dx);

        grady=(FBs(ip,i0,j0+1,k0)-FBs(ip, i0, j0, k0))/(dy);

        aa=0.5*(FBs(ip,i0,j0,k0+1)+FBs(ip,i0,j0+1,k0+1));
        bb=0.5*(FBs(ip,i0,j0,k0-1)+FBs(ip,i0,j0+1,k0-1));
        gradz=(aa-bb)/(2*dz);

        FTEMP(ip, 3)=temp-(hdt)*(cix[ip]*gradx+ciy[ip]*grady+ciz[ip]*gradz);
    }
    gradx=grady=gradz=temp=0;

    // -- Y-direction
    for (int ip=0; ip<npop; ip++) {
        temp=0.5*(FBs(ip,i0, j0, k0)+FBs(ip,i0, j0-1, k0));

        aa=0.5*(FBs(ip,i0+1,j0-1,k0)+FBs(ip,i0+1,j0,k0));
        bb=0.5*(FBs(ip,i0-1,j0-1,k0)+FBs(ip,i0-1,j0,k0));
        gradx=(aa-bb)/(2*dx);

        grady=(FBs(ip,i0,j0,k0)-FBs(ip, i0, j0-1, k0))/(dy);

        aa=0.5*(FBs(ip,i0,j0-1,k0+1)+FBs(ip,i0,j0,k0+1));
        bb=0.5*(FBs(ip,i0,j0-1,k0-1)+FBs(ip,i0,j0,k0-1));
        gradz=(aa-bb)/(2*dz);

        FTEMP(ip, 2)=temp-(hdt)*(cix[ip]*gradx+ciy[ip]*grady+ciz[ip]*gradz);
    }
    gradx=grady=gradz=temp=0;

    // ++ Z-direction
    for (int ip=0; ip<npop; ip++) {
        temp=0.5*(FBs(ip,i0, j0, k0+1)+FBs(ip,i0, j0, k0));

        aa=0.5*(FBs(ip,i0+1,j0,k0)+FBs(ip,i0+1,j0,k0+1));
        bb=0.5*(FBs(ip,i0-1,j0,k0)+FBs(ip,i0-1,j0,k0+1));
        gradx=(aa-bb)/(2*dx);

        aa=0.5*(FBs(ip,i0,j0+1,k0)+FBs(ip,i0,j0+1,k0+1));
        bb=0.5*(FBs(ip,i0,j0-1,k0)+FBs(ip,i0,j0-1,k0+1));
        grady=(aa-bb)/(2*dy);

        gradz=(FBs(ip,i0,j0,k0+1)-FBs(ip, i0, j0, k0))/(dz);

        FTEMP(ip, 5)=temp-(hdt)*(cix[ip]*gradx+ciy[ip]*grady+ciz[ip]*gradz);
    }
    gradx=grady=gradz=temp=0;

    // -- Z-direction
    for (int ip=0; ip<npop; ip++) {
        temp=0.5*(FBs(ip,i0, j0, k0)+FBs(ip,i0, j0, k0-1));

        aa=0.5*(FBs(ip,i0+1,j0,k0)+FBs(ip,i0+1,j0,k0-1));
        bb=0.5*(FBs(ip,i0-1,j0,k0)+FBs(ip,i0-1,j0,k0-1));
        gradx=(aa-bb)/(2*dx);

        aa=0.5*(FBs(ip,i0,j0+1,k0)+FBs(ip,i0,j0+1,k0-1));
        bb=0.5*(FBs(ip,i0,j0-1,k0)+FBs(ip,i0,j0-1,k0-1));
        grady=(aa-bb)/(2*dy);

        gradz=(FBs(ip,i0,j0,k0)-FBs(ip, i0, j0, k0-1))/(dz);

        FTEMP(ip, 4)=temp-(hdt)*(cix[ip]*gradx+ciy[ip]*grady+ciz[ip]*gradz);
    }
    gradx=grady=gradz=temp=0;

    // compute macroscopic quantities and f^eq
    real r, u, v, w;
    const real RT=1./3.;
    real eu, uv;
    real force1[NPOP];

    // X-direction
    //f2ruvw(&r, &u, &v, &w, &FTEMP(0, 0), p);
    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 0), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        FTEMP(ip, 0)=p->fw1*FTEMP(ip, 0)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    //f2ruvw(&r, &u, &v, &w, &FTEMP(0, 1), p);
    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 1), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        //FTEMP(ip, 1)=p->fw1*FTEMP(ip, 1)+p->fw2*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
        FTEMP(ip, 1)=p->fw1*FTEMP(ip, 1)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    // Y-direction
    //f2ruvw(&r, &u, &v, &w, &FTEMP(0, 2), p);
    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 2), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        //FTEMP(ip, 2)=p->fw1*FTEMP(ip, 2)+p->fw2*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
        FTEMP(ip, 2)=p->fw1*FTEMP(ip, 2)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    //f2ruvw(&r, &u, &v, &w, &FTEMP(0, 3), p);
    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 3), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        //FTEMP(ip, 3)=p->fw1*FTEMP(ip, 3)+p->fw2*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
        FTEMP(ip, 3)=p->fw1*FTEMP(ip, 3)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    // Z-direction
    //f2ruvw(&r, &u, &v, &w, &FTEMP(0, 4), p);
    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 4), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        //FTEMP(ip, 4)=p->fw1*FTEMP(ip, 4)+p->fw2*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
        FTEMP(ip, 4)=p->fw1*FTEMP(ip, 4)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    //f2ruvw(&r, &u, &v, &w, &FTEMP(0, 5), p);
    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 5), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        //FTEMP(ip, 5)=p->fw1*FTEMP(ip, 5)+p->fw2*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
        FTEMP(ip, 5)=p->fw1*FTEMP(ip, 5)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    // update final f
    const real cx=(n->dt)/(n->dx);
    const real cy=(n->dt)/(n->dy);
    const real cz=(n->dt)/(n->dz);
    for (int ip=0; ip<npop; ip++) {
        u=FTEMP(ip, 1)-FTEMP(ip, 0);
        v=FTEMP(ip, 3)-FTEMP(ip, 2);
        w=FTEMP(ip, 5)-FTEMP(ip, 4);
        //v=0;
        //w=0;

        r=FBs(ip, i0, j0, k0)*4./3.f-F(ip, i, j, k)*1./3.f;
        F(ip, i, j, k)=r-cx*cix[ip]*u-cy*ciy[ip]*v-cz*ciz[ip]*w;
    }
}

__device__ void f2ruvw(real *r, real *u, real *v, real *w, const real f[], const DugksParameters *p)
{
    const int npop=19;

    // int *cix=&(p->cixyzo[1*npop]);
    // int *ciy=&(p->cixyzo[2*npop]);
    // int *ciz=&(p->cixyzo[3*npop]);
    //         0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
    const int cix[]= {0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
    const int ciy[]= {0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1};
    const int ciz[]= {0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1};

    *r=*u=*v=*w=0;
    for (int ip=0; ip<npop; ip++) {
        *r+=f[ip];
        *u+=cix[ip]*f[ip];
        *v+=ciy[ip]*f[ip];
        *w+=ciz[ip]*f[ip];
    }
    *u/=*r;
    *v/=*r;
    *w/=*r;
}

__device__ void f2ruvwh(real *r, real *u, real *v, real *w, const real f[], real force1[], const DugksParameters *p, const HydroParameters *hp)
{
    const int npop=19;

    // int *cix=&(p->cixyzo[1*npop]);
    // int *ciy=&(p->cixyzo[2*npop]);
    // int *ciz=&(p->cixyzo[3*npop]);
    //         0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
    const int cix[]= {0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
    const int ciy[]= {0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1};
    const int ciz[]= {0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1};

    *r=*u=*v=*w=0;
    for (int ip=0; ip<npop; ip++) {
        *r+=f[ip];
        *u+=cix[ip]*f[ip];
        *v+=ciy[ip]*f[ip];
        *w+=ciz[ip]*f[ip];
    }
    *u/=*r;
    *v/=*r;
    *w/=*r;

    *u+=0.5*p->hdt*hp->ambient_fx;
    *v+=0.5*p->hdt*hp->ambient_fy;
    *w+=0.5*p->hdt*hp->ambient_fz;
    // *u+=0.5*p->hdt*FXX;
    // *v+=0.5*p->hdt*FYY;
    // *w+=0.5*p->hdt*FZZ;

    for (int ip=0; ip<npop; ip++) {
        //force1[ip]=(FXX*(cix[ip]-*u) + FYY*(ciy[ip]-*v) + FZZ*(ciz[ip]-*w))/p->RT;
        force1[ip]=(hp->ambient_fx*(cix[ip]-*u) + hp->ambient_fy*(ciy[ip]-*v) + hp->ambient_fz*(ciz[ip]-*w))/p->RT;
    }
}

//===========================================================================
// Computes fluxes and evolves distribution function. Uses trilinear
// interolation.
//===========================================================================
__global__ void calc_fx_tri(real f[], const real fb[], const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp)
{
    int i=threadIdx.x+blockIdx.x*blockDim.x+1;
    int j=threadIdx.y+blockIdx.y*blockDim.y+1;
    int k=threadIdx.z+blockIdx.z*blockDim.z+1;

    int i0=threadIdx.x+1;
    int j0=threadIdx.y+1;
    int k0=threadIdx.z+1;

    const real& hdt=p->hdt;
    const int npop=NPOP;
    //         0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
    const int cix[]= {0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
    const int ciy[]= {0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1};
    const int ciz[]= {0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1};
    // const int *cix=&(p->cixyzo[1*n->npop]);
    // const int *ciy=&(p->cixyzo[2*n->npop]);
    // const int *ciz=&(p->cixyzo[3*n->npop]);

    const int& lx=n->lx;
    const int& ly=n->ly;
    const int& lz=n->lz;

    const real& dx=n->dx;
    const real& dy=n->dy;
    const real& dz=n->dz;

    extern __shared__ real fbs[];

    #include "block_copy.h"
    __syncthreads();

    // fluxes for all six faces
    real f_temp[6][19];
#define FTEMP(ip, s) f_temp[s][ip]

    int ii, jj, kk;
    real mx0, my0, mz0;
    real mx1, my1, mz1;

    FTEMP(0, 1)=0.5*(FBs(0, i0, j0, k0)+FBs(0,i0+1, j0, k0));
    FTEMP(0, 0)=0.5*(FBs(0, i0, j0, k0)+FBs(0,i0-1, j0, k0));

    FTEMP(0, 3)=0.5*(FBs(0, i0, j0, k0)+FBs(0, i0, j0+1, k0));
    FTEMP(0, 2)=0.5*(FBs(0, i0, j0, k0)+FBs(0, i0, j0-1, k0));

    FTEMP(0, 5)=0.5*(FBs(0, i0, j0, k0)+FBs(0, i0, j0, k0+1));
    FTEMP(0, 4)=0.5*(FBs(0, i0, j0, k0)+FBs(0, i0, j0, k0-1));

    // X-direction
    for (int ip=1; ip<npop; ip++) {
        //! mx, my and mz are fractions used in interpolation
        mx1=(dx/2.-hdt*cix[ip])/dx;
        my1=(-hdt*ciy[ip])/dy;
        mz1=(-hdt*ciz[ip])/dz;

        mx0=1-mx1;
        my0=my1;
        mz0=mz1;

        //! positive x
        ii=1;
        jj=-ciy[ip];
        kk=-ciz[ip];

        FTEMP(ip, 1)=trilinear(
                         FBs(ip,    i0,    j0, k0),
                         FBs(ip, i0+ii,    j0, k0),
                         FBs(ip,    i0, j0+jj, k0),
                         FBs(ip, i0+ii, j0+jj, k0),
                         FBs(ip,    i0,    j0, k0+kk),
                         FBs(ip, i0+ii,    j0, k0+kk),
                         FBs(ip,    i0, j0+jj, k0+kk),
                         FBs(ip, i0+ii, j0+jj, k0+kk),
                         (mx1), fabs(my1), fabs(mz1));

        //! negative x
        ii=-1;

        FTEMP(ip, 0)=trilinear(
                         FBs(ip,    i0,    j0, k0),
                         FBs(ip, i0+ii,    j0, k0),
                         FBs(ip,    i0, j0+jj, k0),
                         FBs(ip, i0+ii, j0+jj, k0),
                         FBs(ip,    i0,    j0, k0+kk),
                         FBs(ip, i0+ii,    j0, k0+kk),
                         FBs(ip,    i0, j0+jj, k0+kk),
                         FBs(ip, i0+ii, j0+jj, k0+kk),
                         (mx0), fabs(my0), fabs(mz0));
    }

    // Wall boundary conditions
    if (n->periods[0]==0) {
        if (i==1) {
            FTEMP(1, 0)=FTEMP(2, 0);
            FTEMP(7, 0)=FTEMP(10,0);
            FTEMP(9, 0)=FTEMP(8, 0);
            FTEMP(11, 0)=FTEMP(14, 0);
            FTEMP(13, 0)=FTEMP(12, 0);
        }
        if (i==lx) {
            FTEMP(2, 1)=FTEMP(1, 1);
            FTEMP(10,1)=FTEMP(7,1);
            FTEMP(8, 1)=FTEMP(9, 1);
            FTEMP(14, 1)=FTEMP(11, 1);
            FTEMP(12, 1)=FTEMP(13, 1);
        }
    }

    // Y-direction
    for (int ip=1; ip<npop; ip++) {
        mx1=(-hdt*cix[ip])/dx;
        my1=(dy/2.-hdt*ciy[ip])/dy;
        mz1=(-hdt*ciz[ip])/dz;

        mx0=mx1;
        my0=1-my1;
        mz0=mz1;

        //! positive y
        ii=-cix[ip];
        jj=1;
        kk=-ciz[ip];

        //! fixme: do we change order of abcdefgh for other panels ?
        FTEMP(ip, 3)=trilinear(
                         FBs(ip,    i0,    j0, k0),
                         FBs(ip, i0+ii,    j0, k0),
                         FBs(ip,    i0, j0+jj, k0),
                         FBs(ip, i0+ii, j0+jj, k0),
                         FBs(ip,    i0,    j0, k0+kk),
                         FBs(ip, i0+ii,    j0, k0+kk),
                         FBs(ip,    i0, j0+jj, k0+kk),
                         FBs(ip, i0+ii, j0+jj, k0+kk),
                         fabs(mx1), (my1), fabs(mz1));

        //! negative y
        jj=-1;

        FTEMP(ip, 2)=trilinear(
                         FBs(ip,    i0,    j0, k0),
                         FBs(ip, i0+ii,    j0, k0),
                         FBs(ip,    i0, j0+jj, k0),
                         FBs(ip, i0+ii, j0+jj, k0),
                         FBs(ip,    i0,    j0, k0+kk),
                         FBs(ip, i0+ii,    j0, k0+kk),
                         FBs(ip,    i0, j0+jj, k0+kk),
                         FBs(ip, i0+ii, j0+jj, k0+kk),
                         fabs(mx0), (my0), fabs(mz0));
    }

    // Z-direction
    for (int ip=1; ip<npop; ip++) {
        mx1=(-hdt*cix[ip])/dx;
        my1=(-hdt*ciy[ip])/dy;
        mz1=(dz/2.-hdt*ciz[ip])/dz;

        mx0=mx1;
        my0=my1;
        mz0=1-mz1;

        //! positive z
        ii=-cix[ip];
        jj=-ciy[ip];
        kk=1;

        FTEMP(ip, 5)=trilinear(
                         FBs(ip,    i0,    j0, k0),
                         FBs(ip, i0+ii,    j0, k0),
                         FBs(ip,    i0, j0+jj, k0),
                         FBs(ip, i0+ii, j0+jj, k0),
                         FBs(ip,    i0,    j0, k0+kk),
                         FBs(ip, i0+ii,    j0, k0+kk),
                         FBs(ip,    i0, j0+jj, k0+kk),
                         FBs(ip, i0+ii, j0+jj, k0+kk),
                         fabs(mx1), fabs(my1), (mz1));

        //! negative z
        kk=-1;

        FTEMP(ip, 4)=trilinear(
                         FBs(ip,    i0,    j0, k0),
                         FBs(ip, i0+ii,    j0, k0),
                         FBs(ip,    i0, j0+jj, k0),
                         FBs(ip, i0+ii, j0+jj, k0),
                         FBs(ip,    i0,    j0, k0+kk),
                         FBs(ip, i0+ii,    j0, k0+kk),
                         FBs(ip,    i0, j0+jj, k0+kk),
                         FBs(ip, i0+ii, j0+jj, k0+kk),
                         fabs(mx0), fabs(my0), (mz0));
    }

    // compute macroscopic quantities and f^eq
    real r, u, v, w;
    real force1[NPOP];

    // X-direction
    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 0), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        FTEMP(ip, 0)=p->fw1*FTEMP(ip, 0)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 1), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        FTEMP(ip, 1)=p->fw1*FTEMP(ip, 1)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    // Y-direction
    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 2), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        FTEMP(ip, 2)=p->fw1*FTEMP(ip, 2)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 3), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        FTEMP(ip, 3)=p->fw1*FTEMP(ip, 3)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    // Z-direction
    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 4), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        FTEMP(ip, 4)=p->fw1*FTEMP(ip, 4)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    f2ruvwh(&r, &u, &v, &w, &FTEMP(0, 5), force1, p, hp);
    for (int ip=0; ip<npop; ip++) {
        FTEMP(ip, 5)=p->fw1*FTEMP(ip, 5)+(p->fw2+p->fw0*force1[ip])*feq_s(r, u, v, w, p->tp[ip], cix[ip], ciy[ip], ciz[ip]);
    }

    // update final f
    const real cx=(n->dt)/(n->dx);
    const real cy=(n->dt)/(n->dy);
    const real cz=(n->dt)/(n->dz);
    for (int ip=0; ip<npop; ip++) {
        u=FTEMP(ip, 1)-FTEMP(ip, 0);
        v=FTEMP(ip, 3)-FTEMP(ip, 2);
        w=FTEMP(ip, 5)-FTEMP(ip, 4);

        r=FBs(ip, i0, j0, k0)*4./3.-F(ip, i, j, k)*1./3.;
        F(ip, i, j, k)=r-cx*cix[ip]*u-cy*ciy[ip]*v-cz*ciz[ip]*w;
    }

    return ;
}

//===========================================================================
// Perform trilinear interpolation
//===========================================================================
__device__ real trilinear(real a, real b, real c, real d, real e, real f, real g, real h, real mx, real my, real mz)
{
    real x00, x01, x10, x11;
    real y00, y11;

    x00=(1-mx)*a+mx*b;
    x01=(1-mx)*c+mx*d;
    x10=(1-mx)*e+mx*f;
    x11=(1-mx)*g+mx*h;

    y00=(1-my)*x00+my*x01;
    y11=(1-my)*x10+my*x11;

    return (1-mz)*y00+mz*y11;
}

//===========================================================================
// Print some parameters to check their proper loading by GPU
//===========================================================================
__global__ void test_params(NumericalParameters *n, DugksParameters *p)
{
    int ii=threadIdx.x+blockIdx.x*blockDim.x;
    int jj=threadIdx.y+blockIdx.y*blockDim.y;
    int kk=threadIdx.z+blockIdx.z*blockDim.z;

    int& npop=n->npop;

    int *cix=&(p->cixyzo[1*npop]);
    int *ciy=&(p->cixyzo[2*npop]);
    int *ciz=&(p->cixyzo[3*npop]);

    if (ii==1 && jj==1 && kk==1) {
        printf("Parameters in thread 1,1,1:\n");
        printf("%f\n", p->RT);
        printf("%f\n", p->visc);
        printf("%f\n", p->hdt);
        printf("%f\n", p->bw0);
        printf("%f\n", p->bw1);
        printf("%f\n", p->bw2);
        printf("%f\n", p->fw0);
        printf("%f\n", p->fw1);
        printf("%f\n", p->fw2);
        printf("%f\n", p->qw1);
        printf("%f\n", p->qw2);
        printf("%f\n", p->gw1);
        printf("%f\n", p->gw2);

        printf("cixyz %d:\n", npop);
        for (int i=0; i<npop; i++) {
            printf("%d\t%d %d %d\n", i, cix[i], ciy[i], ciz[i]);
        }
    }
}

// ==========================================================================
// Here we compute macroscopic quantities from distribution function.
// =========================================================================
__global__ void calc_macro(real r[], real u[], real v[], real w[], const real f[],
                           const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp)
{
    int ii=threadIdx.x+blockIdx.x*blockDim.x;
    int jj=threadIdx.y+blockIdx.y*blockDim.y;
    int kk=threadIdx.z+blockIdx.z*blockDim.z;

    int i=ii+1;
    int j=jj+1;
    int k=kk+1;

    const int *cix=&(p->cixyzo[1*n->npop]);
    const int *ciy=&(p->cixyzo[2*n->npop]);
    const int *ciz=&(p->cixyzo[3*n->npop]);

    const int& npop=n->npop;
    const int& lx=n->lx;
    const int& ly=n->ly;
    const int& lz=n->lz;

    real rr, uu, vv, ww;

    rr=uu=vv=ww=0;
    for (int ip=0; ip<npop; ip++) {
        real ff=F(ip, i, j, k);
        rr+=ff;
        uu+=cix[ip]*ff;
        vv+=ciy[ip]*ff;
        ww+=ciz[ip]*ff;
    }
    uu/=rr;
    vv/=rr;
    ww/=rr;

    //    int ijk=i*(n->ly+2)*(n->lz+2)+j*(n->lz+2)+k;
    int ijk=k*(n->ly+2)*(n->lx+2)+j*(n->lx+2)+i;
    r[ijk]=rr;
    u[ijk]=uu+hp->ambient_fx*p->hdt;
    v[ijk]=vv+hp->ambient_fy*p->hdt;
    w[ijk]=ww+hp->ambient_fz*p->hdt;
    // u[ijk]=uu+FXX*p->hdt;
    // v[ijk]=vv+FYY*p->hdt;
    // w[ijk]=ww+FZZ*p->hdt;
}

__global__ void calc_force(int step, real forcex[], real forcey[], real forcez[],
                           const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp)
{
    int i=threadIdx.x+blockIdx.x*blockDim.x+1;
    int j=threadIdx.y+blockIdx.y*blockDim.y+1;
    int k=threadIdx.z+blockIdx.z*blockDim.z+1;

    const int& lx=n->lx;
    const int& ly=n->ly;
    const int& lz=n->lz;

    FORCEX(i,j,k)=hp->ambient_fx;
    FORCEY(i,j,k)=hp->ambient_fy;
    FORCEZ(i,j,k)=hp->ambient_fz;
    // FORCEX(i,j,k)=FXX;
    // FORCEY(i,j,k)=FYY;
    // FORCEZ(i,j,k)=FZZ;
    return ;
}

__global__ void calc_perturb_force(int step, real forcex[], real forcey[], real forcez[],
                                   const NumericalParameters *n, const DugksParameters *p, const HydroParameters *hp)
{
    int i=threadIdx.x+blockIdx.x*blockDim.x+1;
    int j=threadIdx.y+blockIdx.y*blockDim.y+1;
    int k=threadIdx.z+blockIdx.z*blockDim.z+1;

    const int& lx=n->lx;
    const int& ly=n->ly;
    const int& lz=n->lz;
    const real& xwidth=n->xwidth;
    const real& ywidth=n->ywidth;
    const real& zwidth=n->zwidth;
    const real& visc=p->visc;

    const real beta=3;
    const real Tpd=2000;
    const real gamma=2;
    const real Re_tau=180;
    const real ustar=2*Re_tau*visc/xwidth;

    const real FF=2*pow(ustar,2)/xwidth;
    real A0=5*40.0/n->nx*beta*sin(pi2*step*n->dt/Tpd);

    real xx, yy, zz;

    if (3<i && i<n->nx/4+2) {
        xx=pi2*(i-3)*n->dx/(xwidth/4.);
        yy=pi2*(j-0.5)*n->dy/ywidth;
        zz=pi2*((k-0.5)*n->dz+n->zoffset)/zwidth;

        FORCEX(i,j,k)=0.5*FF*A0*xwidth/4*(1-cos(xx))*cos(beta*yy)*cos(gamma*zz);
        FORCEY(i,j,k)=0.5*FF*A0*ywidth/gamma*sin(xx)*cos(beta*yy)*sin(gamma*zz);
        FORCEZ(i,j,k)=FF*(1-A0*zwidth/beta*sin(xx)*sin(beta*yy)*cos(gamma*zz));
    }
    else {
        FORCEX(i,j,k)=hp->ambient_fx;
        FORCEY(i,j,k)=hp->ambient_fy;
        FORCEZ(i,j,k)=hp->ambient_fz;
    }
}

// ========================================================================
// Calculate f^bar
// ========================================================================
__global__ void calc_fb1(real fb[], real f[], real forcex[], real forcey[], real forcez[], const NumericalParameters *n, const DugksParameters *p)
{
    int i=threadIdx.x+blockIdx.x*blockDim.x+1;
    int j=threadIdx.y+blockIdx.y*blockDim.y+1;
    int k=threadIdx.z+blockIdx.z*blockDim.z+1;

    const int& lx=n->lx;
    const int& ly=n->ly;
    const int& lz=n->lz;
    const real& bw1=p->bw1;
    const real& bw2=p->bw2;
    const real& bw0=p->bw0;
    // const int& npop=n->npop;
    // const int *cix=&(p->cixyzo[1*n->npop]);
    // const int *ciy=&(p->cixyzo[2*n->npop]);
    // const int *ciz=&(p->cixyzo[3*n->npop]);

    // const int lx=256;
    // const int ly=256;
    // const int lz=256;
    const int npop=NPOP;
    //         0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
    const int cix[]= {0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
    const int ciy[]= {0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1};
    const int ciz[]= {0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1};
    const real ww0 = 1.0/3.0;
    const real ww1 = 1.0/18.0;
    const real ww2 = 1.0/36.0;
    const real tp[]= {ww0,
                      ww1, ww1, ww1, ww1, ww1, ww1,
                      ww2, ww2, ww2, ww2, ww2, ww2, ww2, ww2, ww2, ww2, ww2, ww2
                     };

    real ff[NPOP];

    real rr, uu, vv, ww;
    rr=uu=vv=ww=0;
    for (int ip=0; ip<npop; ip++) {
        ff[ip]=F(ip, i, j, k);
        rr+=ff[ip];
        uu+=cix[ip]*ff[ip];
        vv+=ciy[ip]*ff[ip];
        ww+=ciz[ip]*ff[ip];
    }
    uu/=rr;
    vv/=rr;
    ww/=rr;

    uu+=p->hdt*FORCEX(i,j,k);
    vv+=p->hdt*FORCEY(i,j,k);
    ww+=p->hdt*FORCEZ(i,j,k);

    const real RT=1./3.f;
    real feqa[NPOP];
    real eu, uv;
    for (int ip=0; ip<npop; ip++) {
        eu=(cix[ip]*uu+ciy[ip]*vv+ciz[ip]*ww)/RT;
        uv=(uu*uu+vv*vv+ww*ww)/(RT);
        feqa[ip]=tp[ip]*rr*(1.0f+eu+0.5f*(eu*eu-uv));
    }

    for (int ip=0; ip<npop; ip++) {
        real temp = ((cix[ip]-uu)*FORCEX(i,j,k)+(ciy[ip]-vv)*FORCEY(i,j,k)+(ciz[ip]-ww)*FORCEZ(i,j,k))/RT;
        FB(ip, i, j, k)=bw1*ff[ip]+(bw2+bw0*temp)*feqa[ip];
    }
}

// ==========================================================================
// Equilibrium f for single direction
// ==========================================================================
__device__ real feq_s(real r, real u, real v, real w, real tp, int cix, int ciy, int ciz)
{
    const real RT=1./3.f;
    real eu=(cix*u+ciy*v+ciz*w)/RT;
    real uv=(u*u+v*v+w*w)/(RT);
    return tp*r*(1.0f+eu+0.5f*(eu*eu-uv));
}

// ==========================================================================
// Apply BC copy in X-direction
// ==========================================================================
__global__ void extend_fb_x(real *fb, NumericalParameters *n, DugksParameters *p)
{
    int i=threadIdx.x+blockIdx.x*blockDim.x;
    int j=threadIdx.y+blockIdx.y*blockDim.y;
    int k=threadIdx.z+blockIdx.z*blockDim.z;

    int& npop=n->npop;
    int& lx=n->lx;
    int& ly=n->ly;
    int& lz=n->lz;

    if (j<ly+2 && k<lz+2) {
        for (int ip=0; ip<npop; ip++) {
            FB(ip,    0, j, k)=FB(ip, lx, j, k);
            FB(ip, lx+1, j, k)=FB(ip, 1, j, k);
        }
    }
}

// ==========================================================================
// Apply wall-BC
// ==========================================================================
__global__ void extend_fb_x_wall(real *fb, NumericalParameters *n, DugksParameters *p)
{
    int i=threadIdx.x+blockIdx.x*blockDim.x;
    int j=threadIdx.y+blockIdx.y*blockDim.y;
    int k=threadIdx.z+blockIdx.z*blockDim.z;

    int& npop=n->npop;
    int& lx=n->lx;
    int& ly=n->ly;
    int& lz=n->lz;

    if (j<ly+2 && k<lz+2) {
        for (int ip=0; ip<npop; ip++) {
            FB(ip,    0, j, k)=2*FB(ip, 1, j, k)-FB(ip, 2, j, k);
            FB(ip, lx+1, j, k)=2*FB(ip, lx, j, k)-FB(ip, lx-1, j, k);
        }
    }
}

// ==========================================================================
// Apply BC copy in Y-direction
// ==========================================================================
__global__ void extend_fb_y(real *fb, NumericalParameters *n, DugksParameters *p)
{
    int i=threadIdx.x+blockIdx.x*blockDim.x;
    int j=threadIdx.y+blockIdx.y*blockDim.y;
    int k=threadIdx.z+blockIdx.z*blockDim.z;

    int& npop=n->npop;
    int& lx=n->lx;
    int& ly=n->ly;
    int& lz=n->lz;

    if (i<lx+2 && k<lz+2) {
        for (int ip=0; ip<npop; ip++) {
            FB(ip, i,    0, k)=FB(ip, i, ly, k);
            FB(ip, i, ly+1, k)=FB(ip, i, 1, k);
        }
    }
}

// ==========================================================================
// Apply BC copy in Z-direction
// ==========================================================================
__global__ void extend_fb_z(real *fb, NumericalParameters *n, DugksParameters *p)
{
    int i=threadIdx.x+blockIdx.x*blockDim.x;
    int j=threadIdx.y+blockIdx.y*blockDim.y;
    int k=threadIdx.z+blockIdx.z*blockDim.z;

    const int& npop=n->npop;
    const int& lx=n->lx;
    const int& ly=n->ly;
    const int& lz=n->lz;

    if (i<lx+2 && j<ly+2) {
        for (int ip=0; ip<npop; ip++) {
            FB(ip, i, j,    0)=FB(ip, i, j, lz);
            FB(ip, i, j, lz+1)=FB(ip, i, j, 1);
        }
    }
}
