#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <cstring>

#include "precision.h"
#include "numerics.h"
#include "lbm_d3q19.h"
#include "config.h"

// source:https://docs.open-mpi.org/en/v5.0.x/tuning-apps/networking/cuda.html
// see the Section 11.2.6.17
#if !defined(OPEN_MPI) || !OPEN_MPI
#error This source code uses an Open MPI-specific extension
#endif

/* Needed for MPIX_Query_cuda_support(), below */
#include "mpi-ext.h"

using namespace Numerics;

int myrank, nproc;
MPI_Status status;
int ierr;

// Three-dimensional cartesian topology is assumed,
// but for now only 1D domain decomposition in
// Z-direction is used
extern const int ndims=3;
int dims[ndims];
int periods[ndims];
const int reorder=0;
MPI_Comm cartesian;
int coords[ndims];
int neighbors[ndims*ndims*ndims];
int xprev, xnext;
int yprev, ynext;
int zprev, znext;

// Physical subdomain boundaries
real xoffset, xoffset1;
real yoffset, yoffset1;
real zoffset, zoffset1;

double timestart, timeend;
double timestart1, timeend1;

// To send/recv halo layers
real* zbuffer0;
real* zbuffer1;
real* zbuffer00;
real* zbuffer11;

//===========================================================================
// Check am I root thread
//===========================================================================
bool root()
{
    return (myrank==0) ? true : false;
}

void free_parallel(void )
{
    free(zbuffer0);
    free(zbuffer1);
    free(zbuffer00);
    free(zbuffer11);
    
    MPI_Finalize();
}

//===========================================================================
// Initialize 1D domain decomposition in Z-direction
//===========================================================================
void init_parallel(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
    if (root()) printf("This MPI library has CUDA-aware support.\n", MPIX_CUDA_AWARE_SUPPORT);
#elif defined(MPIX_CUDA_AWARE_SUPPORT) && !MPIX_CUDA_AWARE_SUPPORT
    if (root()) printf("This MPI library does not have CUDA-aware support.\n");
#else
    if (root()) printf("This MPI library cannot determine if there is CUDA-aware support.\n");
#endif /* MPIX_CUDA_AWARE_SUPPORT */

    if (myrank==0) {
        std::cout << "Starting the MPI..." << std::endl;
        std::cout << "nproc: " << nproc << ", myrank: " << myrank << std::endl;
    }

    Config config;
    periods[0]=std::stoi(config.get("periods_x"));
    periods[1]=std::stoi(config.get("periods_y"));
    periods[2]=std::stoi(config.get("periods_z"));
    
    dims[0]=1;
    dims[1]=1;
    dims[2]=nproc;

    lx=nx/dims[0];
    ly=ny/dims[1];
    lz=nz/dims[2];

    zbuffer0 = (real *)malloc(19*(lx+2)*(ly+2)*sizeof(real));
    zbuffer1 = (real *)malloc(19*(lx+2)*(ly+2)*sizeof(real));
    zbuffer00 = (real *)malloc(19*(lx+2)*(ly+2)*sizeof(real));
    zbuffer11 = (real *)malloc(19*(lx+2)*(ly+2)*sizeof(real));

    ierr=MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartesian);
    if (ierr!=0) {
        std::cout << "ERROR: Cartesian topology error. " << ierr << std::endl;
    }

    MPI_Cart_coords(cartesian, myrank, ndims, coords);

    MPI_Cart_shift(cartesian, 0, 1, &xprev, &xnext);
    MPI_Cart_shift(cartesian, 1, 1, &yprev, &ynext);
    MPI_Cart_shift(cartesian, 2, 1, &zprev, &znext);

    printf("myrank: %4d, coords: [%4d, %4d, %4d], neighbors: [%4d %4d, %4d %4d, %4d %4d]\n", myrank,
           coords[0], coords[1], coords[2],
           xprev, xnext, yprev, ynext, zprev, znext);

    MPI_Barrier(cartesian);

    xoffset=coords[0]*lx*dx;
    xoffset1=(coords[0]+1)*lx*dx;

    yoffset=coords[1]*ly*dy;
    yoffset1=(coords[1]+1)*ly*dy;

    zoffset=coords[2]*lz*dz;
    zoffset1=(coords[2]+1)*lz*dz;

    // get all neighbors
    int temp_coords[ndims];
    int temp_rank;

    for (int i=-1; i<=1; i++) {
        for (int j=-1; j<=1; j++) {
            for (int k=-1; k<=1; k++) {
                temp_coords[0]=coords[0]+i;
                temp_coords[1]=coords[1]+j;
                temp_coords[2]=coords[2]+k;

                int ijk=(k+1)*ndims*ndims+(j+1)*ndims+(i+1);
                neighbors[ijk]=MPI_PROC_NULL;

                if (periods[0]==0 && (temp_coords[0]<0) || temp_coords[0]>dims[0]-1)
                    continue;

                if (periods[1]==0 && (temp_coords[1]<0) || temp_coords[1]>dims[1]-1)
                    continue;

                if (periods[2]==0 && (temp_coords[2]<0) || temp_coords[2]>dims[2]-1)
                    continue;

                MPI_Cart_rank(cartesian, temp_coords, &temp_rank);
                neighbors[ijk]=temp_rank;
            }
        }
    }

    if (myrank==0)
        for (int i=-1; i<=1; i++) {
            for (int j=-1; j<=1; j++) {
                for (int k=-1; k<=1; k++) {
                    int ijk=(k+1)*ndims*ndims+(j+1)*ndims+(i+1);
                    printf("myrank: %4d [%2d, %2d, %2d] = %4d\n", myrank, i, j, k, neighbors[ijk]);
                }
            }
        }

    MPI_Barrier(cartesian);

    return;
}

// ==========================================================================
// Exchange halo layers with neighbors
// ==========================================================================
void exchange(real* fb)
{
    int length=19*(lx+2)*(ly+2);

#define FB(ip,i,j,k)     fb[(k)*(ly+2)*(lx+2)*NPOP + (j)*(lx+2)*NPOP + (i)*NPOP + ip]

    // Use CUDA aware MPI
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
    MPI_Sendrecv(&FB(0, 0, 0,    1), length, MY_MPI_REAL, zprev, 0,
                 &FB(0, 0, 0, lz+1), length, MY_MPI_REAL, znext, 0,
                 cartesian, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&FB(0, 0, 0, lz), length, MY_MPI_REAL, znext, 0,
                 &FB(0, 0, 0,  0), length, MY_MPI_REAL, zprev, 0,
                 cartesian, MPI_STATUS_IGNORE);
    cudaDeviceSynchronize();
#else

    // Use non CUDA aware MPI
    int lengthb=19*(lx+2)*(ly+2)*sizeof(real);

    cudaMemcpy(zbuffer0, &FB(0, 0, 0, 1), lengthb, cudaMemcpyDeviceToHost);
    cudaMemcpy(zbuffer1, &FB(0, 0, 0, lz), lengthb, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

#ifdef PREC4
    MPI_Sendrecv(zbuffer0, length, MPI_FLOAT, zprev, 0,
                 zbuffer11, length, MPI_FLOAT, znext, 0,
                 cartesian, MPI_STATUS_IGNORE);

    MPI_Sendrecv(zbuffer1, length, MPI_FLOAT, znext, 0,
                 zbuffer00, length, MPI_FLOAT, zprev, 0,
                 cartesian, MPI_STATUS_IGNORE);
#else
    MPI_Sendrecv(zbuffer0, length, MPI_DOUBLE, zprev, 0,
                 zbuffer11, length, MPI_DOUBLE, znext, 0,
                 cartesian, MPI_STATUS_IGNORE);

    MPI_Sendrecv(zbuffer1, length, MPI_DOUBLE, znext, 0,
                 zbuffer00, length, MPI_DOUBLE, zprev, 0,
                 cartesian, MPI_STATUS_IGNORE);

    cudaMemcpy(&FB(0, 0, 0, 0), zbuffer00, lengthb, cudaMemcpyHostToDevice);
    cudaMemcpy(&FB(0, 0, 0, lz+1), zbuffer11, lengthb, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
#endif
#endif

    MPI_Barrier(cartesian);
}

// ==========================================================================
// Collect local arrays into global arrays, used for I/O
// ==========================================================================
void collect_ruvw(real *rhg, real *uxg, real *uyg, real *uzg,
                  real *rho, real *ux, real *uy, real *uz)
{
    real *rho1;
    real *ux1;
    real *uy1;
    real *uz1;

    const int length2=(lx+2)*(ly+2)*(lz+2);
    const int length=(lx+0)*(ly+0)*(lz+0);

    rho1=new real[length];
    ux1 =new real[length];
    uy1 =new real[length];
    uz1 =new real[length];

    // Here we first copy from halo-extended fields to non-halo layer
    // fields, in other words from
    // ux(0:lx+1, 0:ly+1, 0:lz+1) -> ux(1:lx, 1:ly, 1:lz),
    // and of course the last ux has zero based indexing, and
    // spans from 0:lx-1, ...
    for (int k=0; k<lz; k++) {
        for (int j=0; j<ly; j++) {
            for (int i=0; i<lx; i++) {
                int ijk0=k*(ly)*(lx)+j*(lx)+i;
                int ijk2=(k+1)*(ly+2)*(lx+2)+(j+1)*(lx+2)+(i+1);

                rho1[ijk0]=rho[ijk2];
                ux1[ijk0]=ux[ijk2];
                uy1[ijk0]=uy[ijk2];
                uz1[ijk0]=uz[ijk2];
            }
        }
    }

    if (root()) {
        std::memcpy(rhg, rho1, length*sizeof(real));
        std::memcpy(uxg, ux1, length*sizeof(real));
        std::memcpy(uyg, uy1, length*sizeof(real));
        std::memcpy(uzg, uz1, length*sizeof(real));

        // Now collect data from all processors
        MPI_Status status;
        for (int rank=0; rank<nproc; rank++) {
            if (rank!=myrank) {
#ifdef PREC4
                MPI_Recv(rho1, length,MPI_FLOAT, rank, 1, cartesian, &status);
                MPI_Recv(ux1, length, MPI_FLOAT, rank, 2, cartesian, &status);
                MPI_Recv(uy1, length, MPI_FLOAT, rank, 3, cartesian, &status);
                MPI_Recv(uz1, length, MPI_FLOAT, rank, 4, cartesian, &status);
#else
                MPI_Recv(rho1, length,MPI_DOUBLE, rank, 1, cartesian, &status);
                MPI_Recv(ux1, length, MPI_DOUBLE, rank, 2, cartesian, &status);
                MPI_Recv(uy1, length, MPI_DOUBLE, rank, 3, cartesian, &status);
                MPI_Recv(uz1, length, MPI_DOUBLE, rank, 4, cartesian, &status);
#endif
                size_t memoffset=rank*(lz)*(ny)*(nx);
                std::memcpy(&rhg[memoffset], rho1, length*sizeof(real));
                std::memcpy(&uxg[memoffset], ux1, length*sizeof(real));
                std::memcpy(&uyg[memoffset], uy1, length*sizeof(real));
                std::memcpy(&uzg[memoffset], uz1, length*sizeof(real));
            }
        }
    }
    else {
#ifdef PREC4
        MPI_Send(rho1, length, MPI_FLOAT, 0, 1, cartesian);
        MPI_Send(ux1, length, MPI_FLOAT, 0, 2, cartesian);
        MPI_Send(uy1, length, MPI_FLOAT, 0, 3, cartesian);
        MPI_Send(uz1, length, MPI_FLOAT, 0, 4, cartesian);
#else
        MPI_Send(rho1, length,MPI_DOUBLE, 0, 1, cartesian);
        MPI_Send(ux1, length, MPI_DOUBLE, 0, 2, cartesian);
        MPI_Send(uy1, length, MPI_DOUBLE, 0, 3, cartesian);
        MPI_Send(uz1, length, MPI_DOUBLE, 0, 4, cartesian);
#endif
    }

    delete [] rho1;
    delete [] ux1;
    delete [] uy1;
    delete [] uz1;
}

// ==========================================================================
// Combine all f's into single f array. Used to dump state into binary
// checkpoint file.
// ==========================================================================
void collect_f(real *fg, const real *f)
{
    // fg is of the size nx*ny*nz*19
    //  f is of the size (lx+2)*(ly+2)*(lz+2)*19
    real *temp_f0;
    const int length2=(lx+2)*(ly+2)*(lz+2);
    const int length0=(lx+0)*(ly+0)*(lz+0);

    temp_f0=new real[DVM::npop*length0];

    //
    for (int k=0; k<lz; k++) {
        for (int j=0; j<ly; j++) {
            for (int i=0; i<lx; i++) {
                for (int ip=0; ip<DVM::npop; ip++) {
                    int ipijk0=k*ly*lx*DVM::npop+j*lx*DVM::npop+i*DVM::npop+ip;
                    int ipijk2=(k+1)*(ly+2)*(lx+2)*DVM::npop+(j+1)*(lx+2)*DVM::npop+(i+1)*DVM::npop+ip;

                    temp_f0[ipijk0]=f[ipijk2];
                }
            }
        }
    }

    if (!root()) {
        MPI_Send(temp_f0, DVM::npop*length0, MY_MPI_REAL, 0, myrank, cartesian);
    }
    else {
        std::memcpy(fg, temp_f0, DVM::npop*length0*sizeof(real));

        // Now collect data from all processors
        MPI_Status status;
        for (int rank=0; rank<nproc; rank++) {
            if (rank!=myrank) {
                MPI_Recv(temp_f0, DVM::npop*length0, MY_MPI_REAL, rank, rank, cartesian, &status);

                size_t memoffset=rank*(lx*ly*lz)*DVM::npop;
                std::memcpy(&fg[memoffset], temp_f0, DVM::npop*length0*sizeof(real));
            }
        }
    }

    delete [] temp_f0;
}
