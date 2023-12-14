#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <mpi.h>

#include "precision.h"
#include "numerics.h"
#include "dugks.h"
#include "parallel.h"
#include "lbm_d3q19.h"

using namespace Numerics;

void save_tecplot_serial(const real rhog[], const real uxg[], const real uyg[], const real uzg[], const int nx, const int ny, const int nz, const int step);
void statistics_serial(const real rhg[], const real uxg[], const real uyg[], const real uzg[], const int nx, const int ny, const int nz, const int step);

//===========================================================================
// Save/Load data to continue the simulation
//===========================================================================
void dump(int step, const real *f)
{
    real *fg = nullptr;

    if (root()) {
        fg=new real[DVM::npop*nx*ny*nz];
    }

    collect_f(fg, f);

    if (root()) {
        char temp[16];
        sprintf(temp, "ddf%8.8d.dat\0", step);
        std::ofstream ostrm(temp, std::ios::binary);
        ostrm.write((char *)fg, DVM::npop*nx*ny*nz*sizeof(real));
        ostrm.close();

        delete [] fg;
    }
}

//===========================================================================
// Save/Load data to continue the simulation
//===========================================================================
void read_dump(int step, real *f)
{
    real *temp_f = new real[DVM::npop*lx*ly*lz];

    if (root()) {
        char temp[16];
        sprintf(temp, "ddf%8.8d.dat\0", step);
        std::ifstream istrm(temp, std::ios::binary);

        real *fg = new real[DVM::npop*nx*ny*nz];
        istrm.read((char *)fg, DVM::npop*(nx)*(ny)*(nz)*sizeof(real));

        std::memcpy(temp_f, fg, DVM::npop*lx*ly*lz*sizeof(real));

        for (int rank=0; rank<nproc; rank++) {
            if (rank!=myrank) {
                size_t memoffset = rank*DVM::npop*lx*ly*lz;
                MPI_Send(&fg[memoffset], DVM::npop*lx*ly*lz, MY_MPI_REAL, rank, rank, cartesian);
            }
        }

        delete [] fg;
    }
    else {
        MPI_Status status;

        MPI_Recv(temp_f, DVM::npop*lx*ly*lz, MY_MPI_REAL, 0, myrank, cartesian, &status);
    }

    // copy to local array which has halo layer
    for (int k=0; k<lz; k++) {
        for (int j=0; j<ly; j++) {
            for (int i=0; i<lx; i++) {
                for (int ip=0; ip< DVM::npop; ip++) {
                    int ijk0=k*lx*ly*DVM::npop+j*lx*DVM::npop+i*DVM::npop+ip;
                    int ijk2=(k+1)*(lx+2)*(ly+2)*DVM::npop+(j+1)*(lx+2)*DVM::npop+(i+1)*DVM::npop+ip;
                    f[ijk2]=temp_f[ijk0];
                }
            }
        }
    }

    delete [] temp_f;
}

// ==========================================================================
// Output simulation results in 3D ASCII Tecplot file
// ==========================================================================
void output(int step, real *rho, real *ux, real *uy, real *uz)
{
    std::ofstream file;

    file.open("ruvw.tec");

    file << "TITLE = \"RHO-UVW 3D-Volume Data\"" << std::endl;
    file << "VARIABLES = \"X\", \"Y\", \"Z\", \"RHO\", \"UX\", \"UY\", \"UZ\"" << std::endl;
    file << "ZONE I=" << nx << ", J=" << ny << ", K=" << nz << " F=POINT" << std::endl;

    for (int i=1; i<=nx; i++) {
        for (int j=1; j<=ny; j++) {
            for (int k=1; k<=nz; k++) {
                int ijk=i*(ly+2)*(lz+2)+j*(lz+2)+k;
                file << (i-0.5)*dx << " " << (j-0.5)*dy << " " << (k-0.5)*dz << " "
                     << rho[ijk] << " " << ux[ijk] << " " << uy[ijk] << " " << uz[ijk] << std::endl;
            }
        }
    }
}

// ==========================================================================
// Save 2D slices
// ==========================================================================
void output_2d(int step, real *rho, real *ux, real *uy, real *uz)
{
    real *rhg;
    real *uxg;
    real *uyg;
    real *uzg;

    const int nxnynz=nx*ny*nz;

    if (root()) {
        rhg=new real[nxnynz];
        uxg=new real[nxnynz];
        uyg=new real[nxnynz];
        uzg=new real[nxnynz];
    }

    collect_ruvw(rhg, uxg, uyg, uzg, rho, ux, uy, uz);

    if (root()) {
        save_tecplot_serial(rhg, uxg, uyg, uzg, nx, ny, nz, step);
        statistics_serial(rhg, uxg, uyg, uzg, nx, ny, nz, step);
    }

    if (root()) {
        delete [] rhg;
        delete [] uxg;
        delete [] uyg;
        delete [] uzg;
    }
}

// ==========================================================================
// Save 3D arrays to ASCII Tecplot files
// ==========================================================================
void save_tecplot_serial(const real rhg[],
                         const real uxg[], const real uyg[], const real uzg[],
                         const int nx, const int ny, const int nz, const int step)
{
    char temp[16];
    std::ofstream file;

    // x-slice
    sprintf(temp, "s2X%8.8d.tec\0", step);
    file.open(temp);
    file << "TITLE = \"Example: Simple 2D-Slice Data\"" << std::endl;
    file << "VARIABLES = \"X\", \"Y\", \"Z\", \"RHO\", \"UX\", \"UY\", \"UZ\"" << std::endl;
    file << "ZONE I=" << 1 << ", J=" << ny << ", K=" << nz << " F=POINT" << std::endl;

    for (int k=0; k<nz; k++) {
        for (int j=0; j<ny; j++) {
            int i=nx/2;
            int ijk=k*(ny+0)*(nx+0)+j*(nx+0)+i;
            file << (i+0.5)*dx << " " << (j+0.5)*dy << " " << (k+0.5)*dz << " "
                 << rhg[ijk] << " "
                 << uxg[ijk] << " "
                 << uyg[ijk] << " "
                 << uzg[ijk] << std::endl;
        }
    }
    file.close();

    // y-slice
    sprintf(temp, "s2Y%8.8d.tec\0", step);
    file.open(temp);
    file << "TITLE = \"Example: Simple 2D-Slice Data\"" << std::endl;
    file << "VARIABLES = \"X\", \"Y\", \"Z\", \"RHO\", \"UX\", \"UY\", \"UZ\"" << std::endl;
    file << "ZONE I=" << nx << ", J=" << 1 << ", K=" << nz << " F=POINT" << std::endl;

    for (int k=0; k<nz; k++) {
        int j=ny/2;
        for (int i=0; i<nx; i++) {
            int ijk=k*(ny+0)*(nx+0)+j*(nx+0)+i;
            file << (i+0.5)*dx << " " << (j+0.5)*dy << " " << (k+0.5)*dz << " "
                 << rhg[ijk] << " "
                 << uxg[ijk] << " "
                 << uyg[ijk] << " "
                 << uzg[ijk] << std::endl;
        }
    }
    file.close();

    // z-slice
    sprintf(temp, "s2Z%8.8d.tec\0", step);
    file.open(temp);
    file << "TITLE = \"Example: Simple 2D-Slice Data\"" << std::endl;
    file << "VARIABLES = \"X\", \"Y\", \"Z\", \"RHO\", \"UX\", \"UY\", \"UZ\"" << std::endl;
    file << "ZONE I=" << nx << ", J=" << ny << ", K=" << 1 << " F=POINT" << std::endl;

    for (int j=0; j<ny; j++) {
        int k=nz/2;
        for (int i=0; i<nx; i++) {
            int ijk=k*(ny+0)*(nx+0)+j*(nx+0)+i;
            file << (i+0.5)*dx << " " << (j+0.5)*dy << " " << (k+0.5)*dz << " "
                 << rhg[ijk] << " "
                 << uxg[ijk] << " "
                 << uyg[ijk] << " "
                 << uzg[ijk] << std::endl;
        }
    }
    file.close();
}

// ==========================================================================
// Compute basic statistics required for channel flow simulation
// ==========================================================================
void statistics_serial(const real rhg[],
                       const real uxg[], const real uyg[], const real uzg[],
                       const int nx, const int ny, const int nz, const int step)
{
    real rhomean, umean, vmean, wmean;
    real rhomin, uxmin, uymin, uzmin;
    real rhomax, uxmax, uymax, uzmax;

    rhomean=umean=vmean=wmean=0;
    rhomin=uxmin=uymin=uzmin=99999;
    rhomax=uxmax=uymax=uzmax=-99999;

    for (int k=0; k<nz; k++) {
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                int ijk=k*ny*nx+j*nx+i;
                rhomean+=rhg[ijk];
                umean+=uxg[ijk];
                vmean+=uyg[ijk];
                wmean+=uzg[ijk];

                if (rhg[ijk]<rhomin) rhomin=rhg[ijk];
                if (uxg[ijk]<uxmin) uxmin=uxg[ijk];
                if (uyg[ijk]<uymin) uymin=uyg[ijk];
                if (uzg[ijk]<uzmin) uzmin=uzg[ijk];

                if (rhg[ijk]>rhomax) rhomax=rhg[ijk];
                if (uxg[ijk]>uxmax) uxmax=uxg[ijk];
                if (uyg[ijk]>uymax) uymax=uyg[ijk];
                if (uzg[ijk]>uzmax) uzmax=uzg[ijk];
            }
        }
    }
    rhomean/=nx*ny*nz;
    umean/=nx*ny*nz;
    vmean/=nx*ny*nz;
    wmean/=nx*ny*nz;

    FILE *file;
    file=fopen("diag.dat", "a");
    fprintf(file, "%8.8d  % 7.7e % 7.7e % 7.7e % 7.7e  % 7.7e % 7.7e % 7.7e % 7.7e  % 7.7e % 7.7e % 7.7e % 7.7e\n", step, rhomean, umean, vmean, wmean,
            rhomin, uxmin, uymin, uzmin,
            rhomax, uxmax, uymax, uzmax);
    fclose(file);

    // Now compute profiles
    real uave[nx];
    real vave[nx];
    real wave[nx];
    real stress_xz[nx];
    real stress_xy[nx];
    real stress_yz[nx];
    real usq[nx];
    real vsq[nx];
    real wsq[nx];

    for (int i=0; i<nx; i++) {
        uave[i]=vave[i]=wave[i]=0;
        stress_xz[i]=stress_xy[i]=stress_yz[i]=0;
        usq[i]=vsq[i]=wsq[i]=0;

        for (int k=0; k<nz; k++) {
            for (int j=0; j<ny; j++) {
                int ijk=k*ny*nx+j*nx+i;

                uave[i]+=uxg[ijk];
                vave[i]+=uyg[ijk];
                wave[i]+=uzg[ijk];

                stress_xz[i]+=uxg[ijk]*uzg[ijk];
                stress_xy[i]+=uxg[ijk]*uyg[ijk];
                stress_yz[i]+=uyg[ijk]*uzg[ijk];

                usq[i]+=uxg[ijk]*uxg[ijk];
                vsq[i]+=uyg[ijk]*uyg[ijk];
                wsq[i]+=uzg[ijk]*uzg[ijk];
            }
        }

        // average
        uave[i]/=(ny*nz);
        vave[i]/=(ny*nz);
        wave[i]/=(ny*nz);

        stress_xz[i]/=(ny*nz);
        stress_xy[i]/=(ny*nz);
        stress_yz[i]/=(ny*nz);

        usq[i]/=(ny*nz);
        vsq[i]/=(ny*nz);
        wsq[i]/=(ny*nz);

        // final calc.
        usq[i]-=uave[i]*uave[i];
        vsq[i]-=vave[i]*vave[i];
        wsq[i]-=wave[i]*wave[i];

        stress_xz[i]-=uave[i]*wave[i];
        stress_xy[i]-=uave[i]*vave[i];
        stress_yz[i]-=vave[i]*wave[i];
    }

    file=fopen("profiles-2.dat", "a");
    fprintf(file, "%8.8d\n", step);
    for (int i=0; i<nx; i++) {
        fprintf(file, "% 7.7e % 7.7e % 7.7e   % 7.7e % 7.7e % 7.7e   % 7.7e % 7.7e % 7.7e \n",
                uave[i], vave[i], wave[i],
                stress_xz[i], stress_xy[i], stress_yz[i],
                usq[i], vsq[i], wsq[i]);
    }
    fclose(file);
}
