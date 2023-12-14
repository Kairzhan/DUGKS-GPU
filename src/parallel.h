#ifndef PARALLEL_H_
#define PARALLEL_H_

extern int nproc;
extern int myrank;

extern double timestart, timeend;
extern double timestart1, timeend1;

extern real xoffset;
extern real xoffset1;
extern real yoffset;
extern real yoffset1;
extern real zoffset;
extern real zoffset1;

extern int xprev, xnext;
extern int yprev, ynext;
extern int zprev, znext;

extern MPI_Comm cartesian;

extern const int ndims;
extern const int dims[];
extern const int periods[];

void init_parallel(int argc, char* argv[]);
void free_parallel(void );
bool root(void );

void exchange(real*);
void collect_ruvw(real *rhg, real *uxg, real *uyg, real *uzg,
		  real *rho, real *ux, real *uy, real *uz);
void collect_f(real *fg, const real *f);

#endif
