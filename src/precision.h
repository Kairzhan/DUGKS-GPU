#ifndef PRECISION_H_
#define PRECISION_H_

#ifdef PREC4
typedef float real;
#else
typedef double real;
#endif

#ifdef PREC4
#define MY_MPI_REAL MPI_FLOAT
#else
#define MY_MPI_REAL MPI_DOUBLE
#endif

#endif
