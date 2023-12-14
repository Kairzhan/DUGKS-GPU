#include <iostream>

#include "precision.h"

//***************************************************************************
// D3Q19 Discrete Velocity Model (DVM) is used in this code.
//***************************************************************************
namespace DVM {

    extern const int npop=19;
    //         0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
    extern const int cix[]={0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
    extern const int ciy[]={0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1};
    extern const int ciz[]={0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1};
    extern const int opp[]={0,  2,  1,  4,  3,  6,  5, 10,  9,  8,  7, 14, 13, 12, 11, 18, 17, 16, 15};

    extern const real ww0 = 1.0/3.0;
    extern const real ww1 = 1.0/18.0;
    extern const real ww2 = 1.0/36.0;

    extern const real RT=1./3.;
}
