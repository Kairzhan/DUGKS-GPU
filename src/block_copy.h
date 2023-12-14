// copy memory block into shared memory

// internal nodes
for (int ip=0; ip<npop; ip++) {
    FBs(ip, i0, j0, k0)=FB(ip, i, j, k);
 }

// ///////////////////////////
if (i0==1) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0-1, j0, k0)=FB(ip, i-1, j, k);
    }
 }
if (i0==blockDim.x) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0+1, j0, k0)=FB(ip, i+1, j, k);
    }
 }

if (j0==1) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0, j0-1, k0)=FB(ip, i, j-1, k);
    }
 }
if (j0==blockDim.y) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0, j0+1, k0)=FB(ip, i, j+1, k);
    }
 }

if (k0==1) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0, j0, k0-1)=FB(ip, i, j, k-1);
    }
 }
if (k0==blockDim.z) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0, j0, k0+1)=FB(ip, i, j, k+1);
    }
 }

// ///////////////////////////
// ij
if (i0==1 && j0==1) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0-1, j0-1, k0)=FB(ip, i-1, j-1, k);
    }
 }
if (i0==blockDim.x && j0==1) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0+1, j0-1, k0)=FB(ip, i+1, j-1, k);
    }
 }

if (i0==1 && j0==blockDim.y) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0-1, j0+1, k0)=FB(ip, i-1, j+1, k);
    }
 }
if (i0==blockDim.x && j0==blockDim.y) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0+1, j0+1, k0)=FB(ip, i+1, j+1, k);
    }
 }

// ik
if (i0==1 && k0==1) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0-1, j0, k0-1)=FB(ip, i-1, j, k-1);
    }
 }
if (i0==blockDim.x && k0==1) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0+1, j0, k0-1)=FB(ip, i+1, j, k-1);
    }
 }

if (i0==1 && k0==blockDim.z) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0-1, j0, k0+1)=FB(ip, i-1, j, k+1);
    }
 }
if (i0==blockDim.x && k0==blockDim.z) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0+1, j0, k0+1)=FB(ip, i+1, j, k+1);
    }
 }

// jk
if (j0==1 && k0==1) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0, j0-1, k0-1)=FB(ip, i, j-1, k-1);
    }
 }
if (j0==blockDim.y && k0==1) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0, j0+1, k0-1)=FB(ip, i, j+1, k-1);
    }
 }

if (j0==1 && k0==blockDim.z) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0, j0-1, k0+1)=FB(ip, i, j-1, k+1);
    }
 }
if (j0==blockDim.y && k0==blockDim.z) {
    for (int ip=0; ip<npop; ip++) {
        FBs(ip, i0, j0+1, k0+1)=FB(ip, i, j+1, k+1);
    }
 }

// corner points
for (int ip=0; ip<npop; ip++) {
    if (i0==1) {
        if (j0==1 && k0==1)
            FBs(ip, i0-1, j0-1, k0-1)=FB(ip, i-1, j-1, k-1);
        if (j0==1 && k0==blockDim.z)
            FBs(ip, i0-1, j0-1, k0+1)=FB(ip, i-1, j-1, k+1);

        if (j0==blockDim.y && k0==1)
            FBs(ip, i0-1, j0+1, k0-1)=FB(ip, i-1, j+1, k-1);
        if (j0==blockDim.y && k0==blockDim.z)
            FBs(ip, i0-1, j0+1, k0+1)=FB(ip, i-1, j+1, k+1);
    }
    if (i0==blockDim.x) {
        if (j0==1 && k0==1)
            FBs(ip, i0+1, j0-1, k0-1)=FB(ip, i+1, j-1, k-1);
        if (j0==1 && k0==blockDim.z)
            FBs(ip, i0+1, j0-1, k0+1)=FB(ip, i+1, j-1, k+1);

        if (j0==blockDim.y && k0==1)
            FBs(ip, i0+1, j0+1, k0-1)=FB(ip, i+1, j+1, k-1);
        if (j0==blockDim.y && k0==blockDim.z)
            FBs(ip, i0+1, j0+1, k0+1)=FB(ip, i+1, j+1, k+1);
    }
 }
