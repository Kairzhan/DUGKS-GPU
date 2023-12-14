#include <iostream>
#include <mpi.h>

#include "config.h"

int device_id;
cudaDeviceProp prop;

//===========================================================================
// Performs basic configurations for GPU. Single-node GPU configuration is
// assumed, so that each MPI thread will select one GPU.
//===========================================================================
void init_cuda(int myrank)
{
    Config config;
    
    const int gpusPerThread=std::stoi(config.get("gpusPerThread"));;    
    device_id=myrank%gpusPerThread;
    
    cudaSetDevice(device_id);
    cudaGetDeviceProperties(&prop, device_id);

    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    for (int device = 0; device < deviceCount; ++device) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, device);
        printf("Device %d has compute capability %d.%d.\n",
           device, deviceProp.major, deviceProp.minor);
    }
    
    return;
}
