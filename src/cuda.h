extern int device_id;
extern cudaDeviceProp prop;

// Here we define explicitly the max. size
// of the shared memory a block can use.
#ifndef SHMEM_LIMIT
#define SHMEM_LIMIT 49152
#endif

void init_cuda(int myrank);
