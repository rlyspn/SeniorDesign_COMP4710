#include <cuda.h>
#include <curand_kernel.h>
#include <stdio.h>

__global__ void generate( curandState* globalState, double* dev_storage,  int* seed ) 
{
        int idx = threadIdx.x + blockDim.x * blockIdx.x;

        int seed_2 = *seed;
        curand_init (1234, idx, 0, &globalState[idx] );

        //curandState localState = globalState[idx];
        double stable_storage = curand_uniform_double( &globalState[idx] );
        dev_storage = &stable_storage;
        //globalState[idx] = localState; 
}

int main( int argc, char** argv) 
{ 
         int seed = ( int) time(NULL);
        int N = 1;
        curandState* devStates;
         int *device_seed;

        printf("Seed: %d\n", seed);

        cudaMalloc ((void **) &device_seed, sizeof(*device_seed));
        cudaMalloc ((void **) &devStates, N*sizeof( curandState ) );
       
        cudaMemcpy(device_seed, &seed, sizeof(*device_seed), cudaMemcpyHostToDevice);

        double *host_storage; 
        double *dev_storage;

        host_storage = (double *) malloc(sizeof(double));

        *host_storage = -1.0;
        printf("host_storage before = %f\n", *host_storage);
        
        cudaMalloc((void **) &dev_storage, sizeof(double));

        // generate random numbers
        generate <<< 1, 1 >>> ( devStates, dev_storage, device_seed );

        cudaMemcpy(host_storage, dev_storage, sizeof(double), cudaMemcpyDeviceToHost);

        printf("host_storage = %f\n", *host_storage);

        return 0;
}
