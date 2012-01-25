#include <cuda.h>
#include <curand_kernel.h>
#include <stdio.h>

__global__ void setup_kernel ( curandState * state, unsigned long seed )
{
        int id = threadIdx.x;
        curand_init ( seed, id, 0, &state[id] );
} 

__global__ void generate( curandState* globalState, double* dev_storage ) 
{
        int ind = threadIdx.x;
        curandState localState = globalState[ind];
        double stable_storage = curand_uniform_double( &localState );
        dev_storage = &stable_storage;
        globalState[ind] = localState; 
}

int main( int argc, char** argv) 
{ 
        int N = 1;
        curandState* devStates;
        cudaMalloc ( &devStates, N*sizeof( curandState ) );
        
        double *host_storage; 
        double *dev_storage;

        host_storage = (double *) malloc(sizeof(double));

        *host_storage = -1.0;
        cudaMalloc(&dev_storage, sizeof(*dev_storage));

        // setup seeds
        setup_kernel <<< 1, 1 >>> ( devStates, time(NULL) );

        // generate random numbers
        generate <<< 1, 1 >>> ( devStates, dev_storage );

        cudaMemcpy(host_storage, dev_storage, sizeof(*host_storage), cudaMemcpyDeviceToHost);

        printf("host_storage = %f", *host_storage);

        return 0;
}
