#include <cuda.h>
#include <curand_kernel.h>
#include <stdio.h>

__global__ void setup_kernel ( curandState * state, unsigned long seed )
{
        int id = threadIdx.x;
        curand_init ( seed, id, 0, &state[id] );
} 

__global__ void generate( curandState* globalState, float* dev_storage ) 
{
        int ind = threadIdx.x;
        curandState localState = globalState[ind];
        float stable_storage = curand_uniform( &globalState )[ind];
        dev_storage = &stable_storage;
        globalState[ind] = localState; 
}

int main( int argc, char** argv) 
{ 
        int N = 1;
        curandState* devStates;
        cudaMalloc ( &devStates, N*sizeof( curandState ) );
        
        float *host_storage; 
        float *dev_storage;

        host_storage = (float *) malloc(sizeof(float));

        *host_storage = -1.0;
        printf("host_storage before = %f\n", *host_storage);
        cudaMalloc(&dev_storage, sizeof(*dev_storage));

        // setup seeds
        setup_kernel <<< 1, 1 >>> ( devStates, time(NULL) );

        // generate random numbers
        generate <<< 1, 1 >>> ( devStates, dev_storage );

        cudaMemcpy(host_storage, dev_storage, sizeof(*host_storage), cudaMemcpyDeviceToHost);

        printf("host_storage = %f\n", *host_storage);

        return 0;
}
