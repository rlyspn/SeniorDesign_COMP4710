#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>

main(){
    int i, n = 100;
    curandGenerator_t gen;
    double *devData, *hostData;

    hostData = (double *)calloc(n, sizeof(double));

    cudaMalloc((void **)&devData, n * sizeof(double));

    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);

    curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);

    curandGenerateUniformDouble(gen, devData, n);

    cudaMemcpy(hostData, devData, n * sizeof(double), cudaMemcpyDeviceToHost);

    for (i = 0; i < n; i++) {
        printf("%1.4f\n", hostData[i]);
    }

    curandDestroyGenerator(gen);
    cudaFree(devData);
    free(hostData);

    return 0;
}
