#include <cuda.h>
#include <stdio.h>

__global__ void arrayMult(float *a, float *b, float *result, int N){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N){
        result[idx] = a[idx] * b[idx];
    }
}



int main(){
    float *vector_a, *vector_b;
    float *dev_a, *dev_b, *dev_result;
    float *result;

    int N = 16;
    size_t size = N * sizeof(*vector_a);
    vector_a = (float *) malloc(size);
    vector_b = (float *) malloc(size);
    result = (float *) malloc(size);

    cudaMalloc( (void **) &dev_a, size);
    cudaMalloc( (void **) &dev_b, size);
    cudaMalloc( (void **) &dev_result, size);
    int i;

    for (i = 0; i < N; i++){
        vector_a[i] = 2.f;
        vector_b[i] = 4.f;
    }

    cudaMemcpy(dev_a, vector_a, sizeof(float) * N, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, vector_b, sizeof(float) * N, cudaMemcpyHostToDevice);

    int blockSize = 4;
    int nBlocks = 4;

    arrayMult <<<nBlocks, blockSize>>> (dev_a, dev_b, dev_result, N);

    cudaMemcpy(result, dev_result, sizeof(float) * N, cudaMemcpyDeviceToHost);

    float dotProduct = 0;
    for (i = 0; i < N; i++){
        dotProduct = dotProduct + result[i];
    }
    printf("Vector_A: ");
    for (i = 0; i < N; i++){
        printf("%f,", vector_a[i]);
    }
    printf("\nVector_B: ");
    for (i = 0; i < N; i++){
        printf("%f,", vector_b[i]);
        
    }
    
    printf("\nResult = %f\n", dotProduct);
    
}
