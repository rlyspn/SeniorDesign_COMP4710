#include "parallelTest.cuh"

/**
    Wrapper function that will call the global function used to 
    test the functions that calculate indexes in the half array
    used to hold the energies.
*/
void setupGetIndexTest(){
    int numberOfBlocks = 3;
    int threadsPerBlock = 2;
    int totalTests = numberOfBlocks * threadsPerBlock;

    int *xValues;
    int *yValues;
    int *yValues_device;
    int *xValues_device;
    
    size_t xSize = totalTests * sizeof(int);
    
    yValues = (int *) malloc(xSize);
    xValues = (int *)malloc(xSize);
    cudaMalloc((void **) &yValues_device, xSize);
    cudaMalloc((void **) &xValues_device, xSize);
    
    testGetXKernel <<<numberOfBlocks, threadsPerBlock>>>(xValues_device, totalTests);

    cudaMemcpy(xValues, xValues_device, xSize, cudaMemcpyDeviceToHost);

    assert(xValues[0] == 1);
    assert(xValues[1] == 2);
    assert(xValues[2] == 2);
    assert(xValues[3] == 3);
    assert(xValues[4] == 3);
    assert(xValues[5] == 3);

    printf("getXFromIndex Correct\n");

    //test getYFromIndex)
    testGetYKernel <<<numberOfBlocks, threadsPerBlock>>> (xValues_device,
            yValues_device, totalTests);

    cudaMemcpy(yValues, yValues_device, xSize, cudaMemcpyDeviceToHost);

    assert(yValues[0] == 0);
    assert(yValues[1] == 0);
    assert(yValues[2] == 1);
    assert(yValues[3] == 0);
    assert(yValues[4] == 1);
    assert(yValues[5] == 2);

    cudaFree(xValues_device);
    cudaFree(yValues_device);
    free(yValues);
    free(xValues);
}




void setupMakePeriodic(){
    //TODO
}



void setupWrapBox(){
    srand(time(NULL));
    int numberOfTests = 128;
    double box = 10;
    
    double *inputs_host;
    double *inputs_device;
    double *outputs_host;
    size_t inputSize = sizeof(double) * numberOfTests;

    inputs_host = (double *) malloc(inputSize);
    outputs_host = (double *) malloc(inputSize);
    cudaMalloc((void  **) &inputs_device, inputSize);
    
    //generate random numbers
    for(int i = 0; i < numberOfTests; i++){
        inputs_host[i] = ((double) rand());
    }

    //copy data to device
    cudaMemcpy(inputs_device, inputs_host, inputSize, cudaMemcpyHostToDevice);
    
    int threadsPerBlock = numberOfTests / 2;
    int blocks = numberOfTests / threadsPerBlock +
        (numberOfTests % threadsPerBlock == 0 ? 0 : 1);

    //run test on device
    testWrapBoxKernel <<<blocks, threadsPerBlock>>> (inputs_device, &box, numberOfTests);

    cudaMemcpy(outputs_host, inputs_device, inputSize, cudaMemcpyDeviceToHost);

    //check that values are the same as known correct function
    for(int i = 0; i < numberOfTests; i++){
        assert(outputs_host[i] == wrap_into_box(inputs_host[i], box));
    }

    free(inputs_host);
    free(outputs_host);
    cudaFree(inputs_device);

    printf("__device__ wrapBox tested correctly");

}

void setupCalc_lj(){
    //TODO
}


void testGeneratePoints(){
    //init atoms, environment
    int numberOfAtoms = 10;
    Atom *atoms = new Atom[numberOfAtoms];
    for (int i = 0; i < numberOfAtoms; i++){
        atoms[i] = createAtom(i, -1.0, -1.0, -1.0);
    }
    Environment stableEnviro = createEnvironment(5.0, 10.0, 15.0, 1.0, 122.0, numberOfAtoms);

    Environment *enviro = &stableEnviro;

    generatePoints(atoms, enviro);

    printf("Box limits{x:%f, y:%f, z:%f}\n", enviro->x, enviro->y, enviro->z);

    //assert that all atoms positions are in range of the box
    for (int i = 0; i < numberOfAtoms; i++){
        double dim_x = atoms[i].x;
        double dim_y = atoms[i].y;
        double dim_z = atoms[i].z;

        
        printf("Atom %d{x:%f, y:%f, z:%f}\n", i, dim_x, dim_y, dim_z);

        assert(dim_x >= 0.0 && dim_x <= enviro->x &&
               dim_y >= 0.0 && dim_y <= enviro->y &&
               dim_z >= 0.0 && dim_z <= enviro->z);
    }
    printf("testGeneratePoints successful.\n");
}

void testCalcEnergy(){
    //TODO    
}

int main(){
    testGeneratePoints();    
    return 0;
}
