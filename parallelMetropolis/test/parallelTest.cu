#include "parallelTest.cuh"

void setupGetXFromIndex(){
    int numberOfBlocks = 3;
    int threadsPerBlock = 2;

    int *xValues;
    int *xValues_device;
    
    size_t xSize = numberOfBlocks * threadsPerBlock * sizeof(int);
    xValues = (int *)malloc(xSize);
    cudaMalloc((void **) &xValues_device, xSize);
    
    testGetXKernel <<<numberOfBlocks, threadsPerBlock>>>(xValues_device);

    cudaMemcpy(xValues, xValues_device, xSize, cudaMemcpyDeviceToHost);

    assert(xValues[0] == 1);
    assert(xValues[1] == 2);
    assert(xValues[2] == 2);
    assert(xValues[3] == 3);
    assert(xValues[4] == 3);
    assert(xValues[5] == 3);

    printf("getXFromIndex Correct\n");

    cudaFree(xValues_device);
    free(xValues);
}


void setupGetYFromIndex(){
    //TODO
}


void setupMakePeriodic(){
    //TODO
}



void setupWrapBox(){
    //TODO
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

    //assert that all atoms positions are in range of the box
    for (int i = 0; i < numberOfAtoms; i++){
        double dim_x = atoms[i].x;
        double dim_y = atoms[i].y;
        double dim_z = atoms[i].z;
        
        printf("%f, %f, %f\n", dim_x, dim_y, dim_z);

        assert(dim_x >= 0.0 && dim_x <= enviro->x &&
               dim_y >= 0.0 && dim_y <= enviro->y &&
               dim_z >= 0.0 && dim_z <= enviro->z);
    }
    printf("testGeneratePoints successful.");
}

void testCalcEnergy(){
    //TODO    
}

int main(){
    testGeneratePoints();    
    return 0;
}
