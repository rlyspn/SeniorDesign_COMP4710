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

    printf("getYFromIndex Correct.\n");

    cudaFree(xValues_device);
    cudaFree(yValues_device);
    free(yValues);
    free(xValues);
}



/**
  wrapper function for the __device__ makePeriodic function
*/
void setupMakePeriodic(){
    srand(time(NULL));
    int numberOfTests = 128;
    double *box;;
    
    double *inputs_host;
    double *inputs_device;
    double *outputs_host;
    double *dev_box;
    size_t inputSize = sizeof(double) * numberOfTests;

    box = (double *) malloc(sizeof(double));
    *box = 10.0;
    inputs_host = (double *) malloc(inputSize);
    outputs_host = (double *) malloc(inputSize);
    cudaMalloc((void **) &inputs_device, inputSize);
    cudaMalloc((void **) &dev_box, sizeof(double));
    
    //generate random numbers
    for(int i = 0; i < numberOfTests; i++){
        inputs_host[i] = ((double) (rand() % 100));
    }

    //copy data to device
    cudaMemcpy(inputs_device, inputs_host, inputSize, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_box, box, sizeof(double), cudaMemcpyHostToDevice);
    
    int threadsPerBlock = numberOfTests / 2;
    int blocks = numberOfTests / threadsPerBlock +
        (numberOfTests % threadsPerBlock == 0 ? 0 : 1);

    testMakePeriodicKernel <<< blocks, threadsPerBlock >>> (inputs_device,
            dev_box, numberOfTests);

    cudaMemcpy(outputs_host, inputs_device, inputSize, cudaMemcpyDeviceToHost);

    //check that values are the same as known correct function
    for(int i = 0; i < numberOfTests; i++){
        double test_output = make_periodic(inputs_host[i], *box);
        assert(outputs_host[i] == test_output);
    }

    printf("makePeriodic passed Tests\n");

    free(inputs_host);
    free(outputs_host);
    cudaFree(inputs_device);


}



void setupWrapBox(){
 srand(time(NULL));
    int numberOfTests = 128;
    double *box;;
    
    double *inputs_host;
    double *inputs_device;
    double *outputs_host;
    double *dev_box;
    size_t inputSize = sizeof(double) * numberOfTests;

    box = (double *) malloc(sizeof(double));
    *box = 10.0;
    inputs_host = (double *) malloc(inputSize);
    outputs_host = (double *) malloc(inputSize);
    cudaMalloc((void **) &inputs_device, inputSize);
    cudaMalloc((void **) &dev_box, sizeof(double));
    
    //generate random numbers
    for(int i = 0; i < numberOfTests; i++){
        inputs_host[i] = ((double) (rand() % 100));
    }

    //copy data to device
    cudaMemcpy(inputs_device, inputs_host, inputSize, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_box, box, sizeof(double), cudaMemcpyHostToDevice);
    
    int threadsPerBlock = numberOfTests / 2;
    int blocks = numberOfTests / threadsPerBlock +
        (numberOfTests % threadsPerBlock == 0 ? 0 : 1);

    testWrapBoxKernel <<< blocks, threadsPerBlock >>> (inputs_device,
            dev_box, numberOfTests);

    cudaMemcpy(outputs_host, inputs_device, inputSize, cudaMemcpyDeviceToHost);

    //check that values are the same as known correct function
    for(int i = 0; i < numberOfTests; i++){
        double test_output = wrap_into_box(inputs_host[i], *box);
        assert(outputs_host[i] == test_output);
    }

    printf("wrapBox passed Tests\n");

    free(inputs_host);
    free(outputs_host);
    cudaFree(inputs_device);


}

void setupCalc_lj(){
    double kryptonSigma = 3.624;
    double kryptonEpsilon = 0.317;
    int numberOfAtoms = 2;

    Atom *atoms = new Atom[numberOfAtoms];
    double *energy = (double *) malloc(sizeof(double));
    *energy = 1000.f;
    Atom *atoms_device;
    Environment *enviro_device;
    double *energy_device;

    cudaMalloc((void **) &atoms_device, sizeof(Atom) * numberOfAtoms);
    cudaMalloc((void **) &enviro_device, sizeof(Environment));
    cudaMalloc((void **) &energy_device, sizeof(double));

    Environment stableEnviro = createEnvironment(10, 10, 10, .5,
            298.15, numberOfAtoms);

    Environment *enviro = &stableEnviro;
    generatePoints(atoms, enviro);
    atoms[0].sigma = kryptonSigma;
    atoms[0].epsilon = kryptonEpsilon; 

    cudaMemcpy(atoms_device, atoms, sizeof(Atom) * numberOfAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy(enviro_device, enviro, sizeof(Environment), cudaMemcpyHostToDevice);

    testCalcLJ<<<1,1>>>(atoms_device, enviro_device, energy_device);

    cudaMemcpy(energy, energy_device, sizeof(double), cudaMemcpyDeviceToHost);

    double baseEnergy = calculate_energy(atoms, enviro);
    assert((int)(*energy * pow(10.f, 6.f)) == (int)( baseEnergy * pow(10.f,6.f))); 
   
    printf("Calc_lj is correct\n");
    free(atoms);
    free(energy);
    cudaFree(atoms_device);
    cudaFree(enviro_device);
    cudaFree(energy_device);
}


void testGeneratePoints(){
    //init atoms, environment
    int numberOfAtoms = 10;
    Atom *atoms = new Atom[numberOfAtoms];
    for (int i = 0; i < numberOfAtoms; i++){
        atoms[i] = createAtom(i, -1.0, -1.0, -1.0);
    }
    Environment stableEnviro = createEnvironment(5.0, 10.0, 15.0, 1.0, 298.15, numberOfAtoms);

    Environment *enviro = &stableEnviro;

    generatePoints(atoms, enviro);

    //assert that all atoms positions are in range of the box
    for (int i = 0; i < numberOfAtoms; i++){
        double dim_x = atoms[i].x;
        double dim_y = atoms[i].y;
        double dim_z = atoms[i].z;

        assert(dim_x >= 0.0 && dim_x <= enviro->x &&
               dim_y >= 0.0 && dim_y <= enviro->y &&
               dim_z >= 0.0 && dim_z <= enviro->z);
    }
    printf("testGeneratePoints successful.\n");
}

void testCalcEnergy(){
    // the sigma value of krypton used in the LJ simulation
    double kryptonSigma = 3.624;
    // the epsilon value of krypton used in the LJ simulation
    double kryptonEpsilon = 0.317;

    struct timeval le_tvBegin, le_tvEnd, pl_tvBegin, pl_tvEnd;

    //Generate enviorment and atoms
    int numberOfAtoms = 50;
    Environment stableEnviro = createEnvironment(5.0, 10.0, 15.0, 1.0, 298.15, numberOfAtoms);

    Environment *enviro = &stableEnviro;

    Atom *atoms = new Atom[numberOfAtoms];
    for (int i = 0; i < numberOfAtoms; i++){
        atoms[i] = createAtom(i, -1.0, -1.0, -1.0, kryptonSigma, kryptonEpsilon);
    }

    generatePoints(atoms, enviro);
     
    //calculate energy linearly
    gettimeofday(&le_tvBegin,NULL); //start clock for execution time

    double te_linear = calculate_energy(atoms, enviro);

    gettimeofday(&le_tvEnd,NULL); //stop clock for execution time
    long le_runTime = timevaldiff(&le_tvBegin,&le_tvEnd); //get difference in time in milli seconds

    //calculate energy in parallel
    gettimeofday(&pl_tvBegin,NULL); //start clock for execution time

    double te_parallel =  calcEnergyWrapper(atoms, enviro);	 

    gettimeofday(&pl_tvEnd,NULL); //start clock for execution time
    long pl_runTime = timevaldiff(&pl_tvBegin,&pl_tvEnd); //get difference in time in milli seconds


    //Print out Results
    printf("Number of elements: %d\n", numberOfAtoms);
    printf("Linear Total Energy: %f \n", te_linear);
    printf("In %d ms\n", le_runTime);
    printf("Parallel Total Energy: %f \n", te_parallel);
    printf("In %d ms\n", pl_runTime);
    assert((long long) (pow(10, 6) * te_parallel) == (long long) (pow(10, 6) * te_linear));
    printf("testCalcEnergy successful.");


}

int main(){
    setupCalc_lj();
    setupGetIndexTest();
    setupMakePeriodic();
    setupWrapBox();
    testGeneratePoints();
    testCalcEnergy();
    return 0;
}
