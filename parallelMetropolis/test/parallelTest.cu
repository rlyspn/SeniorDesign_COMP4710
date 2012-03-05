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

bool compareDouble(double a, double b, double limit){
    if((a - b) / b < limit)
        return true;
    else
        return false;
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


void testWrapBox(){
 srand(time(NULL));
    int numberOfTests = 128;
    double box;
    
    double *testDoubles;
    size_t inputSize = sizeof(double) * numberOfTests;

    box = 10.f;
    testDoubles = (double *) malloc(inputSize);
    
    //generate random numbers
    for(int i = 0; i < numberOfTests; i++){
        testDoubles[i] = ((double) (rand() % 100));
    }

     //check that values are the same as known correct function
    for(int i = 0; i < numberOfTests; i++){
        double test_output = wrap_into_box(testDoubles[i], box);
        assert(wrapBox(testDoubles[i], box) == test_output);
    }

    printf("wrapBox passed Tests\n");

    free(testDoubles);


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
    printf("\nparallelEnergy = %2.10f\nlinearEnergy = %2.10f\n", *energy, baseEnergy); 
    printf("Calc_lj is correct\n");
    free(atoms);
    free(energy);
    cudaFree(atoms_device);
    cudaFree(enviro_device);
    cudaFree(energy_device);
}


void testGeneratePoints(){
    //init atoms, environment
    int numberOfAtoms = 3000;
    Atom *atoms = (Atom *) malloc(numberOfAtoms * sizeof(Atom));

    for (int i = 0; i < numberOfAtoms; i++){
        atoms[i] = createAtom(i, -1.0, -1.0, -1.0);
    }
    Environment enviro = createEnvironment(5.0, 10.0, 15.0, 1.0, 298.15, numberOfAtoms);

    generatePoints(atoms, &enviro);

    //assert that all atoms positions are in range of the box
    for (int i = 0; i < numberOfAtoms; i++){
        double dim_x = atoms[i].x;
        double dim_y = atoms[i].y;
        double dim_z = atoms[i].z;

        assert(dim_x >= 0.0 && dim_x <= enviro.x &&
               dim_y >= 0.0 && dim_y <= enviro.y &&
               dim_z >= 0.0 && dim_z <= enviro.z);
    }
    printf("testGeneratePoints successful.\n");

    free(atoms);
}

void testCalcEnergy(){
    // the sigma value of krypton used in the LJ simulation
    double kryptonSigma = 3.624;
    // the epsilon value of krypton used in the LJ simulation
    double kryptonEpsilon = 0.317;

    struct timeval le_tvBegin, le_tvEnd, pl_tvBegin, pl_tvEnd;

    //Generate enviorment and atoms
    int numberOfAtoms = 1000;
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
    printf("Linear Total Energy:   %f \n", te_linear);
    printf("In %d ms\n", le_runTime);
    printf("Parallel Total Energy: %f \n", te_parallel);
    printf("In %d ms\n", pl_runTime);
    assert(compareDouble(te_linear, te_parallel, .05));
    printf("testCalcEnergy successful.");

    
}

void testCalcEnergyWithMolecules(){
    // the sigma value of krypton used in the LJ simulation
    double kryptonSigma = 3.624;
    // the epsilon value of krypton used in the LJ simulation
    double kryptonEpsilon = 0.317;

    struct timeval le_tvBegin, le_tvEnd, pl_tvBegin, pl_tvEnd;

    //Generate enviorment and atoms
    int numberOfAtoms = 100;
    Environment stableEnviro = createEnvironment(5.0, 10.0, 15.0, 1.0, 298.15, numberOfAtoms);

    Environment *enviro = &stableEnviro;

    Atom *atoms = new Atom[numberOfAtoms];
    for (int i = 0; i < numberOfAtoms; i++){
        atoms[i] = createAtom(i, -1.0, -1.0, -1.0, kryptonSigma, kryptonEpsilon);
    }
    enviro->numOfMolecules = numberOfAtoms;
    generatePoints(atoms, enviro);
    Molecule *molecules;
    molecules = (Molecule *)malloc(sizeof(Molecule) * numberOfAtoms);
    for(int i = 0; i < numberOfAtoms; i++){
        molecules[i].numOfAtoms = 1;
        molecules[i].atoms = (Atom *)malloc(sizeof(Atom));
        molecules[i].atoms[0] = atoms[i];
    }

    //calculate energy linearly
    gettimeofday(&le_tvBegin,NULL); //start clock for execution time

    double te_linear = calculate_energy(atoms, enviro);

    gettimeofday(&le_tvEnd,NULL); //stop clock for execution time
    long le_runTime = timevaldiff(&le_tvBegin,&le_tvEnd); //get difference in time in milli seconds

    //calculate energy in parallel
    gettimeofday(&pl_tvBegin,NULL); //start clock for execution time

    double te_parallel =  calcEnergyWrapper(molecules, enviro);	 

    gettimeofday(&pl_tvEnd,NULL); //start clock for execution time
    long pl_runTime = timevaldiff(&pl_tvBegin,&pl_tvEnd); //get difference in time in milli seconds


    //Print out Results
    printf("Number of elements: %d\n", numberOfAtoms);
    printf("Linear Total Energy:   %f \n", te_linear);
    printf("In %d ms\n", le_runTime);
    printf("Parallel Total Energy: %f \n", te_parallel);
    printf("In %d ms\n", pl_runTime);
    assert(compareDouble(te_linear, te_parallel, .05));
    printf("testCalcEnergyWithMolecules successful.");

    
}

void testGetMoleculeFromIDWrapper(){
    int numberOfAtoms = 11;
    int numberOfMolecules = 3;
    
    Atom *atoms;
    Molecule *molecules;
    Environment enviro;
    int *answers;

    Atom *atoms_device;
    Molecule *molecules_device;
    int *answers_device;

    enviro.numOfAtoms = numberOfAtoms;
    enviro.numOfMolecules = numberOfMolecules;

    atoms = (Atom *)malloc(sizeof(Atom) * numberOfAtoms);
    molecules = (Molecule *)malloc(sizeof(Molecule) *numberOfMolecules);
    answers = (int *)malloc(sizeof(int) * numberOfAtoms);

    cudaMalloc((void **) &atoms_device, sizeof(Atom) * numberOfAtoms);
    cudaMalloc((void **) &molecules_device, sizeof(Molecule) * numberOfMolecules);
    cudaMalloc((void **) &answers_device, sizeof(int) * numberOfAtoms);

    enviro.numOfAtoms = numberOfAtoms;
    enviro.numOfMolecules = numberOfMolecules;
    
    for(int i = 0; i < numberOfAtoms; i++){
        atoms[i].id = i;
    }
    molecules[0].id = 0;
    molecules[1].id = 2;
    molecules[2].id = 6;


    cudaMemcpy(atoms_device, atoms, sizeof(Atom) * numberOfAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy(molecules_device, molecules, sizeof(Molecule) * numberOfMolecules, cudaMemcpyHostToDevice);

    int numberOfBlocks = 1;
    int threadsPerBlock = 128;
    testGetMoleculeFromID<<<numberOfBlocks,threadsPerBlock>>>(atoms_device,
            molecules_device, enviro, numberOfAtoms, answers_device);
   
    cudaMemcpy(answers, answers_device, sizeof(int) * numberOfAtoms, cudaMemcpyDeviceToHost);

    assert(answers[0] == 0);
    assert(answers[1] == 0);
    assert(answers[2] == 2);
    assert(answers[3] == 2);
    assert(answers[4] == 2);
    assert(answers[5] == 2);
    assert(answers[6] == 6);
    assert(answers[7] == 6);
    assert(answers[8] == 6);
    assert(answers[9] == 6);
    assert(answers[10] == 6);
   
    printf("getMoleculeFromID passed tests\n");

    free(atoms);
    free(molecules);
    free(answers);

    cudaFree(atoms_device);
    cudaFree(molecules_device);
    cudaFree(answers_device);


}


void testCalcBlendingWrapper(){
    double *d1, *d2, *d1_device, *d2_device, *answers, *answers_device;
    int numberOfTests = 5;
    size_t doubleSize = sizeof(double) * numberOfTests;

    d1 = (double *)malloc(doubleSize);
    d2 = (double *)malloc(doubleSize);
    answers = (double *)malloc(doubleSize);
    
    cudaMalloc((void **) &d1_device, doubleSize);
    cudaMalloc((void **) &d2_device, doubleSize);
    cudaMalloc((void **) &answers_device, doubleSize);

    d1[0] = 0.f;
    d2[0] = 0.f;

    d1[1] = 4.5;
    d2[1] = 2.32;
    
    d1[2] = 52.34;
    d2[2] = 5.f;


    d1[3] = 1.f;
    d2[3] = 7.f;

    d1[4] = 34.56;
    d2[4] = 12.7;
    
    cudaMemcpy(d1_device, d1, doubleSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d2_device, d2, doubleSize, cudaMemcpyHostToDevice);

    int blocks = 1;
    int threadsPerBlock = 64;

    testCalcBlending <<<blocks, threadsPerBlock>>>(d1_device, d2_device, answers_device, numberOfTests);

    cudaMemcpy(answers, answers_device, doubleSize, cudaMemcpyDeviceToHost);

    for(int i = 0 ; i < numberOfTests; i++){
        double expected = sqrt(d1[i] * d2[i]);
        assert(answers[i] / sqrt(d1[i] * d2[i]) < 0.01 || answers[i] == expected);
    }

    printf("calcBlending passed tests.\n");
    
    free(d1);
    free(d2);
    free(answers);
    cudaFree(d1_device);
    cudaFree(d2_device);
    cudaFree(answers_device);
}

int main(){
    testCalcBlendingWrapper();
    testGetMoleculeFromIDWrapper();
    testWrapBox();
    setupCalc_lj();
    setupGetIndexTest();
    setupMakePeriodic();
    testGeneratePoints();
    testCalcEnergy();
    testCalcEnergyWithMolecules();
    
    return 0;
}

