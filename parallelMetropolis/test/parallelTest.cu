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

    printf("Testing getting x index\n");
    printf("%d\n", xValues[0]);
    printf("%d\n", xValues[1]);
    printf("%d\n", xValues[2]);
    printf("%d\n", xValues[3]);
    printf("%d\n", xValues[4]);
    printf("%d\n", xValues[5]);

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
    
    printf("Testing getting y index\n");
    printf("%d\n", yValues[0]);
    printf("%d\n", yValues[1]);
    printf("%d\n", yValues[2]);
    printf("%d\n", yValues[3]);
    printf("%d\n", yValues[4]);
    printf("%d\n", yValues[5]);


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
        inputs_host[i] = ((double) rand());
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
        
        printf("inputs_host[%d] = %f | outputs_host[%d] = %f | test_output = %f\n", i, inputs_host[i], i, outputs_host[i], test_output);
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
        double test_output = wrap_into_box(inputs_host[i], box);
        printf("inputs_host[%d] = %f | outputs_host[%d] = %f | test_output = %f\n", i, inputs_host[i], i, outputs_host[i], test_output);
        assert(outputs_host[i] == test_output);
    }

    free(inputs_host);
    free(outputs_host);
    cudaFree(inputs_device);

    printf("__device__ wrapBox tested correctly\n");

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
        
	struct timeval le_tvBegin, le_tvEnd, pl_tvBegin, pl_tvEnd;

	//Generate enviorment and atoms
	 int numberOfAtoms = 10;
	 Environment enviro = createEnvironment(5.0, 10.0, 15.0, 1.0, 122.0, numberOfAtoms);
	 
    Atom *atoms = new Atom[numberOfAtoms];
    for (int i = 0; i < numberOfAtoms; i++){
        atoms[i] = createAtom(i, rand()*enviro.x, rand()*enviro.y, rand()*enviro.z);
    }
	 
	 //make copies of enviornment and atoms for 
	 //parallel portion
	 Environment enviro2 = enviro;
	 
	 Atom *atoms2 = new Atom[numberOfAtoms];
	 memcpy(atoms2,atoms,numberOfAtoms*sizeof(Atom) );
	 
	 	 
	 /*
	 ** Run the Calculation as Linear.
	 */
	 
	  gettimeofday(&le_tvBegin,NULL); //start clock for execution time
	 
	  double te_linear = calculate_energy(atoms, enviro);
	  
	  gettimeofday(&le_tvEnd,NULL); //start clock for execution time
	  long le_runTime = timevaldiff(&le_tvBegin,&le_tvEnd); //get difference in time in milli seconds

	 	 	 
	 
	 /*
	 ** Run the Calculation as Parallel
	 */
	 
	 gettimeofday(&pl_tvBegin,NULL); //start clock for execution time
	  
	 double te_parallel =  calcEnergyWrapper(atoms2, enviro2);	 
	 
	 gettimeofday(&pl_tvEnd,NULL); //start clock for execution time
	 long pl_runTime = timevaldiff(&pl_tvBegin,&pl_tvEnd); //get difference in time in milli seconds

	 
	 /*
	 ** Print out Results
	 */
	 if( te_parallel == te_linear)
	    printf("testCalcEnergy sucessful\n Both total energies equate to the same value.\n");
	 else
	 	 printf("testCalcEnergy failed\n Both total energies equate to different values.\n");
	 
         printf("Number of elements: %d", numberOfAtoms);
	 printf("Linear Total Energy: %f \n", te_linear);
	 printf("In %d ms", le_runTime);
	 printf("Parallel Total Energy: %f \n", te_parallel);
	 printf("In %d ms", pl_runTime);

    
}

int main(){
    setupGetIndexTest();
    setupMakePeriodic();
    setupWrapBox();
    testGeneratePoints();
    return 0;
}
