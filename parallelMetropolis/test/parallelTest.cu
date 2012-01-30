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
__global__ void testGetXKernel(int *xValues){
    int idx =  threadIdx.x + blockIdx.x * blockDim.x;
    
    //xValues[idx] = getXFromIndex(idx); 
}

void setupGetYFromIndex(){
    //TODO
}
__global__ void testGetYKernel(){
    //TODO
}

void setupMakePeriodic(){
    //TODO
}

__global__ void testMakePeriodicKernel(){
    //TODO
}

void setupWrapBox(){
    //TODO
}
__global__ void testWrapBoxKernel(){
    //TODO    
}

void setupCalc_lj(){
    //TODO
}
__global__ void testCalcLJKernel(){
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

    //calculate size of atoms and environemnt structs
    size_t atomsSize = sizeof(*atoms);
    size_t enviroSize = sizeof(*enviro);

    //declare device structs
    Atom *dev_atoms = new Atom[numberOfAtoms];
    Environment *dev_enviro;

    //allocate memory for device structs
    cudaMalloc( (void**) &dev_atoms, atomsSize);
    cudaMalloc( (void**) &dev_enviro, enviroSize);

    //copy local structs to device structs on device
    cudaMemcpy(dev_atoms, atoms, atomsSize, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_enviro, enviro, enviroSize, cudaMemcpyHostToDevice);

    //allocate memory on device for random number generator state
    curandState* devStates;
    cudaMalloc ( &devStates, numberOfAtoms*sizeof( curandState ) );
                
    // setup seeds
    setup_generator <<<5, 2>>> ( devStates, time(NULL) );

    // generate random numbers
    generatePoints <<<5, 2>>> ( devStates, dev_atoms, dev_enviro );

    //copy atoms back to host
    cudaMemcpy(atoms, dev_atoms, atomsSize, cudaMemcpyDeviceToHost);

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
    //testGeneratePoints(); 
	 testCalcEnergy();   
    return 0;
}
