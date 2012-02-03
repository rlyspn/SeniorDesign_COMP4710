#include "metroCudaUtil.cuh"

//const int THREADS_PER_BLOCK = 128;

//calculates X (larger indexed atom) for energy calculation based on index in atom array
__device__ int getXFromIndex(int idx){
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

//calculates Y (smaller indexed atom) for energy calculation based on index in atom array
__device__ int getYFromIndex(int x, int idx){
    return idx - (x * x - x) / 2;
}

//apply periodic boundaries
__device__ double makePeriodic(double x, double box){
    
    while(x < -0.5 * box){
        x += box;
    }

    while(x > 0.5 * box){
        x -= box;
    }

    return x;

}

//keep coordinates with box
__device__ double wrapBox(double x, double box){

    while(x > box){
        x -= box;
    }
    while(x < 0){
        x += box;
    }

    return x;
}

//calculate Lennard-Jones energy between two atoms
__device__ double calc_lj(Atom atom1, Atom atom2, Environment enviro){
    //store LJ constants locally
    double sigma = atom1.sigma;
    double epsilon = atom1.epsilon;
    
    //calculate difference in coordinates
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;

    //magic chemistry
    deltaX = makePeriodic(deltaX, enviro.x);
    deltaY = makePeriodic(deltaY, enviro.y);
    deltaZ = makePeriodic(deltaZ, enviro.z);

    const double r2 = (deltaX * deltaX) +
                      (deltaY * deltaY) + 
                      (deltaZ * deltaZ);

    const double sig2OverR2 = pow(sigma, 2) / r2;
    const double sig6OverR6 = pow(sig2OverR2, 3);
    const double sig12OverR12 = pow(sig6OverR6, 2);
    const double energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
    
    if (r2 == 0){
        return 0;
    }
    else{
        return energy;
    }
}

//asisgn positions for the atoms in parallel
__global__ void assignAtomPositions(double *dev_doubles, Atom *atoms, Environment *enviro){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    //for each randomly generated double...
    if (idx < enviro->numOfAtoms * 3){
        //...assign a single dimension of a single atom...
        int atomIndex = idx / 3;
        int dim_select = idx % 3;

        //...and scale the double to the boundary and map it to the dimension.
        if (dim_select == 0){
            atoms[atomIndex].x = enviro->x * dev_doubles[idx];
        }
        else if (dim_select == 1){
            atoms[atomIndex].y = enviro->y * dev_doubles[idx];
        }
        else{
            atoms[atomIndex].z = enviro->z * dev_doubles[idx];
        }
    }

}

//generate coordinate data for the atoms
void generatePoints(Atom *atoms, Environment *enviro){
    //setup CUDA storage
    int N = enviro->numOfAtoms * 3;
    curandGenerator_t generator;
    double *devDoubles;
    Atom *devAtoms;
    Environment *devEnviro;

    //allocate memory on device
    cudaMalloc((void**)&devDoubles, N * sizeof(double));
    cudaMalloc((void**)&devAtoms, enviro->numOfAtoms * sizeof(Atom));
    cudaMalloc((void**)&devEnviro, sizeof(Environment));

    //copy local data to device
    cudaMemcpy(devAtoms, atoms, enviro->numOfAtoms * sizeof(Atom), cudaMemcpyHostToDevice);
    cudaMemcpy(devEnviro, enviro, sizeof(Environment), cudaMemcpyHostToDevice);

    //generate doubles for all coordinates
    curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(generator, (unsigned int) time(NULL));
    curandGenerateUniformDouble(generator, devDoubles, N);

    //calculate number of blocks required
    int numOfBlocks = N / THREADS_PER_BLOCK + (N % THREADS_PER_BLOCK == 0 ? 0 : 1);

    //assign the doubles to the coordinates
    assignAtomPositions <<< numOfBlocks, THREADS_PER_BLOCK >>> (devDoubles, devAtoms, devEnviro);

    //copy the atoms back to host
    cudaMemcpy(atoms, devAtoms, enviro->numOfAtoms * sizeof(Atom), cudaMemcpyDeviceToHost);

    //cleanup
    curandDestroyGenerator(generator);
    cudaFree(devDoubles);
    cudaFree(devAtoms);
    cudaFree(devEnviro);
}
 
double calcEnergyWrapper(Atom *atoms, Environment *enviro){
    //setup CUDA storage
    double totalEnergy = 0.0;
    Atom *atoms_device;
    double *energySum_device;
    double *energySum_host;
    Environment *enviro_device;

    //calculate CUDA thread mgmt
    int N =(int) ( pow( (float) enviro->numOfAtoms,2)-enviro->numOfAtoms)/2;
    int blocks = N / THREADS_PER_BLOCK + (N % THREADS_PER_BLOCK == 0 ? 0 : 1);

    //The number of bytes of shared memory per block of
    //size_t sharedSize = sizeof(double) * THREADS_PER_BLOCK;
    size_t atomSize = enviro->numOfAtoms * sizeof(Atom);
    size_t energySumSize = blocks * sizeof(double);
    
    //allocate memory on the device
    energySum_host = (double *) malloc(energySumSize);
    cudaMalloc((void **) &atoms_device, atomSize);
    cudaMalloc((void **) &energySum_device, energySumSize);
    cudaMalloc((void **) &enviro_device, sizeof(Environment));

    for(int i = 0; i < blocks; i++){
        energySum_host[i] = 100.f;
    }

    //copy data to the device
    cudaMemcpy(atoms_device, atoms, atomSize, cudaMemcpyHostToDevice);
    cudaMemcpy(enviro_device, enviro, sizeof(Environment), cudaMemcpyHostToDevice);

    calcEnergy <<<blocks, THREADS_PER_BLOCK>>>(atoms_device, enviro_device, energySum_device);
    
    cudaMemcpy(energySum_host, energySum_device, energySumSize, cudaMemcpyDeviceToHost);

    for(int i = 0; i < blocks; i++){
        totalEnergy += energySum_host[i];
        //printf("energySum_host[%d] = %f\n", i, energySum_host[i]);
    }

    //cleanup
    cudaFree(atoms_device);
    cudaFree(energySum_device);
    free(energySum_host);

    return totalEnergy;
}

__global__ void calcEnergy(Atom *atoms, Environment *enviro, double *energySum){

	//need to figure out how many threads per block will be executed
	// must be a power of 2
    __shared__ double cache[THREADS_PER_BLOCK];	
	
	int cacheIndex = threadIdx.x;
    int idx =  blockIdx.x * blockDim.x + threadIdx.x;
	
    double lj_energy;
	
	int N =(int) ( pow( (float) enviro->numOfAtoms,2)-enviro->numOfAtoms)/2;
	
	if(idx < N ){
		//calculate the x and y positions in the Atom array
		int xAtom_pos, yAtom_pos;
		xAtom_pos =  getXFromIndex(idx);
		yAtom_pos =  getYFromIndex(xAtom_pos, idx);
		
		Atom xAtom, yAtom;
		xAtom = atoms[xAtom_pos];
		yAtom = atoms[yAtom_pos];
		
		lj_energy = calc_lj(xAtom,yAtom,*enviro);
    //    cache[cacheIndex] = lj_energy;
	}
	else {
		lj_energy = 0.0;
	}		
	

	 // set the cache values
    cache[cacheIndex] = lj_energy;
	
	 // synchronize threads in this block
    __syncthreads();
	 
	// adds 2 positions together
	int i = blockDim.x/2;
   while (i != 0) {
       if (cacheIndex < i)
           cache[cacheIndex] += cache[cacheIndex + i];
       __syncthreads();
       i /= 2;
   }
	
	// copy this block's sum to the enrgySums array
	// at its block index postition
	if (cacheIndex == 0)
        energySum[blockIdx.x] = cache[0];
		
}


#ifdef DEBUG

//these are all test wrappers for __device__ functions because they cannot be called from an outside source file.

__global__ void testWrapBoxKernel(double *x, double *box, int n){ 
    int idx =  threadIdx.x + blockIdx.x * blockDim.x;

    if (idx < n){
        x[idx] = wrapBox(x[idx], *box);
    }
}

__global__ void testMakePeriodicKernel(double *x, double *box, int n){ 
    int idx =  threadIdx.x + blockIdx.x * blockDim.x;

    if (idx < n){
        x[idx] = makePeriodic(x[idx], *box);
    }   
}

__global__ void testGetYKernel(int *xValues, int *yValues, int n){ 
    int idx =  threadIdx.x + blockIdx.x * blockDim.x;
    
    if (idx < n){
        yValues[idx] = getYFromIndex(xValues[idx], idx);
    }
}

__global__ void testGetXKernel(int *xValues, int n){
    int idx =  threadIdx.x + blockIdx.x * blockDim.x;
    
    if (idx < n){
        xValues[idx] = getXFromIndex(idx); 
    }
}

__global__ void testCalcLJ(Atom *atoms, Environment *enviro, double *energy){
    Atom atom1 = atoms[0];
    Atom atom2 = atoms[1];

    double testEnergy = calc_lj(atom1, atom2, *enviro);
    
    *energy = testEnergy;
}

#endif //DEBUG
