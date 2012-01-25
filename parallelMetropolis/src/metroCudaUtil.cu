#include "metroCudaUtil.cuh"


__device__ int getXFromIndex(int idx){
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv;
}

__device__ int getYFromIndex(int x, int idx){
    return idx - x;
}

__device__ double makePeriodic(double x, const double box){
    
    while(x < -0.5 * box){
        x += box;
    }

    while(x > 0.5 * box){
        x -= box;
    }

    return x;

}

__device__ double wrapBox(double x, double box){

    while(x > box){
        x -= box;
    }
    while(x < 0){
        x += box;
    }

    return x;
}

__device__ double calc_lj(Atom atom1, Atom atom2, Environment enviro){
    double sigma = atom1.sigma;
    double epsilon = atom2.epsilon;
    
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;

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

    return energy;
}

__global__ void setup_generator(curandState *globalState, unsigned long seed)
{
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        curand_init(seed, idx, 0, &globalState[idx]);
} 

__global__ void generatePoints(curandState *globalState, Atom *atoms, Environment *enviro){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curandState localState = globalState[idx];

    atoms[idx].id = idx;
    atoms[idx].x = enviro->x * curand_uniform_double(&localState);
    atoms[idx].y = enviro->y * curand_uniform_double(&localState);
    atoms[idx].z = enviro->z * curand_uniform_double(&localState);
    globalState[idx] = localState; 

}

double calcEnergyWrapper(Atom *atoms, Environment enviro){
    double totalEnergy = 0.0;

    int N =(int) ( pow( (float) enviro.numOfAtoms,2)-enviro.numOfAtoms)/2;
    int threadsPerBlock = 128;
    int blocks = N / threadsPerBlock + (N % threadsPerBlock == 0 ? 0 : 1);

    //The number of bytes of shared memory per block of
    size_t sharedSize = sizeof(double) * threadsPerBlock;

    //allocate memory on the device
    Atom *atoms_device;
    double *energySum_device;
    double *energySum_host;

    size_t atomSize = enviro.numOfAtoms * sizeof(Atom);
    size_t energySumSize = blocks * sizeof(double);

    energySum_host = (double *) malloc(energySumSize);
    cudaMalloc((void **) &atoms_device, atomSize);
    cudaMalloc((void **) &energySum_device, energySumSize);

    cudaMemcpy(atoms_device, atoms, atomSize, cudaMemcpyHostToDevice);

    calcEnergy <<<blocks, threadsPerBlock, sharedSize>>>(atoms_device, enviro, energySum_device);
    
    cudaMemcpy(energySum_host, energySum_device, energySumSize, cudaMemcpyDeviceToHost);

    for(int i = 0; i < blocks; i++){
        totalEnergy += energySum_host[i];
    }

    cudaFree(atoms_device);
    cudaFree(energySum_device);

    return totalEnergy;
}

__global__ void calcEnergy(Atom *atoms, Environment enviro, double *energySum){

	//need to figure out how many threads per block will be executed
	// must be a power of 2
    extern __shared__ double cache[];	
	
	int cacheIndex = threadIdx.x;
	
    
    int idx =  blockIdx.x * blockDim.x + threadIdx.x;
	double lj_energy;
	
	int N =(int) ( pow( (float) enviro.numOfAtoms,2)-enviro.numOfAtoms)/2;
	
	if(idx < N ){
		//calculate the x and y positions in the Atom array
		int xAtom_pos, yAtom_pos;
		xAtom_pos =  getXFromIndex(idx);
		yAtom_pos =  getYFromIndex(xAtom_pos, idx);
		
		Atom xAtom, yAtom;
		xAtom = atoms[xAtom_pos];
		yAtom = atoms[yAtom_pos];
		
		lj_energy = calc_lj(xAtom,yAtom,enviro);	
	}
	else {
		lj_energy = 0;
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
