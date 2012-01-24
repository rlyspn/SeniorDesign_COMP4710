NV=nvcc
SRC=src/
TST=test/

cudaUtil: metroUtil $(SRC)metroCudaUtil.cuh $(SRC)metroCudaUtil.cu
	$(NV) -c $(SRC)metroCudaUtil.cu

metroUtil: $(SRC)metroParallelUtil.h $(SRC)metroParallelUtil.cpp
	$(NV) -c $(SRC)metroParallelUtil.cpp

clean:
	rm -f -R *.o