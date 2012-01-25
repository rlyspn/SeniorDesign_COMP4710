NV=nvcc
SRC=src/
TST=test/
TSTEXE=parallelTest

all: tests

tests: cudaUtil metroUtil baseTest

baseTest: $(TST)baseTests.h $(TST)baseTests.cpp
	$(NV) -c $(TST)baseTests.cpp

cudaUtil: metroUtil $(SRC)metroCudaUtil.cuh $(SRC)metroCudaUtil.cu
	$(NV) -c $(SRC)metroCudaUtil.cu

metroUtil: $(SRC)metroParallelUtil.h $(SRC)metroParallelUtil.cpp
	$(NV) -c $(SRC)metroParallelUtil.cpp
	
clean:
	rm -f -R *.o
	rm -f $(TSTEXE)
