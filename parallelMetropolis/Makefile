NV=nvcc
SRC=src/
UTIL=../Utilities/src/
TST=test/
TSTEXE=parallelTest
FLAGS=-arch sm_13 -lcurand
EXE=kryptonSim
DEMOEXE=parallelDemo

all: tests kryptonSim

demo: cudaUtil metroUtil $(SRC)inClassDemo.cu
	$(NV) $(FLAGS) metroCudaUtil.o metroParallelUtil.o $(SRC)inClassDemo.cu -o $(DEMOEXE)	

kryptonSim: cudaUtil metroUtil $(SRC)kryptonSimulation.cu
	$(NV) $(FLAGS) metroCudaUtil.o metroParallelUtil.o $(SRC)kryptonSimulation.cu -o $(EXE)	

tests: cudaUtil metroUtil baseTest $(TST)parallelTest.cu $(TST)parallelTest.cuh
	$(NV) $(FLAGS) baseTests.o metroCudaUtil.o  metroParallelUtil.o $(TST)parallelTest.cu -o $(TSTEXE)

baseTest: $(TST)baseTests.h $(TST)baseTests.cpp
	$(NV) $(FLAGS) -c $(TST)baseTests.cpp

cudaUtil: metroUtil $(SRC)metroCudaUtil.cuh $(SRC)metroCudaUtil.cu
	$(NV) $(FLAGS) -c $(SRC)metroCudaUtil.cu

metroUtil: $(UTIL)metroUtil.h $(UTIL)metroUtil.cpp
	$(NV) $(FLAGS) -c $(UTIL)metroUtil.cpp
	
clean:
	rm -f -R *.o
	rm -f $(TSTEXE)
	rm -f $(EXE)
	rm -f $(DEMOEXE)
