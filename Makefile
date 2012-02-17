NV=nvcc

SRC=src/
TST=test/
UTIL=Utilities/
PARA=parallelMetropolis/
LIN=LinearMetropolis/

FLAGS=-arch sm_13 -lcurand

TSTEXE=parallelTest
EXE=kryptonSim
DEMOEXE=parallelDemo
UTILTESTEXE=utilTest

all: tests kryptonSim utilTest

demo: cudaUtil metroUtil $(PARA)$(SRC)inClassDemo.cu
	$(NV) $(FLAGS) metroCudaUtil.o metroUtil.o $(PARA)$(SRC)inClassDemo.cu -o $(DEMOEXE)	

kryptonSim: cudaUtil metroUtil $(PARA)$(SRC)kryptonSimulation.cu
	$(NV) $(FLAGS) metroCudaUtil.o metroUtil.o $(PARA)$(SRC)kryptonSimulation.cu -o $(EXE)	

tests: cudaUtil metroUtil baseTest $(PARA)$(TST)parallelTest.cu $(PARA)$(TST)parallelTest.cuh
	$(NV) $(FLAGS) baseTests.o metroCudaUtil.o  metroUtil.o $(PARA)$(TST)parallelTest.cu -o $(TSTEXE)

baseTest: $(PARA)$(TST)baseTests.h $(PARA)$(TST)baseTests.cpp
	$(NV) $(FLAGS) -c $(PARA)$(TST)baseTests.cpp

cudaUtil: metroUtil $(PARA)$(SRC)metroCudaUtil.cuh $(PARA)$(SRC)metroCudaUtil.cu
	$(NV) $(FLAGS) -c $(PARA)$(SRC)metroCudaUtil.cu

metroUtil: OPLSScan zMatrix $(UTIL)$(SRC)metroUtil.h $(UTIL)$(SRC)metroUtil.cpp
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)metroUtil.cpp

zMatrix: $(UTIL)$(SRC)Zmatrix_Scan.cpp $(UTIL)$(SRC)Zmatrix_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)Zmatrix_Scan.cpp

OPLSScan: $(UTIL)$(SRC)Opls_Scan.cpp $(UTIL)$(SRC)Opls_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)Opls_Scan.cpp

utilTest: metroUtil zMatrix OPLSScan
	$(NV) metroUtil.o Zmatrix_Scan.o Opls_Scan.o $(UTIL)$(TST)utilityTests.cpp -o $(UTILTESTEXE)

clean:
	rm -f -R *.o
	rm -f $(TSTEXE)
	rm -f $(EXE)
	rm -f $(DEMOEXE)
	rm -f $(UTILTESTEXE)
