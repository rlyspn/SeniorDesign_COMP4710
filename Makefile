NV=nvcc

SRC=src/
TST=test/
UTIL=Utilities/
PARA=parallelMetropolis/
LIN=LinearMetropolis/
BIN=bin/

FLAGS=-arch sm_13 -lcurand

TSTEXE=parallelTest
EXE=metroSim
KRYPTONEXE=kryptonSim
DEMOEXE=parallelDemo
UTILTESTEXE=utilTest

all: dir tests kryptonSim utilTests metroSim

demo: cudaUtil metroUtil $(PARA)$(SRC)inClassDemo.cu
	$(NV) $(FLAGS) $(BIN)metroCudaUtil.o $(BIN)metroUtil.o $(PARA)$(SRC)inClassDemo.cu -o $(BIN)$(DEMOEXE)	

metroSim: cudaUtil metroUtil zMatrix OPLSScan $(PARA)$(SRC)metropolisSimulation.cu
	#$(NV) $(FLAGS) $(BIN)*.o $(PARA)$(SRC)metropolisSimulation.cu -o $(BIN)$(EXE) 
	$(NV) $(FLAGS) -c $(PARA)$(SRC)metropolisSimulation.cu -o $(BIN)metropolisSimulation.o 
	$(NV) $(FLAGS) $(BIN)metroUtil.o $(BIN)Opls_Scan.o $(BIN)Zmatrix_Scan.o $(BIN)metropolisSimulation.o -o $(BIN)$(EXE) 

kryptonSim: cudaUtil metroUtil $(PARA)$(SRC)kryptonSimulation.cu
	$(NV) $(FLAGS) $(BIN)metroCudaUtil.o $(BIN)metroUtil.o $(PARA)$(SRC)kryptonSimulation.cu -o $(BIN)$(KRYPTONEXE)	

tests: cudaUtil metroUtil baseTest $(PARA)$(TST)parallelTest.cu $(PARA)$(TST)parallelTest.cuh
	$(NV) $(FLAGS) $(BIN)baseTests.o $(BIN)metroCudaUtil.o  $(BIN)metroUtil.o $(PARA)$(TST)parallelTest.cu -o $(BIN)$(TSTEXE)

baseTest: $(PARA)$(TST)baseTests.h $(PARA)$(TST)baseTests.cpp
	$(NV) $(FLAGS) -c $(PARA)$(TST)baseTests.cpp -o $(BIN)baseTests.o

cudaUtil: metroUtil $(PARA)$(SRC)metroCudaUtil.cuh $(PARA)$(SRC)metroCudaUtil.cu
	$(NV) $(FLAGS) -c $(PARA)$(SRC)metroCudaUtil.cu -o $(BIN)metroCudaUtil.o

metroUtil: OPLSScan zMatrix $(UTIL)$(SRC)metroUtil.h $(UTIL)$(SRC)metroUtil.cpp
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)metroUtil.cpp -o $(BIN)metroUtil.o

zMatrix: $(UTIL)$(SRC)Zmatrix_Scan.cpp $(UTIL)$(SRC)Zmatrix_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)Zmatrix_Scan.cpp -o $(BIN)Zmatrix_Scan.o

OPLSScan: $(UTIL)$(SRC)Opls_Scan.cpp $(UTIL)$(SRC)Opls_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)Opls_Scan.cpp -o $(BIN)Opls_Scan.o

utilTests: metroUtil zMatrix OPLSScan
	$(NV) $(BIN)metroUtil.o $(BIN)Zmatrix_Scan.o $(BIN)Opls_Scan.o Utilities/test/utilityTests.cpp -o $(BIN)$(UTILTESTEXE)

dir:
	mkdir -p $(BIN)

clean:
	rm -f $(BIN)*
