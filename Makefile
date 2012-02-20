NV=nvcc

SRC=src/
TST=test/
UTIL=Utilities/
PARA=parallelMetropolis/
LIN=LinearMetropolis/
BIN=bin/

FLAGS=-arch sm_13 -lcurand

TSTEXE=parallelTest
EXE=kryptonSim
DEMOEXE=parallelDemo
UTILTESTEXE=utilTest

all: dir $(BIN)tests $(BIN)kryptonSim $(BIN)utilTests

$(BIN)demo: $(BIN)cudaUtil $(BIN)metroUtil $(PARA)$(SRC)inClassDemo.cu
	$(NV) $(FLAGS) $(BIN)metroCudaUtil.o $(BIN)metroUtil.o $(PARA)$(SRC)inClassDemo.cu -o $(BIN)$(DEMOEXE)	

$(BIN)kryptonSim: $(BIN)cudaUtil $(BIN)metroUtil $(PARA)$(SRC)kryptonSimulation.cu
	$(NV) $(FLAGS) $(BIN)metroCudaUtil.o $(BIN)metroUtil.o $(PARA)$(SRC)kryptonSimulation.cu -o $(BIN)$(EXE)	

$(BIN)tests: $(BIN)cudaUtil $(BIN)metroUtil $(BIN)baseTest $(PARA)$(TST)parallelTest.cu $(PARA)$(TST)parallelTest.cuh
	$(NV) $(FLAGS) $(BIN)baseTests.o $(BIN)metroCudaUtil.o  $(BIN)metroUtil.o $(PARA)$(TST)parallelTest.cu -o $(BIN)$(TSTEXE)

$(BIN)baseTest: $(PARA)$(TST)baseTests.h $(PARA)$(TST)baseTests.cpp
	$(NV) $(FLAGS) -c $(PARA)$(TST)baseTests.cpp -o $(BIN)baseTests.o

$(BIN)cudaUtil: $(BIN)metroUtil $(PARA)$(SRC)metroCudaUtil.cuh $(PARA)$(SRC)metroCudaUtil.cu
	$(NV) $(FLAGS) -c $(PARA)$(SRC)metroCudaUtil.cu -o $(BIN)metroCudaUtil.o

$(BIN)metroUtil: $(BIN)OPLSScan $(BIN)zMatrix $(UTIL)$(SRC)metroUtil.h $(UTIL)$(SRC)metroUtil.cpp
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)metroUtil.cpp -o $(BIN)metroUtil.o

$(BIN)zMatrix: $(UTIL)$(SRC)Zmatrix_Scan.cpp $(UTIL)$(SRC)Zmatrix_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)Zmatrix_Scan.cpp -o $(BIN)Zmatrix_Scan.o

$(BIN)OPLSScan: $(UTIL)$(SRC)Opls_Scan.cpp $(UTIL)$(SRC)Opls_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)Opls_Scan.cpp -o $(BIN)Opls_Scan.o

$(BIN)utilTests: $(BIN)metroUtil $(BIN)zMatrix $(BIN)OPLSScan
	$(NV) $(BIN)metroUtil.o $(BIN)Zmatrix_Scan.o $(BIN)Opls_Scan.o Utilities/test/utilityTests.cpp -o $(BIN)$(UTILTESTEXE)

dir:
	mkdir -p $(BIN)

clean:
	rm -f $(BIN)*
