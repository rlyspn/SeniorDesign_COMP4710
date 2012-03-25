NV=nvcc
#NV=g++

SRC=src/
TST=test/
UTIL=Utilities/
PARA=parallelMetropolis/
LIN=LinearMetropolis/
BIN=bin/

FLAGS=-arch sm_13 -lcurand -g
#FLAGS=-g

TSTEXE=parallelTest
EXE=metroSim
KRYPTONEXE=kryptonSim
DEMOEXE=parallelDemo
UTILTESTEXE=utilTest

all: dir tests kryptonSim utilTests metroSim

demo: cudaUtil metroUtil $(PARA)$(SRC)inClassDemo.cu
	$(NV) $(FLAGS) $(BIN)metroCudaUtil.o $(BIN)metroUtil.o $(PARA)$(SRC)inClassDemo.cu -o $(BIN)$(DEMOEXE)	

metroSim: cudaUtil metroUtil stateScan configScan zMatrix OPLSScan baseTest $(PARA)$(SRC)metropolisSimulation.cu
	$(NV) $(FLAGS) -c $(PARA)$(SRC)metropolisSimulation.cu -o $(BIN)metropolisSimulation.o 
	$(NV) $(FLAGS) $(BIN)geometricUtil.o $(BIN)State_Scan.o $(BIN)metroUtil.o $(BIN)Config_Scan.o $(BIN)metroCudaUtil.o $(BIN)Opls_Scan.o $(BIN)Zmatrix_Scan.o $(BIN)metropolisSimulation.o $(BIN)baseTests.o -o $(BIN)$(EXE) 

kryptonSim: cudaUtil metroUtil baseTest $(PARA)$(SRC)kryptonSimulation.cu
	$(NV) $(FLAGS) $(BIN)metroCudaUtil.o $(BIN)metroUtil.o $(BIN)baseTests.o $(PARA)$(SRC)kryptonSimulation.cu -o $(BIN)$(KRYPTONEXE)	

tests: cudaUtil geoUtil metroUtil baseTest $(PARA)$(TST)parallelTest.cu $(PARA)$(TST)parallelTest.cuh
	$(NV) $(FLAGS) $(BIN)geometricUtil.o $(BIN)baseTests.o $(BIN)metroCudaUtil.o  $(BIN)metroUtil.o $(PARA)$(TST)parallelTest.cu -o $(BIN)$(TSTEXE)

baseTest: $(PARA)$(TST)baseTests.h $(PARA)$(TST)baseTests.cpp
	$(NV) $(FLAGS) -c $(PARA)$(TST)baseTests.cpp -o $(BIN)baseTests.o

cudaUtil: metroUtil baseTest $(PARA)$(SRC)metroCudaUtil.cuh $(PARA)$(SRC)metroCudaUtil.cu
	$(NV) $(FLAGS) -c $(PARA)$(SRC)metroCudaUtil.cu -o $(BIN)metroCudaUtil.o

metroUtil: OPLSScan zMatrix $(UTIL)$(SRC)metroUtil.h $(UTIL)$(SRC)metroUtil.cpp
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)metroUtil.cpp -o $(BIN)metroUtil.o

utilTests: stateScan stateTest configTest metroUtil zMatrix OPLSScan geoTest
	$(NV) $(BIN)geometricUtil.o $(BIN)configurationTest.o $(BIN)Config_Scan.o $(BIN)stateTest.o $(BIN)metroUtil.o $(BIN)State_Scan.o $(BIN)Zmatrix_Scan.o $(BIN)Opls_Scan.o $(BIN)geometricTest.o Utilities/test/utilityTests.cpp -o $(BIN)$(UTILTESTEXE)

zMatrix: geoUtil $(UTIL)$(SRC)Zmatrix_Scan.cpp $(UTIL)$(SRC)Zmatrix_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)Zmatrix_Scan.cpp -o $(BIN)Zmatrix_Scan.o

OPLSScan: $(UTIL)$(SRC)Opls_Scan.cpp $(UTIL)$(SRC)Opls_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)Opls_Scan.cpp -o $(BIN)Opls_Scan.o

geoUtil: $(UTIL)$(SRC)geometricUtil.cpp $(UTIL)$(SRC)geometricUtil.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)geometricUtil.cpp -o $(BIN)geometricUtil.o

stateTest: metroUtil $(UTIL)$(TST)stateTest.cpp $(UTIL)$(TST)stateTest.h 
	$(NV) $(FLAGS) -c $(UTIL)$(TST)stateTest.cpp -o $(BIN)stateTest.o

configScan: $(UTIL)$(SRC)Config_Scan.cpp $(UTIL)$(SRC)Config_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)Config_Scan.cpp -o $(BIN)Config_Scan.o

stateScan: $(UTIL)$(SRC)State_Scan.cpp $(UTIL)$(SRC)State_Scan.h
	$(NV) $(FLAGS) -c $(UTIL)$(SRC)State_Scan.cpp -o $(BIN)State_Scan.o

configTest: configScan $(UTIL)$(TST)configurationTest.h $(UTIL)$(TST)configurationTest.cpp
	$(NV) $(FLAGS) -c $(UTIL)$(TST)configurationTest.cpp -o $(BIN)configurationTest.o

geoTest: geoUtil $(UTIL)$(TST)geometricTest.cpp $(UTIL)$(TST)geometricTest.h
	$(NV) $(FLAGS) -c $(UTIL)$(TST)geometricTest.cpp -o $(BIN)geometricTest.o

dir:
	mkdir -p $(BIN)

clean:
	rm -f $(BIN)*.o
