CC=g++
EXE=metro
TSTEXE=linearTest
SRC=src/
TST=test/

all: runtest metro

runtest: test
	./$(TSTEXE)	

test: util baseTest $(TST)linearTest.cpp
	$(CC) $(TST)linearTest.cpp metroUtil.o baseTests.o -o $(TSTEXE)

metro: util
	$(CC) $(SRC)kryptonMetropolis.cpp metroUtil.o -o $(EXE)

util: $(SRC)metroUtil.h $(SRC)metroUtil.cpp
	$(CC) -c $(SRC)metroUtil.cpp

baseTest: $(TST)baseTests.h $(TST)baseTests.cpp
	$(CC) -c $(TST)baseTests.cpp

clean:
	rm -f -R *.o
	rm -f -R $(EXE)
	rm -f -R $(TSTEXE)
