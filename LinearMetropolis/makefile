CC=g++
EXE=metro
DEMOEXE=LinearDemo
TSTEXE=linearTest
SRC=src/
UTIL=../Utilities/src/
TST=test/

all:  metro metroUtil

#runtest: test
#	./$(TSTEXE)	

#test: util baseTest $(TST)linearTest.cpp
#	$(CC) $(TST)linearTest.cpp metroUtil.o baseTests.o -o $(TSTEXE)

demo: util
	$(CC) $(SRC)inClassLinearDemo.cpp metroUtil.o -o $(DEMOEXE)

#metro: util
#	$(CC) $(SRC)kryptonMetropolis.cpp metroUtil.o -o $(EXE)

metroUtil: $(SRC)metroLinearUtil.cpp $(SRC)metroLinearUtil.h util baseTest
	$(CC) $(SRC)metroLinearUtil.cpp -c

util: $(UTIL)metroUtil.h $(UTIL)metroUtil.cpp
	$(CC) -c $(UTIL)metroUtil.cpp

baseTest: $(TST)baseTests.h $(TST)baseTests.cpp
	$(CC) -c $(TST)baseTests.cpp

clean:
	rm -f -R *.o
	rm -f -R $(EXE)
	rm -f -R $(TSTEXE)
	rm -f -R $(DEMOEXE)
