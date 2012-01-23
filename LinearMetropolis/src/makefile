CC=g++
EXE=metro

all: metro

metro: util
	$(CC) kryptonMetropolis.cpp metroUtil.o -o $(EXE)

util: metroUtil.h metroUtil.cpp
	$(CC) -c metroUtil.cpp

clean:
	rm -f *.o
