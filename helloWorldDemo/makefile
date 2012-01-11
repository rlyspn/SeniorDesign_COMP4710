NV=nvcc
EXE=dotProd
CFLAGS=-o

all: dot

dot: dotProduct.cu
	$(NV) dotProduct.cu $(CFLAGS) $(EXE)

clean:
	rm -f $(EXE)