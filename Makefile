CC=gcc
CFLAGS=-c

FF  = gfortran
AR = ar rcs
CPLUS=g++

DIR=src
Executable=Executable

#ALG = GRAPH_PATH
ALG = ALGEBRA_PATH

#BLAS
BLAS = ./OpenBLAS/lib/libopenblas.a


## Machine Specific optimizations
OPT =   -pthread -D$(ALG) -fwrapv -O2  -Wall -Wstrict-prototypes -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2  -fstack-protector-strong -Wformat -Werror=format-security -fPIC -DUSE_CBLAS
#OPT = -g -fPIC -DUSE_CBLAS
INC =   -I ./$(DIR)  #-I/usr/include/python2.7/ -I/usr/local/lib/python2.7/dist-packages/numpy/core/include/numpy/ -I/usr/include/python2.7

.c.o:
	$(CC) -c  $(OPT) $(INC) $< -o $@



OBJ= $(DIR)/paoloexample.o
OBJ3= $(DIR)/paoloexample_b.o
OBJ2 = $(DIR)/SparseBLAS.o $(DIR)/sorting.o  $(DIR)/parsparsecoo.o $(DIR)/block.o $(DIR)/parsorting.o
vecGen.o: $(DIR)/vecGen.c $(DIR)/vecGen.h
	$(CC) $(CFLAGS) $< -o $@

matGen.o: $(DIR)/matGen.c $(DIR)/matGen.h
	$(CC) $(CFLAGS) $< -o $@

main: $(OBJ) $(DIR)/main.c
	$(CC) $^ -o $(Executable)/$@

par: $(OBJ) $(DIR)/mainPar.c
	$(CC) $^ -o $(Executable)/$@ -g -lpthread

#all: $(OBJ) $(DIR)/main.c
#	$(CC) $^ -o $(Executable)/main


sparsecoo: $(OBJ)
	make lib
	$(CC) $(OPT) $(DIR)/paoloexample.o  -o $(Executable)/stest  -L ./lib -lsparsefastmm  -lpthread  -L/usr/lib/x86_64-linux-gnu -lopenblas

sparsemb: $(OBJ3)
	make lib
	$(CC) $(OPT) $(DIR)/paoloexample_b.o  -o $(Executable)/stest2  -L ./lib -lsparsefastmm  -lpthread  -lopenblas


lib: $(OBJ2) 
	$(AR) ./lib/libsparsefastmm.a $(OBJ2)

clean:

	rm  $(DIR)/*.o lib/* 
