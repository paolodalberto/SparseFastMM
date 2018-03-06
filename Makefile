CC=gcc
CFLAGS=-c

FF  = gfortran
AR = ar rcs

DIR=src
Executable=Executable

## Machine Specific optimizations 
OPT = -O2 -Wall # -pthread -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC

INC = -I ./$(DIR)  #-I/usr/include/python2.7/ -I/usr/local/lib/python2.7/dist-packages/numpy/core/include/numpy/ -I/usr/include/python2.7

.c.o:
	$(CC) -c  $(OPT) $(INC) $< -o $@



OBJ= vecGen.o matGen.o

vecGen.o: $(DIR)/vecGen.c $(DIR)/vecGen.h
	$(CC) $(CFLAGS) $< -o $@

matGen.o: $(DIR)/matGen.c $(DIR)/matGen.h
	$(CC) $(CFLAGS) $< -o $@

all: $(OBJ) $(DIR)/main.c
	$(CC) $^ -o $(Executable)/main


sparsecoo: $(DIR)/SparseBLAS.o $(DIR)/sorting.o  $(DIR)/paoloexample.o $(DIR)/parsparsecoo.o
	$(CC) $(OPT) $(DIR)/paoloexample.o $(DIR)/SparseBLAS.o $(DIR)/sorting.o $(DIR)/parsparsecoo.o -o $(Executable)/stest -lpthread

clean:
	rm *.o $(DIR)/*.o 
