CC=gcc
CFLAGS=-c

DIR=src
OBJ= vecGen.o matGen.o

vecGen.o: $(DIR)/vecGen.c $(DIR)/vecGen.h
	$(CC) $(CFLAGS) $< -o $@

matGen.o: $(DIR)/matGen.c $(DIR)/matGen.h
	$(CC) $(CFLAGS) $< -o $@

main: $(OBJ) $(DIR)/main.c
	$(CC) $(CFLAGS) $^ -o $@

par: $(OBJ) $(DIR)/mainPar.c
	$(CC) $(CFLAGS) $^ -o $@

all: main par
