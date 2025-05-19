
CC=gcc
CFLAGS=-Wall

BIN_PATH=./bin
DIST_PATH=./dist
SRC_PATH=./src

all: $(BIN_PATH)/main

$(BIN_PATH)/main: main.c $(DIST_PATH)/alloc.o $(DIST_PATH)/images.o $(DIST_PATH)/matrix.o $(DIST_PATH)/geom.o 
	$(CC) $(CFLAGS) -I./include -o $@ $^ -lm

$(DIST_PATH)/images.o: $(SRC_PATH)/images.c
	$(CC) $(CFLAGS) -I./include -c -o $@ $^

$(DIST_PATH)/matrix.o: $(SRC_PATH)/matrix.c
	$(CC) $(CFLAGS) -I./include -c -o $@ $^

$(DIST_PATH)/geom.o: $(SRC_PATH)/geom.c
	$(CC) $(CFLAGS) -I./include -c -o $@ $^

$(DIST_PATH)/alloc.o: $(SRC_PATH)/alloc.c
	$(CC) $(CFLAGS) -I./include -c -o $@ $^