CC=g++
CFLAGS=-Isrc
SRC=src/arithmetic.cpp src/io.cpp src/wire.cpp src/main.cpp
BIN=build/ivsolver

.PHONY: clean build test

build:
	@mkdir -p build
	$(CC) $(SRC) -o $(BIN) $(CFLAGS) -std=c++2a -O3

test:
	./build/ivsolver ./res/key_points.in ./res/results

clean:
	@rm -f $(BIN)