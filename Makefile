CXX ?= g++
CXXFLAGS ?= -std=c++2a -O3
CPPFLAGS += -Isrc
LDFLAGS += -pthread
SRC = src/arithmetic.cpp src/io.cpp src/wire.cpp src/field.cpp src/circuit.cpp src/solver.cpp src/main.cpp
HDR = $(wildcard src/*.h)
BIN = build/ivsolver

# WebAssembly build (requires the Emscripten SDK: https://emscripten.org).
EMXX ?= em++
WASM_SRC = src/arithmetic.cpp src/io.cpp src/wire.cpp src/field.cpp src/circuit.cpp src/solver.cpp src/wasm_api.cpp
WASM_OUT = build/wasm/ivsolver.mjs
WASM_FLAGS = -std=c++2a -O3 -Isrc -lembind -fexceptions \
	-sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createIVSolver \
	-sENVIRONMENT=web,worker,node -sALLOW_MEMORY_GROWTH=1 -sSINGLE_FILE=1

.PHONY: build clean test wasm

build: $(BIN)

$(BIN): $(SRC) $(HDR)
	@mkdir -p build
	$(CXX) $(SRC) -o $(BIN) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS)

wasm: $(WASM_OUT)

$(WASM_OUT): $(WASM_SRC) $(HDR)
	@mkdir -p build/wasm
	$(EMXX) $(WASM_SRC) -o $(WASM_OUT) $(WASM_FLAGS)

test: $(BIN)
	python3 tests/run_tests.py $(BIN)

clean:
	@rm -f $(BIN)
	@rm -rf build/wasm
