CPP=g++
CFLAGS=-Wall -Wextra --std=c++11
LFLAFS=
SRC=$(wildcard src/*.cpp)
OBJ=$(patsubst src/%.cpp,build/%.o,$(SRC))

.PHONY: all clean run

all: main

run: main
	./main
	xdg-open ./image.png

main: $(OBJ)
	$(CPP) $(LFLAGS) $^ -o $@

build/%.o: src/%.cpp src/%.hpp build
	$(CPP) $(CFLAGS) -c $< -o $@

build:
	mkdir build

clean:
	rm -rf build

