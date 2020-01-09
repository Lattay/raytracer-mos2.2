CPP=g++
CFLAGS=-Wall -Wextra -Wpedantic --std=c++11
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

build/%.o: src/%.cpp src/%.hpp
	$(CPP) $(CFLAGS) -c $< -o $@

tags: $(SRC)
	ctags --extras=rq $(SRC)

build:
	mkdir build

clean:
	rm -rf build

wipe: clean
	rm -f *.png
	rm -f main
	rm -f tags
