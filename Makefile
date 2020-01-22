CPP=g++
DEBUG=0

ifeq ($(DEBUG),0)
	CFLAGS=-Wall -Wextra -Wpedantic --std=c++11 -fopenmp -O3
else
	CFLAGS=-Wall -Wextra -Wpedantic --std=c++11 -fopenmp -g
endif

LFLAGS=-lgomp
SRC=$(wildcard src/*.cpp)
OBJ=$(patsubst src/%.cpp,build/%.o,$(SRC))

.PHONY: all clean run

all: main

run: main
	./main
	xdg-open ./image.png &

debug: main
	gdb ./main

main: $(OBJ)
	$(CPP) $(LFLAGS) $^ -o $@

build/%.o: src/%.cpp src/%.hpp build/
	$(CPP) $(CFLAGS) -c $< -o $@

tags: $(SRC)
	ctags --extras=rq $(SRC)

build/:
	mkdir -p build

clean:
	rm -rf build

wipe: clean
	rm -f *.png
	rm -f main
	rm -f tags
