CPP=g++
DEBUG=0
OMP=1
CFLAGS=-Wall -Wextra -Wpedantic --std=c++11
ifeq ($(DEBUG),1)
	CFLAGS+=-g -Og
else
	CFLAGS+=-O3 
endif
ifeq ($(OMP),1)
	CFLAGS+=-fopenmp
endif

LFLAGS=-lgomp
SRC=$(wildcard src/*.cpp)
OBJ=$(patsubst src/%.cpp,build/%.o,$(SRC))

.PHONY: all clean run

all: main

run: main
	./main
	/home/theo/.local/bin/xdg-open ./image.png &

debug: main
	gdb ./main

main: $(OBJ)
	$(CPP) $(CFLAGS) $(LFLAGS) $^ -o $@

build/%.o: src/%.cpp src/%.hpp | build/
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
