CC=g++
CFLAGS=-Wall -Wpedantic -Wextras
LFLAFS=

.PHONY: all clean

all: main

main: build/main.o
	$(CC) $(LFLAGS) $^ -o $@

build/%.o: src/%.cpp src/%.hpp
	$(CC) $(CFLAGS) $< -o $@

build:
	mkdir build

clean:
	rm -rf build

