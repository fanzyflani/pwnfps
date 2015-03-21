CC ?= gcc

pwnfps_CFLAGS = \
	-g -O3 -fopenmp -ffast-math -funroll-loops -Wall -Wextra \
	-Wno-unused-parameter -Wno-unused-but-set-variable -Wno-unused-variable \
	`sdl-config --cflags` `pkg-config lua5.1 --cflags`

pwnfps_LIBS = -lm `sdl-config --libs` `pkg-config lua5.1 --libs`


all: pwnfps

pwnfps:
	mkdir -p build
	$(CC) $(pwnfps_CFLAGS) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o build/$@ main.c $(pwnfps_LIBS)
	cp -f game.lua level.txt build

clean:
	rm -rf build
