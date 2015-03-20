CC ?= gcc
CXX ?= g++

pwnfps_CFLAGS = -g -O3 -Wall -Wextra -Wno-unused-parameter -Wno-unused-but-set-variable \
				`sdl-config --cflags` `pkg-config lua5.1 --cflags`
pwnfps_LIBS   = -lm `sdl-config --libs` `pkg-config lua5.1 --libs`


all: pwnfps

pwnfps:
	mkdir build
	$(CC) $(pwnfps_CFLAGS) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o build/$@ main.c $(pwnfps_LIBS)
	cp -f game.lua level.txt build

clean:
	rm -rf build
