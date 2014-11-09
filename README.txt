(pwnfps engine - game needs a title right now)
raytraced portal engine for maximum euclidfuckery
2014, @fanzyflani

to build

    cc -msse2 -fopenmp -g -O2 -o pwnfps main.c `sdl-config --cflags --libs` -lm

to build on windows, find some magical incantations as it's rather annoying

although it does appear to work! here's what i do:

    mingw32-gcc -msse2 -g -O2 -o pwnfps.exe main.c -Iwinlibs/SDL -Lwinlibs -lmingw32 -lSDLmain -lSDL -lm -Wall -Wextra -Wno-unused-parameter

of course this requires a winlibs/ directory with SDL includes in winlibs/SDL/ and SDL libs in winlibs/

openmp highly recommended if you have more than one core
(unless you like melting holes in your laptop)

there's a chance if you're using a 32-bit x86 system that it will crash due to alignment issues

if you're using a non-x86 system, RIP. only x86 is supported, sorry - i use SSE2 intrinsics.

speaking of which, you need a CPU that supports SSE2 at the least - all 64-bit x86 CPUs support it.

want to have fun? play around with level.txt.

