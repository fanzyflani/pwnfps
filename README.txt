(pwnfps engine - game needs a title right now)
raytraced portal engine for maximum euclidfuckery
2014, @fanzyflani

to build

    gcc -g -fopenmp -ffast-math -funroll-loops -O3 -o pwnfps main.c `sdl-config --cflags --libs` `pkg-config lua-5.1 --cflags --libs` -lm -Wall -Wextra -Wno-unused-parameter

if you're using a debian-derived distro you may want to try lua-51 rather than lua-5.1 or whatever the hell it is

to build on windows, find some magical incantations as it's rather annoying

although it does appear to work! here's what i do:

    mingw32-gcc -mfpmath=both -mstackrealign -fopenmp -O3 -funroll-loops -ffast-math -msse2 -g -o pwnfps.exe main.c -Iwinlibs/SDL -Iwinlibs -Lwinlibs -lmingw32 -lSDLmain -lSDL -lm -llua -Wall -Wextra -Wno-unused-parameter

of course this requires a winlibs/ directory with SDL includes in winlibs/SDL/ and SDL libs + lua stuff in winlibs/

a quick test suggests -mfpmath=sse is best for 64-bit and -mfpmath=both is best for 32-bit

openmp highly recommended if you have more than one core
(unless you like melting holes in your laptop)

there's a chance if you're using a 32-bit x86 system that it will crash due to alignment issues

if you're using a non-x86 system, RIP. only x86 is supported, sorry - i use SSE2 intrinsics.

speaking of which, you need a CPU that supports SSE2 at the least - all 64-bit x86 CPUs support it.

want to have fun? play around with level.txt.

