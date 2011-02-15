#
#
# -DLINUX if on Linux/Unix NO -DLINUX on cygwin
# -DUSE_DOUBLE_PRECISION for double precision (the most efficient option)
# -DSSE2 -O3 for pentium4 optimizations
# -DPARALLEL is you have mpi then compile with mpiCC and run with mpirun 
#
# NOTE: I hate makefiles and I do not use them. 
# If you know how to write a better makefile send it to me please!
#

CC     = g++ -I../Libraries 
# Basic compilation flags (-O2 for speed)
CFLAGS = -DLINUX -O3
# Flags for double precision and sse (P4 optiomizations)
#CFLAGS = -DLINUX -DSSE2 -O3
# Flags for parallel (mpi) double precision and sse
#CFLAGS = -DPARALLEL -DLINUX -DUSE_DOUBLE_PRECISION -DSSE2 -O2

all:
	${CC} ${CFLAGS} example01.cpp -o example01.exe
	${CC} ${CFLAGS} example02.cpp -o example02.exe
	${CC} ${CFLAGS} example03.cpp -o example03.exe
	${CC} ${CFLAGS} example04.cpp -o example04.exe
	${CC} ${CFLAGS} example05.cpp -o example05.exe
	${CC} ${CFLAGS} example06.cpp -o example06.exe
	${CC} ${CFLAGS} example07.cpp -o example07.exe
	${CC} ${CFLAGS} example08.cpp -o example08.exe
	${CC} ${CFLAGS} example09.cpp -o example09.exe
	${CC} ${CFLAGS} example10.cpp -o example10.exe
	${CC} ${CFLAGS} example11.cpp -o example11.exe
	${CC} ${CFLAGS} example12.cpp -o example12.exe
	${CC} ${CFLAGS} example13.cpp -o example13.exe
	${CC} ${CFLAGS} example14.cpp -o example14.exe
	${CC} ${CFLAGS} example15.cpp -o example15.exe
	${CC} ${CFLAGS} test_sse2.cpp -o test_sse2.exe
	${CC} ${CFLAGS} test_actions.cpp -o test_actions.exe
	${CC} ${CFLAGS} gauge.cpp        -o gauge.exe
	${CC} ${CFLAGS} meson.cpp        -o meson.exe
	${CC} ${CFLAGS} ks_meson.cpp     -o ks_meson.exe
	${CC} ${CFLAGS} make_meson.cpp   -o make_meson.exe
	${CC} ${CFLAGS} make_meson2.cpp  -o make_meson2.exe

test:
	./test_actions.exe
	./gauge.exe -hot -nc=3 -L=12x4x4x4 -output=gauge00.mdp
	./gauge.exe -input=gauge00.mdp -beta=6.0 -steps=10 -output=gauge01.mdp
	./meson.exe -input=gauge01.mdp -lkappa=0.1 -hkappa=0.1 -ps -clean
	./ks_meson.exe -input=gauge01.mdp -lmass=0.1 -hmass=0.1 -u0=1.0 -meson=5x5


