# example for foo
.SUFFIXES: .o .cc
CXX = g++
FLAGS = -g -ansi
LIB = -lm -L/usr/local/lib -lfftw3
INCLUDE = -I/usr/local/include
OBJT = main.o 
TARG = fftw_test
${TARG}: ${OBJT}
	${CXX} ${FLAGS} -o ${TARG} ${OBJT} ${LIB}
.cc.o:
	${CXX} ${FLAGS} -c $< ${INCLUDE}
clean:
	rm ${OBJT} ${TARG} *~

