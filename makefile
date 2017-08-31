.PHONY: clean
all = rk4
CC = gcc
LIBS = -lm
FC = gfortran
FCFLAGS = -c

rk4: rk4.o
	$(CC) -o rk4 rk4.o $(LIBS)

clean:
	rm -rf *.o rk4
