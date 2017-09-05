.PHONY: clean
all = rk4
CC = gcc
LIBS = -lm
FC = gfortran

rk4.o: rk4.c
	$(CC) -c rk4.c $(LIBS)
radar5.o: radar5.f
	$(FC) -c radar5.f
decsol.o: decsol.f
	$(FC) -c decsol.f
rk4: rk4.o radar5.o decsol.o
	$(FC) -o rk4 rk4.o radar5.o decsol.o
clean:
	rm -rf *.o rk4
