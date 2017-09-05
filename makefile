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
dc-decdel.o: dc_decdel.f
	$(FC) -c dc_decdel.f
contr5.o: contr5.f
	$(FC) -c contr5.f
rk4: rk4.o radar5.o decsol.o dc_decdel.o contr5.o
	$(FC) -o rk4 rk4.o radar5.o decsol.o dc_decdel.o contr5.o
clean:
	rm -rf *.o rk4
