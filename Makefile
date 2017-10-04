all:
	gcc -c libdcl.c
	ar rc libdcl.a libdcl.o
	ranlib libdcl.a
	gcc example.c -ldcl -lm -lcrypto -lgsl -lgslcblas -lamos -lgfortran -L. -o example
clean:
	rm -f *.o
	rm -f *.o
	rm -f example