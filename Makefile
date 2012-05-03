OBJECTS= tadbit.o main.o tadbit_R.o tadbit_R.so tadbit

all: tadbit tadbit_R.so

clean:
	- rm -f $(OBJECTS)

tadbit_R.so: tadbit.o
	R CMD SHLIB tadbit_R.c tadbit.o

tadbit: libefence.a tadbit.o main.o
	- rm -f tadbit
	cc -std=gnu99 -g -Wall -lpthread -lm \
	tadbit.o libefence.a main.o -o tadbit 

tadbit.o: tadbit.c
	cc -std=gnu99 -fPIC -g -c tadbit.c -o tadbit.o

main.o: main.c
	cc -std=gnu99 -g -c main.c -o main.o

libefence.a: efence.o page.o print.o
	- rm -f libefence.a
	ar crv libefence.a efence.o page.o print.o
	- rm -f efence.o page.o print.o

efence.o: efence.c
	cc -g   -c -o efence.o efence.c

page.o: page.c
	cc -g   -c -o page.o page.c

print.o: print.c
	cc -g   -c -o print.o print.c
