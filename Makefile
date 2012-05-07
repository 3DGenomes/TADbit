OBJECTS= tadbit.o main.o tadbit_R.o tadbit_R.so tadbit

all: tadbit tadbit_R.so

clean:
	- rm -f $(OBJECTS)

tadbit_R.so: tadbit.o
	R CMD SHLIB tadbit_R.c tadbit.o

tadbit: tadbit.o main.o
	- rm -f tadbit
	cc -std=gnu99 -g -lpthread -lm \
	tadbit.o main.o -o tadbit -Wall

tadbit.o: tadbit.c
	cc -std=gnu99 -fPIC -g -c tadbit.c -o tadbit.o -O3 -lpthread -Wall

main.o: main.c
	cc -std=gnu99 -g -c main.c -o main.o -O3
