OBJECTS = tadbit.o tadbit_R.o tadbit.so

all: tadbit.so

clean:
	rm -f $(OBJECTS)

tadbit.so: tadbit.o tadbit_R.c tadbit.h
	R CMD SHLIB tadbit_R.c tadbit.o -o tadbit.so -lm

tadbit.o: tadbit.c tadbit.h
	gcc -std=gnu99 -fPIC -g -O3 -c tadbit.c \
		-lpthread -Wall
