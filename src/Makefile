OBJECTS = tadbit.o tadbit_R.o tadbit.so

all: tadbit.so

clean:
	rm -f $(OBJECTS)
	rm -rf build/

tadbit.so: tadbit.o tadbit_R.c
	R CMD SHLIB tadbit_R.c tadbit.o -o tadbit.so

tadbit.o: tadbit.c
	gcc -std=gnu99 -fPIC -g -O3 -c tadbit.c \
		-lpthread -Wall
