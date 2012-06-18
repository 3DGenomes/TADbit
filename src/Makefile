OBJECTS = tadbit.o tadbit_R.o tadbit.so _pytadbit.so pytadbit.pyc

all: tadbit.so _pytadbit.so

clean:
	rm -f $(OBJECTS)
	rm -rf build/

_pytadbit.so: pytadbit_wrap.c
	python setup.py build_ext --inplace

tadbit.so: tadbit.o tadbit_R.c
	R CMD SHLIB tadbit_R.c tadbit.o -o tadbit.so

tadbit.o: tadbit.c
	gcc -std=gnu99 -fPIC -g -O3 -c tadbit.c -lpthread -Wall

pytadbit_wrap.c: pytadbit.i tadbit.c
	swig -python pytadbit.i
