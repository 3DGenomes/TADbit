vpath %.c ..
vpath %.h ..

P= testset
OBJECTS= tadbit.o
CFLAGS= -I.. `pkg-config --cflags glib-2.0` -g -pg -Wall -std=gnu99 \
	          -O0 -fstrict-aliasing -fprofile-arcs -ftest-coverage
LDLIBS= `pkg-config --libs glib-2.0` -lpthread -lm
CC= gcc
$(P): $(OBJECTS)

tadbit.o: ../tadbit.c ../tadbit.h
	gcc -std=gnu99 -fPIC -g -O3 -c ../tadbit.c -lpthread -DNDEBUG -Wall

clean:
	rm -f testset *.o *.gcda *.gcno *.gcov gmon.out analysis.txt \
		callgrind.out.* cache.txt

test: testset
	gtester --verbose --keep-going testset

fulltest: testset
	gtester --verbose --keep-going -m=thorough testset

debug:
	gdb --command=debug.gdb --args testset
