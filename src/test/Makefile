vpath %.c ..
vpath %.h ..

P= testset
OBJECTS= tadbit.o
CFLAGS= -I.. `pkg-config --cflags glib-2.0` -g -pg -Wall -std=gnu99 \
	          -O3 -fstrict-aliasing -fprofile-arcs -ftest-coverage
LDLIBS= `pkg-config --libs glib-2.0` -lpthread -lm
CC= gcc
$(P): $(OBJECTS)

clean:
	rm -f testset *.o *.gcda *.gcno *.gcov gmon.out analysis.txt \
		callgrind.out.* cache.txt

test:
	gtester --verbose --keep-going testset

debug:
	gdb --command=debug.gdb --args testset
