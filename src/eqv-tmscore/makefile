OBJECTST = matrices.o eqv-drmsd-concat-xyz.o
HEADERST = matrices.h

#CC	= g++ # g++-4.1
#CC	= g++-3.3
#CC	= g++-4.2
#CC	= g++-4.3
CC	= g++-4.4
#CC	= g++-4.5
CFLAGS  = -O3


eqv-drmsd-concat-xyz :	$(OBJECTST) 
	$(CC) $(CFLAGS) -o $@ $(OBJECTST) $(LLIBS)

eqv-drmsd-concat-xyz.o : eqv-drmsd-concat-xyz.cc $(HEADERST)
	$(CC) -c $(CFLAGS) $*.cc

matrices.o : matrices.cc $(HEADERST)
	$(CC) -c $(CFLAGS) $*.cc

clean:
	rm *.o #Train
