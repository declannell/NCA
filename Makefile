CC = g++

CFLAGS = -g -Wall #-D_CC_OVERLAP

LDFLAGS = -lm

GFOBJS = parameters.o 

EXECS = a

##i think the problem is that i dont have explicitly the dmft depends on parameters.

all: $(EXECS)

a: main.o  $(GFOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tags:
	etags *.c *.h

.PHONY: clean tags tests

clean:
	$(RM) *.o $(EXECS) $(TESTS) TAGS tags