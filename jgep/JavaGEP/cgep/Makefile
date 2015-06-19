# CFLAGS=-pg -Wall -g
CFLAGS=-O2

all:
	cc $(CFLAGS) -c gep.c
	cc $(CFLAGS) -c driver.c
	cc $(CFLAGS) -c regress.c
	cc $(CFLAGS) gep.o regress.o driver.o -o driver -lm

clean:
	rm -f driver *.o *~
