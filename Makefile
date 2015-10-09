PROG=corsage

CC=gcc
CFLAGS=-Wall -pedantic -std=gnu99 -O3
INCLUDES=-I.
OBJS=kthread.o bseq.o init.o
LIBS=-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

$(PROG):$(OBJS) main.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

# DO NOT DELETE

bseq.o: bseq.h kseq.h
init.o: bseq.h khash.h corsage.h
fill.o: bseq.h khash.h corsage.h
corr.o: bseq.h khash.h corsage.h
main.o: corsage.h
