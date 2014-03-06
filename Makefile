CC=gcc
CFLAGS=-Wall -pedantic -std=gnu99 -O3

hector:hector.c kseq.h khash.h
	$(CC) $(CFLAGS) hector.c -o $@ -lz

clean:
	rm -fr gmon.out *.o ext/*.o a.out hector *~ *.a *.dSYM session*
