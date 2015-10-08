CC=gcc
CFLAGS=-Wall -pedantic -std=gnu99 -O3

corsage:corsage.c kseq.h khash.h
	$(CC) $(CFLAGS) -pthread corsage.c -o $@ -lz

clean:
	rm -fr gmon.out *.o ext/*.o a.out corsage *~ *.a *.dSYM session*
