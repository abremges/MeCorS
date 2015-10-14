#ifndef BSEQ_H
#define BSEQ_H

#include <stdint.h>

struct bseq_file_s;
typedef struct bseq_file_s bseq_file_t;

typedef struct {
	int l_seq, l_qual, l_comment, rid;
	char *name, *seq, *qual, *comment;
} bseq1_t;

bseq_file_t *bseq_open(const char *fn);
void bseq_close(bseq_file_t *fp);
bseq1_t *bseq_read(bseq_file_t *fp, int chunk_size, int *n_, const int with_comment, const int with_qual);
int bseq_eof(bseq_file_t *fp);

#endif /* BSEQ_H */
