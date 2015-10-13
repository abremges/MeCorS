#ifndef CORSAGE_H
#define CORSAGE_H

#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <zlib.h>

#include "bseq.h"
#include "khash.h"
typedef struct next_base {
    volatile uint8_t a, c, g, t;
} next_base_t;
KHASH_MAP_INIT_INT64(SAG, next_base_t)
khash_t(SAG) *h;

typedef struct {
    char *one, *two;
	int k, min_cov;
    int n_threads, batch_size;
	uintmax_t n_total1, n_total2, n_meta;
} corsage_t;

extern int verbose;
extern corsage_t opt;

extern unsigned char seq_fwd_table[128];

int main_init();
int main_fill();

#endif /* CORSAGE_H */
