#ifndef MECORS_H
#define MECORS_H

#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <zlib.h>

#include "bseq.h"
#include "khash.h"
typedef struct next_base {
	volatile uint16_t a, c, g, t;
} next_base_t;
KHASH_MAP_INIT_INT64(SAG, next_base_t)
khash_t(SAG) *h;

typedef struct {
	char *one, *two;
	int k, min_cov;
	int n_threads, batch_size;
	uintmax_t n_init, n_fill, n_corr;
} mecors_t;

typedef struct {
	int batch_size, n_processed, n_threads;
	bseq_file_t *fp;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
	int n_seq;
	bseq1_t *seq;
} step_t;

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
double cputime(void);
double realtime(void);

extern int mecors_verbose;
extern double mecors_real_time;
extern mecors_t opt;

extern unsigned char seq_fwd_table[128];

int main_init();
int main_fill();
int main_corr();

#endif /* MECORS_H */
