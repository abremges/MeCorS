#include "corsage.h"

static inline void process_kmer(const uint64_t kmer, const int base) {
    khiter_t iter = kh_get(SAG, h, kmer);
    if (iter != kh_end(h)) {
        next_base_t *p_next = &(kh_value(h, iter));
        switch (base) {
            case 0:
                if (p_next->a < 60000) __sync_fetch_and_add(&(p_next->a), 1);
                break;
            case 1:
                if (p_next->c < 60000) __sync_fetch_and_add(&(p_next->c), 1);
                break;
            case 2:
                if (p_next->g < 60000) __sync_fetch_and_add(&(p_next->g), 1);
                break;
            case 3:
                if (p_next->t < 60000) __sync_fetch_and_add(&(p_next->t), 1);
                break;
        }
    }
}

static inline void process_read(const bseq1_t read, const int k) {
    uint64_t forward = 0;
    uint64_t reverse = 0;
    for (int i = 0, index = 1; i < read.l_seq; ++i, ++index) { // 1-based index
        int c = seq_fwd_table[(int) read.seq[i]];
        if (c < 4) {
            forward = (((forward << 2) | c) & ((1ULL<<k*2)-1));
            reverse = ((reverse >> 2) | (((uint64_t) (c^3)) << ((k*2)-2)));
        } else {
            index = 0;
        }
        if (index >= k) {
            if (i < read.l_seq-1) process_kmer(forward, seq_fwd_table[(int) read.seq[i+1]]); // OK
            if (i >= k) process_kmer(reverse, seq_fwd_table[(int) read.seq[i-k]]^3); // OK
        }
    }
}

static void worker_for(void *_data, long i, int tid) { // kt_for() callback
    step_t *s = (step_t*)_data;
    process_read(s->seq[i], opt.k);
}

static void *worker_pipeline(void *shared, int step, void *in) {
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read batch of sequences
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = bseq_read(p->fp, p->batch_size, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (int i = 0; i < s->n_seq; ++i) {
				s->seq[i].rid = p->n_processed++;
            }
			return s;
		} else {
            free(s);
        }
    } else if (step == 1) { // step 1: update hash with values
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_seq);
		return in;
    } else if (step == 2) { // step 2: output and clean up
        step_t *s = (step_t*)in;
		for (int i = 0; i < s->n_seq; ++i) {
            free(s->seq[i].name);
			free(s->seq[i].seq);
            free(s->seq[i].qual);
		}
		free(s->seq);
		free(s);
	}
    return 0;
}

int main_fill(const corsage_t opt) {
    fprintf(stderr, "Fill.\n");
    pipeline_t pl;
    memset(&pl, 0, sizeof(pipeline_t));
    pl.fp = bseq_open(opt.two);
    if (pl.fp == 0) return -1; //TODO
    pl.n_threads = opt.n_threads, pl.batch_size = opt.batch_size;
    kt_pipeline(opt.n_threads, worker_pipeline, &pl, 3);
    fprintf(stderr, "Done.\n");
    return 0;
}
