#include "corsage.h"

pthread_mutex_t mutex_khash;

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int batch_size, n_processed, n_threads;
	bseq_file_t *fp;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq;
	bseq1_t *seq;
} step_t;

static inline void process_kmer(const uint64_t kmer) {
	khiter_t kh_n = kh_end(h); // first get number of buckets in khash
    khiter_t iter = kh_get(SAG, h, kmer); // then query k-mer position
    if (iter >= kh_n) { // most likely not in there (or khash resized)
		pthread_mutex_lock(&mutex_khash);
		int ret;
		iter = kh_put(SAG, h, kmer, &ret);
		kh_value(h, iter) = (next_base_t) {0,0,0,0};
		pthread_mutex_unlock(&mutex_khash);
    }
}

static inline void process_read(const bseq1_t *read, const int k) {
    uint64_t forward = 0;
    uint64_t reverse = 0;
    for (int i = 0, index = 1; i < read->l_seq; ++i, ++index) {
        int c = seq_fwd_table[(int) read->seq[i]];
        if (c < 4) {
            forward = (((forward << 2) | c) & ((1ULL<<k*2)-1));
            reverse = ((reverse >> 2) | (((uint64_t) (c^3)) << ((k*2)-2)));
        } else {
            index = 0;
        }
        if (index >= k) {
            process_kmer(forward);
            process_kmer(reverse);
        }
    }
}

static void worker_for(void *_data, long i, int tid) { // kt_for() callback
    step_t *s = (step_t*)_data;
    for (i = 0; i < s->n_seq; ++i) {
        process_read(s->seq, 31);
    }
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
    } else if (step == 1) { // step 1: fill hash with kmers
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

int main_init(const corsage_t opt) {
    pthread_mutex_init(&mutex_khash, NULL);

    fprintf(stderr, "Init.\n");
    pipeline_t pl;
    memset(&pl, 0, sizeof(pipeline_t));
    pl.fp = bseq_open(opt.one);
    if (pl.fp == 0) return -1;
    pl.n_threads = opt.n_threads, pl.batch_size = opt.batch_size;
    kt_pipeline(opt.n_threads, worker_pipeline, &pl, 3);
    fprintf(stderr, "Done.\n");

    return 0;
}
