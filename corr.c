#include "corsage.h"

inline int baseCorrect(const int base, const next_base_t next) {
    // TODO also introduce a qvalue filter? always try to correct/confirm low-q bases
	int min_cov = opt.min_cov;

    // Look for support in the metagenome, if min_cov is reached don't try to correct
    if (base == 0 && next.a >= min_cov) return base;
    if (base == 1 && next.c >= min_cov) return base;
    if (base == 2 && next.g >= min_cov) return base;
    if (base == 3 && next.t >= min_cov) return base;

    // The metagenome doesn't support the next base sufficiently, we try to correct
    int majority = (next.a + next.c + next.g + next.t)/2;
    if (next.a > majority && next.a >= min_cov) return 0;
    if (next.c > majority && next.c >= min_cov) return 1;
    if (next.g > majority && next.g >= min_cov) return 2;
    if (next.t > majority && next.t >= min_cov) return 3;

    // We could not correct the next base, indicate failure by returning -1
    return -1;
}

static inline int process_read(const bseq1_t read, const int k) {
    int failed = 0;
    int corrected = 0;
    uint64_t forward = 0;
    for (int i = 0, index = 1; i < read.l_seq-1; ++i, ++index) { // stops 1 char before end to look ahead!
        int c = seq_fwd_table[(int) read.seq[i]];
        if (c < 4) {
            forward = (((forward << 2) | c) & ((1ULL<<k*2)-1));
        } else {
            index = 0;
            ++failed;
        }
        if (index >= k) {
            khiter_t iter = kh_get(SAG, h, forward);
            if (iter != kh_end(h)) {
                next_base_t next = kh_value(h, iter); // supported 32nd bases
                int base = seq_fwd_table[(int)read.seq[i+1]]; // current base
                int sbase = baseCorrect(base, next);
                if (base != sbase) { // Correction needed
                    switch (sbase) {
                        case 0:
                            read.seq[i+1] = 'A';
                            ++corrected;
                            break;
                        case 1:
                            read.seq[i+1] = 'C';
                            ++corrected;
                            break;
                        case 2:
                            read.seq[i+1] = 'G';
                            ++corrected;
                            break;
                        case 3:
                            read.seq[i+1] = 'T';
                            ++corrected;
                            break;
                        default:
                            index = -1; // Skip uncorrectable base
                            ++failed;
                            break;
                    }
                }
            }
        }
    }
    return failed;
}

//////////
// Taken from: https://github.com/lh3/seqtk/blob/master/seqtk.c
//////////
char comp_tab[] = {
    0, 1,	2,	3,	4, 5,	6,	7,	8, 9, 10,	11,	12, 13, 14,	15,
    16, 17, 18,	19,	20, 21, 22,	23,	24, 25, 26,	27,	28, 29, 30,	31,
    32, 33, 34,	35,	36, 37, 38,	39,	40, 41, 42,	43,	44, 45, 46,	47,
    48, 49, 50,	51,	52, 53, 54,	55,	56, 57, 58,	59,	60, 61, 62,	63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	92, 93, 94,	95,
    64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};
void seq_revcomp(bseq1_t seq) {
    int c0, c1;
    for (int i = 0; i < seq.l_seq>>1; ++i) { // reverse complement sequence
        c0 = comp_tab[(int)seq.seq[i]];
        c1 = comp_tab[(int)seq.seq[seq.l_seq - 1 - i]];
        seq.seq[i] = c1;
        seq.seq[seq.l_seq - 1 - i] = c0;
    }
    if (seq.l_seq & 1) // complement the remaining base
    seq.seq[seq.l_seq>>1] = comp_tab[(int)seq.seq[seq.l_seq>>1]];
    if (seq.l_qual) {
        for (int i = 0; i < seq.l_seq>>1; ++i) // reverse quality
        c0 = seq.qual[i], seq.qual[i] = seq.qual[seq.l_qual - 1 - i], seq.qual[seq.l_qual - 1 - i] = c0;
    }
}

static void worker_for(void *_data, long i, int tid) { // kt_for() callback
    step_t *s = (step_t*)_data;
	if (process_read(s->seq[i], opt.k)) {
		seq_revcomp(s->seq[i]);
		process_read(s->seq[i], opt.k);
		seq_revcomp(s->seq[i]);
	}
}

inline void stk_printseq(const bseq1_t s) {
    fputc(s.l_qual? '@' : '>', stdout);
    fputs(s.name, stdout);
    //TODO Take care of comments
    fputc('\n', stdout);
    fputs(s.seq, stdout);
    fputc('\n', stdout);
    if (s.l_qual) {
        fputs("+\n", stdout);
        fputs(s.qual, stdout);
        fputc('\n', stdout);
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
    } else if (step == 1) { // step 1: correct reads
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_seq);
		return in;
    } else if (step == 2) { // step 2: output and clean up
        step_t *s = (step_t*)in;
		for (int i = 0; i < s->n_seq; ++i) {
			stk_printseq(s->seq[i]);
            free(s->seq[i].name);
			free(s->seq[i].seq);
            free(s->seq[i].qual);
		}
		free(s->seq);
		free(s);
	}
    return 0;
}

int main_corr(const corsage_t opt) {
    fprintf(stderr, "Corr.\n");
    pipeline_t pl;
    memset(&pl, 0, sizeof(pipeline_t));
    pl.fp = bseq_open(opt.one);
    if (pl.fp == 0) return -1; //TODO
    pl.n_threads = opt.n_threads, pl.batch_size = opt.batch_size;
    kt_pipeline(opt.n_threads, worker_pipeline, &pl, 3);
    fprintf(stderr, "Done.\n");
    return 0;
}
