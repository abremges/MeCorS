#include "corsage.h"

#define VERSION "0.3.0-alpha"

// A = 00, C = 01, G = 10, T = 11, i.e. XOR with 3 -> compl.
unsigned char seq_fwd_table[128] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};
unsigned char seq_rev_table[4] = {
    'A', 'C', 'G', 'T'
};

corsage_t opt = {.one = NULL, .two = NULL,    \
    .k = 31, .min_cov = 2,                    \
    .n_threads = 16, .batch_size = 100000000, \
    .n_total1 = 0, .n_total2 = 0, .n_meta = 0};

int verbose = 0;

////////////////////////////////////////////////////////////////////////////////
// main                                                                       //
////////////////////////////////////////////////////////////////////////////////

static int usage() {
    fprintf(stderr, "corsage version %s by Andreas Bremges (andreas@cebitec.uni-bielefeld.de)\n\n", VERSION);

    fprintf(stderr, "Usage: corsage [options] -s <SC.fastq> -m <MG.fastq>\n\n");

    fprintf(stderr, "       -s <SC.fastq>  single cell sequencing reads\n");
    fprintf(stderr, "       -m <MG.fastq>  metagenomic sequencing reads\n");
    fprintf(stderr, "                      (fastq or fasta, can be gzip'ed)\n\n");

    fprintf(stderr, "       -k INT         k-mer size for error correction [%i]\n", opt.k);
    fprintf(stderr, "       -c INT         min. coverage in the metagenome [%i]\n\n", opt.min_cov);

    fprintf(stderr, "       -t INT         number of threads [%i]\n", opt.n_threads);
    fprintf(stderr, "       -B NUM         process ~NUM bp in each batch [100M]\n");
    fprintf(stderr, "       -v             be verbose\n\n");
    return 42;
}

int main(int argc, char *argv[]) {
    int c;
    while((c = getopt(argc, argv, "s:m:k:c:t:B:v")) != -1) {
        switch (c) {
            case 's':
                opt.one = optarg;
                break;
            case 'm':
                opt.two = optarg;
                break;
            case 'k':
                opt.k = atoi(optarg);
                if (opt.k < 15) opt.k = 15;
                if (opt.k > 31) opt.k = 31;
                break;
            case 'c':
                opt.min_cov = atoi(optarg);
                if (opt.min_cov < 1) opt.min_cov = 1;
                break;
            case 't':
                opt.n_threads = atoi(optarg);
                if (opt.n_threads < 15) opt.n_threads = 15;
                if (opt.n_threads > 31) opt.n_threads = 31;
                break;
            case 'B': {
                char *p;
                double x = strtod(optarg, &p);
                if (*p == 'G' || *p == 'g') x *= 1e9;
                else if (*p == 'M' || *p == 'm') x *= 1e6;
                else if (*p == 'K' || *p == 'k') x *= 1e3;
                opt.batch_size = (uint64_t)(x + .499);
                if (opt.batch_size < 1) opt.batch_size = 1;
                break; }
            case 'v':
                verbose = 1;
                break;
            default:
                return usage();
        }
    }
    if(!opt.one) return usage(); // || !opt.two) return usage(opt);
    h = kh_init(SAG);
    main_init(opt);

    kh_destroy(SAG, h);
    return 0;
}

// inline void stk_printseq(const bseq1_t *s) {
//     pthread_mutex_lock(&mutex_stdio);
//     fputc(s->l_qual? '@' : '>', stdout);
//     fputs(s->name, stdout);
//     //TODO Take care of comments
//     fputc('\n', stdout);
//     fputs(s->seq, stdout);
//     fputc('\n', stdout);
//     if (s->l_qual) {
//         fputs("+\n", stdout);
//         fputs(s->qual, stdout);
//         fputc('\n', stdout);
//     }
//     pthread_mutex_unlock(&mutex_stdio);
// }
