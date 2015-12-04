#include "mecors.h"

#define VERSION "0.4.0"

mecors_t opt = {.one = NULL, .two = NULL,    \
    .k = 31, .min_cov = 2,                    \
    .n_threads = 16, .batch_size = 100000000, \
    .n_init = 0, .n_fill = 0, .n_corr = 0};

int mecors_verbose = 0;
double mecors_real_time;

////////////////////////////////////////////////////////////////////////////////
// main                                                                       //
////////////////////////////////////////////////////////////////////////////////

static int usage() {
    fprintf(stderr, "MeCorS version %s by Andreas Bremges (andreas@cebitec.uni-bielefeld.de)\n\n", VERSION);

    fprintf(stderr, "Usage: mecors [options] -s <SC.fastq> -m <MG.fastq>\n\n");

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
                if (opt.n_threads < 1) opt.n_threads = 1;
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
                mecors_verbose = 1;
                break;
                default:
                return usage();
            }
        }
        if(!opt.one || !opt.two) return usage();

        mecors_real_time = realtime();
        if (mecors_verbose) fprintf(stderr, "[%.1f] this is MeCorS, version %s\n", realtime() - mecors_real_time, VERSION);
        h = kh_init(SAG);
        if (main_init(opt)) {
            kh_destroy(SAG, h);
            return 1;
        }
        if (main_fill(opt)) {
            kh_destroy(SAG, h);
            return 1;
        }
        if (main_corr(opt)) {
            kh_destroy(SAG, h);
            return 1;
        }
        kh_destroy(SAG, h);
        if (mecors_verbose) fprintf(stderr, "[%.1f] thank you for using MeCorS!\n", realtime() - mecors_real_time);
        return 0;
    }
