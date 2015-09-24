/* The MIT License

   Copyright (c) 2015 Andreas Bremges <andreas@cebitec.uni-bielefeld.de>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#define VERSION "0.2.2"

#include <inttypes.h> // uint64_t
#include <stdio.h>    // printf
#include <unistd.h>   // getopt
#include <zlib.h>     // gzip

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khash.h"
typedef struct next_base {
    uint8_t a, c, g, t;
} next_base_t;
KHASH_MAP_INIT_INT64(SAG, next_base_t)
khash_t(SAG) *h;

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

int k = 31, min_cov = 2;
uintmax_t n_total1 = 0, n_total2 = 0, n_meta = 0;

////////////////////////////////////////////////////////////////////////////////
// Init the data structure, i.e. insert keys into hash (single-cell k-mers)   //
////////////////////////////////////////////////////////////////////////////////

void kmerInitHash(const uint64_t kmer) {
    khiter_t iter = kh_get(SAG, h, kmer);
    if (iter == kh_end(h)) {
        int ret;
        iter = kh_put(SAG, h, kmer, &ret);
        kh_value(h, iter) = (next_base_t) {0,0,0,0};
    }
}

void readInitHash(const kseq_t *read) {
    uint64_t forward = 0;
    uint64_t reverse = 0;
    for (int i = 0, index = 1; i < read->seq.l; ++i, ++index) {
        int c = seq_fwd_table[(int) read->seq.s[i]];
        if (c < 4) {
            forward = (((forward << 2) | c) & ((1ULL<<k*2)-1));
            reverse = ((reverse >> 2) | (((uint64_t) (c^3)) << ((k*2)-2)));
        } else {
            index = 0;
        }
        if (index >= k) {
            kmerInitHash(forward);
            kmerInitHash(reverse);
        }
    }
    ++n_total1;
    if (!(n_total1 % 10000)) fprintf(stderr, "%" PRIuMAX "\n", n_total1);
}

void fileInitHash(const char *file) {
    gzFile fp = strcmp(file, "-") ? gzopen(file, "r") : gzdopen(fileno(stdin), "r");
    kseq_t *r = kseq_init(fp);
    while (kseq_read(r) >= 0) {
        readInitHash(r);
    }
    kseq_destroy(r);
    gzclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
// Fill hash with information about the next base, based on metagenome reads  //
////////////////////////////////////////////////////////////////////////////////

void kmerPlusPlus(const uint64_t kmer, const int base) {
    khiter_t iter = kh_get(SAG, h, kmer);
    if (iter != kh_end(h)) {
        next_base_t next = kh_value(h, iter);
        switch (base) {
            case 0:
                if (next.a < 255) ++next.a;
                break;
            case 1:
                if (next.c < 255) ++next.c;
                break;
            case 2:
                if (next.g < 255) ++next.g;
                break;
            case 3:
                if (next.t < 255) ++next.t;
                break;
        }
    }
}

void readPlusPlus(const kseq_t *read) {
    uint64_t forward = 0;
    uint64_t reverse = 0;
    for (int i = 0, index = 1; i < read->seq.l; ++i, ++index) { // 1-based index
        int c = seq_fwd_table[(int) read->seq.s[i]];
        if (c < 4) {
            forward = (((forward << 2) | c) & ((1ULL<<k*2)-1));
            reverse = ((reverse >> 2) | (((uint64_t) (c^3)) << ((k*2)-2)));
        } else {
            index = 0;
        }
        if (index >= k) {
            if (i < read->seq.l-1) kmerPlusPlus(forward, seq_fwd_table[(int) read->seq.s[i+1]]); // OK
            if (i >= k) kmerPlusPlus(reverse, seq_fwd_table[(int) read->seq.s[i-k]]^3); // OK
        }
    }
    ++n_meta;
    if (!(n_meta % 10000)) fprintf(stderr, "%" PRIuMAX "\n", n_meta);
}

void filePlusPlus(const char *file) {
    gzFile fp = strcmp(file, "-") ? gzopen(file, "r") : gzdopen(fileno(stdin), "r");
    kseq_t *r = kseq_init(fp);
    while (kseq_read(r) >= 0) {
        readPlusPlus(r);
    }
    kseq_destroy(r);
    gzclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
// Error correction                                                           //
////////////////////////////////////////////////////////////////////////////////

inline int baseCorrect(const int base, const next_base_t next) {

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

int readCorrect(const kseq_t *read) {
    int failed = 0;
    int corrected = 0;
    uint64_t forward = 0;
    for (int i = 0, index = 1; i < read->seq.l-1; ++i, ++index) { // stops 1 char before end to look ahead!
        int c = seq_fwd_table[(int) read->seq.s[i]];
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
                int base = seq_fwd_table[(int)read->seq.s[i+1]]; // current base
                int sbase = baseCorrect(base, next);
                if (base != sbase) { // Correction needed
                    switch (sbase) {
                        case 0:
                            read->seq.s[i+1] = 'A';
                            ++corrected;
                            break;
                        case 1:
                            read->seq.s[i+1] = 'C';
                            ++corrected;
                            break;
                        case 2:
                            read->seq.s[i+1] = 'G';
                            ++corrected;
                            break;
                        case 3:
                            read->seq.s[i+1] = 'T';
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
void seq_revcomp(kseq_t *seq) {
    int c0, c1;
    for (int i = 0; i < seq->seq.l>>1; ++i) { // reverse complement sequence
        c0 = comp_tab[(int)seq->seq.s[i]];
        c1 = comp_tab[(int)seq->seq.s[seq->seq.l - 1 - i]];
        seq->seq.s[i] = c1;
        seq->seq.s[seq->seq.l - 1 - i] = c0;
    }
    if (seq->seq.l & 1) // complement the remaining base
    seq->seq.s[seq->seq.l>>1] = comp_tab[(int)seq->seq.s[seq->seq.l>>1]];
    if (seq->qual.l) {
        for (int i = 0; i < seq->seq.l>>1; ++i) // reverse quality
        c0 = seq->qual.s[i], seq->qual.s[i] = seq->qual.s[seq->qual.l - 1 - i], seq->qual.s[seq->qual.l - 1 - i] = c0;
    }
}
inline void stk_printseq(const kseq_t *s) {
    fputc(s->qual.l? '@' : '>', stdout);
    fputs(s->name.s, stdout);
    if (s->comment.l) {
        fputc(' ', stdout);
        fputs(s->comment.s, stdout);
    }
    fputc('\n', stdout);
    fputs(s->seq.s, stdout);
    fputc('\n', stdout);
    if (s->qual.l) {
        fputs("+\n", stdout);
        fputs(s->qual.s, stdout);
        fputc('\n', stdout);
    }
}

void fileCorrect(const char *file) {
    gzFile fp = strcmp(file, "-") ? gzopen(file, "r") : gzdopen(fileno(stdin), "r");
    kseq_t *r = kseq_init(fp);
    while (kseq_read(r) >= 0) {
        if (readCorrect(r)) {
            seq_revcomp(r);
            readCorrect(r);
            seq_revcomp(r);
        }
        stk_printseq(r);
        
        ++n_total2;
        if (!(n_total2 % 10000)) fprintf(stderr, "%" PRIuMAX " / %" PRIuMAX "\n", n_total2, n_total1);
    }
    kseq_destroy(r);
    gzclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
// main                                                                       //
////////////////////////////////////////////////////////////////////////////////

static int usage() {
    fprintf(stderr, "corsage version %s by Andreas Bremges (andreas@cebitec.uni-bielefeld.de)\n\n", VERSION);

    fprintf(stderr, "Usage: corsage [options] -s <SC.fastq> -m <MG.fastq>\n\n");

    fprintf(stderr, "       -s <SC.fastq>  single cell sequencing reads\n");
    fprintf(stderr, "       -m <MG.fastq>  metagenomic sequencing reads\n");
    fprintf(stderr, "                      (fastq or fasta, can be gzip'ed)\n\n");

    fprintf(stderr, "       -k INT         k-mer size for error correction [%i]\n", k);
    fprintf(stderr, "       -c INT         min. coverage in the metagenome [%i]\n\n", min_cov);
    return 42;
}

int main(int argc, char *argv[]) {
    char *one = 0, *two = 0;
    int c;
    while((c = getopt(argc, argv, "s:m:k:c:")) != -1) {
        switch (c) {
            case 's':
                one = optarg;
                break;
            case 'm':
                two = optarg;
                break;
            case 'k':
                k = atoi(optarg);
                if (k < 15) k = 15;
                if (k > 31) k = 31;
                break;
            case 'c':
                min_cov = atoi(optarg);
                if (min_cov < 1) min_cov = 1;
                break;
            default:
                return usage();
        }
    }
    if(!one || !two) return usage();
    h = kh_init(SAG);
    fprintf(stderr, "Init.\n");
    fileInitHash(one);
    if (!(n_total1 % 10000)) fprintf(stderr, "%" PRIuMAX "\n", n_total1);
    fprintf(stderr, "Fill.\n");
    filePlusPlus(two);
    if (!(n_meta % 10000)) fprintf(stderr, "%" PRIuMAX "\n", n_meta);
    fprintf(stderr, "Corr.\n");
    fileCorrect(one);
    fprintf(stderr, "Done.\n");
    if ((n_total2 % 10000)) fprintf(stderr, "%" PRIuMAX " / %" PRIuMAX "\n", n_total2, n_total1);

    kh_destroy(SAG, h);
    return 0;
}
