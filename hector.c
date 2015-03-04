/* The MIT License

   Copyright (c) 2014 Andreas Bremges <andreas@cebitec.uni-bielefeld.de>

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

#define PROGRAM "hector"
#define DESCR   "Hybrid Error Correction of Single Cell Sequencing Reads"
#define VERSION "0.0.1"
#define AUTHOR  "Andreas Bremges"
#define CONTACT "andreas@cebitec.uni-bielefeld.de"

#include <stdint.h> // uint64_t
#include <stdio.h>  // printf
#include <stdlib.h> // qsort
#include <time.h>   // time
#include <unistd.h> // getopt
#include <zlib.h>   // gzip

#include "kseq.h"
KSEQ_INIT(gzFile, gzread) // TODO Check for alternatives (?)

#include "khash.h"
KHASH_MAP_INIT_INT64(SAG, uint32_t) // TODO Again, check usage in e.g. fermi
khash_t(SAG) *h; // TODO Also, rename from SAG

//unsigned int ctable[32] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648 };

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

int k = 23;

////////////////////////////////////////////////////////////////////////////////
// Init the data structure, i.e. insert keys into hash (single-cell k-mers)   //
////////////////////////////////////////////////////////////////////////////////

void kmerInitHash(const uint64_t kmer) {
    khiter_t iter = kh_get(SAG, h, kmer);
    if (iter == kh_end(h)) {
        int ret;
        iter = kh_put(SAG, h, kmer, &ret);
        kh_value(h, iter) = 0;
    }
}

void readInitHash(const kseq_t *read) {
    int index = 0;
    uint64_t forward = 0;
    uint64_t reverse = 0;
    uint64_t mask = (1ULL<<k*2) - 1;
    int i;
    for (i = 0; i < read->seq.l; ++i) {
        index++;
        int c = seq_fwd_table[(int) read->seq.s[i]];
        if (c < 4) {
            forward = (((forward << 2) | c) & mask); // see below
            //if (mask > 0) forward &= mask; <- only a problem if k = 32
            reverse = ((reverse >> 2) | (((uint64_t) (c^3)) << ((k*2)-2)));
        } else {
            index = 0;
        }
        if (index >= k) {
            //uint64_t z = ((forward < reverse) ? forward : reverse);
            kmerInitHash(forward);
            kmerInitHash(reverse);
        }
    }
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
        uint32_t c = (kh_value(h, iter) >> (base*8)) & 255;
        if (c < 31 && !(rand() % (1U<<c))) {//ctable[c])) {
            kh_value(h, iter) += (1U<<(base*8));
        }
    }
}

void readPlusPlus(const kseq_t *read) {
    int index = 0;
    uint64_t forward = 0;
    uint64_t reverse = 0;
    uint64_t mask = (1ULL<<k*2) - 1;
    int i;
    for (i = 0; i < read->seq.l; ++i) {
        index++;
        int c = seq_fwd_table[(int) read->seq.s[i]];
        if (c < 4) {
            forward = (((forward << 2) | c) & mask); // see below
            //if (mask > 0) forward &= mask; <- only a problem if k = 32
            reverse = ((reverse >> 2) | (((uint64_t) (c^3)) << ((k*2)-2))); // TODO Declare also shift constant?
        } else {
            //index = 0;
        }
        if (index >= k) {
            //uint64_t z = ((forward < reverse) ? forward : reverse); // TODO Store (k-1)mers to avoid the sequence lookups?
            if (i < read->seq.l-1) kmerPlusPlus(forward, seq_fwd_table[(int) read->seq.s[i+1]]); // TODO Check for boundaries
            if (i >= k) kmerPlusPlus(reverse, seq_fwd_table[(int) read->seq.s[i-k]]^3); // Double-check here

            // 			uint64_t z = reverse;
            // 			char kmer[k];
            // 			int j;
            // 			for (j = k-1; j >= 0; --j) {
            // 				kmer[j] = seq_rev_table[z&3];
            // 				z = (z >> 2);
            // 			}
            // 			fprintf(stdout, "%s\t%i\n", kmer, seq_fwd_table[(int) read->seq.s[i-k]]^3);
        }
    }
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

#define MASK 0xFFFFFFFFU
int readCorrect(const kseq_t *read) {
    int failed = 0;
    int index = 0;
    uint64_t forward = 0;
    uint64_t mask = (1ULL<<k*2) - 1;
    int i;
    for (i = 0; i < read->seq.l-1; ++i) { // stops 1 char before end to look ahead!
        index++;
        int c = seq_fwd_table[(int) read->seq.s[i]];
        if (c < 4) {
            forward = (((forward << 2) | c) & mask); // see below
            //if (mask > 0) forward &= mask; <- only a problem if k = 32
            //reverse = ((reverse >> 2) | (((uint64_t) (c^3)) << ((k*2)-2))); // TODO Declare also shift constant?
        } else {
            index = 0;
            ++failed; // TODO What happens here?
        }
        if (index >= k) {
            khiter_t iter = kh_get(SAG, h, forward);
            if (iter != kh_end(h)) {
                uint32_t tmp = kh_value(h, iter);
                if (seq_fwd_table[(int)read->seq.s[i+1]] >= 4 || !(tmp&(255<<(seq_fwd_table[(int)read->seq.s[i+1]]<<3)))) { // no support through MG
                    if (tmp^(MASK&(255<<(seq_fwd_table[(int)read->seq.s[i+1]]<<3)))) { // but there is at least one alternative
                        if (tmp&(255<<(seq_fwd_table['A']<<3)) && !(tmp&(MASK^(255<<(seq_fwd_table['A']<<3))))) {
                            read->seq.s[i+1] = 'A';
                        } else if (tmp&(255<<(seq_fwd_table['C']<<3)) && !(tmp&(MASK^(255<<(seq_fwd_table['C']<<3))))) {
                            read->seq.s[i+1] = 'C';
                        } else if (tmp&(255<<(seq_fwd_table['G']<<3)) && !(tmp&(MASK^(255<<(seq_fwd_table['G']<<3))))) {
                            read->seq.s[i+1] = 'G';
                        } else if (tmp&(255<<(seq_fwd_table['T']<<3)) && !(tmp&(MASK^(255<<(seq_fwd_table['T']<<3))))) { // TODO Remove ugly code redundancy
                            read->seq.s[i+1] = 'T';
                        } else { // cannot decide
                            index = -1; // TODO think about it: 0 or -1?
                            ++failed; // TODO Consider coverage?
                        }

                    } else { // no support at all
                        index = -1;
                        ++failed;
                    }
                }
            } // TODO what else?
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
//////////

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
    // TODO Work on a copy of the kseq_t, to allow to output original read if [below]
    while (kseq_read(r) >= 0) {
        if (readCorrect(r)) {
            seq_revcomp(r);
            readCorrect(r);
            seq_revcomp(r);
        }
        stk_printseq(r); // TODO Check number of corrections, suppress overly aggressiveness
    }
    kseq_destroy(r);
    gzclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
// main                                                                       //
////////////////////////////////////////////////////////////////////////////////

static int usage() {
    fprintf(stderr, "pcount [-k INT] -1 <SC.fastq> -2 <MG.fastq>\n");
    return 42;
}

int main(int argc, char *argv[]) {
    char *one = 0, *two = 0;
    int c;
    while((c = getopt(argc, argv, "1:2:k:")) != -1) {
        switch (c) {
            case '1':
            one = optarg;
            break;
            case '2':
            two = optarg;
            break;
            case 'k':
            k = atoi(optarg);
            //if (k < 15) k = 15;
            if (k > 31) k = 31;
            break;
            default:
            return usage();
        }
    }
    if(!one || !two) return usage();
    srand(time(NULL)); // stupid rand() should be sufficient for our needs

    h = kh_init(SAG);
    fprintf(stderr, "Init.\n");
    fileInitHash(one);
    fprintf(stderr, "Fill.\n");
    filePlusPlus(two);
    fprintf(stderr, "Corr.\n");
    fileCorrect(one);
    fprintf(stderr, "Done.\n");

    // 	khiter_t iter;
    // 	for (iter = kh_begin(h); iter != kh_end(h); ++iter) {
    // 		if (kh_exist(h, iter)) {
    // 			char kmer[k];
    // 			int i;
    // 			uint64_t z = kh_key(h, iter);
    // 			for (i = k-1; i >= 0; --i) {
    // 				kmer[i] = seq_rev_table[z&3];
    // 				z = (z >> 2);
    // 			}
    // 			int a = 0, c = 0, g = 0, t = 0;
    // 			uint32_t tmp = kh_value(h, iter);
    // 			a = (tmp >> (seq_fwd_table['A']<<3)) & 255;
    // 			c = (tmp >> (seq_fwd_table['C']<<3)) & 255;
    // 			g = (tmp >> (seq_fwd_table['G']<<3)) & 255;
    // 			t = (tmp >> (seq_fwd_table['T']<<3)) & 255;
    // 			fprintf(stdout, "%s %i,%i,%i,%i ", kmer, a, c, g, t);
    // 			if (tmp&(255<<(seq_fwd_table['A']<<3))) {
    // 				fprintf(stdout, "A");
    // 			}
    // 			if (tmp&(255<<(seq_fwd_table['C']<<3))) {
    // 				fprintf(stdout, "C");
    // 			}
    // 			if (tmp&(255<<(seq_fwd_table['G']<<3))) {
    // 				fprintf(stdout, "G");
    // 			}
    // 			if (tmp&(255<<(seq_fwd_table['T']<<3))) {
    // 				fprintf(stdout, "T");
    // 			}
    // 			fprintf(stdout, "\n");
    // 		}
    // 	}

    kh_destroy(SAG, h);
    return 0;
}
