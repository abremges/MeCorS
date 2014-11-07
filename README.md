hector / corsage [![Build Status](https://magnum.travis-ci.com/abremges/hector.svg?token=Ebg4XZAcowyajM89NgpH&branch=master)](https://magnum.travis-ci.com/abremges/hector)
======

Hybrid Error Correction of Single Cell Sequencing Reads:

An erratic coverage profile and chimeric reads complicate the error correction of single cell data. Only one tool, BayesHammer (from the SPAdes assembler), works on single cell sequencing reads. However, nobody explored the possibility to use the corresponding metagenomic sequence information for single cell error correction.

Algorithm: For each base in the single cell read, we check if it is supported by this base in the metagenome with the same 23-mer prefix. Otherwise we try to correct the base. This is inspired by the error correction algorithm of Heng Li's fermi. If we detect an error but the metagenome offers several alternatives, we first skip this base and then try to correct it using the reverse complement of the read.

