corsage [![Build Status](https://magnum.travis-ci.com/abremges/corsage.svg?token=Ebg4XZAcowyajM89NgpH&branch=master)](https://magnum.travis-ci.com/abremges/corsage)
=====

Metagenome-enabled error correction of single cell sequencing reads

    Usage: corsage [options] -s <SC.fastq> -m <MG.fastq>
    
           -s <SC.fastq>  single cell sequencing reads
           -m <MG.fastq>  metagenomic sequencing reads
                          (fastq or fasta, can be gzip'ed)
    
           -k INT         k-mer size for error correction [31]
           -c INT         min. coverage in the metagenome [2]
    
           -t INT         number of threads [16]
           -B NUM         process ~NUM bp in each batch [100M]
           -v             be verbose

=====
*Heavily relies on on [klib](https://github.com/attractivechaos/klib).
Source code in part adopted from – or inspired by – [Heng Li](https://github.com/lh3).*
