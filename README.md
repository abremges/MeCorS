MeCorS [![Build Status](https://travis-ci.org/abremges/MeCorS.svg?branch=master)](https://travis-ci.org/abremges/MeCorS)
======

Metagenome-enabled error correction of single cell sequencing reads

## Installation
    git clone https://github.com/metagenomics/MeCorS.git
    cd MeCorS && make

## Usage
    mecors [options] -s <SC.fastq> -m <MG.fastq>

    -s <SC.fastq>  single cell sequencing reads
    -m <MG.fastq>  metagenomic sequencing reads
                   (fastq or fasta, can be gzip'ed)

    -k INT         k-mer size for error correction [31]
    -c INT         min. coverage in the metagenome [2]

    -t INT         number of threads [16]
    -B NUM         process ~NUM bp in each batch [100M]
    -v             be verbose

## Citation
MeCorS: Metagenome-enabled error correction of single cell sequencing reads  
[Andreas Bremges](https://github.com/abremges), Esther Singer, Tanja Woyke, Alexander Sczyrba  
*Bioinformatics* **(2016)** doi:[10.1093/bioinformatics/btw144](http://dx.doi.org/10.1093/bioinformatics/btw144)

======
*Heavily relies on on [klib](https://github.com/attractivechaos/klib).
Source code in part adopted from – or inspired by – [Heng Li](https://github.com/lh3).*
