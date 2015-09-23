# Description

Scripts for targeted assembly of MPET fragments using PET reads. These are wrapper scripts around [BioBloomTools](https://github.com/bcgsc/biobloom).

# Usage

```
Usage: kollector-mpet [options] <mpet_read1.fq> <mpet_read2.fq> \
   <pet_read1.fq> <pet_read2.fq> [<pet_read1.fq> <pet_read2.fq>]...

Do a targeted assembly of MPET fragment(s) using ABySS. The
input files are MPET and PET sequencing reads which must be provided
as FASTA/FASTQ pairs.  The input files may be gzipped.

The MPET read pairs act as anchors for recruiting PET reads for
the assembly. As such, the coverage depth of the MPET pairs is not
important to ensure a good outcome.  The most significant factors
affecting the results are the length of the MPET reads, the
sequencing error rate of the MPET reads, and the coverage depth
of the PET reads.

It is recommended to trim the adapter sequences from the MPET
reads before running this script.

Options:
    -a FILE   align reads to reference genome FILE after
              each iteration (for debugging) [disabled]
    -C        don't clean up intermediate files (for debugging)
    -h        show this help message
    -j N      threads [1]
    -k N      k-mer size for ABySS contig assembly [50]
    -K N      k-mer size for read overlap detection [25]
    -m N      max iterations for recruiting PET reads [10]
    -n N      max k-mers to recruit in total [25000]
    -o FILE   output file prefix ['kollector']
    -r FILE   Bloom filter containing repeat k-mers for
              exclusion from scoring calculations; must match
              k-mer size selected with -K opt [disabled]
    -s N      min score for recruiting seed PETs using seed MPETs;
              floating point number between 0 and 1 [0.5]
    -S N      min score for recruiting PETs from seed PETs and MPETs;
              floating point number between 0 and 1 [0.3]
```
