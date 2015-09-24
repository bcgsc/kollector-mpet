# Description

Scripts for targeted assembly of MPET fragments using PET reads. These are wrapper scripts around [BioBloomTools](https://github.com/bcgsc/biobloom).

# Usage

```
$ kollector-mpet -h
Usage: $PROGRAM [options] <mpet_read1.fq> <mpet_read2.fq> \
   <pet_read1.fq> <pet_read2.fq> [<pet_read1.fq> <pet_read2.fq>]...

Description:

Do a targeted assembly of MPET fragment(s) using ABySS. The
input files are MPET and PET sequencing reads which must be provided
as FASTA/FASTQ pairs.  The input files may be gzipped.

The MPET read pairs act as anchors for recruiting PET reads for
the assembly. As such, the coverage depth of the MPET pairs is not
important to ensure a good outcome.  The most significant factors
affecting the results are the the length of the MPET reads, the
sequencing error rate of the MPET reads, and the coverage depth
of the PET reads.

It is recommended to trim the adapter sequences from the MPET
reads before running this script.

Options:

    -a        align reads to reference genome (specified by -g)
              after each iteration (for debugging) [disabled]
    -C        don't clean up intermediate files (for debugging)
    -e        evaluate assembled MPET fragments by aligning
              to the reference sequence (specified by -g)
              [disabled]
    -g FILE   use FILE for the reference genome (required if
              using -a or -e options) [disabled]
    -h        show this help message
    -j N      threads [1]
    -l N      min match length for recruiting seed PETs using
              seed MPETs as bait; only one of the PET reads
              of each pair needs to match [30]
    -L N      min match length for recruiting PETs using
              previously recruited PETs as bait; both reads
              of a PET pair must match with at least this
              length [15]
    -k N      k-mer size for ABySS contig assembly [50]
    -K N      k-mer size for read overlap detection [25]
    -m N      max iterations for recruiting PET reads [10]
    -n N      max k-mers to recruit in total [25000]
    -o FILE   output file prefix ['kollector']
    -r FILE   Bloom filter containing repeat k-mers for
              exclusion from scoring calculations; must match
              k-mer size selected with -K opt [disabled]
    -R FILE   print summary report to FILE (see description
              of fields below) [disabled]

Summary report:

    (1) 'mpet': num input MPET pairs
    (2) 'mpet_assembled': num MPET frags assembled
    (3) 'percent_mpet_assembled': percent MPET frags assembled
    (4) 'kmers_recruited': num distinct k-mers recruited
    (5) 'N50': N50 assembly metric
    (6) 'sum': total number of bases assembled
    (7) 'wallclock_time': running time in seconds
    (8) 'peak_disk_usage': peak disk space usage in bytes
```
