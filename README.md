Script for targeted assembly of MPET fragments using PET reads. These are wrapper scripts around [BioBloomTools](https://github.com/bcgsc/biobloom).

```
$ kollector-mpet
Usage: kollector-mpet [options] <seed_mpet.fa> <pet_read1.fq> <pet_read2.fq>
Gather PET reads for targeted assembly of MPET fragment(s).

Options:
    -a FILE   align reads to reference genome FILE after
              each iteration [disabled]
    -k N      k-mer size for read overlap detection [25]
    -m N      max iterations for recruiting PET reads [10]
    -n N      max k-mers to recruit in total [25000]
    -o FILE   output file prefix ['kollector']
    -s N      min match score for recruiting reads; floating
              point number between 0 and 1 [0.3]
```
