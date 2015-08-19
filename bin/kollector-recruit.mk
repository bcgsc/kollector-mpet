#!/usr/bin/make -Rrf
SHELL=/bin/bash -o pipefail

#------------------------------------------------------------
# params
#------------------------------------------------------------

# k-mer size
k?=25
# threads
j?=1
# max false positive rate permitted for Bloom filters
max_fpr?=0.001
# min match score when recruiting PET reads
s?=0.75
# expected number of Bloom filter elements
n?=30000

#------------------------------------------------------------
# meta rules
#------------------------------------------------------------

.PRECIOUS: %.bf
.PHONY: check-params recruit

default: recruit

recruit: check-params $(name).fa

check-params:
ifndef name
	$(error missing required param 'name' (output file prefix))
endif
ifndef seed
	$(error missing required param 'seed' (FASTA file))
endif
ifndef pe
	$(error missing required param 'pe' (2 FASTA/FASTQ file(s)))
endif

clean:
	rm -f *.fai *.bf *.txt *.seed*.fa _summary.tsv \
		*-4.fa *.partial

#------------------------------------------------------------
# pipeline rules
#------------------------------------------------------------

# index FASTA file
%.fai: %
	samtools faidx $*

# iteratively add PETs with paired matches to Bloom filter
$(name).bf: $(seed).fai $(pe)
	biobloommaker -k $k -p $(name) -f $(max_fpr) -t $j -n $n \
		-r $s $(seed) $(pe)

# build FASTA for recruited PETs
$(name).fa: $(name).bf
	biobloomcategorizer -d $(name) -f $(name).bf -t $j -s 0.95 -e $(pe) | \
		abyss-tofastq --fasta > $@.partial
	mv $@.partial $@
	mv _summary.tsv $(name).tsv
