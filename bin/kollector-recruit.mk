#!/usr/bin/make -Rrf
ifdef profile
SHELL=/usr/bin/time -f '=> kollector-recruit.mk: %e %C' /bin/bash -o pipefail
else
SHELL=/bin/bash -o pipefail
endif

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
.DELETE_ON_ERROR:
.PHONY: check-params recruit

default: recruit

recruit: check-params $(name).fa.gz

check-name-param:
ifndef name
	$(error missing required param 'name' (output file prefix))
endif

check-params: check-name-param
ifndef seed
	$(error missing required param 'seed' (FASTA file))
endif
ifndef pe
	$(error missing required param 'pe' (2 FASTA/FASTQ file(s)))
endif

clean: check-name-param
	rm -f $(name).{fa,fa.gz,fa.fai,bf,txt,tsv}

#------------------------------------------------------------
# pipeline rules
#------------------------------------------------------------

# index FASTA file
%.fai: %
	samtools faidx $*

# iteratively add PETs with paired matches to Bloom filter
$(name).fa.gz: $(seed).fai $(pe)
	biobloommaker -k $k -p $(name) -f $(max_fpr) -t $j -n $n \
		-r $s -P $(if $(subtract),-s $(subtract)) \
		$(if $h,-g $h) $(seed) $(pe) | \
		seqtk seq -A | \
		gzip > $@
