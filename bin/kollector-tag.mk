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
# min match score when recruiting seed PETs using seed MPETs
s?=0.5
# expected number of Bloom filter elements
n?=30000

#------------------------------------------------------------
# meta rules
#------------------------------------------------------------

.PHONY: check-name-param check-params seed clean

default: seed

seed: check-params $(name).seed.fa

check-name-param:
ifndef name
	$(error missing required param 'name' (output file prefix))
endif

check-params: check-name-param
ifndef seed_mp
	$(error missing required param 'seed_mp' (FASTA/FASTQ file(s)))
endif
ifndef pe
	$(error missing required param 'pe' (2 FASTA/FASTQ file(s)))
endif

clean: check-name-param
	rm -f $(name).seed.fa{,.fai} $(name).seed_mp.{bf,txt} $(name).seed_pe.fa.gz

#------------------------------------------------------------
# pipeline rules
#------------------------------------------------------------

# index FASTA file
%.fai: %
	samtools faidx $*

# build Bloom filter for seed MPET
$(name).seed_mp.bf: $(seed_mp).fai
	biobloommaker -k $k -p $(name).seed_mp -f $(max_fpr) -t $j \
		-n $n $(seed_mp) $(if $(subtract),-s $(subtract))

# get seed PET (PETs with single-end matches to seed MPET)
$(name).seed_pe.fa.gz: $(name).seed_mp.bf $(pe)
	biobloomcategorizer -t $j -d $(name).seed_mp -f $(name).seed_mp.bf -t $j \
		-s $s -e -i $(pe) | gzip > $@.partial
	mv $@.partial $@

# combine seed MPET reads and seed PET reads
$(name).seed.fa: $(seed_mp) $(name).seed_pe.fa.gz
	abyss-tofastq --fasta $^ > $@.partial
	mv $@.partial $@
