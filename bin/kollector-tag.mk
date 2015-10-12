#!/usr/bin/make -Rrf
ifdef profile
SHELL=/usr/bin/time -f '=> kollector-tag.mk: %e %C' /bin/bash -o pipefail
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
ifndef mp
	$(error missing required param 'mp' (MPET FASTA/FASTQ file(s)))
endif
ifndef pe
	$(error missing required param 'pe' (FASTA/FASTQ file pair(s)))
endif

clean: check-name-param
	rm -f $(name).seed.fa{,.fai} $(name).seed_mp.{bf,txt} $(name).seed_pe.fa

#------------------------------------------------------------
# pipeline rules
#------------------------------------------------------------

# index FASTA file
%.fai: %
	samtools faidx $*

# interleave MPET and convert to FASTA
# (BioBloomTools requires a single interleaved FASTA for input
# seed sequences.)
$(name).seed_mp.fa: $(mp)
	seqtk mergepe $(mp) | seqtk seq -A > $@.partial
	mv $@.partial $@

# build Bloom filter for seed MPET
$(name).seed_mp.bf: $(name).seed_mp.fa.fai
	biobloommaker -k $k -p $(name).seed_mp -f $(max_fpr) -t $j \
		-n $n $(name).seed_mp.fa $(if $(subtract),-s $(subtract)) \
		$(if $h,-g $h)

# get seed PETs (PETs with single-end matches to seed MPETs)
$(name).seed_pe.fa: $(name).seed_mp.bf $(pe)
	rm -f $@.partial
	for i in $$(seq 1 2 $(words $(pe))); do \
		file1=$$(echo $(pe) | tr -s ' ' '\t' | cut -f$$i); \
		file2=$$(echo $(pe) | tr -s ' ' '\t' | cut -f$$(($$i+1))); \
		biobloomcategorizer -t $j -d $(name).seed_mp -f $(name).seed_mp.bf -t $j \
			-s $s -e -i $$file1 $$file2 >> $@.partial; \
	done
	mv $@.partial $@

# combine seed MPET reads and seed PET reads
$(name).seed.fa: $(name).seed_mp.fa $(name).seed_pe.fa
	abyss-tofastq --fasta $^ > $@.partial
	mv $@.partial $@
