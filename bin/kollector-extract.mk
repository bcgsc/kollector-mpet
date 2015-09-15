#!/usr/bin/make -Rrf
SHELL=/bin/bash -o pipefail

#------------------------------------------------------------
# params
#------------------------------------------------------------

# min MPET fragment length
l?=500
# filename for length-filtered contigs (targeted assembly)
filtered_contigs?=$(name)-contigs.gt$lbp.fa
# threads
j?=1

#------------------------------------------------------------
# meta rules
#------------------------------------------------------------

.PHONY: check-params
default: $(name).fa.gz

check-params:
ifndef name
	$(error missing required param 'name' (output file prefix))
endif
ifndef contigs
	$(error missing required param 'contigs' (assembled contigs FASTA))
endif
ifndef mpet
	$(error missing required param 'mpet' (2 FASTA/FASTQ file(s)))
endif

#------------------------------------------------------------
# alignment rules
#------------------------------------------------------------

# build bwa index for FASTA file
%.bwt: %
	bwa index $<

# filter out short sequences
$(filtered_contigs): $(contigs)
	seqtk seq -L$l $< > $@.partial
	mv $@.partial $@

$(name).sam.gz: $(mpet) $(filtered_contigs).bwt
	bwa mem -a -t$j $(filtered_contigs) $(mpet) | \
		gzip > $@.partial
	mv $@.partial $@

$(name)-ids.txt: $(name).sam.gz
	zcat $< | \
		bioawk -c sam '$$tlen>0 && $$mapq>0 {print $$rname}' | \
		sort | \
		uniq \
		> $@.partial
	mv $@.partial $@

$(name).fa.gz: $(filtered_contigs) $(name)-ids.txt
	seqtk subseq $^ | gzip > $@.partial
	mv $@.partial $@
