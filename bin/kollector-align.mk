#!/usr/bin/make -Rrf
SHELL=/bin/bash -o pipefail

#------------------------------------------------------------
# meta rules
#------------------------------------------------------------

.PRECIOUS: %.sam.gz %.sorted.bam %.sorted.bam.bai
.PHONY: align

default: align

align: check-params $(query).sorted.bam.bai

check-params:
ifndef query
	$(error missing required param 'query' (interleaved FASTA))
endif
ifndef ref
	$(error missing required param 'ref' (ref FASTA/FASTQ file(s)))
endif

#------------------------------------------------------------
# alignment rules
#------------------------------------------------------------

%.sam.gz: %
	bwa mem $(ref) <(abyss-tofastq --fasta $* | awk '(NR-1)%4<=1') \
		<(abyss-tofastq --fasta $* | awk '(NR-1)%4>=2') | \
		gzip > $@.partial
	mv $@.partial $@

%.bam: %.sam.gz
	zcat $< | samtools view -bSo $@ -

%.sorted.bam: %.bam
	samtools sort $< $*.sorted

%.bam.bai: %.bam
	samtools index $<
