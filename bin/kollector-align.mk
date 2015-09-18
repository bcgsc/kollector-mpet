#!/usr/bin/make -Rrf
SHELL=/bin/bash -o pipefail

#------------------------------------------------------------
# parameters
#------------------------------------------------------------

# threads
j?=1

#------------------------------------------------------------
# meta rules
#------------------------------------------------------------

.PRECIOUS: %.sam.gz %.sorted.bam %.sorted.bam.bai
.PHONY: align

default: align

align: check-params $(name).sorted.bam.bai

check-params:
ifndef name
	$(error missing required param 'name' (output file prefix))
endif
ifndef query
	$(error missing required param 'query' (interleaved FASTA))
endif
ifndef ref
	$(error missing required param 'ref' (ref FASTA/FASTQ file(s)))
endif

#------------------------------------------------------------
# alignment rules
#------------------------------------------------------------

$(name).sam.gz: $(query)
	bwa mem -t $j $(ref) <(abyss-tofastq --fasta $^ | awk '(NR-1)%4<=1') \
		<(abyss-tofastq --fasta $^ | awk '(NR-1)%4>=2') | \
		gzip > $@.partial
	mv $@.partial $@

%.bam: %.sam.gz
	zcat $< | samtools view -bSo $@ -

%.sorted.bam: %.bam
	samtools sort $< $*.sorted

%.bam.bai: %.bam
	samtools index $<
