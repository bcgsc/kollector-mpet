#!/usr/bin/make -Rrf
SHELL=/bin/bash -o pipefail

#------------------------------------------------------------
# params
#------------------------------------------------------------

# threads
j?=1

#------------------------------------------------------------
# meta rules
#------------------------------------------------------------

.DELETE_ON_ERROR:
.PHONY:
default: $(name).misassembled-ids.txt $(name).sam2coord.tsv

check-params:
ifndef ref
	$(error missing required param 'ref' (reference genome))
endif
ifndef assembly
	$(error missing required param 'assembly' (assembled MPET fragments))
endif
ifndef name
	$(error missing required param 'name' (output file prefix))
endif

#------------------------------------------------------------
# main rules
#------------------------------------------------------------

# align assembled MPET frags to reference
$(name).sam.gz: $(assembly) $(ref)
	bwa mem -t$j $(ref) $(assembly) | gzip > $@

# identify misassembled MPET fragments
# (sequences with chimeric alignments)
$(name).misassembled-ids.txt: $(name).sam.gz
	zcat $< | awk '!/^@/ && and($$2,2048) {print $$1}' | \
		sort | uniq > $@

# get alignment coords and edit distances for each MPET fragment
$(name).sam2coord.tsv: $(name).misassembled-ids.txt $(name).sam.gz
	zcat $(name).sam.gz | \
		filter-ids $< - | \
		sam2coord --rlen --header > $@
