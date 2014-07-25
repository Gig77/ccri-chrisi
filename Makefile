export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=~/chrisi
TRIM_BEFORE_BASE=8

# download BAM files from IMP
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18376_ACATTA_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18377_GGTGAG_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18378_CGAAGG_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18379_AAGACA_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18380_TAATCG_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18381_CGCAAC_C4993ACXX_5_20140509B_20140509.bam &

# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18370_AATAGC_C4E7NACXX_8_20140603B_20140603.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18371_TTAACT_C4E7NACXX_8_20140603B_20140603.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18372_AATGAA_C4E7NACXX_8_20140603B_20140603.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18373_GATTGT_C4E7NACXX_8_20140603B_20140603.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18374_ATAAGA_C4E7NACXX_8_20140603B_20140603.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18375_GCCACA_C4E7NACXX_8_20140603B_20140603.bam &

# PRADA web site: http://bioinformatics.mdanderson.org/main/PRADA:Overview

# download PRADA reference (6 Gb)
#~/chrisi/data/prada/PRADA-reference.hg19.20130828.tar.gz: 
#	curl bioinformatics.mdanderson.org/Software/PRADA/PRADA-reference.hg19.20130828.tar.gz -o $@.part
#	mv $@.part $@

SAMPLES=18376_ACATTA_C4993ACXX_5_20140509B_20140509 \
		18378_CGAAGG_C4993ACXX_5_20140509B_20140509 \
		18380_TAATCG_C4993ACXX_5_20140509B_20140509 \
		18377_GGTGAG_C4993ACXX_5_20140509B_20140509 \
		18379_AAGACA_C4993ACXX_5_20140509B_20140509 \
		18381_CGCAAC_C4993ACXX_5_20140509B_20140509 \
		18370_AATAGC_C4E7NACXX_8_20140603B_20140603 \
		18371_TTAACT_C4E7NACXX_8_20140603B_20140603 \
		18372_AATGAA_C4E7NACXX_8_20140603B_20140603 \
		18373_GATTGT_C4E7NACXX_8_20140603B_20140603 \
		18374_ATAAGA_C4E7NACXX_8_20140603B_20140603 \
		18375_GCCACA_C4E7NACXX_8_20140603B_20140603 \
		test

all: gsnap htseq qc blast deseq

include ~/generic/scripts/rna-seq/gsnap.mk
include ~/generic/scripts/rna-seq/htseq.mk
include ~/generic/scripts/rna-seq/qc.mk
include ~/generic/scripts/rna-seq/blast.mk


.PHONY: deseq
deseq: deseq/strobl-dox-empty-vs-etv6.deseq2.chipseq-annotated.tsv deseq/zuber-etv6-nodox-vs-dox.deseq2.chipseq-annotated.tsv deseq/zuber-dox-empty-vs-etv6.deseq2.chipseq-annotated.tsv

deseq/strobl-dox-empty-vs-etv6.deseq2.tsv: ~/generic/scripts/rna-seq/diff-exp.R htseq/18378_CGAAGG_C4993ACXX_5_20140509B_20140509.count htseq/18379_AAGACA_C4993ACXX_5_20140509B_20140509.count htseq/18380_TAATCG_C4993ACXX_5_20140509B_20140509.count htseq/18381_CGCAAC_C4993ACXX_5_20140509B_20140509.count
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--control 18378:18378_CGAAGG_C4993ACXX_5_20140509B_20140509.count,18379:18379_AAGACA_C4993ACXX_5_20140509B_20140509.count \
		--experiment 18380:18380_TAATCG_C4993ACXX_5_20140509B_20140509.count,18381:18381_CGCAAC_C4993ACXX_5_20140509B_20140509.count \
		--output-tsv $@.part \
		--output-xls $(firstword $(subst .tsv, ,$@)).xlsx
	mv $@.part $@
	
deseq/zuber-etv6-nodox-vs-dox.deseq2.tsv: ~/generic/scripts/rna-seq/diff-exp.R htseq/18370_AATAGC_C4E7NACXX_8_20140603B_20140603.count htseq/18372_AATGAA_C4E7NACXX_8_20140603B_20140603.count htseq/18371_TTAACT_C4E7NACXX_8_20140603B_20140603.count htseq/18373_GATTGT_C4E7NACXX_8_20140603B_20140603.count
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--control 18370:18370_AATAGC_C4E7NACXX_8_20140603B_20140603.count,18372:18372_AATGAA_C4E7NACXX_8_20140603B_20140603.count \
		--experiment 18371:18371_TTAACT_C4E7NACXX_8_20140603B_20140603.count,18373:18373_GATTGT_C4E7NACXX_8_20140603B_20140603.count \
		--output-tsv $@.part \
		--output-xls $(firstword $(subst .tsv, ,$@)).xlsx
	mv $@.part $@

deseq/zuber-dox-empty-vs-etv6.deseq2.tsv: ~/generic/scripts/rna-seq/diff-exp.R htseq/18375_GCCACA_C4E7NACXX_8_20140603B_20140603.count htseq/18377_GGTGAG_C4993ACXX_5_20140509B_20140509.count htseq/18371_TTAACT_C4E7NACXX_8_20140603B_20140603.count htseq/18373_GATTGT_C4E7NACXX_8_20140603B_20140603.count 
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--control 18375:18375_GCCACA_C4E7NACXX_8_20140603B_20140603.count,18377:18377_GGTGAG_C4993ACXX_5_20140509B_20140509.count \
		--experiment 18371:18371_TTAACT_C4E7NACXX_8_20140603B_20140603.count,18373:18373_GATTGT_C4E7NACXX_8_20140603B_20140603.count \
		--output-tsv $@.part \
		--output-xls $(firstword $(subst .tsv, ,$@)).xlsx
	mv $@.part $@

deseq/%.chipseq-annotated.tsv: deseq/%.tsv ~/chrisi/scripts/annotate-runx1-chipseq-targets.R
	Rscript ~/chrisi/scripts/annotate-runx1-chipseq-targets.R --input $< --output $@.part
	mv $@.part $@

	