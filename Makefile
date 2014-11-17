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

# call snps from rna-seq
# java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2snp <(~/tools/samtools-0.1.19/samtools view -b -u -q 1 18378_CGAAGG_C4993ACXX_5_20140509B_20140509.gsnap.filtered.bam | ~/tools/samtools-0.1.19/samtools mpileup -f ~/chrisi/data/gsnap/human_g1k_v37_etv6runx1.fasta -) --min-avg-qual 20 --min-coverage 10 --min-reads2 2 --p-value 1 --strand-filter 0 --min-var-freq 0.01 | tee 18378.varscan.snp
# cat 18378.varscan.snp | perl -ne 'print "$1\t$2\n" if (/:(\d+):\d:\d:(\d+)\%/);' > 18378.varscan.snp.freq

# check for hygromycine
#java -jar ~/tools/picard-tools-1.114/SamToFastq.jar VALIDATION_STRINGENCY=SILENT INPUT=<(~/tools/samtools-0.1.19/samtools view -f 4 gsnap/18372_AATGAA_C4E7NACXX_8_20140603B_20140603.gsnap.bam) FASTQ=18372.unmapped.fastq
#~/tools/bwa-0.7.9/bwa aln -t 10 hygromycine.fa 18372.unmapped.fastq > 18372.unmapped.sai
#~/tools/bwa-0.7.9/bwa samse hygromycine.fa 18372.unmapped.sai 18372.unmapped.fastq > 18372.unmapped.sam
#~/tools/samtools-0.1.19/samtools view -SF 772 18372.unmapped.sam

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

all: gsnap htseq qc blast deseq fastqc

include ~/generic/scripts/rna-seq/gsnap.mk
include ~/generic/scripts/rna-seq/htseq.mk
include ~/generic/scripts/rna-seq/fastqc.mk
include ~/generic/scripts/rna-seq/qc.mk
include ~/generic/scripts/rna-seq/blast.mk


.PHONY: deseq
deseq: deseq/strobl-dox-empty-vs-etv6.deseq2.chipseq-annotated.tsv deseq/zuber-etv6-nodox-vs-dox.deseq2.chipseq-annotated.tsv deseq/zuber-dox-empty-vs-etv6.deseq2.chipseq-annotated.tsv deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.chipseq-annotated.fuka-ross-boer.tsv

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

deseq/zuber-nodox-empty-vs-etv6.deseq2.tsv: ~/generic/scripts/rna-seq/diff-exp.R htseq/18374_ATAAGA_C4E7NACXX_8_20140603B_20140603.count htseq/18376_ACATTA_C4993ACXX_5_20140509B_20140509.count htseq/18370_AATAGC_C4E7NACXX_8_20140603B_20140603.count htseq/18372_AATGAA_C4E7NACXX_8_20140603B_20140603.count
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--control 18374:18374_ATAAGA_C4E7NACXX_8_20140603B_20140603.count,18376:18376_ACATTA_C4993ACXX_5_20140509B_20140509.count \
		--experiment 18370:18370_AATAGC_C4E7NACXX_8_20140603B_20140603.count,18372:18372_AATGAA_C4E7NACXX_8_20140603B_20140603.count \
		--output-tsv $@.part
	mv $@.part $@

deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.tsv: ~/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	Rscript ~/generic/scripts/rna-seq/diff-exp.R \
		--control 18370:18370_AATAGC_C4E7NACXX_8_20140603B_20140603.count,18372:18372_AATGAA_C4E7NACXX_8_20140603B_20140603.count,18374:18374_ATAAGA_C4E7NACXX_8_20140603B_20140603.count,18375:18375_GCCACA_C4E7NACXX_8_20140603B_20140603.count,18376:18376_ACATTA_C4993ACXX_5_20140509B_20140509.count,18377:18377_GGTGAG_C4993ACXX_5_20140509B_20140509.count,18378:18378_CGAAGG_C4993ACXX_5_20140509B_20140509.count,18379:18379_AAGACA_C4993ACXX_5_20140509B_20140509.count \
		--experiment 18373:18373_GATTGT_C4E7NACXX_8_20140603B_20140603.count,18371:18371_TTAACT_C4E7NACXX_8_20140603B_20140603.count,18380:18380_TAATCG_C4993ACXX_5_20140509B_20140509.count,18381:18381_CGCAAC_C4993ACXX_5_20140509B_20140509.count \
		--output-tsv $@.part
	mv $@.part $@

deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.chipseq-annotated.fuka-ross-boer.tsv: deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.chipseq-annotated.tsv
	Rscript ~/chrisi/scripts/annotate-fuka-ross-boer.R
	
deseq/%.chipseq-annotated.tsv: deseq/%.tsv ~/chrisi/scripts/annotate-runx1-chipseq-targets.R
	Rscript ~/chrisi/scripts/annotate-runx1-chipseq-targets.R --input $< --output $@.part
	mv $@.part $@
		
