export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

# download BAM files from IMP
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18376_ACATTA_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18377_GGTGAG_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18378_CGAAGG_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18379_AAGACA_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18380_TAATCG_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/18381_CGCAAC_C4993ACXX_5_20140509B_20140509.bam &
# wget -c --quiet --no-check-certificate --auth-no-challenge --user 'Christine.Portsmouth' --password 'RS5YC4LE8E' http://ngs.csf.ac.at/data/C4993ACXX_5_20140509B_20140512.bam &

# PRADA web site: http://bioinformatics.mdanderson.org/main/PRADA:Overview

# download PRADA reference (6 Gb)
#~/chrisi/data/prada/PRADA-reference.hg19.20130828.tar.gz: 
#	curl bioinformatics.mdanderson.org/Software/PRADA/PRADA-reference.hg19.20130828.tar.gz -o $@.part
#	mv $@.part $@
	
# build GMAP database
# cd ~/chrisi/data/gsnap
# ~/tools/gmap-2014-05-15/util/gmap_build -d g1k_v37 ~/generic/data/broad/human_g1k_v37.fasta
# cat ~/generic/data/broad/human_g1k_v37.fasta ~/chrisi/results/etv6-runx1.breakpoint.fa > ~/chrisi/data/gsnap/human_g1k_v37_etv6runx1.fasta 
# ~/tools/gmap-2014-05-15/bin/gmap_build --dir=/data/christian/chrisi/data/current/gsnap --db=g1k_v37_etv6runx1 ~/chrisi/data/gsnap/human_g1k_v37_etv6runx1.fasta
# rm ~/chrisi/data/gsnap/human_g1k_v37_etv6runx1.fasta
#
# gunzip -c snp138.txt.gz | ~/tools/gmap-2014-05-15/bin/dbsnp_iit -w 1 > g1k_v37.snp138.txt
# cd g1k_v37_etv6runx1
# cat g1k_v37.snp138.txt | ~/tools/gmap-2014-05-15/bin/iit_store -o g1k_v37.snp138
# ~/tools/gmap-2014-05-15/bin/snpindex -d g1k_v37_etv6runx1 -D g1k_v37_etv6runx1 -v snpfile138
# mv g1k_v37_etv6runx1/g1k_v37.snp138.iit g1k_v37_etv6runx1/g1k_v37_etv6runx1.maps

# for splicesites
# ?? gunzip -c refGene.txt.gz | psl_splicesites -s 1 > mm10splicesites #ok
# cat g1k_v37.splicesites | ~/tools/gmap-2014-05-15/bin/iit_store -o g1k_v37_etv6runx1.maps/g1k_v37.splicesites




#%.fastq: %.bam
#	java -jar ~/tools/picard-tools-1.114/SamToFastq.jar INPUT=$< FASTQ=$@.part
#	mv $@.part $@

SAMPLES=18376_ACATTA_C4993ACXX_5_20140509B_20140509 \
		18378_CGAAGG_C4993ACXX_5_20140509B_20140509 \
		18380_TAATCG_C4993ACXX_5_20140509B_20140509 \
		18377_GGTGAG_C4993ACXX_5_20140509B_20140509 \
		18379_AAGACA_C4993ACXX_5_20140509B_20140509 \
		18381_CGCAAC_C4993ACXX_5_20140509B_20140509 \
		test

all: $(foreach S, $(SAMPLES), htseq/$S.count flagstat/$S.samtools.flagstat)

#----------------------------------------------
# BAM to FASTQ, align with GSNAP, sort & index
#----------------------------------------------
gsnap/%.gsnap.bam: ~/chrisi/data/bam/%.bam
	~/tools/gmap-2014-05-15/src/gsnap \
			--db=g1k_v37_etv6runx1 \
			--dir=/data/christian/chrisi/data/current/gsnap/g1k_v37_etv6runx1 \
			--format=sam \
			--npaths=1 \
			--quiet-if-excessive \
			--nofails \
			--batch=4  \
			--quality-protocol=sanger \
			--print-snps \
			--nthreads=24 \
			--input-buffer-size=5000 \
			--use-splicing=g1k_v37.splicesites \
			--use-snps=g1k_v37.snp138 \
			--genome-unk-mismatch=0 \
			<(java -jar ~/tools/picard-tools-1.114/SamToFastq.jar INPUT=$< FASTQ=/dev/stdout) \
		| ~/tools/samtools-0.1.19/samtools view -Shb - \
		| ~/tools/samtools-0.1.19/samtools sort -m 1000000000 - $@ \
		2>&1 | $(LOG)
	mv $@.bam $@
	~/tools/samtools-0.1.19/samtools index $@ 2>&1 | $(LOG)

flagstat/%.samtools.flagstat: gsnap/%.gsnap.bam
	samtools flagstat $< 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

htseq/%.count: gsnap/%.gsnap.bam ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz
	~/tools/HTSeq-0.6.1/scripts/htseq-count -f bam -t exon -s no $< ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

htseq/combined.count: $(foreach S, $(SAMPLES), htseq/$S.count) ~/chrisi/scripts/combine-counts.R
	Rscript ~/chrisi/scripts/combine-counts.R $(foreach S, $(SAMPLES), htseq/$S.count) 2>&1 1>$@.part | $(LOG)
	mv $@.part $@
	
	
