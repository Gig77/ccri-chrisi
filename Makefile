export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work

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


#%.fastq: %.bam
#	java -jar ~/tools/picard-tools-1.114/SamToFastq.jar INPUT=$< FASTQ=$@.part
#	mv $@.part $@

SAMPLES=18376_ACATTA_C4993ACXX_5_20140509B_20140509 \
		18378_CGAAGG_C4993ACXX_5_20140509B_20140509 \
		18380_TAATCG_C4993ACXX_5_20140509B_20140509 \
		C4993ACXX_5_20140509B_20140512 \
		18377_GGTGAG_C4993ACXX_5_20140509B_20140509 \
		18379_AAGACA_C4993ACXX_5_20140509B_20140509 \
		18381_CGCAAC_C4993ACXX_5_20140509B_20140509 \
		test

all: $(foreach S, $(SAMPLES), htseq/$S.count)

#----------------------------------------------
# BAM to FASTQ, align with GSNAP, sort & index
#----------------------------------------------
gsnap/%.gsnap.bam: ~/chrisi/data/bam/%.bam
	~/tools/gmap-2014-05-15/src/gsnap \
		--db=g1k_v37 \
		--dir=/data/christian/chrisi/data/current/gsnap/g1k_v37 \
		--format=sam \
		--max-mismatches=1 \
		--npaths=1 \
		--quiet-if-excessive \
		--nofails \
		--batch=4  \
		--quality-protocol=sanger \
		--print-snps \
		--nthreads=20 \
		--input-buffer-size=5000 \
		--use-splicing=g1k_v37.splicesites \
		--use-snps=g1k_v37.snp138 \
		<(java -jar ~/tools/picard-tools-1.114/SamToFastq.jar INPUT=$< FASTQ=/dev/stdout) \
	| ~/tools/samtools-0.1.19/samtools view -Sb - \
	| ~/tools/samtools-0.1.19/samtools sort -m 1000000000 - $@
	mv $@.bam $@
	~/tools/samtools-0.1.19/samtools reheader header.sam $@ > $@.reheader
	mv $@.reheader $@
	~/tools/samtools-0.1.19/samtools index $@
	
htseq/%.count: gsnap/%.gsnap.bam 
	~/tools/HTSeq-0.6.1/scripts/htseq-count -f bam -t exon -s no $< ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.gtf.gz > $@.part
	mv $@.part $@

htseq/combined.count: $(foreach S, $(SAMPLES), htseq/$S.count) ~/chrisi/scripts/combine-counts.R
	Rscript ~/chrisi/scripts/combine-counts.R $(foreach S, $(SAMPLES), htseq/$S.count) > $@.part
	mv $@.part $@
	
	
