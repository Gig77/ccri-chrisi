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

# FASTA reference genome
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
	
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
		18370_AATAGC_C4E7NACXX_8_20140603B_20140603 \
		18371_TTAACT_C4E7NACXX_8_20140603B_20140603 \
		18372_AATGAA_C4E7NACXX_8_20140603B_20140603 \
		18373_GATTGT_C4E7NACXX_8_20140603B_20140603 \
		18374_ATAAGA_C4E7NACXX_8_20140603B_20140603 \
		18375_GCCACA_C4E7NACXX_8_20140603B_20140603 \
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

gsnap/%.gsnap.unpaired_uniq.bam: ~/chrisi/data/bam/%.bam
	java -jar ~/tools/picard-tools-1.114/SamToFastq.jar INPUT=$< FASTQ=/dev/stdout \
		| ~/tools/gmap-2014-05-15/src/gsnap \
			--db=g1k_v37_etv6runx1 \
			--dir=/data/christian/chrisi/data/current/gsnap/g1k_v37_etv6runx1 \
			--format=sam \
			--npaths=1000 \
			--quiet-if-excessive \
			--nofails \
			--batch=4 \
			--quality-protocol=sanger \
			--print-snps \
			--nthreads=20 \
			--input-buffer-size=5000 \
			--use-splicing=g1k_v37.splicesites \
			--use-snps=g1k_v37.snp138 \
			--genome-unk-mismatch=0 \
			--split-output=gsnap/$*.gsnap \
			/dev/stdin
			
	# samtools view threw error if no reads in SAM file; required patching file bam_index.c; see http://sourceforge.net/p/samtools/mailman/message/28264249/
	~/tools/samtools-0.1.19/samtools view -Shb gsnap/$*.gsnap.unpaired_uniq | ~/tools/samtools-0.1.19/samtools sort -m 1000000000 - gsnap/$*.gsnap.unpaired_uniq ; rm gsnap/$*.gsnap.unpaired_uniq
	~/tools/samtools-0.1.19/samtools view -Shb gsnap/$*.gsnap.unpaired_mult | ~/tools/samtools-0.1.19/samtools sort -m 1000000000 - gsnap/$*.gsnap.unpaired_mult ; rm gsnap/$*.gsnap.unpaired_mult
	~/tools/samtools-0.1.19/samtools view -Shb gsnap/$*.gsnap.unpaired_mult_xs | ~/tools/samtools-0.1.19/samtools sort -m 1000000000 - gsnap/$*.gsnap.unpaired_mult_xs ; rm gsnap/$*.gsnap.unpaired_mult_xs
	~/tools/samtools-0.1.19/samtools view -Shb gsnap/$*.gsnap.nomapping | ~/tools/samtools-0.1.19/samtools sort -m 1000000000 - gsnap/$*.gsnap.nomapping ; rm gsnap/$*.gsnap.nomapping
	~/tools/samtools-0.1.19/samtools view -Shb gsnap/$*.gsnap.unpaired_transloc | ~/tools/samtools-0.1.19/samtools sort -m 1000000000 - gsnap/$*.gsnap.unpaired_transloc ; rm gsnap/$*.gsnap.unpaired_transloc
	~/tools/samtools-0.1.19/samtools view -Shb gsnap/$*.gsnap.unpaired_circular | ~/tools/samtools-0.1.19/samtools sort -m 1000000000 - gsnap/$*.gsnap.unpaired_circular ; rm gsnap/$*.gsnap.unpaired_circular
	~/tools/samtools-0.1.19/samtools index $@
		
flagstat/%.samtools.flagstat: gsnap/%.gsnap.bam
	samtools flagstat $< 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

htseq/%.count: gsnap/%.gsnap.bam ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.no-rRNA.gtf.gz
	~/tools/HTSeq-0.6.1/scripts/htseq-count -f bam -t exon -s no $< ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.no-rRNA.gtf.gz 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

#----
# quality control
#----

# generate BED files from Ensembl GTF
# zcat ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz | perl -ne 'print "$1\t$2\t$3\t$4\n" if (/^(.*?)\t.*?\texon\t(.*?)\t(.*?)\t.*gene_name \"(.*?)\"/)' > ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.all.exons.bed
# zcat ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz | perl -ne 'print "$1\t$2\t$3\t$4\n" if (/^(.*?)\trRNA\texon\t(.*?)\t(.*?)\t.*gene_name \"(.*?)\"/)' > ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.rRNA.exons.bed
# zcat ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz | perl -ne 'print "$1\t$2\t$3\t$4\n" if (/^(.*?)\tlincRNA\texon\t(.*?)\t(.*?)\t.*gene_name \"(.*?)\"/)' > ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.lincRNA.exons.bed
# zcat ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz | perl -ne 'print "$1\t$2\t$3\t$4\n" if (/^(.*?)\tmiRNA\texon\t(.*?)\t(.*?)\t.*gene_name \"(.*?)\"/)' > ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.miRNA.exons.bed
# zcat ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz | perl -ne 'print "$1\t$3\t$4\t$5\n" if (/^(.*?)\t(Mt_rRNA|Mt_tRNA)\texon\t(.*?)\t(.*?)\t.*gene_name \"(.*?)\"/)' > ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.mtRNA.exons.bed
# zcat ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz | perl -ne 'print "$1\t$3\t$4\t$5\n" if (/^(.*?)\t(misc_RNA|snoRNA|snRNA)\texon\t(.*?)\t(.*?)\t.*gene_name \"(.*?)\"/)' > ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.miscRNA.exons.bed
# zcat ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz | perl -ne 'print "$1\t$3\t$4\t$5\n" if (/^(.*?)\t(nonsense_mediated_decay|non_stop_decay|polymorphic_pseudogene|processed_pseudogene|pseudogene|transcribed_processed_pseudogene|transcribed_unprocessed_pseudogene|translated_processed_pseudogene|unitary_pseudogene|unprocessed_pseudogene)\texon\t(.*?)\t(.*?)\t.*gene_name \"(.*?)\"/)' > ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.pseudogene.exons.bed
# zcat ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz | perl -ne 'print "$1\t$3\t$4\t$5\n" if (/^(.*?)\t(protein_coding)\texon\t(.*?)\t(.*?)\t.*gene_name \"(.*?)\"/)' > ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.protein_coding.exons.bed
# ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.other.exons.bed was created by subtracting the above features from ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.all.exons.bed using '~/tools/bedtools-2.17.0/bin/subtractBed -A'

quality: $(foreach S, $(SAMPLES), rseqc/$S.saturation.pdf) htseq/coverage-saturation-curve.pdf rseqc/allpatients.rRNA.count rseqc/allpatients.read-distribution.txt rseqc/allpatients.splice_junction.pdf rseqc/allpatients.splicing_events.pdf rseqc/allpatients.geneBodyCoverage.pdf rseqc/allpatients.DupRate_plot.pdf rseqc/allpatients.infer_experiment.txt

rseqc/%.stats.txt: gsnap/%.gsnap.unpaired_uniq.bam
	echo -ne $*"\t" > $@.part 
	echo -ne `samtools view -c ~/chrisi/data/bam/$*.bam`"\t" >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.nomapping.bam`"\t" >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.unpaired_mult.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.all.exons.bed`"\n"\
		     `samtools view -c gsnap/$*.gsnap.unpaired_uniq.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.all.exons.bed` | perl -ne '$$sum += $$_; END {print "$$sum\t"}' >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.unpaired_mult.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.rRNA.exons.bed`"\n"\
		     `samtools view -c gsnap/$*.gsnap.unpaired_uniq.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.rRNA.exons.bed` | perl -ne '$$sum += $$_; END {print "$$sum\t"}' >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.unpaired_mult.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.lincRNA.exons.bed`"\n"\
		     `samtools view -c gsnap/$*.gsnap.unpaired_uniq.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.lincRNA.exons.bed` | perl -ne '$$sum += $$_; END {print "$$sum\t"}' >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.unpaired_mult.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.miRNA.exons.bed`"\n"\
		     `samtools view -c gsnap/$*.gsnap.unpaired_uniq.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.miRNA.exons.bed` | perl -ne '$$sum += $$_; END {print "$$sum\t"}' >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.unpaired_mult.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.mtRNA.exons.bed`"\n"\
		     `samtools view -c gsnap/$*.gsnap.unpaired_uniq.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.mtRNA.exons.bed` | perl -ne '$$sum += $$_; END {print "$$sum\t"}' >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.unpaired_mult.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.miscRNA.exons.bed`"\n"\
		     `samtools view -c gsnap/$*.gsnap.unpaired_uniq.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.miscRNA.exons.bed` | perl -ne '$$sum += $$_; END {print "$$sum\t"}' >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.unpaired_mult.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.other.exons.bed`"\n"\
		     `samtools view -c gsnap/$*.gsnap.unpaired_uniq.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.other.exons.bed` | perl -ne '$$sum += $$_; END {print "$$sum\t"}' >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.unpaired_mult.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.protein_coding.exons.bed`"\t" >> $@.part 
	echo -ne `samtools view -c gsnap/$*.gsnap.unpaired_uniq.bam -L ~/generic/data/ensembl/Homo_sapiens.GRCh37.75.protein_coding.exons.bed`"\n" >> $@.part
	mv $@.part $@ 

rseqc/allpatients.stats.txt: $(foreach S, $(SAMPLES), rseqc/$S.stats.txt)
	echo -e "sample\ttotal\tunmapped\texonic\trRNA\tlincRNA\tmiRNA\tmtRNA\tmiscRNA\tother, non-protein-coding\tprotein-coding, multi-mapped\tprotein-coding, unique\n" > $@.part
	cat $^ >> $@.part
	mv $@.part $@ 

rseqc/allpatients.geneBodyCoverage.pdf: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.geneBodyCoverage.pdf)
	gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@

rseqc/%.rseqc.geneBodyCoverage.pdf: gsnap/%.gsnap.bam ~/generic/data/rseqc/hg19_Ensembl.bed
	~/tools/RSeQC-2.3.9/bin/geneBody_coverage.py -i $< -r ~/generic/data/rseqc/hg19_Ensembl.bed -o rseqc/$*.rseqc

rseqc/allpatients.DupRate_plot.pdf: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.DupRate_plot.pdf)
	gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@

rseqc/%.rseqc.DupRate_plot.pdf: gsnap/%.gsnap.bam
	~/tools/RSeQC-2.3.9/bin/read_duplication.py -i $< -o rseqc/$*.rseqc

rseqc/allpatients.read-distribution.txt: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.read-distribution.txt)
	rm -f $@.part
	for S in $^ ; do echo $$S >> $@.part; cat $$S >> $@.part ; done
	mv $@.part $@

rseqc/%.rseqc.read-distribution.txt: gsnap/%.gsnap.bam ~/generic/data/rseqc/hg19_Ensembl.bed
	~/tools/RSeQC-2.3.9/bin/read_distribution.py -i $< -r ~/generic/data/rseqc/hg19_Ensembl.bed > $@.part
	mv $@.part $@

rseqc/allpatients.splice_junction.pdf: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.splice_junction.pdf)
	gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@

rseqc/allpatients.splicing_events.pdf: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.splice_events.pdf)
	gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@

rseqc/%.rseqc.splice_events.pdf rseqc/%.rseqc.splice_junction.pdf: gsnap/%.gsnap.bam ~/generic/data/rseqc/hg19_Ensembl.bed
	~/tools/RSeQC-2.3.9/bin/junction_annotation.py -i $< -r ~/generic/data/rseqc/hg19_Ensembl.bed -o rseqc/$*.rseqc

rseqc/allpatients.rRNA.count: $(foreach S, $(SAMPLES), rseqc/$S.rRNA.count)
	rm -f $@.part
	for S in $^ ; do echo $$S >> $@.part; cat $$S >> $@.part ; done
	mv $@.part $@

rseqc/%.rRNA.count: gsnap/%.gsnap.bam flagstat/%.samtools.flagstat
	echo "Total number of reads:" `head -1 flagstat/$*.samtools.flagstat | cut -f 1 -d ' '` > $@.part
	echo "Number of reads mapping to rRNA genes: " `samtools view $< -L ~/generic/data/ensembl/rRNA.ensembl.biomart.GRCh37.p13.bed | wc -l` >> $@.part
	mv $@.part $@
	
rseqc/allpatients.infer_experiment.txt: $(foreach S, $(SAMPLES), rseqc/$S.infer_experiment.txt)
	rm -f $@.part
	for S in $^ ; do echo $$S >> $@.part; cat $$S >> $@.part ; done
	mv $@.part $@

rseqc/%.infer_experiment.txt: gsnap/%.gsnap.bam ~/generic/data/rseqc/hg19_Ensembl.bed
	~/tools/RSeQC-2.3.9/bin/infer_experiment.py -i $< -r ~/generic/data/rseqc/hg19_Ensembl.bed > $@.part
	mv $@.part $@
	
rseqc/%.saturation.pdf: gsnap/%.gsnap.bam ~/generic/data/rseqc/hg19_Ensembl.bed
	~/tools/RSeQC-2.3.9/bin/RPKM_saturation.py -i $< -r ~/generic/data/rseqc/hg19_Ensembl.bed -o rseqc/$*
	
htseq/coverage-saturation-curve.pdf: $(foreach S, $(SAMPLES), htseq/$S.subsamples.count) ~/generic/scripts/plot_saturation_curve.R
	Rscript ~/generic/scripts/plot_saturation_curve.R --input-dir htseq --min-reads 3 --y-max 30000 --output $@.part
	mv $@.part $@ 

htseq/%.subsamples.count: gsnap/%.gsnap.bam ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.no-rRNA.gtf.gz htseq/%.count
	~/tools/samtools-0.1.19/samtools view -s 22.1 $< | ~/tools/HTSeq-0.6.1/scripts/htseq-count -f sam -t exon -s no - ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.no-rRNA.gtf.gz > $@.part ; \
	for FRACT in 2 3 4 5 6 7 8 9 ; do \
		~/tools/samtools-0.1.19/samtools view -s 22.$$FRACT $< | ~/tools/HTSeq-0.6.1/scripts/htseq-count -f sam -t exon -s no - ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.no-rRNA.gtf.gz > $@.part2 ; \
		paste $@.part <(cut -f 2 $@.part2) > $@.part3 ; mv $@.part3 $@.part ; rm $@.part2 ; \
	done
	echo -e "gene\t10%\t20%\t30%\t40%\t50%\t60%\t70%\t80%\t90%\t100%" | cat - <(paste $@.part <(cut -f 2 htseq/$*.count)) > $@
	rm $@.part
	
	