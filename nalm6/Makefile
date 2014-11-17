export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=~/chrisi

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
		18375_GCCACA_C4E7NACXX_8_20140603B_20140603

all: allsamples.filtered-variants.tsv

#-----------	
# DEDUPLICATE AND MERGE BAMS
#-----------

%.rmdup.bam: ../gsnap/%.bam
	~/tools/samtools-0.1.19/samtools rmdup -s $< $@.part
	mv $@.part $@
	
allsamples.gsnap.filtered.rmdup.bam: $(foreach S, $(SAMPLES), $S.gsnap.filtered.rmdup.bam)
	~/tools/samtools-0.1.19/samtools merge -r $@.part $^
	mv $@.part $@
	~/tools/samtools-0.1.19/samtools index $@

#-----------	
# VARSCAN
#-----------

allsamples.varscan.vcf: allsamples.gsnap.filtered.rmdup.bam allsamples.gsnap.filtered.rmdup.bam.bai
	# -q 10 : remove non-uniquely mapping reads and reads with low mapping quality
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns \
		<(~/tools/samtools-0.1.19/samtools view -u -q 10 $(word 1,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/broad/human_g1k_v37.fasta -) \
		--min-coverage 2 \
		--min-reads2 2 \
		--min-avg-qual 20 \
		--min-var-freq 0.2 \
		--min-freq-for-hom 0.75 \
		--p-value 1 \
		--strand-filter 0 \
		--output-vcf 1 \
		--variants 1 \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

#-----------	
# SNPEFF
#-----------

allsamples.varscan.dbsnp.vcf: allsamples.varscan.vcf ~/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar SnpSift.jar annotate \
		-v ~/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf \
		<(cat $(PWD)/$< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

allsamples.varscan.dbsnp.snpeff.vcf: allsamples.varscan.dbsnp.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx4g -jar snpEff.jar -v -lof GRCh37.75 -stats $(PWD)/allsamples.varscan.dbsnp.snpeff.summary.html $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

allsamples.varscan.dbsnp.snpeff.dbNSFP.vcf: allsamples.varscan.dbsnp.snpeff.vcf ~/generic/data/dbNSFP-2.6/dbNSFP2.6_variant.tsv.gz
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar SnpSift.jar dbnsfp \
		-v ~/generic/data/dbNSFP-2.6/dbNSFP2.6_variant.tsv.gz \
		-collapse \
		-f SIFT_pred,SIFT_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,SiPhy_29way_logOdds,LRT_pred,LRT_score,MutationTaster_pred,MutationTaster_score,MutationAssessor_pred,MutationAssessor_score,FATHMM_pred,FATHMM_score,RadialSVM_pred,RadialSVM_score,GERP++_RS,1000Gp1_AF,1000Gp1_AFR_AF,1000Gp1_EUR_AF,1000Gp1_AMR_AF,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF,Uniprot_acc,Interpro_domain, \
		$(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

#-----------	
# FILTERED VARIANTS LIST
#-----------	

allsamples.filtered-variants.tsv: allsamples.varscan.dbsnp.snpeff.dbNSFP.vcf ~/chrisi/scripts/nalm6/filter-variants.pl
	perl ~/chrisi/scripts/nalm6/filter-variants.pl \
		$< \
		--patient nalm6 \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro ~/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--remission-variants-file ~/hdall/results/remission-variants.tsv.gz \
		--cosmic-mutation-file ~/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv \
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		--clinvar-file ~/generic/data/clinvar/clinvar_20140929.vcf.gz \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@
	
