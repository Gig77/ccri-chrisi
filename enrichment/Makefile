export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

#------
# GSEA
#------

gsea/diff_exp_genes.rnk: deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.tsv /mnt/projects/chrisi/scripts/get-gsea-rnk-file.R
	Rscript /mnt/projects/chrisi/scripts/get-gsea-rnk-file.R
	
gsea/all: gsea/diff_exp_genes.rnk
	java -cp ~/tools/gsea-2.0.13/gsea2-2.0.13.jar -Xmx3048m xtools.gsea.GseaPreranked \
		-rpt_label etv6runx1 \
		-rnk gsea/diff_exp_genes.rnk \
		-gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c1.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c3.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c4.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c6.all.v4.0.symbols.gmt,gseaftp.broadinstitute.org://pub/gsea/gene_sets/c7.all.v4.0.symbols.gmt \
		-collapse false -mode Max_probe -norm meandiv -nperm 100 -scoring_scheme weighted -rpt_label my_analysis -include_only_symbols true -make_sets true \
		-rnd_seed 149 \
		-plot_top_x 300 \
		-set_max 500 \
		-set_min 15 \
		-zip_report false \
		-gui false \
		-out gsea

gsea/chrisi_vs_fuka2011: gsea/diff_exp_genes.rnk
	java -cp /home/STANNANET/christian.frech/tools/gsea-2.0.13/gsea2-2.0.13.jar -Xmx3048m xtools.gsea.GseaPreranked \
		-rpt_label chrisi_vs_fuka2011 \
		-rnk /mnt/projects/chrisi/results/enrichment/gsea/diff_exp_genes.rnk \
		-gmx /mnt/projects/chrisi/results/fuka/tableS1.gmx#fuka2011_top_200,/mnt/projects/chrisi/results/fuka/tableS1.gmx#fuka2011_down_top200,/mnt/projects/chrisi/results/fuka/tableS1.gmx#fuka2011_up_top200 \
		-out /mnt/projects/chrisi/results/enrichment/gsea \
		-collapse false -mode Max_probe -norm meandiv -nperm 100 -scoring_scheme weighted -include_only_symbols true -make_sets true \
		-rnd_seed 149 \
		-plot_top_x 300 \
		-set_max 500 \
		-set_min 15 \
		-zip_report false \
		-gui false 

gsea/hsc: gsea/diff_exp_genes.rnk
	java -cp /home/STANNANET/christian.frech/tools/gsea-2.0.13/gsea2-2.0.13.jar -Xmx3048m xtools.gsea.GseaPreranked \
		-rpt_label chrisi_vs_hsc_signatures \
		-rnk /mnt/projects/chrisi/results/enrichment/gsea/diff_exp_genes.rnk \
		-gmx /mnt/projects/chrisi/results/enrichment/gsea/JAATINEN_HEMATOPOIETIC_STEM_CELL_UP.gmx,/home/STANNANET/christian.frech/tel_aml_kd/data/eppert_2011/LSC_HSC_gene_sets.gmx#EPPERT2011_HSC_R,/home/STANNANET/christian.frech/tel_aml_kd/data/eppert_2011/LSC_HSC_gene_sets.gmx#EPPERT2011_LSC_R \
		-out /mnt/projects/chrisi/results/enrichment/gsea \
		-collapse false -mode Max_probe -norm meandiv -nperm 100 -scoring_scheme weighted -include_only_symbols true -make_sets true \
		-rnd_seed 149 \
		-plot_top_x 300 \
		-set_max 500 \
		-set_min 5 \
		-zip_report false \
		-gui false
		
gsea/primTA-vs-preB: gsea/diff_exp_genes.rnk
	java -cp /home/STANNANET/christian.frech/tools/gsea-2.0.13/gsea2-2.0.13.jar -Xmx3048m xtools.gsea.GseaPreranked \
		-rpt_label telaml_ross_2003_vs_Bcells \
		-rnk /mnt/projects/chrisi/results/enrichment/gsea/diff_exp_genes.rnk \
		-gmx /mnt/projects/chrisi/results/enrichment/gsea/gmttab.TELAMLROSSvsBcells.gmt \
		-out /mnt/projects/chrisi/results/enrichment/gsea \
		-collapse false -mode Max_probe -norm meandiv -nperm 100 -scoring_scheme weighted -include_only_symbols true -make_sets true \
		-rnd_seed 149 \
		-plot_top_x 300 \
		-set_max 500 \
		-set_min 5 \
		-zip_report false \
		-gui false

#------
# DAVID
#------

david/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.up.txt: integrated-expressed-vs-notexpressed.single-factor.deseq2.tsv /mnt/projects/chrisi/scripts/enrichment/david.pl
	perl /mnt/projects/chrisi/scripts/enrichment/david.pl \
		--deseq-file $< \
		--direction up \
		--foldchange 1 \
		--min-padj 0.1 \
		> $@.part
	mv $@.part $@
	
david/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.down.txt: integrated-expressed-vs-notexpressed.single-factor.deseq2.tsv /mnt/projects/chrisi/scripts/enrichment/david.pl
	perl /mnt/projects/chrisi/scripts/enrichment/david.pl \
		--deseq-file $< \
		--direction down \
		--foldchange 1 \
		--min-padj 0.1 \
		> $@.part
	mv $@.part $@

david/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.down.txt.heatmap.numRow.100.largeCatexcl.1e+05.pdf: david/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.up.txt david/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.down.txt /mnt/projects/chrisi/scripts/enrichment/gen-david-heatmap.R
	Rscript /mnt/projects/chrisi/scripts/enrichment/gen-david-heatmap.R
	mv david/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.up.txt.heatmap.numRow.100.largeCatexcl.1e+05.pdf.part david/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.up.txt.heatmap.numRow.100.largeCatexcl.1e+05.pdf
	mv david/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.down.txt.heatmap.numRow.100.largeCatexcl.1e+05.pdf.part david/integrated-expressed-vs-notexpressed.single-factor.deseq2.david.down.txt.heatmap.numRow.100.largeCatexcl.1e+05.pdf
	 