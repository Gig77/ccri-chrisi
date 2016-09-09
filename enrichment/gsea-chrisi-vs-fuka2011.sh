cd /mnt/projects/chrisi/
java -cp /data_synology/software/gsea-2.0.13/gsea2-2.0.13.jar -Xmx3048m xtools.gsea.GseaPreranked \
	-rpt_label chrisi_vs_fuka2011 \
	-rnk /mnt/projects/chrisi/results/gsea/diff_exp_genes.rnk \
	-gmx /mnt/projects/chrisi/results/fuka/tableS1.gmx#fuka2011_top_200,/mnt/projects/chrisi/results/fuka/tableS1.gmx#fuka2011_down_top200,/mnt/projects/chrisi/results/fuka/tableS1.gmx#fuka2011_up_top200 \
	-out /mnt/projects/chrisi/results/gsea \
	-collapse false -mode Max_probe -norm meandiv -nperm 100 -scoring_scheme weighted -include_only_symbols true -make_sets true \
	-rnd_seed 149 \
	-plot_top_x 300 \
	-set_max 500 \
	-set_min 15 \
	-zip_report false \
	-gui false 
