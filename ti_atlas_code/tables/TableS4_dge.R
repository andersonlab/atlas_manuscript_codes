#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/tables/")

library(dplyr)
#f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/dge/mast_random_effect_final/disease_status-run_8-mast-cp10k_greater_1-filter_mt_ribo_ig/disease_status_merged-de_results_all.tsv"
#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"
#f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-all_samples/disease_status_merged-de_results.tsv"

# f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/disease_status_merged-de_results.tsv"
# f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"
# f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-all_samples/disease_status_merged-de_results.tsv"

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/disease_status_merged-de_results.tsv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/disease_status_merged-de_results.tsv"
f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Full/disease_status_merged-de_results.tsv"


df=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df2=read.table(f3, sep="\t", header=TRUE, stringsAsFactors = FALSE)

df$cohort <- "Discovery"
df1$cohort <- "Replication"
df2$cohort <- "Full"

all=rbind(df, df1, df2)

all=all[!is.na(all$qvalue_bh_allcelltypes),]
all$significant <- all$qvalue_bh_allcelltypes < 0.05
all$significant <- factor(all$significant, levels = c(TRUE, FALSE), labels = c("True", "False"))

df_sig=all[all$significant=="True",]


df_sig=df_sig[,c("label","category", "gene_symbol", "log2fc", "qvalue_bh_allcelltypes", "test_statistic",  "mean_counts", "n_cells", "cohort", "significant")]
names(df_sig)[1] <- "cell_type"
names(df_sig)[2] <- "major_cell_population"
names(df_sig)[3] <- "gene_symbols"

df_sig=df_sig[df_sig$cohort!="Full",]

write.table(file="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/tables/TableS4_dge.tsv", df_sig,  sep="\t", row.names = FALSE)



