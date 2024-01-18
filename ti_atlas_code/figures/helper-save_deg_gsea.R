##### DGE

#p="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/"
#p="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/"
#p="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-all_samples/"


p="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/"
p="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/"
p="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Full/"

deg="disease_status_merged-de_results.tsv.gz"
deg="disease_status_merged-gsea_results.tsv.gz"
#deg="disease_status_merged-de_results.tsv"

choose.data=deg

f1=paste0(p, choose.data)
df1=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)

df1=df1[df1$de_method=="mast::singlecell::glmer",]
names(df1)[names(df1) %in% "cell_label"] <- "label__machine_retired"

f1='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
df=read.table(f1, sep=",", header=TRUE, stringsAsFactors = FALSE)

df1=merge(df, df1, by = "label__machine_retired", all.x = TRUE)

#names(df1)[2:4] <- c("celltype_category", "celltype_label", "celltype_category_retired")
#names(df2)[2:4] <- c("major_cell_type", "cell_label", "cell_label_retired")

write.table(df1, file = paste0(p, gsub(".gz", "", choose.data)), sep="\t", row.names = FALSE)

