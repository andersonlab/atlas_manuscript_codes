
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/cluster_markers.method=wilcoxon/adata-normalized_pca-bbknn-umap-clustered-cluster_markers.tsv.gz"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/cluster_markers.method=cellex/adata-normalized_pca-bbknn-umap-clustered-esmu-top_markers.tsv.gz"

#df=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep="\t", header=TRUE, stringsAsFactors = FALSE)

#df$method <- "wilcoxon"
df1$method <- "cellex"

#df=df[c("gene_symbols", "cluster","method", "logfoldchanges", "pvals_adj","scores")]
df1=df1[c("gene_symbols", "cluster","method", "value")]

#table(df$pvals_adj < 0.05)
f3='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
df3=read.table(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)
keepme=c("label", "cluster", "category")
df3=df3[,keepme]
df1=merge(df1, df3, by = "cluster", all.x = TRUE)
df1$cluster <- NULL
table(is.na(df1$value))

names(df1)[names(df1) == "value"] <- "specificity"
names(df1)[names(df1) == "label"] <- "cell_type"
names(df1)[names(df1) == "category"] <- "major_cell_population"


# #Add average expression
# l1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/tables/mean_log1p_cp10k.csv"
# mean_counts=read.csv(l1, sep=",", header=TRUE, stringsAsFactors = FALSE)
# df2 = df1 %>% left_join(mean_counts, by=c("gene_symbols","cell_type"))

write.table(df1, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/tables/TableS3_seg.tsv", sep="\t", row.names=FALSE)



#df2 <- df1 %>% group_by(label) %>% slice_max(value, n =100)
#df2
