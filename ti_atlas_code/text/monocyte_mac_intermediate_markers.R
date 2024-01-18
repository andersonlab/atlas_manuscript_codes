
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

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

df1=df1[df1$label=="Monocytes",]
monocyte_markers=df1[df1$gene_symbols %in% c("FCGR3B", "FCGR1A", "CSF3R"),]
monocyte_markers
round(mean(monocyte_markers$value),2)

monocyte_markers=df1[df1$gene_symbols %in% c("FCER1G", "TYROBP"),]
monocyte_markers
round(mean(monocyte_markers$value),2)
#df2 <- df1 %>% group_by(label) %>% slice_max(value, n =100)
#df2
