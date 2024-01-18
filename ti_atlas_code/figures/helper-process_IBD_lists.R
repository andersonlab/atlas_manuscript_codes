
library(stringr)


#p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/IBD_genes_Carl.csv"
#p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/Monogenic-IBD-list-ESPGHAN-2020-2.csv"

p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/ICDA-IBD-GoldStandard-IBDGenes-CAA.csv"
p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/list_monogenic_ibd_list1.tsv"

#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
library(data.table)
df1=fread(p1, sep=",", header=TRUE, stringsAsFactors = FALSE)
df4=read.table(p4, sep="\t", header=TRUE, stringsAsFactors = FALSE)

df1=df1[!is.na(df1$confidence),]
df1=df1[df1$confidence == 1000,]
df1=data.frame("gene_symbol"=df1$gene)
df1$list <- "IBD_genes"

df4=data.frame("gene_symbol"=df4$Gene)
df4$list <- "Monogenic_IBD_genes"

df=rbind(df1, df4) 

write.table(df, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/IBD_genes_and_Monogenic.tsv", sep="\t")




