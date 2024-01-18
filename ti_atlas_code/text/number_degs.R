library(dplyr)
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

#f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/disease_status_merged-de_results.tsv"
#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"
#f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-all_samples/disease_status_merged-de_results.tsv"

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/disease_status_merged-de_results.tsv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/disease_status_merged-de_results.tsv"


df=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df$freeze <- "Discovery"
df1$freeze <- "Replication"
df$comparison <- "Discovery"
df1$comparison <- "Replication"


df_final=rbind(df, df1)

plot(sort(df_final$test_statistic[df_final$test_statistic>-50]))
summary(df_final$test_statistic[df_final$test_statistic>-50])

# df_final$hit_discovery <- df_final$qvalue_bh_allcelltypes < 0.05 & df_final$freeze == "Discovery"
# df_final$hit_replication <- df_final$qvalue_bh_allcelltypes < 0.05 & df_final$freeze == "Replication"

df_final=df_final[!is.na(df_final$qvalue_bh_allcelltypes),]

#df_final[df_final$qvalue_bh_allcelltypes < 0.05 & df_final$freeze == "Discovery", "hit"] <- "hit_discovery"
#df_final[df_final$qvalue_bh_allcelltypes < 0.05 & df_final$freeze == "Replication", "hit"] <- "hit_replication"


#keep=c("gene_symbol", "label","category", "freeze", "log2fc", "hit_discovery", "hit_replication")
#df_final=df_final[,keep]
# 
# df_final_melted= reshape2::melt(df_final, id.vars=c("gene_symbol", "label", "category", "freeze", "log2fc"), measure.vars=c("hit_discovery", "hit_replication" ))
# 
# df_final_melted$hit_discovery <- df_final_melted$variable == "hit_discovery" & df_final_melted$value == TRUE 
# df_final_melted$hit_replication <- df_final_melted$variable == "hit_replication" & df_final_melted$value == TRUE 
# 
# df_final_melted$value <- NULL
# df_final_melted$variable <- NULL
#df_final_melted[df_final_melted$value == FALSE, "variable"] <- NA

df_plt <- reshape2::dcast(df_final, gene_symbol + label + category ~ freeze , value.var = "log2fc")
df_plt1 <- reshape2::dcast(df_final, gene_symbol + label + category ~ freeze , value.var = "qvalue_bh_allcelltypes")
df_plt2 <- reshape2::dcast(df_final, gene_symbol + label + category ~ freeze , value.var = "n_cells")

df_plt$hit_discovery <- df_plt1$Discovery < 0.05
df_plt$hit_replication <- df_plt1$Replication < 0.05
df_plt$n_cells_discovery <- df_plt2$Discovery 
df_plt$n_cells_replication <- df_plt2$Replication 

#for(i in df_plt$label){
#  hh=df_plt[df_plt$label==i,]
#  print(length(unique(hh$gene_symbol)))
#}

# Genes in common
dis=df_plt[!is.na(df_plt$hit_both),]
dis=dis[dis$hit_both==TRUE,]
dis=dis[!duplicated(dis$gene_symbol),]
length(unique(dis[,"gene_symbol"]))


df_plt$hit_both <- (df_plt$hit_discovery & df_plt$hit_replication)
df_plt$hit_type <- "Neither"
df_plt$hit_type[(df_plt$hit_discovery & !df_plt$hit_replication)] <- "Discovery"
df_plt$hit_type[(!df_plt$hit_discovery & df_plt$hit_replication)] <- "Replication"
df_plt$hit_type[(df_plt$hit_discovery & df_plt$hit_replication)] <- "Discovery & Replication"
df_plt$temp_asdfs_v1 <- df_plt$Discovery
df_plt$temp_asdfs_v2 <- df_plt$Replication
# drop points with missing data
df_plt <- subset(df_plt, !is.na(hit_both))

df_plt$category <- factor(df_plt$category , levels=c("Stem cells", "Enterocyte", "Secretory","B Cell plasma",  "B Cell", "Myeloid","T Cell", "Mesenchymal")) 
df_plt$hit_type <- factor(df_plt$hit_type, levels=c("Discovery & Replication", "Discovery", "Replication", "Neither"))

df_plt=df_plt[df_plt$hit_type!="Neither",]

check <- df_plt[,c("gene_symbol", "label", "hit_type", "Discovery", "Replication")]
check <- check[check$hit_type %in% c("Discovery", "Replication", "Discovery & Replication"),]


length(unique(check[check$hit_type %in% c("Discovery", "Discovery & Replication"),]$gene_symbol))
length(unique(check[check$hit_type %in% c("Replication", "Discovery & Replication"),]$gene_symbol))


# Non unique
length(check[check$hit_type %in% c("Discovery", "Discovery & Replication"),]$gene_symbol)
length(check[check$hit_type %in% c("Replication", "Discovery & Replication"),]$gene_symbol)
