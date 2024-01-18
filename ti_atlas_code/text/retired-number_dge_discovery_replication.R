#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

library(dplyr)
# f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/dge/mast_random_effect_final/disease_status-run_8-mast-cp10k_greater_1-filter_mt_ribo_ig/disease_status_merged-de_results_all.tsv"
# f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"
# f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-all_samples/disease_status_merged-de_results.tsv"
# 
# f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/disease_status_merged-de_results.tsv"
# f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"
# f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-all_samples/disease_status_merged-de_results.tsv"

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/disease_status_merged-de_results.tsv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/disease_status_merged-de_results.tsv"
f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Full/disease_status_merged-de_results.tsv"

df=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df2=read.table(f3, sep="\t", header=TRUE, stringsAsFactors = FALSE)

#df[is.na(df$category),]
#df_sig[is.na(df_sig$category),]

prepare_df <- function(df){
  
  #df=df
  df=df[!is.na(df$qvalue_bh_allcelltypes),]
  df$significant <- df$qvalue_bh_allcelltypes < 0.05
  df$significant <- factor(df$significant, levels = c(TRUE, FALSE), labels = c("True", "False"))
  df$direction_regulation <- ifelse(df$log2fc >= 0, "up-regulated", "down-regulated")
  df$direction_regulation <- factor(df$direction_regulation, levels=c("up-regulated", "down-regulated"))
  
  df_sig=df[df$significant=="True",]
  
  
  #df_sig$category <- factor(df_sig$category, levels=rev(order))
  df_sig[df_sig$category %in% c("T Cell", "B Cell", "B Cell plasma",  "Myeloid"), "compartment"] <- "Immune"
  df_sig[df_sig$category %in% c("Secretory", "Enterocyte", "Stem cells"), "compartment"] <- "Epithelial"
  df_sig[df_sig$category %in% c("Mesenchymal"), "compartment"] <- "Mesenchymal"
  df_sig=df_sig[,c("label", "qvalue_bh_allcelltypes", "category", "gene_symbol", "compartment", "direction_regulation")] #"direction_regulation"
  df3 <- df_sig %>% group_by(label, direction_regulation) %>% dplyr::mutate(count = n_distinct(gene_symbol))
  
  check <- df3 %>% group_by(label, direction_regulation) %>% arrange(count)
  ha <- check[check$direction_regulation == "up-regulated",]
  ha <- distinct(ha[c("label", "direction_regulation", "count")])
  
  missing = unique(pull(df3[!(df3$label %in% ha$label),"label"]))
  print(missing)
  df3$label <- factor(df3$label, levels=c(missing, ha$label))
  
  
  return(df3)
}

df0=prepare_df(df)
df11=prepare_df(df1)
df22=prepare_df(df2)

df0$freeze <- "Discovery"
df11$freeze <- "Replication"
df22$freeze <- "Full"

length(unique(df0$gene_symbol))
length(table(df0$label))
length(table(df11$label))

# Number of unique DEGs in discovery and replication

#all=rbind(df0, df11, df22)
all=rbind(df0, df11, df22)

length(unique(df0$gene_symbol))
length(unique(df11$gene_symbol))
length(unique(df22$gene_symbol))



cell_types=c( "Stem cell LGR5+" ,  "Stem cell MKI67+ (1)" ,  "Stem cell MKI67+ (2)" ,  "Enterocyte precursor crypt OLFM4+ KRT20++" ,
              "Enterocyte progenitor crypt OLFM4++ KRT20+ (1)",  "Enterocyte progenitor crypt OLFM4++ KRT20+ (2)" , "Enterocyte middle villus (1)" ,
              "Enterocyte middle villus (2)" , "Enterocyte top villus" , "Enterocytes BEST4"  ,  "Paneth cell"  ,  "Goblet cell crypt MKI67+" ,  "Goblet cell middle villus" ,
              "Goblet cell top villus" ,    "Endocrine cell" , "Tuft cell" ,    "Fibroblast/Myofibroblasts"  ,
              "Endothelial cell"  ,  "Pericytes",  "Smooth muscle cell" ,  "Mac resident IL10RA+"  ,  "Mac resident IL10RA-" , "Dendritic cell"  ,
              "MoMac IL10RA+", "MoMac IL10RA-" ,  "Monocytes",   "Mast"  ,  "ILC1 CD3D- NCAM1+" ,
              "ILC3 CD3D- IL23R+" ,  "T cell CD4 CD40LG+ (1)" ,
              "T cell CD4 CD40LG+ (2)" , "T cell CD4 CD40LG+ (3)" , "T cell CD4 Treg"  , "T cell CD4 naïve" ,  "T cell CD4 proliferating"  ,
              "T cell CD4- CD8-", "T cell CD8 (1)",  "T cell CD8 (2)" ,  "T cell CD8 (3)" , "T cell gd",
              "B cell" ,"B cell naïve" ,"B cell activated",
              "B cell germinal centre/plasmablasts" ,  "B cell memory (1)"  ,  "B cell memory (2)"  ,  "B cell plasma IgA CD38+"  ,
              "B cell plasma IgA CD38++" ,  "B cell plasma IgA CD38+++" )

#all$label <- factor(all$label , levels=rev(cell_types))

all$category <- factor(all$category , levels=c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal", "Myeloid", "T Cell", "B Cell", "B Cell plasma")) 

#all$freeze <- factor(all$freeze, levels=c("Discovery", "Replication", "discovery_and_Replication"), 
#                     labels=c("Discovery auto-annotated dataset", "Replication auto-annotated dataset", "Full auto-annotated dataset"))

all$freeze <- factor(all$freeze, levels=c("Discovery", "Replication", "Full"), 
                     labels=c("Discovery", "Replication", "Full cohort"))


all$direction_regulation <- factor(all$direction_regulation, levels=c("down-regulated", "up-regulated"), labels=c("down", "up"))



#####################

length(unique(all[all$freeze=="Discovery",]$gene_symbol))
length(unique(all[all$freeze=="Replication",]$gene_symbol))

length(unique(all[all$freeze %in% c("Discovery", "Discovery & Replication"),]$gene_symbol))
length(unique(all[all$freeze %in% c("Replication", "Discovery & Replication"),]$gene_symbol))


check =all[all$freeze=="Discovery",]
all$category

#length(unique(all[(all$freeze=="Discovery" & all$category=="B Cell plasma"),]$gene_symbol))
#length(unique(all[(all$freeze=="Replication" & all$category=="B Cell plasma"),]$gene_symbol))

