#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")


library(dplyr)

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/disease_status_merged-de_results.tsv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"


df=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep="\t", header=TRUE, stringsAsFactors = FALSE)

prepare_df <- function(df){
  
  df$significant <- df$qvalue_bh_allcelltypes < 0.05
  df$significant <- factor(df$significant, levels = c(TRUE, FALSE), labels = c("True", "False"))
  df_sig=df[df$significant=="True",]
  df_sig$direction_regulation <- ifelse(df_sig$log2fc >= 0, "up-regulated", "down-regulated")
  df1=df_sig[,c("gene_symbol", "category", "label", "qvalue_bh_allcelltypes")]
  df2 <- df1 %>% group_by(label) %>%  dplyr::summarize(n_degs = n())
  df1=df[,c("label", "gene_symbol")]
  df3 <- df1 %>% group_by(label) %>%  dplyr::summarize(n_genes_tested = n())
  df1=distinct(df[,c("label", "n_cells")])
  df4=merge(df2, df3, by="label")
  df5=merge(df4, df1, by="label")
  df1=distinct(df[,c("label", "category")])
  df6=merge(df5, df1, by="label")
}


df0=prepare_df(df)
df11=prepare_df(df1)

df0$freeze <- "Discovery"
df11$freeze <- "Replication"

all=rbind(df0, df11)


# library(ggrepel)
# plt <- ggplot(all, aes(x=n_cells, y=n_degs, colour=category)) + geom_point() + facet_wrap(vars(freeze))#+ annotation_custom(grob1) #+ geom_smooth(method=lm, se=FALSE) + geom_point() 
# plt <- plt + theme_bw() + xlab("Number of cells") + ylab("Number of DEG (fdr<0.05)")#Number of significant differentially expressed genes + xlab("Number of genes tested") + ylab("Number of genes tested")
# plt <- plt + labs(color="Major cell type") + scale_colour_manual(values=palette) + geom_text_repel(aes(label=label), size=3)
# plt <- plt + theme(text = element_text(size=15)) 
# plt
# ggsave(filename="dge_power.png", plot=plt, width=15, height = 6)



#all$freeze <-  factor(all$freeze, levels=c("Discovery", "Replication", "discovery_and_Replication"),
#                      labels=c("Discovery dataset", "Replication dataset", "Discovery and replication"))

all1=all[c("label","category", "n_degs", "freeze")]
all1=all1[all1$freeze != "Replication",]
#all1=all1[all1$freeze != "Discovery",]
try=reshape2::dcast(all1, label + category ~ freeze, value.var = "n_degs")


all$freeze <- factor(all$freeze, levels=c("Discovery", "Replication"),
                     labels=c("Discovery", "Replication"))

all$frac_dge = (all$n_degs / all$n_genes_tested) * 100
all$category <- factor(all$category , levels=c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal", "Myeloid", "T Cell", "B Cell", "B Cell plasma"))

all[all$category %in% c("T Cell", "B Cell", "B Cell plasma",  "Myeloid"), "compartment"] <- "Immune"
all[all$category %in% c("Secretory", "Enterocyte", "Stem cells"), "compartment"] <- "Epithelial"
all[all$category %in% c("Mesenchymal"), "compartment"] <- "Mesenchymal"

all1=ddply(all, .(freeze, compartment), summarize, Mean=mean(n_degs))
all1

frz="Discovery"
all1_mean=all1[all1$freeze==frz,]

Y=all1_mean[all1_mean$compartment=="Epithelial",]$Mean 
X1=all1_mean[all1_mean$compartment=="Immune",]$Mean
X2=all1_mean[all1_mean$compartment=="Mesenchymal",]$Mean

fold_change_epi_immune=(Y-X1)/X1
fold_change_epi_mesch=(Y-X2)/X2
fold_change_epi_immune
fold_change_epi_mesch

frz="Replication"
all1_mean=all1[all1$freeze==frz,]

Y=all1_mean[all1_mean$compartment=="Epithelial",]$Mean 
X1=all1_mean[all1_mean$compartment=="Immune",]$Mean
X2=all1_mean[all1_mean$compartment=="Mesenchymal",]$Mean

fold_change_epi_immune=(Y-X1)/X1
fold_change_epi_mesch=(Y-X2)/X2
fold_change_epi_immune
fold_change_epi_mesch



