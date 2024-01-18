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

# Number of unique DEGs in discovery and replication

#all=rbind(df0, df11, df22)
all=rbind(df0, df11, df22)

length(unique(df0$gene_symbol))
length(unique(df11$gene_symbol))
length(unique(df22$gene_symbol))


table(df0$label)
table(df11$label)
table(df22$label)

# 
# cell_types=c( "Stem cell LGR5+" ,  "Stem cell MKI67+ (1)" ,  "Stem cell MKI67+ (2)" ,  "Enterocyte precursor crypt OLFM4+ KRT20++" ,
#               "Enterocyte progenitor crypt OLFM4++ KRT20+ (1)",  "Enterocyte progenitor crypt OLFM4++ KRT20+ (2)" , "Enterocyte middle villus (1)" ,
#               "Enterocyte middle villus (2)" , "Enterocyte top villus" , "Enterocytes BEST4"  ,  "Paneth cell"  ,  "Goblet cell crypt MKI67+" ,  "Goblet cell middle villus" ,
#               "Goblet cell top villus" ,    "Endocrine cell" , "Tuft cell" ,    "Fibroblast/Myofibroblasts"  ,
#               "Endothelial cell"  ,  "Pericytes",  "Smooth muscle cell" ,  "Mac resident IL10RA+"  ,  "Mac resident IL10RA-" , "Dendritic cell"  ,
#               "MoMac IL10RA+", "MoMac IL10RA-" ,  "Monocytes",   "Mast"  ,  "ILC1 CD3D- NCAM1+" ,
#               "ILC3 CD3D- IL23R+" ,  "T cell CD4 CD40LG+ (1)" ,
#               "T cell CD4 CD40LG+ (2)" , "T cell CD4 CD40LG+ (3)" , "T cell CD4 Treg"  , "T cell CD4 naïve" ,  "T cell CD4 proliferating"  ,
#               "T cell CD4- CD8-", "T cell CD8 (1)",  "T cell CD8 (2)" ,  "T cell CD8 (3)" , "T cell gd",
#               "B cell" ,"B cell naïve" ,"B cell activated",
#               "B cell germinal centre/plasmablasts" ,  "B cell memory (1)"  ,  "B cell memory (2)"  ,  "B cell plasma IgA CD38+"  ,
#               "B cell plasma IgA CD38++" ,  "B cell plasma IgA CD38+++" )
# 
# #all$label <- factor(all$label , levels=rev(cell_types))

all$category <- factor(all$category , levels=c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal", "Myeloid", "T Cell", "B Cell", "B Cell plasma")) 

#all$freeze <- factor(all$freeze, levels=c("Discovery", "Replication", "discovery_and_Replication"), 
#                     labels=c("Discovery auto-annotated dataset", "Replication auto-annotated dataset", "Full auto-annotated dataset"))

all$freeze <- factor(all$freeze, levels=c("Discovery", "Replication", "Full"), 
                     labels=c("Discovery", "Replication", "Full cohort"))


all$direction_regulation <- factor(all$direction_regulation, levels=c("down-regulated", "up-regulated"), labels=c("down", "up"))

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
all$label <- factor(all$label , levels=rev(oo$label))   

palette=c("skyblue3", "tomato")
plt <- ggplot(all, aes(x=label, fill=direction_regulation)) + geom_bar(stat="count", position = position_stack(vjust = .5)) 
plt <- plt + facet_grid(category ~ freeze, scales = "free", space = "free") + scale_fill_manual(values=palette)
plt <- plt + coord_flip() + ylab("Number of dysregulated genes") + labs(fill="Gene regulation") + ggplot2::theme_bw()
plt <- plt + theme(legend.title = element_blank(), 
                   axis.title.y=element_blank(), 
                   axis.title.x = element_text(size=15), 
                   strip.text.y = element_text(angle=360, size=12, hjust=0),
                   strip.text.x = element_text(size = 15), 
                   legend.text = element_text(size=15),
                   strip.background = element_blank()) 
plt <- plt + scale_y_continuous(limits=c(0, 700), breaks=c(0, 100,300, 400, 700))
plt

setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")
#ggsave(plt, filename=paste("dge_total_number.png"), width=12, height =10, dpi = 300)
#ggsave(plt, filename=paste("dge_total_number.pdf"), width=10, height =10, device = cairo_pdf)
ggsave(plot=plt, filename="dge_total_number.png", width = 10,height = 10) 
ggsave(plot=plt, filename="dge_total_number.pdf", width = 10,height = 10, device = cairo_pdf, dpi=200)

















#####################

length(unique(all[all$freeze=="Discovery",]$gene_symbol))
length(unique(all[all$freeze=="Replication",]$gene_symbol))









keep=c("category", "gene_symbol", "direction_regulation", "compartment", "freeze")
all1=all[,keep]
all2=all1[duplicated(all1),]

ord_cat=c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal", "Myeloid", "T Cell", "B Cell", "B Cell plasma")
all2$category <- factor(all2$category , levels=rev(ord_cat)) 

plt <- ggplot(all2, aes(x=label, fill=direction_regulation)) + geom_bar(stat="count", position = position_stack(vjust = .5)) 
plt <- plt + facet_grid(compartment ~ freeze, scales = "free", space = "free") + scale_fill_manual(values=palette)
plt <- plt + coord_flip() + ylab("Number of dysregulated genes") + labs(fill="Gene regulation") + ggplot2::theme_bw()
plt <- plt + theme(legend.title = element_blank(), 
                   axis.title.y=element_blank(), 
                   axis.title.x = element_text(size=15), 
                   strip.text.y = element_text(angle=360, size=13, hjust=0),
                   strip.text.x = element_text(size = 15), 
                   legend.text = element_text(size=15),
                   axis.text.y=element_text(size=15),
                   axis.text.x=element_text(size=10),
                   strip.background = element_blank()) 
plt <- plt + scale_y_continuous(limits=c(0, 1200), breaks=c(0, 300,600, 1200))
plt

ggsave(plot=plt, filename="dge_total_number.png", width = 9,height = 5) 















#plt <- plt + geom_text_repel(aes(label=label), size=3)
#plt

# library(ggplot2)
# 
# g <- ggplot_gtable(ggplot_build(plt))
# stripr <- which(grepl('strip-r', g$layout$name))
# fills <-c("#B5BD61", #Stem
#           "#E377C2", #Enterocyte
#           "#279E68", #Secretory
#           "#AEC7E8", #Tufn
#           "#AA40FC", #Mesenchymal
#           "#8C564B", #Myeloid
#           "#D62728", #mast
#           "#17BECF", #T cell
#           "#1F77B4", #B cell
#           "#FF7F0E" #plasma
# ) # "#17BECF",
# #"#AEC7E8"
# 
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.draw(g)

#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

#ggsave(plt, filename=paste("dge_total_number.png"), width=12, height =10, dpi = 300)
ggsave(plt, filename=paste("dge_total_number.pdf"), width=12, height =10, device = cairo_pdf)


