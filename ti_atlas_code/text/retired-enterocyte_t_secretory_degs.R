#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

library(dplyr)
#f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/dge/mast_random_effect_final/disease_status-run_8-mast-cp10k_greater_1-filter_mt_ribo_ig/disease_status_merged-de_results_all.tsv"
#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"
#f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-all_samples/disease_status_merged-de_results.tsv"

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/disease_status_merged-de_results.tsv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"
f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-all_samples/disease_status_merged-de_results.tsv"


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
  df_sig=df_sig[,c("label", "qvalue_bh_allcelltypes", "category", "gene_symbol", "compartment", "direction_regulation", "n_cells")] #"direction_regulation"
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

all$freeze <- factor(all$freeze, levels=c("Discovery", "Replication", "Full"), 
                     labels=c("Discovery", "Replication", "Full cohort"))

keep=c("freeze", "category",  "label", "n_cells")
all1=all[,keep]
all2=unique(all1)
x1 <- all2 %>% group_by(freeze, category) %>% mutate(sum_n_cells=sum(n_cells))
keep=c("freeze", "category","sum_n_cells")
x1=x1[,keep]
x1=unique(x1)


keep=c("category", "gene_symbol", "freeze")
all1=all[,keep]
#all2=all1[duplicated(all1),]
#all2=all1[unique(all1),]
x2 <- all1 %>% group_by(freeze, category) %>%  mutate(union_degs = n())
keep=c("freeze", "category","union_degs")
x2=x2[,keep]
x2=unique(x2)

#ord_cat=c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal", "Myeloid", "T Cell", "B Cell", "B Cell plasma")

x3 <- left_join(x1, x2, by=c("freeze", "category"))


file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/palette.tsv"
pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
pal=unique(pal[,c("category", "palette_category")])
palette=pal$palette_category
names(palette) <- pal$category


#####
haa <- x3 %>% group_by(freeze) %>% do(model = lm(sum_n_cells ~ union_degs, data = .))

dff=list()
for(i in 1:length(unique(x3$freeze))){
  dff[[i]]=sqrt(summary(haa$model[[i]])$r.squared)
}
dff=data.frame(unlist(dff))
dff$freeze <- as.character(haa$freeze)
names(dff) <- c("R", "freeze")

x3=merge(x3, dff, by = "freeze", all.x = TRUE)
x3$freeze_title <- paste0(x3$freeze,": R=", round(x3$R, 2))

#T cell, enterocyte, secretory

x5 = x3[c("freeze","category", "union_degs")]
x5=x5[x5$freeze=="Discovery",]
p1 <- x5[x5$category == 'T Cell',]$union_degs / x5[x5$category == 'Mesenchymal',]$union_degs
p2 <- x5[x5$category == 'Enterocyte',]$union_degs / x5[x5$category == 'Mesenchymal',]$union_degs
p3 <- x5[x5$category == 'Secretory',]$union_degs / x5[x5$category == 'Mesenchymal',]$union_degs
min(p1,p2,p3)


x5 = x3[c("freeze","category", "union_degs")]
x5=x5[x5$freeze=="Replication",]
p1 <- x5[x5$category == 'T Cell',]$union_degs / x5[x5$category == 'Mesenchymal',]$union_degs
p2 <- x5[x5$category == 'Enterocyte',]$union_degs / x5[x5$category == 'Mesenchymal',]$union_degs
p3 <- x5[x5$category == 'Secretory',]$union_degs / x5[x5$category == 'Mesenchymal',]$union_degs
min(p1,p2,p3)

#B plasma

x5 = x3[c("freeze","category", "union_degs", "sum_n_cells")]
x5=x5[x5$freeze=="Discovery",]
p1 <- x5[x5$category == 'B Cell plasma',]$union_degs 
p11 <- x5[x5$category == 'B Cell plasma',]$sum_n_cells
p1
p11


x5 = x3[c("freeze","category", "union_degs", "sum_n_cells")]
x5=x5[x5$freeze=="Replication",]
p1 <- x5[x5$category == 'B Cell plasma',]$union_degs 
p11 <- x5[x5$category == 'B Cell plasma',]$sum_n_cells
p1
p11





# ADD RATIO OF CELLS
file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/ratio_of_cells.csv"
ratio=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
ratio=ratio[ratio$cohort=="Discovery",]

fold_change=(ratio[ratio$disease_status=="cd",]$n_cells-ratio[ratio$disease_status=="healthy",]$n_cells)/ratio[ratio$disease_status=="healthy",]$n_cells

ratio.df=data.frame(ratio[!duplicated(ratio$celltype_label__machine),]$celltype_label__machine,fold_change)





############
#df3 <- reshape2::dcast(df1, category + label + freeze  ~ disease_status , value.var = "percentage")
#df4=df[,c("freeze", "label", "disease_status", "sanger_sample_id", "nr_cells")]
#df4=unique(df4)
#df5 <- reshape2::dcast(df4, freeze + label +  sanger_sample_id ~ disease_status, value.var = "nr_cells")
df7 <- all %>% group_by(freeze, label) %>% dplyr::summarize(fc=(n_degs[disease_status=="cd"] + n_degs[disease_status=="healthy"])/n_degs[disease_status=="healthy"])









# library(ggrepel)
# plt <- ggplot(all, aes(x=n_cells, y=n_degs, colour=category)) + geom_point() + facet_wrap(vars(freeze))#+ facet_wrap(vars(freeze))#+ annotation_custom(grob1) #+ geom_smooth(method=lm, se=FALSE) + geom_point() 
# plt <- plt + theme_bw() + xlab("Number of cells") + ylab("Number of DEG (fdr<0.05)")#Number of significant differentially expressed genes + xlab("Number of genes tested") + ylab("Number of genes tested")
# plt <- plt + labs(color="Major cell type") + scale_colour_manual(values=palette) + geom_text_repel(aes(label=label), size=3)
# plt <- plt + theme(text = element_text(size=15)) 
# plt
# 
# 
# 
# ggsave(filename="dge_power_dge.png", plot=plt, width=10, height = 8)
# 
# 
# all1=all[c("label","category", "n_cells", "freeze")]
# all1=all1[all1$freeze != "Replication",]
# #all1=all1[all1$freeze != "Discovery",]
# try=reshape2::dcast(all1, label + category ~ freeze, value.var = "n_cells")
# 
# library(ggrepel)
# plt <- ggplot(try, aes(x=Discovery, y=discovery_and_Replication, colour=category)) + geom_point() + facet_wrap(vars(freeze))#+ annotation_custom(grob1) #+ geom_smooth(method=lm, se=FALSE) + geom_point() 
# plt <- plt + theme_bw() + xlab("Discovery dataset (n cells)") + ylab("Discovery and replication dataset (n cells)")#Number of significant differentially expressed genes + xlab("Number of genes tested") + ylab("Number of genes tested")
# plt <- plt + labs(color="Major cell type") + scale_colour_manual(values=palette) + geom_text_repel(aes(label=label), size=3)
# plt <- plt + theme(text = element_text(size=15)) + scale_y_continuous(limits=c(-100, 70000), breaks=c(0, 10000,30000, 70000)) + scale_x_continuous(limits=c(-100, 70000), breaks=c(0, 10000,30000, 70000))
# plt
# ggsave(filename="dge_power1.png", plot=plt, width=10, height = 8)



# library(ggrepel)
# plt <- ggplot(try, aes(x=Replication, y=discovery_and_Replication, colour=category)) + geom_point() #+ facet_wrap(vars(freeze))#+ annotation_custom(grob1) #+ geom_smooth(method=lm, se=FALSE) + geom_point() 
# plt <- plt + theme_bw() + xlab("Replication dataset (n cells)") + ylab("Discovery and replication dataset (n cells)")#Number of significant differentially expressed genes + xlab("Number of genes tested") + ylab("Number of genes tested")
# plt <- plt + labs(color="Major cell type") + scale_colour_manual(values=palette) + geom_text_repel(aes(label=label), size=3)
# plt <- plt + theme(text = element_text(size=15)) + scale_y_continuous(limits=c(-100, 70000), breaks=c(0, 10000,30000, 70000)) + scale_x_continuous(limits=c(-100, 70000), breaks=c(0, 10000,30000, 70000))
# plt
# ggsave(filename="dge_power2.png", plot=plt, width=10, height = 8)

