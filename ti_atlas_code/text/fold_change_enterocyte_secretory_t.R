#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

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

check <- x1 %>% group_by(freeze) %>% dplyr::summarize(sum_n_cells_freeze=sum(sum_n_cells))
check

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

x4=x3[!(x3$category %in% c("B Cell plasma")),]
haa <- x4 %>% group_by(freeze) %>% do(model = lm(sum_n_cells ~ union_degs, data = .))

dff=list()
for(i in 1:length(unique(x3$freeze))){
  dff[[i]]=sqrt(summary(haa$model[[i]])$r.squared)
}
dff=data.frame(unlist(dff))
dff$freeze <- as.character(haa$freeze)
names(dff) <- c("R", "freeze")

x3=merge(x3, dff, by = "freeze", all.x = TRUE)
x3$freeze_title <- paste0(x3$freeze,": R=", round(x3$R, 2))


#make order of freeze
#ll=unique(x3[,c("R", "freeze_title")])
#x3$freeze_title <- factor(x3$freeze_title , levels=ll[order(ll$R),]$freeze_title)
ll=unique(x3[,c("freeze_title","freeze")])
ll <- ll %>% arrange(factor(freeze, levels = c("Discovery","Replication", "Full cohort"))) 
x3$freeze_title <- factor(x3$freeze_title , levels=ll$freeze_title)

# Add ratio information
f1='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/ratio_n_cells.tsv'
oo=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
x3=merge(x3, oo, by = c("freeze", "category"), all.x = TRUE)


#Get data to repel labels
all5=unique(x3[,c("freeze","freeze_title", "category", "union_degs", "sum_n_cells", "fc_sum_n_cells")])
my_order <- all5 %>% group_by(freeze) %>% top_n(8, union_degs)
cell_types_label = my_order[my_order$freeze == "Replication",]$category
all_repel=all5[all5$category %in% cell_types_label,]

# Number of DEGS B plasma
x3[(x3$category=="B Cell plasma") & (x3$freeze=="Discovery"),]$union_degs
x3[(x3$category=="B Cell plasma") & (x3$freeze=="Replication"),]$union_degs

x3[(x3$category=="B Cell plasma") & (x3$freeze=="Discovery"),]$sum_n_cells
x3[(x3$category=="B Cell plasma") & (x3$freeze=="Replication"),]$sum_n_cells

# Fold change of DEGs epithelial vs B plasma and mesenchymal
all1[all1$category %in% c("Enterocyte", "Secretory", "T Cell"),"compartment"] <- "Enterocyte_Secretory_T"
#all1[all1$category %in% c("B Cell plasma", "Mesenchymal"),"compartment"] <- "Plasma_Mesenchymal"
all1[all1$category %in% c("B Cell plasma"),"compartment"] <- "Plasma"

#all6=all1[all1$compartment %in% c("Enterocyte_Secretory_T", "Plasma_Mesenchymal"),]
all6=all1[all1$compartment %in% c("Enterocyte_Secretory_T", "Plasma"),]

x2 <- all6 %>% group_by(freeze, compartment) %>%  mutate(union_degs = n())
keep=c("freeze", "compartment","union_degs")
x2=x2[,keep]
x2=unique(x2)

#x2[x2$freeze=="Discovery" & x2$compartment=="Enterocyte_Secretory_T",]$union_degs/x2[x2$freeze=="Discovery" & x2$compartment=="Plasma_Mesenchymal",]$union_degs
#x2[x2$freeze=="Replication" & x2$compartment=="Enterocyte_Secretory_T",]$union_degs/x2[x2$freeze=="Replication" & x2$compartment=="Plasma_Mesenchymal",]$union_degs

fc_d=x2[x2$freeze=="Discovery" & x2$compartment=="Enterocyte_Secretory_T",]$union_degs/x2[x2$freeze=="Discovery" & x2$compartment=="Plasma",]$union_degs
fc_r=x2[x2$freeze=="Replication" & x2$compartment=="Enterocyte_Secretory_T",]$union_degs/x2[x2$freeze=="Replication" & x2$compartment=="Plasma",]$union_degs

fc_d,fc_r
