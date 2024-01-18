#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

library(dplyr)
library(ggplot2)
library(plyr)
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

write.table(x3, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/dge_power.tsv",
            sep="\t", row.names = FALSE)


file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/palette.tsv"
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
f1='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/ratio_n_cells.tsv'
oo=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
x3=merge(x3, oo, by = c("freeze", "category"), all.x = TRUE)


#Get data to repel labels
all5=unique(x3[,c("freeze","freeze_title", "category", "union_degs", "sum_n_cells", "fc_sum_n_cells")])
my_order <- all5 %>% group_by(freeze) %>% top_n(8, union_degs)
cell_types_label = my_order[my_order$freeze == "Replication",]$category
all_repel=all5[all5$category %in% cell_types_label,]


ppalette <- as.character(rev(c("navy", "steelblue", "steelblue", "steelblue", "navy", "steelblue", "steelblue", "steelblue")))
library(paletteer)
paletter=paletteer_c("ggthemes::Red-Blue Diverging", 10)
#paletter=paletteer_c("ggthemes::Red-Blue-Biscuit Diverging", 50)
paletter=as.vector(paletter)
min_scale=min(x3$fc_sum_n_cells)
max_scale=max(x3$fc_sum_n_cells)

library(ggrepel)
plt <- ggplot(x3, aes(x=sum_n_cells, y=union_degs, colour=fc_sum_n_cells)) + geom_point(size=3) + facet_grid(.~freeze_title)#+ facet_wrap(vars(freeze))#+ annotation_custom(grob1) #+ geom_smooth(method=lm, se=FALSE) + geom_point() 
plt <- plt + theme_bw() + xlab("Total number of cells") + ylab("Total number of dysregulated genes")#Number of significant differentially expressed genes + xlab("Number of genes tested") + ylab("Number of genes tested")
plt <- plt + geom_text_repel(data=all_repel, aes(label=category), size=4)  #+ scale_colour_manual(values=palette) 
plt <- plt + theme(text = element_text(size=15), 
                   strip.background = element_blank(),
                   axis.title.y=element_text(size=17), 
                   axis.title.x = element_text(size=17),
                   axis.text.x = element_text(size=12, angle=45, vjust=0.5),
                   axis.text.y = element_text(size=12),
                   strip.text.x = element_text(size = 17),
                   legend.text = element_text(size=15) #,legend.position = "none"
) + scale_x_continuous(labels = scales::comma)
#plt <- plt + scale_colour_gradient(low = "#357fb9", high = "#ca4b41") 
plt <- plt + scale_colour_gradientn(colours = rev(paletter), na.value = "white", limit = c(min_scale,max_scale))
plt <- plt + guides(size = "none") + labs(colour = "Ratio")
plt

ggsave(plot=plt, filename="dge_power.png", width = 12,height = 5) 
ggsave(plot=plt, filename="dge_power.pdf", width = 12,height = 5, device = cairo_pdf, dpi=200)








##########
# Power depends on ratio ?
x3[x3$fc_sum_n_cells >2, "balanced_cd_h_ratio"] <- "yes"
x3[x3$fc_sum_n_cells <=2, "balanced_cd_h_ratio"] <- "no"

all_repel[all_repel$fc_sum_n_cells >2, "balanced_cd_h_ratio"] <- "yes"
all_repel[all_repel$fc_sum_n_cells <=2, "balanced_cd_h_ratio"] <- "no"

library(ggrepel)
plt <- ggplot(x3, aes(x=fc_sum_n_cells, y=union_degs, colour=balanced_cd_h_ratio)) + geom_point(size=3) + facet_grid(.~freeze_title)#+ facet_wrap(vars(freeze))#+ annotation_custom(grob1) #+ geom_smooth(method=lm, se=FALSE) + geom_point() 
plt <- plt + theme_bw() + xlab("Ratio CD vs Healthy") + ylab("Total number of dysregulated genes")#Number of significant differentially expressed genes + xlab("Number of genes tested") + ylab("Number of genes tested")
plt <- plt + geom_text_repel(data=all_repel, aes(label=category), size=4)  #+ scale_colour_manual(values=palette) 
plt <- plt + theme(text = element_text(size=15), 
                   strip.background = element_blank(),
                   axis.title.y=element_text(size=17), 
                   axis.title.x = element_text(size=17),
                   axis.text.x = element_text(size=12, angle=45, vjust=0.5),
                   axis.text.y = element_text(size=12),
                   strip.text.x = element_text(size = 17),
                   legend.text = element_text(size=15) #,legend.position = "none"
) + scale_x_continuous(labels = scales::comma)
#plt <- plt + scale_colour_gradient(low = "#357fb9", high = "#ca4b41") 
#plt <- plt + scale_colour_gradientn(colours = rev(paletter), na.value = "white", limit = c(min_scale,max_scale))
plt <- plt + guides(size = "none") + labs(colour = "Unbalanced ratio")
plt

#ggsave(plot=plt, filename="dge_power.png", width = 12,height = 5) 



#x3$freeze_title <- factor(x3$freeze_title, levels=c("Discovery: R=0.93", "Replication: R=0.86", "Full cohort: R=0.94"))
#ord_cat = c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal", "Myeloid", "T Cell", "B Cell", "B Cell plasma")




library(ggrepel)
plt <- ggplot(x3, aes(x=sum_n_cells, y=union_degs, colour=fc_sum_n_cells)) + geom_point(size=3) + facet_grid(.~freeze_title)#+ facet_wrap(vars(freeze))#+ annotation_custom(grob1) #+ geom_smooth(method=lm, se=FALSE) + geom_point() 
plt <- plt + theme_bw() + xlab("Total number of cells") + ylab("Total number of dysregulated genes")#Number of significant differentially expressed genes + xlab("Number of genes tested") + ylab("Number of genes tested")
plt <- plt + geom_text_repel(data=all_repel, aes(label=category), size=4) #+ scale_colour_manual(values=palette) 
plt <- plt + theme(text = element_text(size=15), 
                   strip.background = element_blank(),
                   axis.title.y=element_text(size=17), 
                   axis.title.x = element_text(size=17),
                   axis.text.x = element_text(size=12, angle=45, vjust=0.5),
                   axis.text.y = element_text(size=12),
                   strip.text.x = element_text(size = 17),
                   legend.text = element_text(size=15),
                   legend.position = "none") + scale_x_continuous(labels = scales::comma)
plt <- plt #+ scale_y_continuous(limits=c(0, 1500), breaks=c(0, 500,1000, 1500))
#plt <- plt + ggplot2::geom_smooth(method = "lm",se = FALSE,color = "black",linetype = 'longdash',alpha = 0.5)
plt


library(ggrepel)
plt <- ggplot(x3, aes(x=sum_n_cells, y=union_degs, colour=category)) + geom_point(aes(size=fc_sum_n_cells)) + facet_grid(.~freeze_title)#+ facet_wrap(vars(freeze))#+ annotation_custom(grob1) #+ geom_smooth(method=lm, se=FALSE) + geom_point() 
plt <- plt + theme_bw() + xlab("Total number of cells") + ylab("Total number of dysregulated genes")#Number of significant differentially expressed genes + xlab("Number of genes tested") + ylab("Number of genes tested")
plt <- plt + geom_text_repel(data=all_repel, aes(label=category), size=4) #+ scale_colour_manual(values=palette) 
plt <- plt + theme(text = element_text(size=15), 
                   strip.background = element_blank(),
                   axis.title.y=element_text(size=17), 
                   axis.title.x = element_text(size=17),
                   axis.text.x = element_text(size=12, angle=45, vjust=0.5),
                   axis.text.y = element_text(size=12),
                   strip.text.x = element_text(size = 17),
                   legend.text = element_text(size=15),
                   legend.position = "none") + scale_x_continuous(labels = scales::comma)
plt <- plt #+ scale_y_continuous(limits=c(0, 1500), breaks=c(0, 500,1000, 1500))
#plt <- plt + ggplot2::geom_smooth(method = "lm",se = FALSE,color = "black",linetype = 'longdash',alpha = 0.5)
plt


ggsave(plot=plt, filename="dge_power_col.category.png", width = 12,height = 5) 



#+ scale_colour_gradientn(colours = c("wnite","red"),                                   values = c(0,1,2,3, 4))

#+ scale_colour_gradient(low = "steelblue3",
#                                   high = "red")#+ scale_y_continuous(limits=c(0, 1500), breaks=c(0, 500,1000, 1500))
#plt <- plt + ggplot2::geom_smooth(method = "lm",se = FALSE,color = "black",linetype = 'longdash',alpha = 0.5)
plt




x4=x3[x3$category %in% c("T Cell", "Enterocyte", "Secretory"),]











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

