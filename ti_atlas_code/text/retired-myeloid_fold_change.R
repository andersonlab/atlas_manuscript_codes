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


file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures/palette.tsv"
pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
palette=pal$palette_category
names(palette) <- pal$category

my_order <- all %>% group_by(freeze) %>% top_n(10, frac_dge)
cell_types_label = my_order[my_order$freeze == "Replication",]$label



haa <- all %>% group_by(freeze) %>% do(model = lm(n_cells ~ frac_dge, data = .))

dff=list()
for(i in 1:length(unique(all$freeze))){
  dff[[i]]=sqrt(summary(haa$model[[i]])$r.squared)
}
dff=data.frame(unlist(dff))
dff$freeze <- as.character(haa$freeze)
names(dff) <- c("R", "freeze")

all=merge(all, dff, by = "freeze", all.x = TRUE)
all$freeze_title <- paste0(all$freeze,": R=", round(all$R, 2))

# ADD RATIO OF CELLS
file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/ratio_of_cells.csv"
ratio=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
ratio=ratio[ratio$cohort=="Discovery",]
fold_change=(ratio[ratio$disease_status=="cd",]$n_cells-ratio[ratio$disease_status=="healthy",]$n_cells)/ratio[ratio$disease_status=="healthy",]$n_cells
ratio.df=data.frame(ratio[!duplicated(ratio$celltype_label__machine),]$celltype_label__machine,fold_change)
names(ratio.df)[1]="label__machine_retired"

f1='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
df=read.table(f1, sep=",", header=TRUE, stringsAsFactors = FALSE)
keep=c("label__machine_retired", "label")
df=df[,keep]
ratio.df=merge(ratio.df, df, by = "label__machine_retired", all.x = TRUE)
keep=c("label", "fold_change")
ratio.df=ratio.df[,keep]

# Add ratio 
all1=merge(all, ratio.df, by="label")


all_repel=all[all$label %in% cell_types_label,]






library(ggrepel)
plt <- ggplot(all1, aes(x=n_cells, y=frac_dge, colour=fold_change)) + geom_point() + facet_wrap(vars(freeze_title))#+ facet_wrap(vars(freeze))#+ annotation_custom(grob1) #+ geom_smooth(method=lm, se=FALSE) + geom_point() 
plt <- plt + theme_bw() + xlab("Number of cells") + ylab("Percent of DEG") + scale_colour_gradient2(low = "white", mid = "bisque", high = "red", midpoint = 2)#Number of significant differentially expressed genes + xlab("Number of genes tested") + ylab("Number of genes tested")
plt
plt <- plt + labs(color="Major cell type")  + geom_text_repel(data=all_repel, aes(label=label), size=3) #+ scale_colour_manual(values=palette)
plt <- plt + theme(text = element_text(size=15), strip.background = element_blank(), strip.text.x = element_text(size = 17)) 
#plt <- plt + ggplot2::geom_smooth(method = "lm",se = FALSE,color = "black",linetype = 'longdash',alpha = 0.5)
plt










