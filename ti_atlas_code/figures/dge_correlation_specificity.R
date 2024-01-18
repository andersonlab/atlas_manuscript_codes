
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

#f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/cellex/freeze03-ti-cd_healthy-discovery_labeled.predicted_celltype.esmu.csv.gz"
#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/cellex/freeze03-ti-cd_healthy-replication_labeled.predicted_celltype.esmu.csv.gz"

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/cellex/freeze03-ti-cd_healthy-disc_MT100.predicted_celltype.esmu.csv.gz"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/cellex/freeze03-ti-cd_healthy-repl_MT100.predicted_celltype.esmu.csv.gz"


#f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-all_samples/disease_status_merged-de_results.tsv"

df=read.table(f1, sep=",", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep=",", header=TRUE, stringsAsFactors = FALSE)

df$freeze <- "Discovery"
df1$freeze <- "Replication"

library(reshape2)
df_melt=melt(df)
df1_melt=melt(df1)
all_melt=rbind(df_melt, df1_melt)
names(all_melt) <- c("gene", "freeze", "cell_label", "specificity")


f3='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
df3=read.table(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)
keepme=c("label", "label__machine_retired")
df3=df3[,keepme]
names(df3)[names(df3)=="label__machine_retired"] <- "cell_label"

all_melt=merge(all_melt, df3, by = "cell_label", all.x = TRUE)
all_melt$cell_label <- NULL


l1=all_melt[all_melt$freeze=="Discovery",]
l1$freeze <- NULL
l1=reshape2::dcast(l1, gene ~ label , value.var = "specificity")
rownames(l1) <- l1$gene
l1$gene <- NULL
l1=as.matrix(l1)

l2=all_melt[all_melt$freeze=="Replication",]
l2$freeze <- NULL
l2=reshape2::dcast(l2, gene ~ label , value.var = "specificity")
rownames(l2) <- l2$gene
l2$gene <- NULL
l2=as.matrix(l2)
# 
# all=rbind(df, df1)
# 
# rownames(df) <- df$gene
# rownames(df1) <- df1$gene
# df$gene <- NULL
# df1$gene <- NULL

# rownames(df) <- df$gene
# rownames(df1) <- df1$gene
# df$gene <- NULL
# df1$gene <- NULL
# df$freeze<- NULL
# df1$freeze <- NULL

cr1=cor(l1, method="pearson", use = "pairwise.complete.obs")
cr2=cor(l2, method="pearson", use = "pairwise.complete.obs")


get_upper_tri <- function(cor_mat){
  cor_mat[lower.tri(cor_mat)]<- NA
  return(cor_mat)
}

#cr1 <- get_upper_tri(cr1)
#cr2 <- get_upper_tri(cr2)

cr1_melt <- reshape2::melt(cr1, na.rm = TRUE)
cr2_melt <- reshape2::melt(cr2, na.rm = TRUE)

cr1_melt$freeze <- "Discovery"
cr2_melt$freeze <- "Replication"

hh=rbind(cr1_melt, cr2_melt)

palette <-c("#B5BD61", #Stem
            "#E377C2", #Enterocyte
            "#279E68", #Secretory
            "#AEC7E8", #Tufn
            "#AA40FC", #Mesenchymal
            "#8C564B", #Myeloid
            "#D62728", #mast
            "#17BECF", #T cell
            "#1F77B4", #B cell
            "#FF7F0E") #plasma
fair_cols <- c("steelblue", "white", "red")
names(fair_cols) <- c(-1,0,1) 
lower=-1
upper=1
mid=0
#hh$freeze <- factor(hh$freeze, levels=c("discovery_dataset", "replication_dataset"),
#                               labels=c("Discovery", "Replication"))

mypalette=c("yellowgreen", "steelblue3")

names(hh)[3] <- "Correlation"

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
hh$Var1 <- factor(hh$Var1, levels=rev(oo$label))
hh$Var2 <- factor(hh$Var2, levels=oo$label)

write.table(hh, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/specificity.tsv", sep="\t", row.names = FALSE)

#mypalette=c("yellowgreen", "steelblue3")#e24e1b
library(ggplot2)
plt <- ggplot(data = hh, aes(Var2, Var1, fill = Correlation)) + geom_tile() + facet_grid(~ freeze, scales="free", space="free") + geom_tile(colour="white", size=0.25)
plt <- plt + scale_fill_gradient2(low = "#124e78", high = "coral3", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") 
plt <- plt + theme(strip.text.x = element_text(size = 25), 
                   axis.title.x = element_blank(), 
                   axis.title.y = element_blank(),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, size = 11, hjust = 1),
                   axis.text.y = element_text(size = 11),
                   strip.background = element_blank())
#plt <- plt + geom_dendro(cr1h)#+ scale_fill_manual(values=mypalette) + scale_alpha(range = c(0, 1))
plt <- plt #+ guides(fill= "none")
#plt <- plt + theme_minimal()  #+ coord_fixed()
plt
#ggsave(plt, filename=paste("dge_correlation_specificity.pdf"), width=20, height =12, device = cairo_pdf)
#ggsave(plt, filename=paste("dge_correlation_specificity.png"), width=20, height =12)





#f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/gene_mapping.tsv"
#x=read.table(f3, sep="\t", header=TRUE, stringsAsFactors = FALSE)

