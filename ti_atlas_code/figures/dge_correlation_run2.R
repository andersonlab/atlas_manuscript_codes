
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

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

l1=l1[,oo$label]
l2=l2[,oo$label]

cr1=cor(l1, method="pearson", use = "pairwise.complete.obs")
cr2=cor(l2, method="pearson", use = "pairwise.complete.obs")

# getUpper.tri<-function(mat){
#   lt<-mat
#   lt[lower.tri(mat)]<- "NA"
#   mat<-as.data.frame(lt)
#   mat
# }
# 
# getLower.tri<-function(mat){
#   upper<-mat
#   upper[upper.tri(mat)] <- "NA"
#   #upper[diag(as.matrix(mat))]<-""
#   mat<-as.data.frame(upper)
#   mat
# }

#cr1 <- getUpper.tri(cr1)
#cr2 <- getUpper.tri(cr2)

# Add fold change
cr11 <- read.table("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/correlation_fc_discovery.tsv", sep="\t")
cr22 <- read.table("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/correlation_fc_replication.tsv", sep="\t")

names(cr11) <- row.names(cr1)
names(cr22) <- row.names(cr2)
row.names(cr11) <- row.names(cr1)
row.names(cr22) <- row.names(cr2)


# Discovery
ll1 <- matrix(NA, nrow = 49, ncol = 49)
ll1[upper.tri(ll1)] <- cr1[upper.tri(cr1)]
#ll1[upper.tri(ll1)] <- cr11[upper.tri(cr11)]
ll1[lower.tri(ll1)] <- cr11[lower.tri(cr11)]


# Replication
ll2 <- matrix(NA, nrow = 49, ncol = 49)
ll2[upper.tri(ll2)] <- cr2[upper.tri(cr2)]
#ll2[upper.tri(ll2)] <- cr22[upper.tri(cr22)]
ll2[lower.tri(ll2)] <- cr22[lower.tri(cr22)]


#cr11 <- getLower.tri(cr11)
#cr22 <- getLower.tri(cr22)

#cr1 <- as.matrix(cr1)
#cr2 <- as.matrix(cr2)
#cr11 <- as.matrix(cr11)
#cr22 <- as.matrix(cr22)

diag(ll1) <- 1
diag(ll2) <- 1

#ll1 <- cr1 + cr11
#ll2 <- cr2 + cr22

ll1=as.data.frame(ll1)
ll2=as.data.frame(ll2)

names(ll1) <- row.names(cr1)
row.names(ll1) <- row.names(cr1)
names(ll2) <- row.names(cr2)
row.names(ll2) <- row.names(cr2)

cr1_melt <- reshape2::melt(as.matrix(ll1), na.rm = TRUE)
cr2_melt <- reshape2::melt(as.matrix(ll2), na.rm = TRUE)

cr1_melt$freeze <- "Discovery"
cr2_melt$freeze <- "Replication"

hh=rbind(cr1_melt, cr2_melt)
names(hh)[3] <- "Correlation"

#fair_cols <- c("steelblue", "white", "red")
#names(fair_cols) <- c(-1,0,1) 
#lower=-1
#upper=1
#mid=0
#hh$freeze <- factor(hh$freeze, levels=c("discovery_dataset", "replication_dataset"),
#                               labels=c("Discovery", "Replication"))

mypalette=c("yellowgreen", "steelblue3")

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

#hh$Var1 <- as.character(hh$Var1)
#hh$Var2 <- as.character(hh$Var2)


hh1 <- hh
hh1$Var1 <- factor(hh1$Var1, levels=rev(oo$label))
hh1$Var2 <- factor(hh1$Var2, levels=oo$label)

#write.table(hh, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/specificity.tsv", sep="\t", row.names = FALSE)


#ppp=data.frame(cr1["Dendritic cell",], ll1["Dendritic cell",])

# hh1=hh[hh$freeze=="Discovery",]
# hh2=hh[hh$freeze=="Replication",]
# 
# 
# 
# hh1$Var1 <- factor(hh1$Var1, levels=rev(oo$label))
# hh1$Var2 <- factor(hh1$Var2, levels=oo$label)
# 
# hh2$Var1 <- factor(hh2$Var1, levels=rev(oo$label))
# hh2$Var2 <- factor(hh2$Var2, levels=oo$label)
# 
# hh3=rbind(hh1, hh2)
# 
# hh3$Var1
# 
# hh$freeze=factor(hh$freeze, levels=c("Discovery", "Replication"))


keep=c("category", "label")
oo=oo[,keep]
names(oo)[2] <- "Var2"
#names(hh1)[1] <- "Var1"
hh1=merge(hh1, oo, by="Var2")
hh1$category <- factor(hh1$category, levels=unique(oo$category))
#hh1$label <- factor(hh1$label, levels=rev(oo$label))
#extra <- 0.2
#max_catn <- max(hh1$Var1) + 1 + extra

hh1[hh1$Var2==hh1$Var1, "Correlation"] <- NA

hh1_sub=hh1[hh1$Correlation>0.70,]
hh1_sub=hh1_sub[!is.na(hh1_sub$Var2),]

library(paletteer)
paletter=paletteer_c("ggthemes::Red-Blue-White Diverging", 50)
paletter=as.vector(paletter)
#hh1$category <- factor(hh1$category, labels=c("","","","","","","",""))
#hh1_sub$category <- factor(hh1_sub$category, labels=c("","","","","","","",""))
#mypalette=c("yellowgreen", "steelblue3")#e24e1b
library(ggplot2)
library(ggh4x)
min_cor=min(hh1$Correlation, na.rm=TRUE)
max_cor=max(hh1$Correlation, na.rm=TRUE)
library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))



plt <- ggplot(data = hh1, aes(Var2,Var1,  fill = Correlation)) + geom_tile() + geom_tile(colour="white", size=0.25) 
plt <- plt + facet_nested(~freeze, scales="free", space="free") 
#plt <- plt + geom_point(data = hh1_sub, aes(x=Var2, y=Var1), shape=8, size=2)
#plt <- plt + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limit = c(min_cor,1), space = "Lab", name="Pearson\nCorrelation") 
#plt <- plt + scale_fill_gradientn(colours = rev(paletter),space = "Lab", na.value = "grey", limit = c(min_cor,max_cor))
plt <- plt + scale_fill_gradientn(colours = myPalette(100), space = "Lab", na.value = "grey", limit = c(min_cor,max_cor))
#plt <- plt + scale_fill_gradientn(colours = rev(paletter),space = "Lab", na.value = "grey", limit = c(-1,1))

  #scale_fill_gradient2(low = "#124e78", high = "red", mid = "white", midpoint = 0, na.value = "white", limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") 
plt <- plt + theme(strip.text.x = element_text(size = 20), 
                   axis.title.x = element_blank(), 
                   axis.title.y = element_blank(),
                   axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5,  hjust = 1),
                   axis.text.y = element_text(size = 13),
                   strip.background = element_blank())
#plt <- plt + geom_hline(yintercept = max_catn - 0.5 - (extra/2), colour='red')
#plt <- plt + geom_dendro(cr1h)#+ scale_fill_manual(values=mypalette) + scale_alpha(range = c(0, 1))
#plt <- plt + geom_vline(xintercept = c(3.5, 11.5, 16.5, 20.5, 27.5, 40.5, 46.5), linetype="dashed")
#plt <- plt + geom_rect(data = hh1_sub, mapping = aes(xmin = 0.5, xmax = 16.5, ymin = 33.5, ymax = 49.5), col = "black", fill = NA, linetype = 2)
#plt <- plt + geom_rect(data = hh1_sub, mapping = aes(xmin = 0.5, xmax = 24.5, ymin = 23.5, ymax = 25.5), col = "black", fill = NA, linetype = 2)
#plt <- plt + geom_rect(data = hh1_sub, mapping = aes(xmin = 24.5, xmax = 26.5, ymin = 0.5, ymax = 23.5), col = "black", fill = NA, linetype = 2)

#plt <- plt + geom_rect(data = hh1_sub, mapping = aes(xmin = 0.5, xmax = 49.5, ymin = 2.5, ymax = 3.5), col = "black", fill = NA, linetype = 2)
#plt <- plt + geom_rect(data = hh1_sub, mapping = aes(xmin = 27.5, xmax = 40.5, ymin = 9.5, ymax = 22.5), col = "black", fill = NA, linetype = 2)
#+ guides(fill= "none")
#plt <- plt + theme_minimal()  #+ coord_fixed()
#plt + ggforce::geom_mark_rect(aes(group = category))
plt

#ggsave(plt, filename=paste("dge_correlation_specificity.pdf"), width=20, height =12, device = cairo_pdf)
ggsave(plt, filename=paste("dge_correlation_specificity_fold_change.png"), width=20, height =12)
ggsave(plt, filename=paste("dge_correlation_specificity_fold_change.pdf"), width=20, height =12, device = cairo_pdf, dpi=200)




#f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/gene_mapping.tsv"
#x=read.table(f3, sep="\t", header=TRUE, stringsAsFactors = FALSE)

