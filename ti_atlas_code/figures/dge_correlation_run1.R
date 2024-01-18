setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")



#f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/disease_status_merged-de_results.tsv"
#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/disease_status_merged-de_results.tsv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/disease_status_merged-de_results.tsv"
f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Full/disease_status_merged-de_results.tsv"

df=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df2=read.table(f3, sep="\t", header=TRUE, stringsAsFactors = FALSE)


df$freeze <- "Discovery"
df1$freeze <- "Replication"
df2$freeze <- "Full"

all=rbind(df, df1, df2)

all$significant <- all$qvalue_bh_allcelltypes < 0.05
all$significant <- factor(all$significant, levels = c(TRUE, FALSE), labels = c("True", "False"))
all$direction_regulation <- ifelse(all$log2fc >= 0, "up-regulated", "down-regulated")
all$direction_regulation <- factor(all$direction_regulation, levels=c("up-regulated", "down-regulated"))

all[all$category %in% c("T Cell", "B Cell", "B Cell plasma", "Mast", "Myeloid"), "compartment"] <- "Immune"
all[all$category %in% c("Secretory", "Enterocyte", "Stem cells", "Tuft cell"), "compartment"] <- "Epithelial"
all[all$category %in% c("Mesenchymal"), "compartment"] <- "Mesenchymal"

table(all$significant)

all=all[,c("label", "gene_symbol", "log2fc", "freeze")] #"direction_regulation"

all1=all[all$freeze=="Discovery",]
all2=all[all$freeze=="Replication",]
all3=all[all$freeze=="Full",]

df3 <- reshape2::dcast(all, freeze + label ~ gene_symbol, value.var= "log2fc")



l1=df3[df3$freeze=="Discovery",]
rownames(l1) <- l1$label
l1$freeze <- NULL
l1$label <- NULL
l1=as.matrix(l1)

l2=df3[df3$freeze=="Replication",]
rownames(l2) <- l2$label
l2$freeze <- NULL
l2$label <- NULL
l2=as.matrix(l2)

l3=df3[df3$freeze=="Full",]
rownames(l3) <- l3$label
l3$freeze <- NULL
l3$label <- NULL
l3=as.matrix(l3)

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

l1=l1[oo$label,]
l2=l2[oo$label,]
l3=l3[oo$label,]

#Remove immunoglobuling
ff1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/repo/sc_nextflow/studies/gut-freeze003/data-variable_gene_filter.tsv"
toremove=read.table(ff1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

toremove = c(toremove$gene_symbol) # mhc, 

`%ni%` <- Negate(`%in%`)
l1=l1[,colnames(l1) %ni% toremove]
l2=l2[,colnames(l2) %ni% toremove]
l3=l3[,colnames(l3) %ni% toremove]


cr1=cor(t(l1), method="pearson", use = "pairwise.complete.obs")
cr2=cor(t(l2), method="pearson", use = "pairwise.complete.obs")
cr3=cor(t(l3), method="pearson", use = "pairwise.complete.obs")

#cr1=cor(na.omit(t(l1)))
#cr2=cor(na.omit(t(l2)))




# Save for alter
# getLower.tri<-function(mat){
#   upper<-mat
#   upper[upper.tri(mat)]<- 0
#   mat<-as.data.frame(upper)
#   mat
# }


write.table(cr1, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/correlation_fc_discovery.tsv", sep="\t", row.names = FALSE,  col.names = FALSE)
write.table(cr2, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/correlation_fc_replication.tsv", sep="\t", row.names = FALSE,  col.names = FALSE)
write.table(cr3, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/correlation_fc_full.tsv", sep="\t", row.names = FALSE,  col.names = FALSE)


#Manuscript text

cr1_melt <- reshape2::melt(cr1, na.rm = TRUE)
cr2_melt <- reshape2::melt(cr2, na.rm = TRUE)

cr1_melt$freeze <- "Discovery"
cr2_melt$freeze <- "Replication"

hh=rbind(cr1_melt, cr2_melt)
names(hh)[3] <- "Correlation"

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

hh$Var1 <- factor(hh$Var1, levels=rev(oo$label))
hh$Var2 <- factor(hh$Var2, levels=oo$label)

# Add category
oo1=oo[,c("label", "category")]
oo2=oo[,c("label", "category")]
names(oo1)[1] <- "Var1"
names(oo2)[1] <- "Var2"

hh=merge(hh, oo1, by="Var1")
names(hh)[names(hh)=="category"] <- "category_Var1"

hh=merge(hh, oo2, by="Var2")
names(hh)[names(hh)=="category"] <- "category_Var2"

ll1=hh[hh$freeze=="Discovery",]
ll2=hh[hh$freeze=="Replication",]

ll11=ll1[ll1$Correlation !=1,]
ll22=ll2[ll2$Correlation !=1,]

ll11=ll11[order(ll11$Correlation, decreasing=TRUE),]
ll22=ll22[order(ll22$Correlation, decreasing=TRUE),]


# Check correlation

# Discovery
a1=ll11[ll11$category_Var1=="Stem cells" & ll11$category_Var2=="Stem cells",]
a1=a1[order(a1$Correlation, decreasing=TRUE),]
#a1=a1[duplicated(a1$Correlation),]
#round(max(a1$Correlation),2)
round(mean(a1$Correlation),2)

# Replication
a2=ll22[ll22$category_Var1=="Stem cells" & ll22$category_Var2=="Stem cells",]
a2=a2[order(a2$Correlation, decreasing=TRUE),]
#a2=a2[duplicated(a2$Correlation),]
#round(max(a2$Correlation),2)
round(mean(a2$Correlation),2)


# Stem cell and Enterocyte
c1=ll11[ll11$category_Var1=="Stem cells" & ll11$category_Var2=="Enterocyte",]
c1=c1[order(c1$Correlation, decreasing=TRUE),]
c1=c1[c(grep("precursor", c1$Var2), grep("progenitor", c1$Var2)),]
#c1=c1[c(grep("precursor", c1$Var2), grep("progenitor", c1$Var2)),]
round(mean(c1$Correlation),2)

# Stem cell and Enterocyte
c1=ll22[ll22$category_Var1=="Stem cells" & ll22$category_Var2=="Enterocyte",]
c1=c1[order(c1$Correlation, decreasing=TRUE),]
c1=c1[c(grep("precursor", c1$Var2), grep("progenitor", c1$Var2)),]
#c1=c1[c(grep("precursor", c1$Var2), grep("progenitor", c1$Var2)),]
round(mean(c1$Correlation),2)




# # Enterocyte
# b1=ll11[ll11$category_Var1=="Enterocyte" & ll11$category_Var2=="Enterocyte",]
# b1=b1[order(b1$Correlation, decreasing=TRUE),]
# b1=b1[c(grep("precursor", b1$Var1), grep("progenitor", b1$Var1)),]
# b1=b1[c(grep("precursor", b1$Var2), grep("progenitor", b1$Var2)),]
# #b1=b1[duplicated(b1$Correlation),]
# round(mean(b1$Correlation),2)
# 
# b2=ll22[ll22$category_Var1=="Enterocyte" & ll22$category_Var2=="Enterocyte",]
# b2=b2[order(b2$Correlation, decreasing=TRUE),]
# b2=b2[c(grep("precursor", b2$Var1), grep("progenitor", b2$Var1)),]
# b2=b2[c(grep("precursor", b2$Var2), grep("progenitor", b2$Var2)),]
# round(mean(b2$Correlation),2)
# 
# 
# 
# # B plasma and T
# c1=ll11[ll11$category_Var1=="B Cell plasma" & ll11$category_Var2=="T Cell",]
# c1=c1[order(c1$Correlation, decreasing=TRUE),]
# c1=c1[c1$Var1=="B cell plasma IgA (1)",]
# round(mean(c1$Correlation),2)
# 
# c1=ll22[ll22$category_Var1=="B Cell plasma" & ll22$category_Var2=="T Cell",]
# c1=c1[order(c1$Correlation, decreasing=TRUE),]
# c1=c1[c1$Var1=="B cell plasma IgA (1)",]
# round(mean(c1$Correlation),2)


 #plot(l1[47,], l1[1,])
# l1[47,]
# # Show correlation
# #cr11=melt(cr1)
# #l11=as.data.frame(l1)
# allmisscols <- colnames(l1)[apply(l1,2, function(x) all(is.na(x)))]  
# #allmisscols <- colnames(l2)[apply(l2,2, function(x) all(is.na(x)))]  
# length(allmisscols)
# `%ni%` <- Negate(`%in%`)
# l11=l1[,colnames(l1) %ni% allmisscols]
# #l11=l2[,colnames(l2) %ni% allmisscols]
# 
# #l11=l1[!is.na(l1)]
# # PCA
# 
# #l11[is.na(l11)] <- 0
# 
# 
# pc <- prcomp(cr2,
#               center = TRUE,
#               scale. = TRUE)
#  
# df <- as.data.frame(pc$x[,1:2])
# df$label <- rownames(df)
# 
#  f3='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
#  df3=read.table(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)
#  keepme=c("label", "category")
#  df3=df3[,keepme]
# names(df)[3] <- "label"
# df=merge(df, df3, by="label") 
#  
# 
# file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/palette.tsv"
# pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
# pal=unique(pal[,c("palette_category", "category")])
# palette=pal$palette_category
# names(palette) <- pal$category
# df$category <- factor(df$category , levels=c("Stem cells", "Enterocyte", "Secretory","B Cell plasma",  "B Cell", "Myeloid","T Cell", "Mesenchymal")) 
# 
# 
# #ggplot method
# library(ggrepel) 
# 
# p1 <- ggplot(df, aes(PC1, PC2, colour = category)) +
#    geom_point(size = 3) + scale_colour_manual(values=palette) + theme_classic()  +
#    stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
#                 data = df[df$V3 == "1" | df$V3 == "2",], size = 1) + geom_text_repel(df, mapping=aes(x=PC1, y=PC2, label=label), size=4)
#  
# p1
 
 
# library(devtools)
# #install_github("vqv/ggbiplot")
# library(ggbiplot)
# g <- ggbiplot(pc,
#               obs.scale = 1,
#               var.scale = 1,
#               groups = training$Species,
#               ellipse = TRUE,
#               circle = TRUE,
#               ellipse.prob = 0.68)
# g <- g + scale_color_discrete(name = '')
# g <- g + theme(legend.direction = 'horizontal',
#                legend.position = 'top')
# print(g)




# Uncomment to see Figure


cr1_melt <- reshape2::melt(cr1, na.rm = TRUE)
cr2_melt <- reshape2::melt(cr2, na.rm = TRUE)
cr3_melt <- reshape2::melt(cr3, na.rm = TRUE)

cr1_melt$freeze <- "Discovery"
cr2_melt$freeze <- "Replication"
cr3_melt$freeze <- "Full"


hh=rbind(cr1_melt, cr2_melt, cr3_melt)

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

hh$freeze <- factor(hh$freeze, levels=c("Discovery", "Replication", "Full"))
#hh=hh[hh$freeze=="Full",]
#mypalette=c("yellowgreen", "steelblue3")#e24e1b
library(ggplot2)
library(RColorBrewer)
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))

min_cor=min(hh$Correlation, na.rm=TRUE)
max_cor=max(hh$Correlation, na.rm=TRUE)

plt <- ggplot(data = hh, aes(Var2, Var1, fill = Correlation)) + geom_tile() + facet_grid(~ freeze, scales="free", space="free") + geom_tile(colour="white", size=0.25)
plt <- plt + scale_fill_gradientn(colours = myPalette(100), space = "Lab", na.value = "grey", limit = c(min_cor,max_cor))
#plt <- plt + scale_fill_gradient2(low = "#124e78", high = "coral3", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation")
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












#ggsave(plt, filename="dge_correlation_fold_change.png", width=20, height =12, dpi = 300)
#ggsave(plt, filename=paste("dge_correlation_fold_change.pdf"), width=20, height =12, device = cairo_pdf)
#ggsave(plt, filename=paste("dge_correlation_fold_change.png"), width=20, height =12)


# # Check highest correlation
# 
# ll=hh[hh$Var1 != hh$Var2,]
# ll=ll[ll$Correlation > 0.8,]
# 
# file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
# oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
# oo=oo[,c("category", "label")]
# oo1=oo
# 
# names(oo)[2] <- "Var1"
# names(oo1)[2] <- "Var2"
# 
# ll1=merge(ll, oo, by="Var1")
# ll2=merge(ll1, oo1, by="Var2")






hh1=hh[hh$Correlation > 0.9 & hh$Correlation !=1,]
summary(hh$Correlation)



# melted_cormat$Correlation_threshold <- ifelse(melted_cormat$Correlation >= 0.4, 1, 0)
# library(ggplot2)
# plt <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = freeze)) + geom_tile(aes(alpha = Correlation_threshold)) + facet_grid(~ freeze, scales="free", space="free")
# plt <- plt #+ scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") 
# plt <- plt + theme(strip.text.x = element_text(size = 19), 
#                    axis.title.x = element_blank(), 
#                    axis.title.y = element_blank(),
#                    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 11, hjust = 1),
#                    axis.text.y = element_text(size = 11),
#                    strip.background = element_blank())
# plt <- plt + scale_fill_manual(values=mypalette) + scale_alpha(range = c(0, 1))
# plt <- plt + guides(fill= "none")
# #plt <- plt + theme_minimal()  #+ coord_fixed()
# plt

# papa=as.matrix(dist(t(dist_mat)))
# 
# 
# hclust_avg <- hclust(dist_mat, method = 'average')
# dhc <- as.dendrogram(hclust_avg)
# mydendro <- ggdendro::ggdendrogram(hclust_avg, rotate = TRUE, size = 2) 
# mydendro
