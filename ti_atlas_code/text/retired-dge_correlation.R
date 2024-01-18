setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")



#f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/disease_status_merged-de_results.tsv"
#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/disease_status_merged-de_results.tsv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/disease_status_merged-de_results.tsv"


df=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep="\t", header=TRUE, stringsAsFactors = FALSE)

df$freeze <- "Discovery"
df1$freeze <- "Replication"
all=rbind(df, df1)

all$significant <- all$qvalue_bh_allcelltypes < 0.05
all$significant <- factor(all$significant, levels = c(TRUE, FALSE), labels = c("True", "False"))
all$direction_regulation <- ifelse(all$log2fc >= 0, "up-regulated", "down-regulated")
all$direction_regulation <- factor(all$direction_regulation, levels=c("up-regulated", "down-regulated"))

all[all$category %in% c("T Cell", "B Cell", "B Cell plasma", "Mast", "Myeloid"), "compartment"] <- "Immune"
all[all$category %in% c("Secretory", "Enterocyte", "Stem cells", "Tuft cell"), "compartment"] <- "Epithelial"
all[all$category %in% c("Mesenchymal"), "compartment"] <- "Mesenchymal"


all=all[,c("label", "gene_symbol", "log2fc", "freeze")] #"direction_regulation"

all1=all[all$freeze=="Discovery",]
all2=all[all$freeze=="Replication",]

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

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

l1=l1[oo$label,]
l2=l2[oo$label,]

cr1=cor(t(l1), method="pearson", use = "pairwise.complete.obs")
cr2=cor(t(l2), method="pearson", use = "pairwise.complete.obs")

#cr1=cor(na.omit(t(l1)))
#cr2=cor(na.omit(t(l2)))

# Save for alter
# getLower.tri<-function(mat){
#   upper<-mat
#   upper[upper.tri(mat)]<- 0
#   mat<-as.data.frame(upper)
#   mat
# }


#write.table(cr1, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/correlation_fc_discovery.tsv", sep="\t", row.names = FALSE,  col.names = FALSE)
#write.table(cr2, "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/correlation_fc_replication.tsv", sep="\t", row.names = FALSE,  col.names = FALSE)


# Continue 


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


# Stem
a1=ll11[ll11$category_Var1=="Stem cells" & ll11$category_Var2=="Stem cells",]
a1=a1[order(a1$Correlation, decreasing=TRUE),]
#a1=a1[duplicated(a1$Correlation),]
#round(max(a1$Correlation),2)
round(mean(a1$Correlation),2)

a2=ll22[ll22$category_Var1=="Stem cells" & ll22$category_Var2=="Stem cells",]
a2=a2[order(a2$Correlation, decreasing=TRUE),]
#a2=a2[duplicated(a2$Correlation),]
#round(max(a2$Correlation),2)
round(mean(a2$Correlation),2)


# Enterocyte
b1=ll11[ll11$category_Var1=="Enterocyte" & ll11$category_Var2=="Enterocyte",]
b1=b1[order(b1$Correlation, decreasing=TRUE),]
b1=b1[c(grep("precursor", b1$Var1), grep("progenitor", b1$Var1)),]
b1=b1[c(grep("precursor", b1$Var2), grep("progenitor", b1$Var2)),]
#b1=b1[duplicated(b1$Correlation),]
round(mean(b1$Correlation),2)

b2=ll22[ll22$category_Var1=="Enterocyte" & ll22$category_Var2=="Enterocyte",]
b2=b2[order(b2$Correlation, decreasing=TRUE),]
b2=b2[c(grep("precursor", b2$Var1), grep("progenitor", b2$Var1)),]
b2=b2[c(grep("precursor", b2$Var2), grep("progenitor", b2$Var2)),]
round(mean(b2$Correlation),2)

# Stem cell and Enterocyte
c1=ll11[ll11$category_Var1=="Stem cells" & ll11$category_Var2=="Enterocyte",]
c1=c1[order(c1$Correlation, decreasing=TRUE),]
#c1=c1[c(grep("precursor", c1$Var1), grep("progenitor", c1$Var1)),]
#c1=c1[c(grep("precursor", c1$Var2), grep("progenitor", c1$Var2)),]
round(mean(b1$Correlation),2)

# B plasma and T
c1=ll11[ll11$category_Var1=="B Cell plasma" & ll11$category_Var2=="T Cell",]
c1=c1[order(c1$Correlation, decreasing=TRUE),]
c1=c1[c1$Var1=="B cell plasma IgA (1)",]
round(mean(c1$Correlation),2)

c1=ll22[ll22$category_Var1=="B Cell plasma" & ll22$category_Var2=="T Cell",]
c1=c1[order(c1$Correlation, decreasing=TRUE),]
c1=c1[c1$Var1=="B cell plasma IgA (1)",]
round(mean(c1$Correlation),2)

