
library(stringr)

# 
# f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/cellex/freeze03-ti-cd_healthy-disc_MT100.predicted_celltype.esmu.csv.gz"
# f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/cellex/freeze03-ti-cd_healthy-repl_MT100.predicted_celltype.esmu.csv.gz"
# 
# df=read.table(f1, sep=",", header=TRUE, stringsAsFactors = FALSE)
# df1=read.table(f2, sep=",", header=TRUE, stringsAsFactors = FALSE)
# 
# df$freeze <- "Discovery"
# df1$freeze <- "Replication"
# 
# library(reshape2)
# df_melt=melt(df)
# df1_melt=melt(df1)
# all_melt=rbind(df_melt, df1_melt)
# names(all_melt) <- c("gene", "freeze", "cell_label", "specificity")
# all_melt=all_melt[!is.na(all_melt$specificity),]



f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/cellex/freeze03-ti-cd_healthy-disc_MT100.predicted_celltype.esmu.csv.gz"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/cellex/freeze03-ti-cd_healthy-repl_MT100.predicted_celltype.esmu.csv.gz"

df=read.table(f1, sep=",", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep=",", header=TRUE, stringsAsFactors = FALSE)

df$cohort <- "Discovery"
df1$cohort <- "Replication"

library(reshape2)
df_melt=melt(df)
df1_melt=melt(df1)
all_melt=rbind(df_melt, df1_melt)
names(all_melt) <- c("gene", "cohort", "cell_label", "specificity")
all_melt=all_melt[!is.na(all_melt$specificity),]


f3='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
df3=read.table(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)
keepme=c("label", "label__machine_retired", "category")
df3=df3[,keepme]
names(df3)[names(df3)=="label__machine_retired"] <- "cell_label"
all_melt=merge(all_melt, df3, by = "cell_label", all.x = TRUE)
all_melt$cell_label <- NULL

df_plt=all_melt
df_plt$category <- factor(df_plt$category , levels=c("Stem cells", "Enterocyte", "Secretory","Myeloid",  "B Cell", "T Cell", "Mesenchymal","B Cell plasma")) 

#df_plt <- reshape2::dcast(df_plt, gene + label + category ~ freeze , value.var = "specificity")

names(df_plt)[1] <- "gene_ensemble"

p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/gene_mapping.tsv"
mapping=read.table(p1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df_plt=merge(df_plt, mapping, by="gene_ensemble")
df_plt$gene_ensemble <- NULL

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/IBD_genes_and_Monogenic.tsv"
mylist=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
overlap=intersect(mylist[mylist$list=="Monogenic_IBD_genes",]$gene_symbol,mylist[mylist$list=="IBD_genes",]$gene_symbol) 
#mylist=mylist[!mylist[mylist$list=="IBD_genes",]$gene_symbol %in% overlap,]
mylistg=unique(mylist$gene_symbol)

df_plt=df_plt[df_plt$gene_symbol %in% mylistg,]
df_plt=merge(df_plt, mylist, by="gene_symbol")

# mylistg=c("PRDM1", "IL2RA", "CARD9", "IL23R", "PTPN22", "IKZF1", "INPP5E", 'EP300', 
#          "EBF1", "IFIH1", "NFKB1", "SMAD3", "LRRK2", "GPR35", "JAK2", "SLC22A5", "HNF4A", "TYK2", "NOD2", "NKX2-3")
# df_plt=df_plt[df_plt$gene_symbol %in% mylistg,]


# library(paletteer)
# paletter=paletteer_c("ggthemes::Red-Blue-White Diverging", 50)
# paletter=as.vector(paletter)
# max_scale=max(df_plt$specificity)
# min_scale=min(df_plt$specificity)

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
ord=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
df_plt$label <- factor(df_plt$label, levels= rev(ord$label))

# Select IBD genes
 select_list="IBD_genes"
#select_list="Monogenic_IBD_genes"
 df_plt=df_plt[df_plt$list==select_list,]


# Cluster
xx=df_plt[,c("gene_symbol", "label", "specificity")]
xx1=reshape2::dcast(df_plt, gene_symbol~label, value.var = "specificity", fun.aggregate=mean)
rownames(xx1) <- xx1$gene_symbol
xx1$gene_symbol <- NULL
dd=dist(xx1)
#hc <- hclust(dd, "ave")
hc <- hclust(dd, "complete")
#hc <- hclust(dd, "ward.D")
gene_order=hc$labels[hc$order]
gene_order

gene_order=c("RNF186","MST1"  ,  "FUT2"  ,  "HNF4A" ,
 "PDLIM5" , "NR5A2" , "SMAD3"  , "SHARPIN" ,"ATG16L1", "SDF2L1" , "TMEM258",
"SLC39A8","PPIF"  ,   "PLCG2"  , "ATG4C" ,  "TYK2"  ,  "IFIH1"  , "RELA"  ,  "STAT3" ,  
"JAK2"  , "ITGAV"  , "OSMR"  ,"DOK2"  ,  "IL10RA"  ,"ITGA4" ,  "ADCY7" ,  "IL2RA" ,  "IL10"  ,"TAGAP" ,  "TNFAIP3",  "IL23R"  , "IL18RAP" ,"PTPN22" , 
"SP140","CARMIL2" ,"CCR7"  ,  "BACH2"  , "LRRK2" ,  "NCF4"  ,"CARD9"  , "SLAMF8" ,   "HGFAC"  , "PTAFR" ,  "NOD2"   )

length(unique(gene_order))

df_plt$gene_symbol <- factor(df_plt$gene_symbol, levels=gene_order)
#df_plt$gene_symbol <- factor(df_plt$gene_symbol, levels=mylistg)

# Add DGE
f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Full/disease_status_merged-de_results.tsv"
df2=read.table(f3, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df2=df2[!is.na(df2$qvalue_bh_allcelltypes),]
df2_sig <- df2[df2$qvalue_bh_allcelltypes < 0.05,]
df2_sig=df2_sig[df2_sig$gene_symbol %in% mylistg,]
df2_sig=df2_sig[c("gene_symbol", "label")]
df2_sig$significant<-"True"

df_plt_sig=left_join(df_plt, df2_sig, by=c("gene_symbol", "label"))
df_plt_sig=df_plt_sig[!is.na(df_plt_sig$significant),]

ppalette <- as.character(rev(c("navy", "steelblue", "steelblue", "steelblue", "navy", "steelblue", "steelblue", "steelblue")))
library(paletteer)
paletter=paletteer_c("ggthemes::Red-Blue Diverging", 10)
paletter=as.vector(paletter)

ggthemes$excel$themes
max(df_plt[df_plt$cohort=="Discovery",]$specificity)
max(df_plt[df_plt$cohort=="Replication",]$specificity)


df_plt1 <- df_plt %>% group_by(gene_symbol, label, category, list) %>% dplyr::summarize(specificity_mean=mean(specificity))
colnames(df_plt)
  
library(ggplot2)
plt <- ggplot2::ggplot(df_plt1, ggplot2::aes(x=gene_symbol, y=label, fill=specificity_mean)) 
#plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x=gene_symbol, y=label, fill=specificity)) 
#plt <- plt + scale_fill_gradient2(low = "blue",  high = "red") #, breaks = c(-3,-1,0,1,3)
plt <- plt + scale_fill_gradientn(colours = rev(paletter), na.value = "white", limit = c(0,1))
#plt <- plt + scale_fill_gradientn(low = "blue", high = "red") #, na.value = "white", limit = c(-max_scale,max_scale))
#plt <- plt + scale_fill_gradientn(colours="heat") 
plt <- plt + geom_tile(data = df_plt1, aes(x=gene_symbol, y=label)) 
#plt <- plt + geom_bar(data = df_plt, aes(x=gene_symbol, y=label)) 
#plt
#plt <- plt + geom_tile(data = df_plt_sig, aes(x=gene_symbol, y=label ), colour = 'black', size=0.7) 
#plt <- plt + geom_point(data = df_plt_sig, aes(x=gene_symbol, y=label), shape=8, size=1)

#plt <- ggplot2::ggplot(df_sig, ggplot2::aes(x=gene_symbol, y=label)) 
#plt <- plt + ggplot2::geom_point(aes(size=neg_log10, colour=log2fc)) + scale_colour_gradient2(low = "bisque", mid = "darkgoldenrod1", high = "red", midpoint = 0.4) #, breaks = c(-3,-1,0,1,3)

#plt <- plt + facet_grid(vars(category) , scales="free", space="free")  # shape=significant  #, shape=significant 
plt <- plt #+ geom_tile(data = df_sig, aes(x=gene_symbol, y=label ), colour = 'black')
plt <- plt + facet_grid(vars(category), scales="free", space="free") # + scale_x_discrete(drop = TRUE)
#plt <- plt + ggtitle("GWAS associated IBD genes")
#plt <- plt + ggtitle("Monogenic IBD genes")
plt <- plt + theme( strip.placement = "outside",  
                    axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1),
                    axis.text.y = element_text(size=10),
                    axis.title.x = element_blank(), axis.title.y = element_blank(),
                    strip.text.y = element_text(size = 19, angle=360, hjust=0),
                    plot.title = element_text(size = 24),
                    #strip.text.x = element_text(size = 19),
                    strip.background = element_blank(),
                    legend.key.size = unit(1, 'cm'), #change legend key size
                    legend.key.height = unit(1, 'cm'), #change legend key height
                    legend.key.width = unit(1, 'cm'), #change legend key width
                    legend.title = element_text(size=17), #change legend title font size
                    legend.text = element_text(size=15))

#+ ggplot2::scale_colour_gradientn(colours = terrain.colors(10))
#plt <- plt + scale_colour_brewer(palette="Set1", direction=-1) #+ scale_shape_manual(values=c(17, 16)) #+ scale_size_manual(values=c(3, 1.7)) #+ scale_colour_discrete(na.translate = F) 
#plt <- plt + scale_y_discrete(drop = FALSE) #+ ggtitle("Epithelial")
#plt <- plt + coord_flip()
#plt <- plt + scale_fill_manual(values=mypalette)
plt <- plt  #+ scale_size_continuous(range=c(0, 6), breaks=c(0,1, 3, 6)) #breaks = c(-1,-0.5,0,0.5,1)
#plt <- plt + scale_colour_gradient2(breaks = c(-3,-1,0,1,3))
plt <- plt + scale_shape_manual(values=c(17, 16)) #+ guides(fill= "none")
plt <- plt + labs(fill="Specificity")
plt

ggsave(plot=plt, filename="dge_heatmap_IBD_genes.png", width = 15,height = 10) 
ggsave(plot=plt, filename="dge_heatmap_IBD_genes.pdf", width = 15,height = 10, device = cairo_pdf, dpi=200)

#ggsave(plot=plt, filename="dge_heatmap_IBD_genes_monogenic.png", width = 20,height = 10) 
#ggsave(plot=plt, filename="dge_heatmap_IBD_genes_ramnik.png", width = 15,height = 10) 
