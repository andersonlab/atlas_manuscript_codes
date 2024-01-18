
library(stringr)


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

# f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/IBD_genes_and_Monogenic.tsv"
# x=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
# overlap=intersect(x[x$list=="Monogenic_IBD_genes",]$gene_symbol,x[x$list=="IBD_genes",]$gene_symbol) 
# x=x[!x$list=="IBD_genes",]
# xg=unique(x$gene_symbol)
#df_plt=df_plt[df_plt$gene_symbol %in% mylistg,]
#df_plt=merge(df_plt, mylist, by="gene_symbol")


f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/monogenic_ibd.csv"
y=read.table(f1, sep=",", header=TRUE, stringsAsFactors = FALSE)
yg=unique(y$gene_symbol)
#not=xg[!xg %in% yg]
df_plt=df_plt[df_plt$gene_symbol %in% yg,]

length(yg)
length(unique(df_plt$gene_symbol))
yg[!yg %in% unique(df_plt$gene_symbol)]

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
ord=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
df_plt$label <- factor(df_plt$label, levels= rev(ord$label))


# Cluster order
# xx=df_plt[,c("gene_symbol", "label", "specificity")]
# xx1=reshape2::dcast(df_plt, gene_symbol~label, value.var = "specificity", fun.aggregate=mean)
# rownames(xx1) <- xx1$gene_symbol
# xx1$gene_symbol <- NULL
# dd=dist(xx1)
# hc <- hclust(dd, "complete")
# gene_order=hc$labels[hc$order]

#Manual order
gene_order= c("GUCY2C" ,  "SLC37A4" ,  "AGR2"   ,  "FERMT1" ,  "SLCO2A1"  ,"ALPI" ,   
              "SLC26A3" , "SLC9A3","IL37" ,  "ADA2" , "NLRC4" ,"NCF2" , "TYMP" , "LACC1" ,"SYK" ,  "CYBB" , "NCF4" , "NCF1" , "BTK" ,  "BACH2" , "IL10"   , "IL2RB"  ,  "IL2RG"   ,
              "CD3G"  ,   "ZAP70"  ,  "CTLA4" ,   "FOXP3"  ,  "IL21"   ,  "CD40LG" , "ICOS"  ,  "IL2RA"  , "ITGB2"   , "IL10RA" ,  "ADA"   ,   "TNFAIP3" ,  "ARPC1B" , "CARD11" ,  "CARMIL2" , "IKZF1"  ,  "DOCK2"  ,  "WAS"   ,  "TGFB1"   , "TRIM22" , 
              "CARD8"  ,  "PIK3CD"  ,  "WIPF1",  "TGFBR2" ,     "DCLRE1C" , "CYBC1" ,   "HPS3"  ,   "CYBA"   ,  "MALT1"   ,  "SIRT1"   , "IKBKG"  ,  "PTPN2"  , 
 "PIK3R1" ,  "HSPA1L" ,    "SAMD9" ,   "CASP8"  ,  "STIM1"  , "STAT1"  ,  "TGFBR1" ,    "G6PC3" ,   "IRF2BP2",   "RBCK1" ,   "MVK"  ,    "LRBA"  ,   "TTC37"  ,   "CDC42"  ,  "NOP10" ,   "HPS1"  , 
"BCL10" ,   "RIPK1" ,   "SKIV2L"  , "STXBP3" ,  "ITCH"   ,  "XIAP"   ,  "COL7A1"  , "WNT2B"  ,  "SCGN"  ,   "MASP2"  ,  "ANKZF1" ,  "HPS4"   ,  "TRNT1"   , "PRKCD" ,  
 "ELF4" ,    "IL10RB",   "RAG1" ,    "TTC7A" ,   "NPC1"   , "LIG4"  ,   "PI4KA"  ,   "CD55"  ,   "NFKBIA" ,  "NFAT5"  ,  "RELA" ,   
 "STAT3" ,      "DKC1"  ,   "POLA1", "FMNL2")   


#unique(df_plt$gene_symbol)[!unique(df_plt$gene_symbol) %in% gene_order]
#yg[!yg %in% gene_order]

length(unique(gene_order))
#`%ni%` <- Negate(`%in%`)
#mylistg[mylistg %ni% unique(df_plt$gene_symbol)]

df_plt$gene_symbol <- factor(df_plt$gene_symbol, levels=gene_order)
#df_plt$gene_symbol <- factor(df_plt$gene_symbol, levels=mylistg)

# Add DGE
# f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Full/disease_status_merged-de_results.tsv"
# df2=read.table(f3, sep="\t", header=TRUE, stringsAsFactors = FALSE)
# df2=df2[!is.na(df2$qvalue_bh_allcelltypes),]
# df2_sig <- df2[df2$qvalue_bh_allcelltypes < 0.05,]
# df2_sig=df2_sig[df2_sig$gene_symbol %in% yg,]
# df2_sig=df2_sig[c("gene_symbol", "label")]
# df2_sig$significant<-"True"
# length(unique(df2_sig$gene_symbol))
# 
# df_plt_sig=left_join(df_plt, df2_sig, by=c("gene_symbol", "label"))
# df_plt_sig=df_plt_sig[!is.na(df_plt_sig$significant),]

ppalette <- as.character(rev(c("navy", "steelblue", "steelblue", "steelblue", "navy", "steelblue", "steelblue", "steelblue")))
library(paletteer)
paletter=paletteer_c("ggthemes::Red-Blue Diverging", 10)
paletter=as.vector(paletter)

#max(df_plt[df_plt$cohort=="Discovery",]$specificity)
#max(df_plt[df_plt$cohort=="Replication",]$specificity)

df_plt1 <- df_plt %>% group_by(gene_symbol, label, category) %>% dplyr::summarize(specificity_mean=mean(specificity))


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
#plt <- plt + ggtitle("Curated IBD genes")
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

ggsave(plot=plt, filename="dge_heatmap_IBD_genes_monogenic.png", width = 20,height = 10) 
ggsave(plot=plt, filename="dge_heatmap_IBD_genes_monogenic.pdf", width = 20,height = 10, device = cairo_pdf, dpi=200)

#ggsave(plot=plt, filename="dge_heatmap_IBD_genes_ramnik.png", width = 15,height = 10) 
