#p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/dge/mast_random_effect_final/disease_status-run_8-mast-cp10k_greater_1-filter_mt_ribo_ig/"
#p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-old_samples_only/cell_prob_filter_onestd/"
#p5="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-new_samples_only/cell_prob_filter_onestd/"

library(stringr)

#p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/"
#p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/"

p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/"
p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/"
p5="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Full/"


setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")

f1=paste0(p1, "disease_status_merged-de_results.tsv")
f4=paste0(p4, "disease_status_merged-de_results.tsv")
f5=paste0(p5, "disease_status_merged-de_results.tsv")

df1=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df4=read.table(f4, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df5=read.table(f5, sep="\t", header=TRUE, stringsAsFactors = FALSE)


df1$comparison <- "Discovery"
df4$comparison <- "Replication"
df5$comparison <- "Full"

df1=df1[, c("gene_symbol", "label", "qvalue_bh_allcelltypes",  "n_cells", "log2fc", "mean_cp10k", "mean_counts",  "comparison", "category")] #"n_cells",
df4=df4[, c("gene_symbol", "label", "qvalue_bh_allcelltypes", "n_cells", "log2fc", "mean_cp10k", "mean_counts",  "comparison", "category")] #"n_cells"
df5=df5[, c("gene_symbol", "label", "qvalue_bh_allcelltypes", "n_cells", "log2fc", "mean_cp10k", "mean_counts",  "comparison", "category")] #"n_cells"


df=rbind(df1, df4, df5) 
df$neg_log10 <- -log10(df$qvalue_bh_allcelltypes)
df$significant <- df$qvalue_bh_allcelltypes < 0.05
df$significant <- factor(df$significant, levels = c(TRUE, FALSE), labels = c("True", "False"))
df$direction_regulation <- ifelse(df$log2fc >= 0, "up-regulated", "down-regulated")

df=df[df$comparison=="Full",]


######

#psm=unique(df_sig[grep("^PSM", as.character(df_sig$gene_symbol)), "gene_symbol"])
#hla=unique(df_sig[grep("HLA", as.character(df_sig$gene_symbol)), "gene_symbol"])
#hla=hla[!hla=="HHLA2"]

# file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/Reactome/FOLDING_ASSEMBLY_PEPTIDE_LOADING_OF_CLASS_I_MHC.tsv"
# oo1=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
# pl1=data.frame(annot_id="FOLDING_ASSEMBLY_PEPTIDE_LOADING_OF_CLASS_I_MHC",
#                gene_symbol= unique(sapply(str_split(oo1$MoleculeName, " "), "[[", 2)), stringsAsFactors=FALSE)
# pl1=pl1[!pl1$gene_symbol=="8",]
# #pl1$gene_symbol
# 

#mhc1_non_class= c("MICA", "MICB", "ULBP1", "ULBP2", "CD1A","CD1B","CD1C","CD1D","CD1E", "HFE", "MR1", "PROCR", "AZGP1", "UL18", "UL142", "UL37")
#mhc=c("HLA-DRA","HLA-DRB1",  "HLA-DPA1","HLA-DMA","HLA-DMB","HLA-DPB1","HLA-DQA1", "HLA-DRB5", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
#      "HLA-DOA", "HLA-DOB", "CD74", "CTSS")
mhc=c("IFNGR1")
#       "HLA-DQA1", "HLA-DRB5", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
#"HLA-DOA", "HLA-DOB"


# Pathway specific 

# filenames <- list.files("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/Reactome", pattern="*.tsv", full.names=TRUE)
# filenames
# ldf <- lapply(filenames, read.table)
# #res <- lapply(ldf, summary)
# names(ldf) <- gsub(".tsv", "", sapply(strsplit(filenames, "/" ), "[[", 8))
# library(plyr)
# dff <- plyr::ldply(ldf, data.frame)
# names(dff) <- c("annot_id", "gene_symbol")
# dff=dff[!dff$gene_symbol=="8",]
# dff <- dff %>% dplyr::distinct( annot_id, gene_symbol)
# unique(dff$annot_id)
# dff=dff[!(dff$annot_id %in% c("INTERFERON_ALPHA_BETA_SIGNALING", 
#                               "INTERFERON_GAMMA_SIGNALING", 
#                               "MHC_CLASS_II_ANTIGEN_PRESENTATION",
#                               "CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION")),] #"CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION"
# unique(dff$annot_id)


#################


#pl2=pl2[!pl2$gene_symbol=="8",]

# pl1$gene_symbol %in% pl2$gene_symbol
# pl2$gene_symbol %in% pl1$gene_symbol
# 
# pl=data.frame("annot_id_combined"= "CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION_FOLDING_ASSEMBLY_PEPTIDE_LOADING",
#               "gene_symbol"=unique(pl1$gene_symbol, pl2$gene_symbol), stringsAsFactors = FALSE)

#############
genes=mhc
#genes=mhc

df=df[df$gene_symbol %in% genes,]
df=df[!is.na(df$gene_symbol),]
df=df[!is.na(df$log2fc),]
#df=df[df$category %in% c("Stem cells", "Enterocyte", "Secretory"),]

df$category <- factor(df$category , levels=c("Stem cells", "Enterocyte", "Secretory",  "Myeloid",  "B Cell","T Cell","Mesenchymal", "B Cell plasma"))

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

df$label= as.character(df$label)
df$label <- factor(df$label, levels=oo$label)                            

# df$gene_symbol <- factor(df$gene_symbol, levels=c("HLA-A" ,"HLA-B","HLA-C","HLA-E", "HLA-F", "HLA-G","HLA-H", 
#                                                   "B2M", "TAP1"  , "TAP2" ,"TAPBP", 
#                                                   "HSPA5","CANX","CALR" ,"PDIA3",
#                                                   "SAR1B" , "SEC13", "SEC23A","SEC31A", "SEC24A","SEC24B","SEC24C", "SEC24D", 
#                                                   "ERAP2" , "ERAP1",
#                                                   "BCAP31" ,"BECN1",
#                                                   "PIK3R4","PIK3C3", 
#                                                   "ATG14" ))


#Only take genes upregulated in >2 epithelial cell types
# help=df[(df$category %in% c("Stem cells" , "Enterocyte", "Secretory")) & df$qvalue_bh_allcelltypes < 0.05,]
# help1=rev(sort(table(help$gene_symbol))>2)
# top_genes=names(help1[which(help1==TRUE)])
# 
# # Order genes
# order1=c("HLA-A" ,"HLA-B","HLA-C","HLA-E", "HLA-F", "HLA-G","HLA-H", "B2M", "TAP1","TAP2", "TAPBP", "CANX","CALR")
# order2=top_genes[grep("PSM", top_genes)]
# order2=names(table(help$gene_symbol)[order2])
# new_order=c(order1,order2, top_genes[!(top_genes %in% c(order1, order2))])

new_order=mhc
top_genes=mhc

df9=df[df$gene_symbol %in% top_genes,]
df9$gene_symbol <- factor(df9$gene_symbol, levels=new_order)
df9$label <- factor(df9$label, levels=rev(oo[oo$label %in% df9$label,]$label))

df_sig=df9[df9$significant=="True",]
df_sig$label <- factor(df_sig$label, levels=oo[oo$label %in% df_sig$label,]$label)



#df_sig[df_sig$category %in% c("Stem cells", "Enterocyte", "Secretory"), "compartment"] <- "Epithelial"
#df_sig[df_sig$category %in% c("Myeloid", "T Cell", "B Cell", "B Cell plasma"), "compartment"] <- "Immune"
#df_sig[df_sig$category %in% c("Mesenchymal"), "compartment"] <- "Mesenchymal"



mypalette=c("yellowgreen", "steelblue3")

###########Cap to Epi and Myeloid
#df_sig=df_sig[df_sig$category %in% c("Stem cells", "Enterocyte", "Secretory"),]


#df_sig=df_sig[df_sig$category %in% c("Stem cells", "Enterocyte", "Secretory",  "Myeloid", "B Cell"),]

library(paletteer)
paletter=paletteer_c("ggthemes::Red-Blue-White Diverging", 50)
paletter=as.vector(paletter)
max_scale=max(abs(min(df9$log2fc)), max(df9$log2fc))

library(ggplot2)
plt <- ggplot2::ggplot(df, ggplot2::aes(x=gene_symbol, y=label, fill=log2fc)) 
plt <- plt + geom_tile() 
#plt <- plt + scale_fill_gradient2(low = "navyblue", mid = "bisque", high = "red", midpoint = 0) #, breaks = c(-3,-1,0,1,3)
plt <- plt + scale_fill_gradientn(colours = rev(paletter), na.value = "white", limit = c(-max_scale,max_scale))
plt <- plt + geom_tile(data = df_sig, aes(x=gene_symbol, y=label ), colour = 'black', size=0.7) + geom_point(data = df_sig, aes(x=gene_symbol, y=label), shape=8, size=1)
#plt <- ggplot2::ggplot(df_sig, ggplot2::aes(x=gene_symbol, y=label)) 
#plt <- plt + ggplot2::geom_point(aes(size=neg_log10, colour=log2fc)) + scale_colour_gradient2(low = "bisque", mid = "darkgoldenrod1", high = "red", midpoint = 0.4) #, breaks = c(-3,-1,0,1,3)

plt <- plt + facet_grid(vars(category) , scales="free", space="free")  # shape=significant  #, shape=significant 
plt <- plt #+ geom_tile(data = df_sig, aes(x=gene_symbol, y=label ), colour = 'black')
#plt <- plt + facet_wrap(vars(type)) + scale_x_discrete(drop = TRUE)
#plt <- plt + ggtitle("CLASS I MHC MEDIATED ANTIGEN PROCESSING PRESENTATION")
plt <- plt + theme( strip.placement = "outside",  
                    axis.text.x = element_text(size=17, angle = 90, vjust = 0.5, hjust=1),
                    axis.text.y = element_text(size=17),
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
plt <- plt + labs(fill="Full cohort\nlog2FC")
plt


setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

ggsave(plt, filename=paste("dge_ifng_receptor.png"), width=7, height =5)
ggsave(plt, filename=paste("dge_ifng_receptor.pdf"), width=7, height =5, device = cairo_pdf, dpi=200)
#ggsave(plt, filename=paste("dge_mhc12.png"), width=13, height =12)
#ggsave(plt, filename=paste("dge_antigen_presentation.pdf"), width=20, height =14, device = cairo_pdf)











# ###
# 
# dff2 <- df_sig %>% group_by(gene_symbol, category, comparison) %>% dplyr::summarize(mean_log2fc = mean(log2fc))
# 
# 
# library(ggplot2)
# plt <- ggplot2::ggplot(dff2, ggplot2::aes(x=gene_symbol, y=category, fill=mean_log2fc)) 
# plt <- plt + geom_tile()
# plt <- plt + scale_fill_gradient2(low = "steelblue", mid = "bisque", high = "red", midpoint = 0) #, breaks = c(-3,-1,0,1,3)
# plt <- plt + geom_tile(data = dff2, aes(x=gene_symbol, y=category), colour = 'black')
# plt <- plt + facet_grid( ~ comparison , scales="free", space="free")  # shape=significant  #, shape=significant 
# plt <- plt #+ geom_tile(data = df_sig, aes(x=gene_symbol, y=label ), colour = 'black')
# plt <- plt + theme( strip.placement = "outside",  
#                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#                     axis.title.x = element_blank(), axis.title.y = element_blank(),
#                     strip.text.x = element_text(size = 12), 
#                     strip.text.y = element_text(size = 12, angle=360, hjust=0),
#                     strip.background = element_blank()) #+ ggplot2::scale_colour_gradientn(colours = terrain.colors(10))
# plt <- plt + scale_shape_manual(values=c(17, 16)) #+ guides(fill= "none")
# plt <- plt + labs(fill="log2FC")
# plt
# 
# 
# 
# 
# 
# 
