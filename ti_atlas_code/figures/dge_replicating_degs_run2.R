#p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/dge/mast_random_effect_final/disease_status-run_8-mast-cp10k_greater_1-filter_mt_ribo_ig/"
#p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-old_samples_only/cell_prob_filter_onestd/"
#p5="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-new_samples_only/cell_prob_filter_onestd/"

library(stringr)

#p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/"
#p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/"

p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/"
p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/"
p5="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Full/"


f1=paste0(p1, "disease_status_merged-de_results.tsv")
f4=paste0(p4, "disease_status_merged-de_results.tsv")
f5=paste0(p5, "disease_status_merged-de_results.tsv")

df1=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df4=read.table(f4, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df5=read.table(f5, sep="\t", header=TRUE, stringsAsFactors = FALSE)


df1$cohort <- "Discovery"
df4$cohort <- "Replication"
df5$cohort <- "Full"

df1=df1[, c("gene_symbol", "label", "qvalue_bh_allcelltypes",  "n_cells", "log2fc", "mean_cp10k", "mean_counts",  "cohort", "category")] #"n_cells",
df4=df4[, c("gene_symbol", "label", "qvalue_bh_allcelltypes", "n_cells", "log2fc", "mean_cp10k", "mean_counts",  "cohort", "category")] #"n_cells"
df5=df5[, c("gene_symbol", "label", "qvalue_bh_allcelltypes", "n_cells", "log2fc", "mean_cp10k", "mean_counts",  "cohort", "category")] #"n_cells"


df=rbind(df1, df4, df5) 
df$neg_log10 <- -log10(df$qvalue_bh_allcelltypes)
df$significant <- df$qvalue_bh_allcelltypes < 0.05
df$significant <- factor(df$significant, levels = c(TRUE, FALSE), labels = c("True", "False"))
df$direction_regulation <- ifelse(df$log2fc >= 0, "up-regulated", "down-regulated")

df=df[!df$cohort=="Full",]
df=df[df$category %in% c("Stem cells", "Enterocyte", "Secretory"),]

######

#genes=c("PIGR", "IFI27", "LCN2", "B2M", "REG1A", "REG1B")


genes=c("PIGR", "IFI27", "LCN2", "B2M", "REG1A", "REG1B")
df=df[df$gene_symbol %in% genes,]
df=df[!is.na(df$gene_symbol),]
df=df[!is.na(df$log2fc),]
df$category <- factor(df$category , levels=c("Stem cells", "Enterocyte", "Secretory"))

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

df$label= as.character(df$label)
df$label <- factor(df$label, levels=oo$label)                            

#df[df$category %in% c("Stem cells", "Enterocyte", "Secretory"), "compartment"] <- "Epithelial"
#df[df$category %in% c("Myeloid"), "compartment"] <- "Myeloid"

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
#df$label <- factor(df$label, levels=rev(oo$label))
df$label <- factor(df$label, levels=oo[oo$label %in% df$label,]$label)


#df9=df
#df9$gene_symbol <- factor(df9$gene_symbol, levels=rev(genes))
#df9$label <- factor(df9$label, levels=rev(oo[oo$label %in% df9$label,]$label))

df_sig=df[df$significant=="True",]


df$gene_symbol <- factor(df$gene_symbol, levels=rev(genes))
df_sig$gene_symbol <- factor(df_sig$gene_symbol, levels=rev(genes))
#df_sig$label <- factor(df_sig$label, levels=oo[oo$label %in% df_sig$label,]$label)
#df_sig[df_sig$category %in% c("Mesenchymal"), "compartment"] <- "Mesenchymal"

mypalette=c("yellowgreen", "steelblue3")

###########Cap to Epi and Myeloid
#df_sig=df_sig[df_sig$category %in% c("Stem cells", "Enterocyte", "Secretory"),]


#df_sig=df_sig[df_sig$category %in% c("Stem cells", "Enterocyte", "Secretory",  "Myeloid", "B Cell"),]

# library(paletteer)
# paletter=paletteer_c("ggthemes::Red-Blue-White Diverging", 50)
# paletter=as.vector(paletter)
# max_scale=max(abs(min(df9$log2fc)), max(df9$log2fc))
# min_scale=min(df9$log2fc)
palette=c("skyblue3", "tomato")

df$direction_regulation <- factor(df$direction_regulation, levels=c("down-regulated", "up-regulated"), labels=c("down", "up"))
df_sig$direction_regulation <- factor(df_sig$direction_regulation, levels=c("down-regulated", "up-regulated"), labels=c("down", "up"))


file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/replicating_degs.tsv"
check1=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
check1=check1[check1$gene_symbol %in% genes,]
check1$Discovery <- NULL
check1$Replication <- NULL
check1$category <- NULL
#check1[c("gene_symbol", "")]


df4=merge(df, check1, by=c("gene_symbol", "label"), all=TRUE)
df4_rep=df4 %>% subset(hit_type=="Discovery & Replication")
df4[is.na(df4$hit_type),]$hit_type <- "none"
#df_sig=merge(df_sig, rep, by="gene_symbol")



library(paletteer)
paletter=paletteer_c("ggthemes::Red-Blue Diverging", 50)
paletter=as.vector(paletter)
max_scale=max(abs(min(df$log2fc)), max(df$log2fc))
min_scale=min(df$log2fc)

#pick_genes=c("PIGR", "IFI27","LCN2", "B2M", 'REG1A', "REG1B")
#a <- ifelse(genes %in% pick_genes, "red", "black")
#rep$a <- a



library(ggplot2)
plt <- ggplot2::ggplot(df4, ggplot2::aes(y=gene_symbol, x=label, fill=log2fc)) 
#plt <- plt + geom_tile() 3
plt <- plt + scale_fill_gradient2(low = "navyblue", mid = "bisque", high = "red", midpoint = 0.3, limit = c(min_scale,max_scale)) #, breaks = c(-3,-1,0,1,3)
#plt <- plt + scale_fill_manual(values=palette)
#plt <- plt + scale_fill_gradient2(low = "navyblue", mid = "bisque", high = "red", midpoint = 0, limit = c(min_scale,max_scale)) #, breaks = c(-3,-1,0,1,3)
#plt <- plt + scale_fill_gradientn(colours = rev(paletter), na.value = "white") #, limit = c(min_scale,max_scale)# limit = c(-max_scale,max_scale))
plt <- plt + geom_tile(data = df_sig, aes(y=gene_symbol, x=label ), colour = 'black', size=0.7) 
plt <- plt + geom_point(data = df4_rep, aes(y=gene_symbol, x=label), shape=8, size=1)
#plt <- ggplot2::ggplot(df_sig, ggplot2::aes(x=gene_symbol, y=label)) 
#plt <- plt + ggplot2::geom_point(aes(size=neg_log10, colour=log2fc)) + scale_colour_gradient2(low = "bisque", mid = "darkgoldenrod1", high = "red", midpoint = 0.4) #, breaks = c(-3,-1,0,1,3)

plt <- plt + facet_grid(~cohort , scales="free", space="free")  # shape=significant  #, shape=significant 
plt <- plt #+ geom_tile(data = df_sig, aes(x=gene_symbol, y=label ), colour = 'black')
plt <- plt + facet_nested(~cohort+category, scales = 'free', space = 'free') 
#plt <- plt + ggtitle("CLASS I MHC MEDIATED ANTIGEN PROCESSING PRESENTATION")
plt <- plt + theme( strip.placement = "outside",  
                    axis.text.x = element_text(size=13, angle = 45, vjust = 1, hjust=1),
                    #axis.text.y = element_text(size=10,  colour=rev(a)),
                    axis.text.y = element_text(size=17),
                    axis.title.x = element_blank(), axis.title.y = element_blank(),
                    #strip.text.y = element_text(size = 32, angle=360, hjust=0),
                    #df_sig$label <- factor(df_sig$label, levels=oo[oo$label %in% df_sig$label,]$label)
                    #df_sig$label <- factor(df_sig$label, levels=oo[oo$label %in% df_sig$label,]$label)
                    strip.text.x = element_text(size = 18),
                    plot.title = element_text(size = 26),
                    #strip.text.x = element_text(size = 19),
                    strip.background = element_blank(),
                    #axis.text.x = ggtext::element_markdown(colour = a),
                    #legend.title = element_blank(), 
                    #legend.key.size = unit(1, 'cm'), #change legend key size
                    #legend.key.height = unit(1, 'cm'), #change legend key height
                    #legend.key.width = unit(1, 'cm'), #change legend key width
                    #legend.title = element_text(size=17), #change legend title font size
                    legend.text = element_text(size=15),
                    )#legend.position = "bottom"

#+ ggplot2::scale_colour_gradientn(colours = terrain.colors(10))
#plt <- plt + scale_colour_brewer(palette="Set1", direction=-1) #+ scale_shape_manual(values=c(17, 16)) #+ scale_size_manual(values=c(3, 1.7)) #+ scale_colour_discrete(na.translate = F) 
#plt <- plt + scale_y_discrete(drop = FALSE) #+ ggtitle("Epithelial")
#plt <- plt + coord_flip()
#plt <- plt + scale_fill_manual(values=mypalette)
plt <- plt  #+ scale_size_continuous(range=c(0, 6), breaks=c(0,1, 3, 6)) #breaks = c(-1,-0.5,0,0.5,1)
#plt <- plt + scale_colour_gradient2(breaks = c(-3,-1,0,1,3))
plt <- plt + scale_shape_manual(values=c(17, 16)) #+ guides(fill= "none")
plt <- plt + labs(fill="Log2FC")
plt

ggsave(plt, filename=paste("dge_top_replicating.png"), width=16, height =6)
ggsave(plt, filename=paste("dge_top_replicating.pdf"), width=16, height =6, device = cairo_pdf, dpi=200)


# rep$gene_symbol <- factor(rep$gene_symbol , levels=rev(rep$gene_symbol))
# 
# max_scale=max(abs(min(rep$mean_fc_r)), max(rep$mean_fc_r))
# min_scale=min(rep$mean_fc_r)
# 
# pick_genes=c("PIGR")
# 
# a <- ifelse(rep$gene_symbol %in% pick_genes, "red", "blue")
# rep$a <- a
# 
# gplt <- ggplot(rep, aes(x=gene_symbol, y=n_replicating, size=mean_fc_d, colour=mean_fc_r)) + geom_point() 
# gplt <- gplt + ylab("Number of cell types replicating") 
# gplt <- gplt + scale_colour_gradient2(low = "skyblue3", mid = "bisque", high = "tomato", midpoint = 0.5, limit = c(min_scale,max_scale)) #+ labs(fill="Gene regulation") #+ ggplot2::theme_bw()
# gplt <- gplt + theme(legend.title = element_blank(), 
#                     axis.title.y = element_blank(), 
#                     axis.title.x = element_blank(), 
#                     axis.text.x = ggtext::element_markdown(colour = a),
#                     strip.text.y = element_text(angle=360, size=12, hjust=0),
#                     strip.text.x = element_text(size = 15), 
#                     legend.text = element_text(size=15),
#                     strip.background = element_blank(), 
#                     ) #+ scale_x_discrete(aes(labels=a))
# #gplt <- gplt + scale_fill_gradient2(low = "navyblue", mid = "bisque", high = "red", midpoint = 0, limit = c(min_scale,max_scale)) 
# gplt <- gplt + theme_classic() + coord_flip() + labs(colour="Replication mean log2FC", size="Discovery mean log2FC") + xlab("")
# #plt <- plt + scale_y_continuous(limits=c(0, 700), breaks=c(0, 100,300, 400, 700))
# gplt
# 
# 
# #all <- ggpubr::ggarrange(plt, gplt, heights = c(3, 0.3), widths = c(1, 0.3), nrow=1)
# #all
# 
# all <- grid::grid.draw(egg::ggarrange(plots=list(plt, gplt), widths = c(2,0.5), nrow=1))
# all
# #ggsave(plt, filename=paste("dge_antigen_presentation.png"), width=20, height =14)
# #ggsave(plt, filename=paste("dge_mhc12.png"), width=15, height =14)
# #ggsave(plt, filename=paste("dge_antigen_presentation.pdf"), width=20, height =14, device = cairo_pdf)
# 
# 
# 
# 
# 
# setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")
# 
# ggsave(all, filename=paste("dge_hla_combined.png"), width=24, height =12)
# #ggsave(all, filename=paste("absolute_expression_mhc1.png"), width=24, height =12)
# 
# 
# 
# 
# 
# 
# 
# 
# # ###
# # 
# # dff2 <- df_sig %>% group_by(gene_symbol, category, comparison) %>% dplyr::summarize(mean_log2fc = mean(log2fc))
# # 
# # 
# # library(ggplot2)
# # plt <- ggplot2::ggplot(dff2, ggplot2::aes(x=gene_symbol, y=category, fill=mean_log2fc)) 
# # plt <- plt + geom_tile()
# # plt <- plt + scale_fill_gradient2(low = "steelblue", mid = "bisque", high = "red", midpoint = 0) #, breaks = c(-3,-1,0,1,3)
# # plt <- plt + geom_tile(data = dff2, aes(x=gene_symbol, y=category), colour = 'black')
# # plt <- plt + facet_grid( ~ comparison , scales="free", space="free")  # shape=significant  #, shape=significant 
# # plt <- plt #+ geom_tile(data = df_sig, aes(x=gene_symbol, y=label ), colour = 'black')
# # plt <- plt + theme( strip.placement = "outside",  
# #                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
# #                     axis.title.x = element_blank(), axis.title.y = element_blank(),
# #                     strip.text.x = element_text(size = 12), 
# #                     strip.text.y = element_text(size = 12, angle=360, hjust=0),
# #                     strip.background = element_blank()) #+ ggplot2::scale_colour_gradientn(colours = terrain.colors(10))
# # plt <- plt + scale_shape_manual(values=c(17, 16)) #+ guides(fill= "none")
# # plt <- plt + labs(fill="log2FC")
# # plt
# # 
# # 
# # 
# # 
# # 
# # 
