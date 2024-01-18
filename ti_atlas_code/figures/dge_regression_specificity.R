
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

#f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/cellex/freeze03-ti-cd_healthy-discovery_labeled.predicted_celltype.esmu.csv.gz"
#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/cellex/freeze03-ti-cd_healthy-replication_labeled.predicted_celltype.esmu.csv.gz"

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/cellex/freeze03-ti-cd_healthy-disc_MT100.predicted_celltype.esmu.csv.gz"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/cellex/freeze03-ti-cd_healthy-repl_MT100.predicted_celltype.esmu.csv.gz"

df=read.table(f1, sep=",", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep=",", header=TRUE, stringsAsFactors = FALSE)

df$freeze <- "Discovery"
df1$freeze <- "Replication"

library(reshape2)
df_melt=melt(df)
df1_melt=melt(df1)
all_melt=rbind(df_melt, df1_melt)
names(all_melt) <- c("gene", "freeze", "cell_label", "specificity")
all_melt=all_melt[!is.na(all_melt$specificity),]


f3='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
df3=read.table(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)
keepme=c("label", "label__machine_retired", "category")
df3=df3[,keepme]
names(df3)[names(df3)=="label__machine_retired"] <- "cell_label"
all_melt=merge(all_melt, df3, by = "cell_label", all.x = TRUE)
all_melt$cell_label <- NULL

df_plt=all_melt
df_plt$category <- factor(df_plt$category , levels=c("Stem cells", "Enterocyte", "Secretory","B Cell plasma",  "B Cell", "Myeloid","T Cell", "Mesenchymal")) 

df_plt <- reshape2::dcast(df_plt, gene + label + category ~ freeze , value.var = "specificity")


# haa <- df_plt %>% group_by(category) %>% do(model = lm(Discovery ~ Replication, data = .))
# dff=list()
# for(i in 1:length(haa$category)){
#   dff[[i]]=sqrt(summary(haa$model[[i]])$r.squared)
# }
# dff=data.frame(unlist(dff))
# dff$category <- as.character(haa$category)
# names(dff) <- c("R", "category")
# df_plt=merge(df_plt, dff, by = "category", all.x = TRUE)
# df_plt$category_title <- paste0(df_plt$category,": R=", round(df_plt$R, 2))

#####

haa <- df_plt %>% group_by(label) %>% do(model = lm(Discovery ~ Replication, data = .))
dff=list()
for(i in 1:length(haa$label)){
  dff[[i]]=sqrt(summary(haa$model[[i]])$r.squared)
}
dff=data.frame(unlist(dff))
dff$label <- as.character(haa$label)
names(dff) <- c("R", "label")
df_plt=merge(df_plt, dff, by = "label", all.x = TRUE)

#median=summary(df_plt$R)[3]

df_plt$R_high <- df_plt$R >= 0.6
df_plt$R_high = factor(df_plt$R_high, levels=c(TRUE, FALSE), labels=c("Pearson's R > 0.6", "Pearson's R < 0.6"))

####Remember: lm(y~x) not same as lm(x~y)
palette <-c("#B5BD61", #Stem
            "#E377C2", #Enterocyte
            "#279E68", #Secretory
            "#AEC7E8", #Tufn
            "#AA40FC", #Mesenchymal
            "#8C564B", #Myeloid
            "#D62728", #mast
            "#17BECF", #T cell
            "#1F77B4", #B cell
            "#FF7F0E" #plasma
) # "#17BECF",
keepmee=c("label", "R")
df_plt_repel=df_plt[,keepmee]
df_plt_repel=distinct(df_plt_repel)
df_plt_repel[df_plt_repel$R < 0.6,]
#df_plt$category_title <- paste0(df_plt$category,": R=", round(df_plt$R, 2))
#mypalette=c("red",  "yellowgreen", "steelblue3","grey")
#names(mypalette) <- c("Discovery & Replication", "Discovery", "Replication", "Neither")

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/palette.tsv"
pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
palette=pal$palette_label
names(palette) <- pal$label

df_plt[df_plt$category %in% c("Stem cells", "Enterocyte", "Secretory"),"compartment"]  <- "Epithelial"
df_plt[df_plt$category %in% c("B Cell plasma", "B Cell", "Myeloid", "T Cell", "Mesenchymal"),"compartment"] <- "Immune and Mesenchymal"
#df_plt[df_plt$category %in% c("Mesenchymal"),"compartment"] <- "Mesenchymal"

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
ord=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

df_plt$label <- factor(df_plt$label, levels= ord$label)
#palette=palette[ord$label]


check=unique(df_plt[,c("label","category", "R")])
check <- check %>% group_by(category) %>% dplyr::summarise(R_mean=mean(R))
df_plt=merge(df_plt, check, by="category")
df_plt$category_title <- paste0(df_plt$category,": mean R=", round(df_plt$R_mean, 2))


# Plot the merged data ####################################################
#mypalette = c("deepskyblue3", "darkseagreen3") #, "coral3"
plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = Discovery, y =  Replication, group=label, colour=label)) #x = temp_asdfs_v2, y =  temp_asdfs_v1
plt <- plt + ggplot2::theme_bw(base_size = 12) + xlim(-0.05,1) + ylim(-0.05,1) + scale_colour_manual(values=palette) #+ geom_point(colour="steelblue")
plt <- plt + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'dashed', alpha = 0.7, fullrange=TRUE) #+ stat_poly_eq() #, colour="steelblue"
plt <- plt + ggplot2::geom_abline(slope = 1, intercept = 0)
plt <- plt + ggplot2::geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5)
plt <- plt + ggplot2::geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.5)
#plt <- plt + ggplot2::scale_color_manual(values=mypalette)
plt <- plt + ggplot2::labs(x = "Discovery specificity",y = "Replication specificity",color = "") #FDR < 0.05
plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, desc(R)) , scales = "free", ncol=4)
#plt <- plt + ggplot2::facet_wrap(~ label, scales = "free", ncol=5)
plt <- plt + theme(strip.background = element_blank(), strip.text=element_text(size=11), legend.position = "none") 
#plt <- plt + ggplot2::scale_color_brewer(palette = "Dark2")
plt
ggsave(plt, filename=paste("dge_regression_specificity.png"), width=10, height =5)
ggsave(plt, filename=paste("dge_regression_specificity.pdf"), width=10, height =5, device = cairo_pdf, dpi=200)



# Check differences between specificity in two cohorts
file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/tsv/gene_mapping.tsv"
mapping=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)


names(df_plt)[names(df_plt)=="gene"] <- "gene_ensemble"
df_plt=merge(df_plt, mapping, by="gene_ensemble")

df_plt1=df_plt[df_plt$label %in% c("T cell CD4- CD8-"),]
df_plt1=df_plt1[!is.na(df_plt1$Discovery),]
df_plt1=df_plt1[!is.na(df_plt1$Replication),]
df_plt2=df_plt[df_plt$label %in% c( "T cell CD8+ FGFBP2+ effector"),]
df_plt2=df_plt2[!is.na(df_plt2$Discovery),]
df_plt2=df_plt2[!is.na(df_plt2$Replication),]

df_plt11=df_plt1[df_plt1$Discovery > df_plt1$Replication,]
df_plt12=df_plt1[df_plt1$Discovery < df_plt1$Replication,]
unique(df_plt11$gene_symbol)


# Check specificity of HLA genes
df_hla=df_plt[grep("HLA-", df_plt$gene_symbol),]

plt <- ggplot2::ggplot(df_hla, ggplot2::aes(x = category, y =  Discovery, fill=gene_symbol)) #x = temp_asdfs_v2, y =  temp_asdfs_v1
plt <- plt + ggplot2::theme_bw(base_size = 12) + geom_boxplot()#+ xlim(-0.05,1) + ylim(-0.05,1) #+ scale_colour_manual(values=palette) #+ geom_point(colour="steelblue")
#plt <- plt + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'dashed', alpha = 0.7, fullrange=TRUE) #+ stat_poly_eq() #, colour="steelblue"
#plt <- plt + ggplot2::geom_abline(slope = 1, intercept = 0)
#plt <- plt + ggplot2::geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5)
#plt <- plt + ggplot2::geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.5)
#plt <- plt + ggplot2::scale_color_manual(values=mypalette)
plt <- plt + ggplot2::labs(x = "Discovery specificity",y = "Replication specificity",color = "") #FDR < 0.05
#plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, desc(R)) , scales = "free", ncol=4)
#plt <- plt + ggplot2::facet_wrap(~ label, scales = "free", ncol=5)
plt <- plt + theme(strip.background = element_blank(), strip.text=element_text(size=11)) 
#plt <- plt + ggplot2::scale_color_brewer(palette = "Dark2")
plt








# leg1 <-
#   df_plt %>%
#   filter(label %in% ord[ord$category=="Stem cells","label"]) %>%
#   ggplot(ggplot2::aes(x = Discovery, y =  Replication, group=label, color=label)) + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5) + 
#   scale_colour_manual(values = palette[ord[ord$category=="Stem cells","label"]], name = "Stem cells") #+ theme(legend.position="bottom")
# 
# leg2 <-
#   df_plt %>%
#   filter(label %in% ord[ord$category=="Enterocyte","label"]) %>%
#   ggplot(ggplot2::aes(x = Discovery, y =  Replication, group=label, color=label)) + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5) + 
#   scale_colour_manual(values = palette[ord[ord$category=="Enterocyte","label"]], name = "Enterocyte") #+ theme(legend.position="bottom")
# 
# leg3 <-
#   df_plt %>%
#   filter(label %in% ord[ord$category=="Secretory","label"]) %>%
#   ggplot(ggplot2::aes(x = Discovery, y =  Replication, group=label, color=label)) + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5) + 
#   scale_colour_manual(values = palette[ord[ord$category=="Secretory","label"]], name = "Secretory") #+ theme(legend.position="bottom")
# 
# leg4 <-
#   df_plt %>%
#   filter(label %in% ord[ord$category=="B Cell plasma","label"]) %>%
#   ggplot(ggplot2::aes(x = Discovery, y =  Replication, group=label, color=label)) + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5) + 
#   scale_colour_manual(values = palette[ord[ord$category=="B Cell plasma","label"]], name = "B Cell plasma") #+ theme(legend.position="bottom")
# 
# leg5 <-
#   df_plt %>%
#   filter(label %in% ord[ord$category=="B Cell","label"]) %>%
#   ggplot(ggplot2::aes(x = Discovery, y =  Replication, group=label, color=label)) + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5) + 
#   scale_colour_manual(values = palette[ord[ord$category=="B Cell","label"]], name = "B Cell") #+ theme(legend.position="bottom")
# 
# leg6 <-
#   df_plt %>%
#   filter(label %in% ord[ord$category=="Myeloid","label"]) %>%
#   ggplot(ggplot2::aes(x = Discovery, y =  Replication, group=label, color=label)) + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5) + 
#   scale_colour_manual(values = palette[ord[ord$category=="Myeloid","label"]], name = "Myeloid") #+ theme(legend.position="bottom")
# 
# leg7 <-
#   df_plt %>%
#   filter(label %in% ord[ord$category=="T Cell","label"]) %>%
#   ggplot(ggplot2::aes(x = Discovery, y =  Replication, group=label, color=label)) + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5) + 
#   scale_colour_manual(values = palette[ord[ord$category=="T Cell","label"]], name = "T Cell") #+ theme(legend.position="bottom")
# 
# leg8 <-
#   df_plt %>%
#   filter(label %in% ord[ord$category=="Mesenchymal","label"]) %>%
#   ggplot(ggplot2::aes(x = Discovery, y =  Replication, group=label, color=label)) + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5) + 
#   scale_colour_manual(values = palette[ord[ord$category=="Mesenchymal","label"]], name = "Mesenchymal") #+ theme(legend.position="bottom")
# 
# 
# library(cowplot)
# plt1 <- plot_grid(
#   plt
#   , plot_grid(get_legend(leg1 ), 
#               get_legend(leg2 ), 
#               get_legend(leg3 ),
#               get_legend(leg4 ),
#               get_legend(leg5 ),
#               get_legend(leg6 ),
#               get_legend(leg7 ),
#               get_legend(leg8 ), 
#               nrow=1), nrow=2)
#     
#   
# 
# #ggsave(plt, filename=paste("dge_fc_compare.png"), width=12, height =5, dpi = 300)
# ggsave(plt, filename=paste("dge_regression_specificity.png"), width=10, height =5)
# 
# 
# 
# plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = Discovery, y =  Replication, group=label, color=label)) #x = temp_asdfs_v2, y =  temp_asdfs_v1
# plt <- plt + ggplot2::theme_bw(base_size = 12) + scale_colour_manual(values=palette) + theme(legend.title = element_blank())
# plt <- plt + theme(strip.background = element_blank())
# #plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, desc(R)) , scales = "free", ncol=4)
# #plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Neither'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
# #plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
# #plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
# #plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery & Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
# plt <- plt + ggplot2::geom_smooth( method = "lm", se = FALSE,  linetype = 'longdash',size = 1) #+ stat_poly_eq()
# plt <- plt + ggplot2::labs(x = "Discovery log2 fold change",y = "Replication log2 fold change",color = "") #FDR < 0.05
# 
# #plt <- plt + ggrepel::geom_label_repel( aes(label=label))
# 
# plt
# ggsave(plt, filename=paste("dge_fc_compare_lm_single.png"), width=8, height =6)
# 
# 
# 
# # Plot the merged data ####################################################
# plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = Discovery, y =  Replication, group=label, colour=compartment)) #x = temp_asdfs_v2, y =  temp_asdfs_v1
# plt <- plt + ggplot2::theme_bw(base_size = 12) + scale_colour_manual(values=mypalette) #+ theme(legend.position='none')
# plt <- plt + theme(strip.background = element_blank(), legend.key.height = unit(1,"cm"))
# #plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, desc(R)) , scales = "free", ncol=4)
# #plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Neither'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
# #plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
# #plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
# #plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery & Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
# plt <- plt + ggplot2::geom_smooth( method = "lm", se = FALSE,  linetype = 'longdash',size = 0.5) #+ stat_poly_eq()
# #plt <- plt + ggplot2::geom_smooth(df_plt, ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2), method = "lm", se = FALSE, colour="black",  linetype = 'longdash',size = 0.5) #+ stat_poly_eq()
# 
# plt <- plt + ggplot2::labs(x = "Discovery log2 fold change",y = "Replication log2 fold change",color = "") #FDR < 0.05
# 
# plt
# ggsave(plt, filename=paste("dge_fc_compare_lm.png"), width=8, height =6)
# 
# 
# 
# 
# 
# 
# 
# 
# 
