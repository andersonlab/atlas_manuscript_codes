library(dplyr)
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

#f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/disease_status_merged-de_results.tsv"
#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/disease_status_merged-de_results.tsv"
#f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-all_samples/disease_status_merged-de_results.tsv"

f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/disease_status_merged-de_results.tsv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/disease_status_merged-de_results.tsv"


df=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df1=read.table(f2, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df$freeze <- "Discovery"
df1$freeze <- "Replication"
df$comparison <- "Discovery"
df1$comparison <- "Replication"


df_final=rbind(df, df1)

plot(sort(df_final$test_statistic[df_final$test_statistic>-50]))
summary(df_final$test_statistic[df_final$test_statistic>-50])

# df_final$hit_discovery <- df_final$qvalue_bh_allcelltypes < 0.05 & df_final$freeze == "Discovery"
# df_final$hit_replication <- df_final$qvalue_bh_allcelltypes < 0.05 & df_final$freeze == "Replication"

df_final=df_final[!is.na(df_final$qvalue_bh_allcelltypes),]
table(df_final$qvalue_bh_allcelltypes < 0.05)


#df_final[df_final$qvalue_bh_allcelltypes < 0.05 & df_final$freeze == "Discovery", "hit"] <- "hit_discovery"
#df_final[df_final$qvalue_bh_allcelltypes < 0.05 & df_final$freeze == "Replication", "hit"] <- "hit_replication"


#keep=c("gene_symbol", "label","category", "freeze", "log2fc", "hit_discovery", "hit_replication")
#df_final=df_final[,keep]
# 
# df_final_melted= reshape2::melt(df_final, id.vars=c("gene_symbol", "label", "category", "freeze", "log2fc"), measure.vars=c("hit_discovery", "hit_replication" ))
# 
# df_final_melted$hit_discovery <- df_final_melted$variable == "hit_discovery" & df_final_melted$value == TRUE 
# df_final_melted$hit_replication <- df_final_melted$variable == "hit_replication" & df_final_melted$value == TRUE 
# 
# df_final_melted$value <- NULL
# df_final_melted$variable <- NULL
#df_final_melted[df_final_melted$value == FALSE, "variable"] <- NA

df_plt <- reshape2::dcast(df_final, gene_symbol + label + category ~ freeze , value.var = "log2fc")
df_plt1 <- reshape2::dcast(df_final, gene_symbol + label + category ~ freeze , value.var = "qvalue_bh_allcelltypes")
df_plt2 <- reshape2::dcast(df_final, gene_symbol + label + category ~ freeze , value.var = "n_cells")


#table(df_plt1$Discovery < 0.05)
#table(df_plt1$Replication< 0.05)

df_plt$hit_discovery <- df_plt1$Discovery < 0.05
df_plt$hit_replication <- df_plt1$Replication < 0.05
df_plt$n_cells_discovery <- df_plt2$Discovery 
df_plt$n_cells_replication <- df_plt2$Replication 

#for(i in df_plt$label){
#  hh=df_plt[df_plt$label==i,]
#  print(length(unique(hh$gene_symbol)))
#}

# Genes in common
dis=df_plt[!is.na(df_plt$hit_both),]
dis=dis[dis$hit_both==TRUE,]
dis=dis[!duplicated(dis$gene_symbol),]
length(unique(dis[,"gene_symbol"]))


df_plt$hit_both <- (df_plt$hit_discovery & df_plt$hit_replication)
df_plt$hit_type <- "Neither"
df_plt$hit_type[(df_plt$hit_discovery & !df_plt$hit_replication)] <- "Discovery"
df_plt$hit_type[(!df_plt$hit_discovery & df_plt$hit_replication)] <- "Replication"
df_plt$hit_type[(df_plt$hit_discovery & df_plt$hit_replication)] <- "Discovery & Replication"
df_plt$temp_asdfs_v1 <- df_plt$Discovery
df_plt$temp_asdfs_v2 <- df_plt$Replication

# drop points with missing data
df_plt <- subset(df_plt, !is.na(hit_both))

df_plt$category <- factor(df_plt$category , levels=c("Stem cells", "Enterocyte", "Secretory","B Cell plasma",  "B Cell", "Myeloid","T Cell", "Mesenchymal")) 
df_plt$hit_type <- factor(df_plt$hit_type, levels=c("Discovery & Replication", "Discovery", "Replication", "Neither"))


#df_plt=df_plt[df_plt$category %in%  c("Stem cells", "Enterocyte", "Secretory","Tuft cell") ,]

# check=df_plt[df_plt$category %in% c("Stem cells", "Enterocyte", "Secretory","Tuft cell"),]
# length(unique(check[check$hit_type == "Discovery & Replication",]$gene_symbol))
# length(unique(check[check$hit_type == "Discovery",]$gene_symbol))
# length(unique(check[check$hit_type == "Replication",]$gene_symbol))

# Check how highly variable genes replicate
# ff1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/repo/sc_nextflow/studies/gut-freeze003/data-variable_gene_filter.tsv"
# toremove=read.table(ff1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
# #toremove$gene_symbol
# #df_plt=df_plt[!(df_plt$gene_symbol %in% toremove$gene_symbol),]
# df_plt=df_plt[(df_plt$gene_symbol %in% toremove$gene_symbol),]
# df_plt=merge(df_plt, toremove, by="gene_symbol")
# 
# toremove=c("PIGR", "REG1A", "REG1B")
# df_plt=df_plt[(df_plt$gene_symbol %in% toremove),]

#####3

haa <- df_plt %>% group_by(label) %>% do(model = lm(temp_asdfs_v1 ~ temp_asdfs_v2, data = .))
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


#df_plt$category_title <- paste0(df_plt$category,": R=", round(df_plt$R, 2))
#mypalette=c("red",  "yellowgreen", "steelblue3","grey")
#names(mypalette) <- c("Discovery & Replication", "Discovery", "Replication", "Neither")


keepmee=c("label", "R")
df_plt_repel=df_plt[,keepmee]
df_plt_repel=distinct(df_plt_repel)


#mypalette=c("red",  "yellowgreen", "steelblue3","grey")
#names(mypalette) <- c("Discovery & Replication", "Discovery", "Replication", "Neither")

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/palette.tsv"
pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
palette=pal$palette_label
names(palette) <- pal$label


df_plt[df_plt$category %in% c("Stem cells", "Enterocyte", "Secretory"),"compartment"]  <- "Epithelial"
df_plt[df_plt$category %in% c("B Cell plasma", "B Cell", "Myeloid", "T Cell", "Mesenchymal"),"compartment"] <- "Immune and Mesenchymal"

#df_plt$n_cells_low <- (df_plt$n_cells_discovery < 900 & df_plt$n_cells_replication < 900)
#df_plt$n_cells_low = factor(df_plt$n_cells_low, levels=c(TRUE, FALSE), labels=c("n cells < 800", "n cells > 800"))
#df_plt9=df_plt[(df_plt$n_cells_discovery < 800 | df_plt$n_cells_replication < 800),]
#df_plt9=df_plt9[!is.na(df_plt9$label),]
#table(df_plt9$n_cells_low, df_plt9$label)

check=unique(df_plt[,c("label","category", "R")])
check <- check %>% group_by(category) %>% dplyr::summarise(R_mean=mean(R))
df_plt=merge(df_plt, check, by="category")
df_plt$category_title <- paste0(df_plt$category,": mean R=", round(df_plt$R_mean, 2))


library(ggplot2)
#mypalette = c("red", "steelblue")
#mypalette = c("deepskyblue3", "darkseagreen3")
plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, group=label, colour=label)) #, color=n_cells_low#x = temp_asdfs_v2, y =  temp_asdfs_v1
#plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, group=label, colour=gene_type)) #, color=n_cells_low#x = temp_asdfs_v2, y =  temp_asdfs_v1
plt <- plt + ggplot2::theme_bw(base_size = 12) + xlim(-1,1) + ylim(-1,1) 
plt <- plt + scale_colour_manual(values=palette)#+ ggplot2::scale_color_manual(values=mypalette)  #+ geom_point(colour="steelblue")
plt <- plt + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5, fullrange=TRUE) #+ stat_poly_eq()
#plt <- plt + ggplot2::geom_point()
plt <- plt + ggplot2::geom_abline(slope = 1, intercept = 0)
plt <- plt + ggplot2::geom_hline(yintercept = 0, linetype = 'dotted', alpha = 0.5)
plt <- plt + ggplot2::geom_vline(xintercept = 0, linetype = 'dotted', alpha = 0.5)
plt <- plt + theme(legend.position="none", legend.key.size = unit(2, 'cm'),
                                                                   legend.title = element_text(size=15),
                                                                   legend.text = element_text(size=15))
plt <- plt + ggplot2::labs(x = "Discovery log2 FC",y = "Replication log2 FC",color = "") #FDR < 0.05
plt <- plt + ggplot2::theme(strip.background = element_blank(), strip.text=element_text(size=11))
plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, dplyr::desc(R_mean)) , scales = "free", ncol=4)
#plt <- plt + ggplot2::facet_wrap(~ category , scales = "free", ncol=4)

#plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, desc(R_mean)) , scales = "free", ncol=4)
#plt <- plt + ggplot2::facet_wrap(~ label, scales = "free", ncol=5)
#plt <- plt + ggplot2::scale_color_brewer(palette = "Dark2")
plt

ggsave(plt, filename=paste("dge_regression_fold_change.png"), width=10, height =6)
ggsave(plt, filename=paste("dge_regression_fold_change.pdf"), width=10, height =6, device = cairo_pdf, dpi=200)


# Check only highly varible genes

library(ggplot2)
#mypalette = c("red", "steelblue")
#mypalette = c("deepskyblue3", "darkseagreen3")
#plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, group=label, colour=label)) #, color=n_cells_low#x = temp_asdfs_v2, y =  temp_asdfs_v1
plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, group=label, colour=gene_type)) #, color=n_cells_low#x = temp_asdfs_v2, y =  temp_asdfs_v1
plt <- plt + ggplot2::theme_bw(base_size = 12) + xlim(-1,1) + ylim(-1,1) 
#plt <- plt + scale_colour_manual(values=palette)#+ ggplot2::scale_color_manual(values=mypalette)  #+ geom_point(colour="steelblue")
plt <- plt + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5, fullrange=TRUE) #+ stat_poly_eq()
plt <- plt + ggplot2::geom_point()
plt <- plt + ggplot2::geom_abline(slope = 1, intercept = 0)
plt <- plt + ggplot2::geom_hline(yintercept = 0, linetype = 'dotted', alpha = 0.5)
plt <- plt + ggplot2::geom_vline(xintercept = 0, linetype = 'dotted', alpha = 0.5)
plt <- plt + theme(legend.position="bottom", legend.key.size = unit(2, 'cm'),
                   legend.title = element_text(size=15),
                   legend.text = element_text(size=15))
plt <- plt + ggplot2::labs(x = "Discovery log2 FC",y = "Replication log2 FC",color = "") #FDR < 0.05
plt <- plt + ggplot2::theme(strip.background = element_blank(), strip.text=element_text(size=11))
plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, dplyr::desc(R_mean)) , scales = "free", ncol=4)
#plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, desc(R_mean)) , scales = "free", ncol=4)
#plt <- plt + ggplot2::facet_wrap(~ label, scales = "free", ncol=5)
#plt <- plt + ggplot2::scale_color_brewer(palette = "Dark2")
plt

# Check differences between cell types


# Check differences between specificity in two cohorts

#df_plt11=df_plt[df_plt$Discovery > df_plt$Replication,]
df_plt11=df_plt
#df_plt11=df_plt11[df_plt11$category %in% c("T Cell"),]
#df_plt11=df_plt11[df_plt11$category %in% c("Myeloid", "Mesenchymal"),]
df_plt11$diff=abs(df_plt11$Discovery - df_plt11$Replication)
df_plt11=df_plt11[!is.na(df_plt11$Discovery),]
df_plt11=df_plt11[!is.na(df_plt11$Replication),]

check=df_plt11 %>% arrange(category, desc(diff)) %>% group_by(category) %>% top_n(10) %>% arrange(category, desc(diff))
check1=check[check$category %in% c("B Cell", "T Cell,", "Myeloid", "Mesenchymal"),]
check2=check[check$category %in% c("Enterocyte"),]





df_plt1=df_plt1[!is.na(df_plt1$Discovery),]
df_plt1=df_plt1[!is.na(df_plt1$Replication),]


df_plt11=df_plt1[df_plt1$Discovery > df_plt1$Replication,]
df_plt11$diff=df_plt11$Discovery - df_plt11$Replication

df_plt12=df_plt1[df_plt1$Discovery < df_plt1$Replication,]
unique(df_plt11$gene_symbol)










# g <- ggplot_gtable(ggplot_build(plt))
# strip_both <- which(grepl('strip-', g$layout$name))
# #colors = unique(pal$palette_category)
# colors <- c("#1F77B4", "#1F77B4", "#1F77B4", "#1F77B4", "#1F77B4", "#1F77B4", "#1F77B4", "#1F77B4")
# colors = c("deepskyblue3","deepskyblue3","deepskyblue3", "darkseagreen3",
#            "darkseagreen3","darkseagreen3","darkseagreen3","darkseagreen3")
# k <- 1
# for (i in strip_both) {
#   print(i)
#   j <- which(grepl("text", g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- colors[k]
#   k <- k+1
# }
# grid.draw(g)
# g <- ggplot_gtable(ggplot_build(plt))
# stripr <- which(grepl('strip-r', g$layout$name))
# fills <- unique(pal$palette_category)
# #names(palette) <- pal$category
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.draw(g)


# Plot the merged data ####################################################


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
plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = n_cells_discovery, y =  n_cells_replication, colour=category)) #x = temp_asdfs_v2, y =  temp_asdfs_v1
plt <- plt + geom_point(size=3) + scale_colour_manual(values=palette)
#plt <- plt + ggplot2::geom_smooth(method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5) #+ stat_poly_eq()
plt

# Plot the merged data ####################################################
mypalette = c("steelblue", "lightgrey")
plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, group=label, colour=compartment)) #x = temp_asdfs_v2, y =  temp_asdfs_v1
plt <- plt + ggplot2::theme_bw(base_size = 12) + scale_colour_manual(values=mypalette) #+ theme(legend.position='none')
plt <- plt + theme(strip.background = element_blank(), legend.key.height = unit(1,"cm"))
#plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, desc(R)) , scales = "free", ncol=4)
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Neither'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery & Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
plt <- plt + ggplot2::geom_smooth( method = "lm", se = FALSE,  linetype = 'longdash',size = 0.5) #+ stat_poly_eq()
#plt <- plt + ggplot2::geom_smooth(df_plt, ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2), method = "lm", se = FALSE, colour="black",  linetype = 'longdash',size = 0.5) #+ stat_poly_eq()

plt <- plt + ggplot2::labs(x = "Discovery log2 fold change",y = "Replication log2 fold change",color = "") #FDR < 0.05

plt
ggsave(plt, filename=paste("dge_fc_compare_lm.png"), width=8, height =6)


plt <- ggplot2::ggplot(df_plt, ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, group=label, color=label)) #x = temp_asdfs_v2, y =  temp_asdfs_v1
plt <- plt + ggplot2::theme_bw(base_size = 12) + scale_colour_manual(values=palette) + theme(legend.title = element_blank())
plt <- plt + theme(strip.background = element_blank())
#plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, desc(R)) , scales = "free", ncol=4)
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Neither'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery & Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
plt <- plt + ggplot2::geom_smooth( method = "lm", se = FALSE,  linetype = 'longdash',size = 1) #+ stat_poly_eq()
plt <- plt + ggplot2::labs(x = "Discovery log2 fold change",y = "Replication log2 fold change",color = "") #FDR < 0.05

#plt <- plt + ggrepel::geom_label_repel( aes(label=label))

plt
ggsave(plt, filename=paste("dge_fc_compare_lm_single.png"), width=8, height =6)




#plt <- plt + geom_text(data = df_plt_repel,
#plt <- plt + ggrepel::geom_label_repel(data = df_plt_repel, aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, label=as.factor(label)))

plt

plt <- plt + ggplot2::geom_abline(slope = 1, intercept = 0)
plt <- plt + ggplot2::geom_hline(yintercept = 0, linetype = 'dotted', alpha = 0.5)
plt <- plt + ggplot2::geom_vline(xintercept = 0, linetype = 'dotted', alpha = 0.5)
plt <- plt + ggplot2::scale_color_manual(values=mypalette)
plt <- plt + ggplot2::labs(x = "Discovery log2 fold change",y = "Replication log2 fold change",color = "") #FDR < 0.05
#plt <- plt + ggplot2::facet_wrap(~ reorder(category_title, desc(R)) , scales = "free", ncol=4)
#plt <- plt + ggplot2::facet_wrap(~ label, scales = "free", ncol=5)
#plt <- plt + ggplot2::scale_color_brewer(palette = "Dark2")
#plt <- plt + ggplot2::geom_point( ggplot2::aes(color = hit_type)) 
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Neither'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
#plt <- plt + geom_point(data = subset(df_plt, hit_type == 'Discovery & Replication'), ggplot2::aes(x = temp_asdfs_v1, y =  temp_asdfs_v2, color = hit_type))
plt

#mypalette=c("red",  "yellowgreen", "steelblue3","grey")
#mypalette=c("red",  "steelblue3", "yellowgreen", "grey")
a

plt
#ggsave(plt, filename=paste("dge_fc_compare.png"), width=12, height =5, dpi = 300)
ggsave(plt, filename=paste("dge_fc_compare.pdf"), width=12, height =6, device = cairo_pdf)













