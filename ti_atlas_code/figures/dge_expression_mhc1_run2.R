
#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/00055-study_overview/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")


library(ggplot2)
f="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/absolute_expression_mhc12.csv"
df6=read.csv(f, sep=",", header=TRUE, stringsAsFactors = FALSE)
df6$category <- factor(df6$category , levels=c("Stem cells", "Enterocyte", "Secretory",  "Myeloid", "B Cell", "T Cell", "Mesenchymal", "B Cell plasma")) 
df6$disease_status <- factor(df6$disease_status, levels=c("cd", "healthy"), labels = c("Crohn's disease", "Healthy"))

df7=df6
df7=df7[!df7$cohort=="Full",]
#df7=df7[df7$disease_status=="Crohn's disease",]


mhc1_classical=c("HLA-A", "HLA-B", "HLA-C")
#mhc1_classical=c("HLA-A",  "HLA-C")
#mhc1_classical=c("HLA-B") #"HLA-A", "HLA-B", 
mhc1_nonclassical=c("HLA-E", "HLA-F","HLA-G")
#mhc1_non_class= c("MICA", "MICB", "ULBP1", "ULBP2", "CD1A","CD1B","CD1C","CD1D","CD1E", "HFE", "MR1", "PROCR", "AZGP1", "UL18", "UL142", "UL37")
mhc2=c("HLA-DRA","HLA-DRB1",  "HLA-DPA1","HLA-DMA","HLA-DMB", "HLA-DQA1", "HLA-DRB5", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
       "HLA-DOA", "HLA-DOB", "HLA-DPB1")



#df7=df7[(df7$gene_symbol %in% mhc2),]


# df66 <- df7 %>% group_by(disease_status, gene_symbol, cohort, category,label, sex, age_imputed, inflammation_status, sanger_sample_id) %>% dplyr::summarize(mean_log1p_cp10k=mean(log1p_cp10k), 
#                                                                                                                                                             mean_apoptosis_score=mean(apoptosis_score),
#                                                                                                                                                             mean_interferon_gamma_score=mean(interferon_gamma_score),
#                                                                                                                                                             mean_interferon_alpha_beta_score=mean(interferon_alpha_beta_score),
#                                                                                                                                                             mean_costim_receptor_mean=mean(costim_receptor_mean),
#                                                                                                                                                             mean_costim_ligand_mean=mean(costim_ligand_mean),
#                                                                                                                                                             mean_mitochondrial_mean=mean(mitochondrial_mean))
# 


df7=df7[,c("disease_status", "sanger_sample_id", "gene_symbol","log1p_cp10k", "label", "category", "cohort")]
df66 <- df7 %>% group_by(disease_status, gene_symbol, cohort, category, label, sanger_sample_id) %>% dplyr::summarize(mean_log1p_cp10k=mean(log1p_cp10k))

df66[df66$gene_symbol %in% mhc1_classical, "type"] <- "mhc1_classical"
df66[df66$gene_symbol %in% mhc1_nonclassical, "type"] <- "mhc1_nonclassical"
df66[df66$gene_symbol %in% mhc2, "type"] <- "mhc2"


df7=df66[(df66$gene_symbol %in% mhc1_classical | df66$gene_symbol %in% mhc1_nonclassical),] # | df66$gene_symbol %in% mhc2

df7=df7[df7$disease_status=="Healthy",]
df7=df7[df7$category %in% c("Stem cells", "Enterocyte", "Secretory", "Myeloid"),]

#Add significance level
p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/"
p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/"
f1=paste0(p1, "disease_status_merged-de_results.tsv")
f4=paste0(p4, "disease_status_merged-de_results.tsv")
df1=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df4=read.table(f4, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df1$cohort <- "Discovery"
df4$cohort <- "Replication"
df1=df1[, c("gene_symbol", "label", "qvalue_bh_allcelltypes",  "n_cells", "log2fc", "mean_cp10k", "mean_counts",  "cohort", "category")] #"n_cells",
df4=df4[, c("gene_symbol", "label", "qvalue_bh_allcelltypes", "n_cells", "log2fc", "mean_cp10k", "mean_counts",  "cohort", "category")] #"n_cells"
df=rbind(df1, df4) 
df=df[df$category %in% c("Stem cells", "Enterocyte", "Secretory",  "Myeloid"),]
df$neg_log10 <- -log10(df$qvalue_bh_allcelltypes)
df$significant <- df$qvalue_bh_allcelltypes < 0.05
df$significant <- factor(df$significant, levels = c(TRUE, FALSE), labels = c("*", ""))
df=df[!df$cohort=="Full",]
df_sig=df
#df_sig=df[df$significant=="True",]
#df_sig=df_sig[grep("HLA-",df_sig$gene_symbol),]
df_sig=df_sig[df_sig$gene_symbol %in% c("HLA-A", "HLA-B", "HLA-C","HLA-E", "HLA-F","HLA-G"),]
df_sig=df_sig[,c("gene_symbol", "label", "cohort", "significant")]


df8=merge(df7, df_sig, by=c("gene_symbol","label", "cohort"), all=TRUE)
#df8[is.na(df8$significant),] <- "not tested"
#ll=df8[,c("gene_symbol", "label", "cohort", "category", "significant")]
df9=unique(df8[,c("gene_symbol", "label", "cohort", "category", "significant")])

df7[df7$category %in% c("Stem cells", "Enterocyte", "Secretory"), "compartment"] <- "Epithelial"
df7[df7$category %in% c("Myeloid"), "compartment"] <- "Myeloid"

my_comparisons <- list( c("mild", "moderate"), c("mild", "severe"))
palette=c("#003049", "#669bbc")


file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
df$label= as.character(df$label)
df=df[df$category %in% c("Stem cells", "Enterocyte", "Secretory", "Myeloid"),]
df7$label <- factor(df7$label, levels=oo$label)  
df9$label <- factor(df9$label, levels=oo$label)  

df7[df7$cohort =="Discovery","cohort"] <- "Discovery (Healthy individuals)"
df9[df9$cohort =="Discovery","cohort"] <- "Discovery (Healthy individuals)"

df7[df7$cohort =="Replication","cohort"] <- "Replication (Healthy individuals)"
df9[df9$cohort =="Replication","cohort"] <- "Replication (Healthy individuals)"


#df7=df7[df7$cohort=="Discovery",]
#df9=df9[df9$cohort=="Discovery (Healthy patients)",]
#df8=df8[df8$significant=="*",]
#df10=unique(df8[,c("gene_symbol","category", "cohort", "significant")])
#df10=df10[!is.na(df10$gene_symbol),]

gplt = ggplot(df7, aes(x=category, y=mean_log1p_cp10k, fill=gene_symbol))  + facet_nested(~cohort+compartment, scales = 'free', space = 'free') 
#gplt = gplt + geom_text(data = df10, aes(x=category, y=5, label = significant)) # position = position_dodge(width = .75)
#+ geom_jitter(cex=0.3)
gplt = gplt + theme_bw() + scale_fill_brewer() + ylab("Mean Expression") + theme() + geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75))
gplt = gplt 
#gplt = gplt + facet_grid(~compartment, scales = 'free')
#gplt = gplt + geom_smooth(aes(group=gene_symbol), colour="black", method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5, fullrange=TRUE) #+ stat_regline_equation(label.y = 4, aes(label = ..rr.label..), colour='black')  #+ stat_poly_eq()#+ geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) #+ scale_fill_manual(values=palette) 
gplt = gplt + theme( #axis.title.x = element_text(size=10),
  strip.text.x = element_text(size = 15),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 13), 
  #strip.text.x = element_text(size = 17), 
  axis.text.y = element_text(size = 13), 
  axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust=1),
  strip.background = element_blank(),
  legend.title = element_blank(),
  plot.title = element_text(size = 17)) #+ theme(axis_text_x=element_text(angle=90))
#gplt = gplt + ggtitle("Discovery") #+ stat_compare_means(aes(label = ..p.signif..), method = "t.test",comparisons = my_comparisons) #+ ggtitle("MHC-I classical (HLA-A,B,C)")
gplt

ggsave(plot=gplt, filename="absolute_expression_mhc1.png", width = 10,height = 3) 
ggsave(plot=gplt, filename="absolute_expression_mhc1.pdf", width = 10,height = 3, device = cairo_pdf, dpi=200)





# 
# 
# 
# 
# 
# df8=df66[(df66$gene_symbol %in% mhc2),] # | df66$gene_symbol %in% mhc2
# df8=df8[df8$gene_symbol %in% c("HLA-DRA", "HLA-DRB1","HLA-DPA1"),]
# df8$gene_symbol <- factor(df8$gene_symbol, levels=c("HLA-DRA", "HLA-DRB1","HLA-DPA1"))
# df8=df8[df8$category %in% c("Enterocyte","Secretory","Stem cells", "Myeloid", "B Cell"),]
# 
# my_comparisons <- list( c("mild", "moderate"), c("mild", "severe"))
# palette=c("#003049", "#669bbc")
# gplt = ggplot(df8, aes(x=reorder(category, desc(mean_log1p_cp10k)), y=mean_log1p_cp10k, fill=gene_symbol)) 
# #+ geom_jitter(cex=0.3)
# gplt = gplt + theme_bw() + scale_fill_brewer() + xlab("Mean Apoptosis score") + ylab("Mean Expression") + theme() + geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75))
# #gplt = gplt + geom_smooth(aes(group=gene_symbol), colour="black", method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5, fullrange=TRUE) #+ stat_regline_equation(label.y = 4, aes(label = ..rr.label..), colour='black')  #+ stat_poly_eq()#+ geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) #+ scale_fill_manual(values=palette) 
# gplt = gplt + theme( #axis.title.x = element_text(size=10),
#   strip.text.x = element_text(size = 15),
#   axis.title.x = element_blank(),
#   axis.title.y = element_text(size = 13), 
#   #strip.text.x = element_text(size = 17), 
#   axis.text.y = element_text(size = 13), 
#   axis.text.x = element_text(size = 13, angle = 90, vjust = 1, hjust=1),
#   strip.background = element_blank(),
#   legend.title = element_blank(),
#   plot.title = element_text(size = 17)) + facet_wrap(~disease_status, scales = 'free_x', nrow=1) #+ theme(axis_text_x=element_text(angle=90))
# gplt = gplt
# 
# #+ ggtitle("MHC-I classical and non-classical") #+ stat_compare_means(aes(label = ..p.signif..), method = "t.test",comparisons = my_comparisons) #+ ggtitle("MHC-I classical (HLA-A,B,C)")
# gplt
# 
# 
# 
# 
# 
# # Correlation with apoptosis
# my_comparisons <- list( c("mild", "moderate"), c("mild", "severe"))
# palette=c("#003049", "#669bbc")
# gplt = ggplot(df66, aes(y=mean_log1p_cp10k, x=mean_apoptosis_score, colour=gene_symbol)) #+ geom_jitter(cex=0.3)
# gplt = gplt + theme_bw() + scale_fill_brewer() + xlab("Mean Apoptosis score") + ylab("Mean Expression") + theme() + geom_point() 
# gplt = gplt + geom_smooth(aes(group=gene_symbol), colour="black", method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5, fullrange=TRUE) #+ stat_regline_equation(label.y = 4, aes(label = ..rr.label..), colour='black')  #+ stat_poly_eq()#+ geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) #+ scale_fill_manual(values=palette) 
# gplt = gplt + theme( axis.title.x = element_text(size=10),
#                      strip.text.x = element_text(size = 15),
#                      #axis.title.x = element_blank(),
#                      #strip.text.x = element_text(size = 17), 
#                      axis.text.y = element_text(size = 15), 
#                      strip.background = element_blank(),
#                      #legend.position = "none",
#                      plot.title = element_text(size = 17)) + facet_wrap(~disease_status+category, scales = 'free_x', nrow=2) #+ theme(axis_text_x=element_text(angle=90))
# gplt = gplt #+ ggtitle("MHC-I classical and non-classical") #+ stat_compare_means(aes(label = ..p.signif..), method = "t.test",comparisons = my_comparisons) #+ ggtitle("MHC-I classical (HLA-A,B,C)")
# gplt
# 
# 
# 
# # Correlation with mitochondrial
# my_comparisons <- list( c("mild", "moderate"), c("mild", "severe"))
# palette=c("#003049", "#669bbc")
# gplt = ggplot(df7, aes(y=mean_log1p_cp10k, x=mean_mitochondrial_mean, colour=gene_symbol)) #+ geom_jitter(cex=0.3)
# gplt = gplt + theme_bw() + scale_fill_brewer() + xlab("Mean mitochondrial score") + ylab("Mean Expression") + theme() + geom_point() 
# gplt = gplt + geom_smooth(aes(group=gene_symbol), colour="black", method = "lm",se = FALSE, linetype = 'longdash',alpha = 0.5, fullrange=TRUE)
# gplt = gplt + ggplot2::geom_abline(slope = 1, intercept = 0, colour="red", linetype = "dashed")#+ stat_regline_equation(label.y = 4, aes(label = ..rr.label..), colour='black')  #+ stat_poly_eq()#+ geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) #+ scale_fill_manual(values=palette) 
# gplt = gplt + theme( axis.title.x = element_text(size=10),
#                      strip.text.x = element_text(size = 15),
#                      #axis.title.x = element_blank(),
#                      #strip.text.x = element_text(size = 17), 
#                      axis.text.y = element_text(size = 15), 
#                      strip.background = element_blank(),
#                      #legend.position = "none",
#                      plot.title = element_text(size = 17)) + facet_wrap(~disease_status+category, scales = 'free_x', nrow=2) #+ theme(axis_text_x=element_text(angle=90))
# gplt = gplt #+ ggtitle("MHC-I classical and non-classical") #+ stat_compare_means(aes(label = ..p.signif..), method = "t.test",comparisons = my_comparisons) #+ ggtitle("MHC-I classical (HLA-A,B,C)")
# gplt
# 




