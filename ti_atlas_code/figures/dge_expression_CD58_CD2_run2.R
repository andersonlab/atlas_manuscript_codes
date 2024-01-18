setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/absolute_expression_SLCO2A1.csv"
#f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/absolute_expression_CD58_CD2.csv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/umap_cd58_cd2.csv"
df2=read.csv(f2,  sep=",")
names(df2)
#df2=df2[,c("celltype_category", "celltype_label" ,"CD58_log1p_cp10k", "CD2_log1p_cp10k" )]

df3=melt(df2, id.vars=c("UMAP1", "UMAP2", "celltype_category", "celltype_label"))

df4=df3[df3$variable=="CD58_log1p_cp10k",]
df4_more=df4[df4$value >= 1,]
df4_less=df4[df4$value < 1,]

df5=df3[df3$variable=="CD2_log1p_cp10k",]
df5_more=df5[df5$value >= 1,]
df5_less=df5[df5$value < 1,]


labels4=df4 %>% group_by(celltype_category) %>% dplyr::summarise(UMAP1=mean(UMAP1), UMAP2=mean(UMAP2))
labels5=df5 %>% group_by(celltype_category) %>% dplyr::summarise(UMAP1=mean(UMAP1), UMAP2=mean(UMAP2))
#subset_me=c("Endothelial cell", "Enterocyte middle villus (1)","Enterocyte middle villus (2)", "Enterocyte top villus","Goblet cell top villus")
#labels2=labels2[labels2$celltype_label %in%  subset_me,]


p = ggplot(df4_less, aes(x=UMAP1, y=UMAP2, color=value))
p = p + scale_colour_gradient2(low = "snow", mid = "snow3", high = "red", midpoint=0) #, breaks = c(-3,-1,0,1,3)
p = p + geom_point(size = 0.1) #+ plt9.scale_color_manual(values=palette)
p = p + theme_void()
p = p + geom_point(data=df4_more, aes(x=UMAP1, y=UMAP2, color=value), size = 0.8) #+ plt9.scale_color_manual(values=palette)
p = p + labs(color="CD58")
p = p + geom_text(data=labels4, mapping=aes(x=UMAP1 , y=UMAP2, label=celltype_category ) , color='black', size=7)
#p = p + geom_text(data=labels2, mapping=aes(x=UMAP1 , y=UMAP2-1, label=celltype_label ) , color='black', size=4)
p


p1 = ggplot(df5_less, aes(x=UMAP1, y=UMAP2, color=value))
p1 = p1 + scale_colour_gradient2(low = "snow", mid = "snow3", high = "red", midpoint=0) #, breaks = c(-3,-1,0,1,3)
p1 = p1 + geom_point(size = 0.1) #+ plt9.scale_color_manual(values=palette)
p1 = p1 + theme_void()
p1 = p1 + geom_point(data=df5_more, aes(x=UMAP1, y=UMAP2, color=value), size = 0.8) #+ plt9.scale_color_manual(values=palette)
p1 = p1 + labs(color="CD2")
p1 = p1 + geom_text(data=labels5, mapping=aes(x=UMAP1 , y=UMAP2, label=celltype_category ) , color='black', size=7)
#p = p + geom_text(data=labels2, mapping=aes(x=UMAP1 , y=UMAP2-1, label=celltype_label ) , color='black', size=4)
p1


file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
#oo=oo[oo$category %in% c("Stem cells", "Secretory", "Enterocyte", "T Cell"),]
names(oo)[3] <- "celltype_label"


df3_more=df3[df3$value >= 1,]
df3_less=df3[df3$value < 1,]

check1 <- df3_more %>% group_by(celltype_label, variable) %>% dplyr::summarize(sum=length(variable), mean=mean(value))
check1=merge(check1, oo, by="celltype_label")
check1$celltype_label <- factor(check1$celltype_label, levels=rev(oo$celltype_label))
check1$category <- factor(check1$category, levels=c("Stem cells", "Enterocyte", "Secretory",  "Myeloid",  "B Cell","T Cell","Mesenchymal", "B Cell plasma"))
min_use=min(check1$mean)
check1$variable <- factor(check1$variable , levels=c("CD58_log1p_cp10k", "CD2_log1p_cp10k"), label=c("CD58", "CD2"))

p2 = ggplot(check1, aes(x=variable, y=celltype_label, colour=mean, size=sum))
p2 = p2 + scale_colour_gradient2(low = "snow", mid = "snow3", high = "red", midpoint=1) + theme_bw()  #, breaks = c(-3,-1,0,1,3)
p2 = p2 + facet_grid(vars(category), scales="free", space="free")
p2 = p2 + theme(axis.title.x = element_blank(), 
                axis.title.y = element_blank(),
                axis.text.x=element_text(size=10, angle=90),
                axis.text.y=element_text(size=10),
                strip.text.y = element_text(size = 19, angle=360, hjust=0),
                strip.background = element_blank(),
                strip.text = element_text(size=12))
p2 = p2 + labs(colour="Mean expression", size="Number of cells")
p2 = p2 + geom_point() #+ plt9.scale_color_manual(values=palette)
p2
#check1$variable

all <- egg::ggarrange(plots=list(p,
                                 ggplot() + theme_void(), 
                                 p1, 
                                 ggplot() + theme_void(),
                                 p2), 
                      widths = c(3, 0.2, 3, 0.2, 0.3), 
                      heigths=c(1, 0.2,1, 0.2, 0.1),
                      ncol=5, nrow=1)

all
ggsave(plot=all, filename=paste("umap_CD2_CD58.png"), width=22, height =8, dpi=200)
#ggsave(plot=all, filename=paste("umap_CD2_CD58.pdf"), width=22, height =8, dpi=300, device = cairo_pdf)

#df7=df7[,c("disease_status", "sanger_sample_id", "gene_symbol","log1p_cp10k", "label", "category", "cohort")]
#df66 <- df7 %>% group_by(disease_status, gene_symbol, cohort, category, label, sanger_sample_id) %>% dplyr::summarize(mean_log1p_cp10k=mean(log1p_cp10k))





gplt = ggplot(df2, aes(x=category, y=log1p_cp10k))  + facet_nested(~gene_symbol, scales = 'free', space = 'free') 
#gplt = gplt + geom_text(data = df10, aes(x=category, y=5, label = significant)) # position = position_dodge(width = .75)
#+ geom_jitter(cex=0.3)
gplt = gplt + theme_bw() + scale_fill_brewer() + xlab("Mean Apoptosis score") + ylab("Mean Expression") + theme() + geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75))
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







