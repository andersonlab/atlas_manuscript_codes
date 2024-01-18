f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/umap_gene.csv"
df1=read.csv(f1,  sep=",")


df1=df1[df1$celltype_category %in% c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal"),]
labels=df1 %>% group_by(celltype_category) %>% dplyr::summarise(UMAP1=mean(UMAP1), UMAP2=mean(UMAP2))
labels2=df1 %>% group_by(celltype_label) %>% dplyr::summarise(UMAP1=mean(UMAP1), UMAP2=mean(UMAP2)) %>% ungroup()
subset_me=c("Endothelial cell", "Enterocyte middle villus (1)","Enterocyte middle villus (2)", "Enterocyte top villus","Goblet cell top villus")
labels2=labels2[labels2$celltype_label %in%  subset_me,]


max(df1$UMAP1)
min(df1$UMAP1)
max(df1$UMAP2)

table(df1$UMAP2 > 8)
table(df1$UMAP1 < -6)

df1=df1[!df1$UMAP2 > 6,]
df1=df1[!df1$UMAP1 < -6,]

df1_more=df1[df1$SLCO2A1_log1p_cp10k >= 1,]
df1_less=df1[df1$SLCO2A1_log1p_cp10k < 1,]

max(df1$SLCO2A1_log1p_cp10k)


p = ggplot(df1_less, aes(x=UMAP1, y=UMAP2, color=SLCO2A1_log1p_cp10k))
p = p + scale_colour_gradient2(low = "snow", mid = "snow3", high = "red", midpoint=0) #, breaks = c(-3,-1,0,1,3)
p = p + geom_point(size = 0.1) #+ plt9.scale_color_manual(values=palette)
p = p + theme_void()
p = p + geom_point(data=df1_more, aes(x=UMAP1, y=UMAP2, color=SLCO2A1_log1p_cp10k), size = 0.8) #+ plt9.scale_color_manual(values=palette)
p = p + labs(color="SLCO2A1")
p = p + geom_text(data=labels, mapping=aes(x=UMAP1 , y=UMAP2, label=celltype_category ) , color='black', size=7)
p = p + geom_text(data=labels2, mapping=aes(x=UMAP1 , y=UMAP2-1, label=celltype_label ) , color='black', size=4)
p

#ggsave(p, filename=paste("umap_gene.png"), width=9, height =6)


check1 <- df1_more %>% group_by(celltype_label) %>% dplyr::summarize(sum=length(celltype_label), mean=mean(SLCO2A1_log1p_cp10k))
check1$gene_symbols="SLCO2A1"

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
oo=oo[oo$category %in% c("Stem cells", "Secretory", "Enterocyte", "Mesenchymal"),]
names(oo)[3] <- "celltype_label"
check1=merge(check1, oo, by="celltype_label")

min_use=min(check1$mean)
check1$category__machine <- NULL
check1$label__machine <- NULL
check1$cluster <- NULL

ll=data.frame(celltype_label=oo[!oo$celltype_label %in% check1$celltype_label,]$celltype_label,
              category=oo[!oo$celltype_label %in% check1$celltype_label,]$category)
ll$sum=0
ll$mean=0
ll$gene_symbols="SLCO2A1"          

check1=rbind(check1, ll[,names(check1)])   
check1$celltype_label <- factor(check1$celltype_label, levels=rev(oo$celltype_label))
check1$category <- factor(check1$category, levels=c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal"))

p1 = ggplot(check1, aes(x=gene_symbols, y=celltype_label, colour=mean, size=sum))
p1 = p1 + scale_colour_gradient2(low = "snow", mid = "snow3", high = "red", midpoint=0) + theme_bw()  #, breaks = c(-3,-1,0,1,3)
p1 = p1 + facet_grid(vars(category), scales="free", space="free")
p1 = p1 + theme(axis.title.x = element_blank(), 
                axis.title.y = element_blank(),
                axis.text.x=element_text(size=10),
                axis.text.y=element_text(size=14),
                strip.background = element_blank(),
                strip.text = element_text(size=12))
p1 = p1 + labs(colour="Mean expression", size="Number of cells")
p1 = p1 + geom_point() #+ plt9.scale_color_manual(values=palette)
p1


all <- egg::ggarrange(plots=list(p,ggplot() + theme_void(), p1), 
                      widths = c(3, 0.2, 0.3), 
                      heigths=c(1, 0.2, 0.1),
                      ncol=3, nrow=1)

ggsave(plot=all, filename=paste("umap_SLCO2A1.png"), width=15, height =6, dpi=200)
#ggsave(plot=all, filename=paste("umap_SLCO2A1.pdf"), width=15, height =6, device = cairo_pdf, dpi=200)


# #all = grid::grid.draw()
# #grid::grid.draw(all)
# #class(all)
# 
# #f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/absolute_expression_SLCO2A1.csv"
# f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/tables_helper/mean_expression-both_cohorts.csv"
# df2=read.csv(f2,  sep=",")
# 
# f3='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
# df3=read.table(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)
# keepme=c("label", "label__machine_retired", "category")
# df3=df3[,keepme]
# names(df3)[names(df3)=="label__machine_retired"] <- "cell_type"
# df2=merge(df2, df3, by = "cell_type", all.x = TRUE)
# df2$cell_type <- NULL
# 
# df2=df2[df2$gene_symbols=="SLCO2A1",]
# df2=df2[df2$category %in% c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal"),]
# 
# max(df2$mean_log1p_cp10k)
# 
# df4=df1
# df4$gene_symbols="SLCO2A1"
# 
# 
# p1 = ggplot(df2, aes(x=gene_symbols, y=label, colour=mean_log1p_cp10k, size=nonzero_counts))
# p1 = p1 + scale_colour_gradient2(low = "snow3", mid = "snow3", high = "red", midpoint=0.1) #, breaks = c(-3,-1,0,1,3)
# p1 = p1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
# p1 = p1 + labs(fill="log1p_cp10k", size="nonzero_counts")
# p1 = p1 + geom_point() #+ plt9.scale_color_manual(values=palette)
# p1
# 
# 





