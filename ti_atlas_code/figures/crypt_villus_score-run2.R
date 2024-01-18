
file="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/crypt_villus_score.csv"
#df=read.csv(file, sep=",", header=TRUE, stringsAsFactors = FALSE)
df=read.table(file, sep="\t", header=TRUE)

df=df[!df$celltype_label %in% c("Paneth cell", "Tuft cell", "Endocrine cell"),]

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/palette.tsv"
pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)



# col1=RColorBrewer::brewer.pal(length(unique(df[df['celltype_category']=="Enterocyte", "celltype_label"])), "Set1")
# names(col1)= unique(df[df['celltype_category']=="Enterocyte", "celltype_label"])
# 
# palette=c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
# 
# col1=palette[1:length(unique(df[df['celltype_category']=="Enterocyte", "celltype_label"]))]
# names(col1)= unique(df[df['celltype_category']=="Enterocyte", "celltype_label"])
# 
# #palette=c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
# 
# #col1=palette[2:4]
# col2=palette[1:length(unique(df[df['celltype_category']=="Secretory", "celltype_label"]))]
# names(col2)=sort(unique(df[df['celltype_category']=="Secretory", "celltype_label"]))
# 
# #col2=palette[1:6]
# #names(col2)=sort(unique(df[df['celltype_category']=="Enterocyte", "cell_type_label"]))

palette=pal$palette_label
names(palette) <- pal$label

#palette=c(col1, col2)

#col2=palette[1:length(unique(pull(df["celltype_label"])))]
#names(col2)=sort(unique(pull(df["celltype_label"])))
df1=df[df$celltype_category %in% c("Stem cells", "Enterocyte"),]
df1=df1[!(df1$celltype_label %in% c("Enterocytes BEST4")),]

gplt = ggplot(df1, aes(x=reorder(celltype_label, crypt_villus_score_moor_top_entero), y=crypt_villus_score_moor_top_entero, colour=celltype_label))
gplt = gplt + geom_jitter(cex=0.5) + ylab("Crypt - villus score") #+ facet_wrap(vars(celltype_category), scales='free')
gplt = gplt + geom_violin(alpha = 0.8)  #+ facet_wrap(vars(celltype_category), scales='free')
gplt = gplt + theme_bw() + coord_flip() + scale_colour_manual(values=palette) + ggtitle("Enterocyte")
gplt = gplt + theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), legend.position = "none", 
                    strip.background = element_blank(), text = element_text(size = 20)) + scale_y_continuous(limits=c(-0.8, 3), breaks=c(0,1,2,3))
gplt


df2=df[df$celltype_category %in% c("Secretory"),]

gplt1 = ggplot(df2, aes(x=reorder(celltype_label, crypt_villus_score_moor_top_goblet), y=crypt_villus_score_moor_top_entero, colour=celltype_label))
gplt1 = gplt1 + geom_jitter(cex=0.5) + ylab("Crypt - villus score") #+ facet_wrap(vars(celltype_category), scales='free')
gplt1 = gplt1 + geom_violin(alpha = 0.8)
gplt1 = gplt1 + theme_bw() + coord_flip() + scale_colour_manual(values=palette) + ggtitle("Secretory")
gplt1 = gplt1 + theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(),  legend.position = "none", 
                    strip.background = element_blank(), text = element_text(size = 20)) + scale_y_continuous(limits=c(-0.8, 2), breaks=c(0,1,2))
gplt1

g <- cowplot::plot_grid(gplt, gplt1, align = "hv")
g

library(sjPlot)          
#g <- gridExtra::grid.arrange(gplt, gplt1)
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/00055-study_overview/")
ggsave(plot=g, filename="crypt_villus_score.png", width = 16,height = 5, dpi=150) 
#ggsave(plot=g, filename="crypt_villus_score.svg", width = 16,height = 5, dpi=150)
#save_plot(plot=g, filename="crypt_villus_score.svg", width = 16,height = 5, dpi=150)



gplt = ggplot(df1, aes(x=reorder(celltype_label, crypt_villus_score_moor_top_entero), y=crypt_villus_score_moor_top_entero, colour=celltype_label))
gplt = gplt + geom_jitter(cex=0.5) + ylab("Crypt - villus score") #+ facet_wrap(vars(celltype_category), scales='free')
gplt = gplt + theme_bw() 
gplt = gplt + coord_flip() + scale_colour_manual(values=palette) + ggtitle("Enterocyte")
gplt = gplt + theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), legend.position = "none", 
                    strip.background = element_blank(), text = element_text(size = 20)) + scale_y_continuous(limits=c(-0.8, 3), breaks=c(0,1,2,3))
gplt = gplt #+ scale_y_discrete(position = "right")
gplt

#g <- gridExtra::grid.arrange(gplt, gplt1)
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/00055-study_overview/")

ggsave(plot=gplt, filename="crypt_villus_score_ichg.png", width = 16,height = 5) 




