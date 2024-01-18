
#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/00055-study_overview/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/00055-study_overview/")


library(ggplot2)
f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/cell_proportions_discovery.csv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/cell_proportions_replication.csv"
f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/cell_proportions_full.csv"


df11=read.csv(f1, sep=",", header=TRUE, stringsAsFactors = FALSE)
df22=read.csv(f2, sep=",", header=TRUE, stringsAsFactors = FALSE)
df33=read.csv(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)

names(df11)[3] <- "label__machine_retired"
names(df22)[3] <- "label__machine_retired"
names(df33)[3] <- "label__machine_retired"

df11$freeze <- "Discovery"
df22$freeze <- "Replication"
df33$freeze <- "Full"

df=rbind(df11, df22, df33)

f3='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
ww=read.table(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)
df=merge(df, ww, by = "label__machine_retired", all.x = TRUE)

#nr_cells_cd=sum(df[df$disease_status=="cd",]$nr_cells)
#nr_cells_h=sum(df[df$disease_status=="healthy",]$nr_cells)

#library(dplyr)
#df <- df %>% group_by(freeze, disease_status) %>% mutate(sum_nr_cells = sum(nr_cells)) %>% ungroup() 
#df$percentage = (df$nr_cells/df$sum_nr_cells)*100


df1 <- df %>% group_by(freeze, disease_status, label, category) %>% dplyr::summarize(sd = sd(nr_cells), mean=mean(nr_cells), percentage=mean(nr_cells)/sum(nr_cells))
df1 <- df %>% group_by(freeze, disease_status, label, category) %>% dplyr::summarize(sd = sd(nr_cells), mean=mean(nr_cells), percentage=mean(nr_cells)/sum(nr_cells))


# labels=c( "Stem cell LGR5+" ,  "Stem cell MKI67+ (1)" ,  "Stem cell MKI67+ (2)" ,  "Enterocyte precursor crypt OLFM4+ KRT20++" ,
#               "Enterocyte progenitor crypt OLFM4++ KRT20+ (1)",  "Enterocyte progenitor crypt OLFM4++ KRT20+ (2)" , "Enterocyte middle villus (1)" ,
#               "Enterocyte middle villus (2)" , "Enterocyte top villus" , "Enterocytes BEST4"  ,  "Paneth cell"  ,  "Goblet cell crypt MKI67+" ,  "Goblet cell middle villus" ,
#               "Goblet cell top villus" ,    "Endocrine cell" , "Tuft cell" ,    "Fibroblast/Myofibroblasts"  ,
#               "Endothelial cell"  ,  "Pericytes",  "Smooth muscle cell" ,  "Mac resident IL10RA+"  ,  "Mac resident IL10RA-" , "Dendritic cell"  ,
#               "MoMac IL10RA+", "MoMac IL10RA-" ,  "Monocytes",   "Mast"  ,  "ILC1 CD3D- NCAM1+" ,
#               "ILC3 CD3D- IL23R+" ,  "T cell CD4 CD40LG+ (1)" ,
#               "T cell CD4 CD40LG+ (2)" , "T cell CD4 CD40LG+ (3)" , "T cell CD4 Treg"  , "T cell CD4 naïve" ,  "T cell CD4 proliferating"  ,
#               "T cell CD4- CD8-", "T cell CD8 (1)",  "T cell CD8 (2)" ,  "T cell CD8 (3)" , "T cell gd",
#               "B cell" ,"B cell naïve" ,"B cell activated",
#               "B cell germinal centre/plasmablasts" ,  "B cell memory (1)"  ,  "B cell memory (2)"  ,  "B cell plasma IgA CD38+"  ,
#               "B cell plasma IgA CD38++" ,  "B cell plasma IgA CD38+++" )


file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

df1$label <- factor(df1$label , levels=rev(oo$label))   
df1$category <- factor(df1$category , levels=c("Stem cells", "Enterocyte", "Secretory",  "Tuft cell", "Mesenchymal",  "Myeloid",  "Mast", "T Cell", "B Cell", "B Cell plasma")) 
df1$disease_status <- factor(df1$disease_status, levels=c("cd", "healthy"), labels=c("Crohn's disease", "Healthy"))
df1$freeze <- factor(df1$freeze, levels=c("Discovery", "Replication", "Full"), 
                     labels=c("Discovery", "Replication", "Full cohort"))

palette=c("#003049", "#669bbc")


gplt = ggplot(df1, aes(x=label, y=mean,fill = disease_status)) + ylab("Number of cells")
gplt = gplt + ggplot2::theme_bw(base_size = 12)
gplt = gplt + theme(strip.background = element_blank())
gplt = gplt + facet_grid(category~freeze, scales = "free", space = "free")
gplt = gplt + geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=palette)#+ scale_y_reverse(breaks=c(-300, -200, -100, 0, 100,200,300), labels=c("300", "200","100", "0", "100", "200", "300"))
gplt = gplt + geom_errorbar(aes(x=label, ymin=mean, ymax=mean+sd),  width=.2, position=position_dodge(.9)) 
gplt = gplt + coord_flip() 
gplt = gplt + theme(axis.title.y = element_blank(), 
                    axis.title.x = element_text(size=15), 
                    strip.text.y = element_text(angle=360, size=12, hjust=0), 
                    strip.text.x = element_text(size = 15), 
                    legend.text=element_text(size=15), 
                    legend.title = element_blank())
gplt = gplt + scale_y_continuous(limits=c(0, 1300), breaks=c(0, 500,1000, 1300))
gplt


ggsave(plot=gplt, filename="cell_proportions_cohorts.png", width = 12,height = 12) 
ggsave(plot=gplt, filename="cell_proportions_cohorts.pdf", width = 12,height = 12, device = cairo_pdf, dpi=200)












plot(df$nr_cells)
########
library(rstatix)
library(ggsignif)
library(ggpubr)


library(rstatix)
stat.test <- df %>% group_by(freeze, label) %>% t_test(nr_cells ~ disease_status) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
stat.test <- stat.test %>% add_xy_position(x = "label", dodge = 0.1)
stat.test=stat.test[!stat.test$p.adj.signif=="ns",]

stat.test=stat.test[,c("freeze", "label", "group1", "group2", "p.adj.signif")]
stat.test.melted= reshape2::melt(stat.test, id.vars=c("freeze", "label", "p.adj.signif"), measure.vars=c("group1", "group2"))
stat.test.melted$variable <- NULL
names(stat.test.melted)[4] <- "p.adj.signif"

df3=merge(df2, stat.test, by=c("freeze", "label"))
df3=df3[df3$disease_status == "Healthy",]

gplt = ggplot(df2, aes(x=label, y=percentage,fill = disease_status)) + ylab("Mean percentage of total cells (%)")
gplt = gplt + ggplot2::theme_bw(base_size = 12)
gplt = gplt + theme(strip.background = element_blank())
gplt = gplt + facet_wrap(vars(freeze)) 
gplt = gplt + geom_bar(stat="identity", position="identity") + scale_fill_manual(values=palette)#+ scale_y_reverse(breaks=c(-300, -200, -100, 0, 100,200,300), labels=c("300", "200","100", "0", "100", "200", "300"))
gplt = gplt + geom_errorbar(aes(x=label, ymin=percentage, ymax=percentage+sd),  width=.3) 
gplt = gplt + coord_flip() + theme(axis.title.y = element_blank()) + theme(strip.text.x = element_text(size = 15), legend.title = element_blank())
gplt = gplt + geom_text(df3, mapping=aes(x=label, y=1.3, label = p.adj.signif, vjust=-0.3))
gplt

ggsave(plot=gplt, filename="cell_proportions_cohorts_test.png", width = 10,height = 8) 

#########

df3 <- reshape2::dcast(df1, category + label + freeze  ~ disease_status , value.var = "percentage")
df4=df[,c("freeze", "label","category", "disease_status", "sanger_sample_id", "nr_cells")]
df4=unique(df4)
df5 <- reshape2::dcast(df4, freeze + label +  sanger_sample_id ~ disease_status, value.var = "nr_cells")
df6 <- df4 %>% group_by(freeze, category, label, disease_status) %>% dplyr::summarize(sum_cells=sum(nr_cells))
df7 <- df6 %>% group_by(freeze, category,  label) %>% dplyr::summarize(fc=(sum_cells[disease_status=="cd"] + sum_cells[disease_status=="healthy"])/sum_cells[disease_status=="healthy"])

df7=df7[df7$category %in% c("Myeloid"),]

x=df7[df7$freeze=="Discovery",]
x=x[order(x$fc, decreasing = TRUE),]

df7$label<- factor(df7$label, levels=rev(x$label))

gplt1 = ggplot(df7, aes(x=label, y=fc)) 
gplt1 = gplt1 + facet_wrap(vars(freeze))
gplt1 = gplt1 + ggplot2::theme_bw(base_size = 12)
gplt1 = gplt1 + theme(strip.background = element_blank())
gplt1 = gplt1 + geom_bar(stat="identity", position="identity") + scale_fill_manual(values=palette)#+ scale_y_reverse(breaks=c(-300, -200, -100, 0, 100,200,300), labels=c("300", "200","100", "0", "100", "200", "300"))
gplt1 = gplt1 + coord_flip() + theme(axis.title.y = element_blank()) + theme(strip.text.x = element_text(size = 15),legend.title = element_blank()) + ylab("Log 2 fold change of mean percentage (%)")
gplt1


##########
df3 <- reshape2::dcast(df1, category + label + freeze  ~ disease_status , value.var = "percentage")
df3$fc=log2(df3$cd-df3$healthy/df3$healthy)
palette=c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')

gplt1 = ggplot(df3, aes(x=label, y=fc, fill=category)) 
gplt1 = gplt1 + facet_wrap(vars(freeze))
gplt1 = gplt1 + ggplot2::theme_bw(base_size = 12)
gplt1 = gplt1 + theme(strip.background = element_blank())
gplt1 = gplt1 + geom_bar(stat="identity", position="identity") + scale_fill_manual(values=palette)#+ scale_y_reverse(breaks=c(-300, -200, -100, 0, 100,200,300), labels=c("300", "200","100", "0", "100", "200", "300"))
gplt1 = gplt1 + coord_flip() + theme(axis.title.y = element_blank()) + theme(strip.text.x = element_text(size = 15),legend.title = element_blank()) + ylab("Log 2 fold change of mean percentage (%)")
gplt1
ggsave(plot=gplt1, filename="cell_proportions_cohorts_imbalance.png", width = 10,height = 8) 

#ggsave(plot=gplt, filename="cell_proportions_fc.png", width = 8,height = 6) 
#ggsave(plot=gplt, filename="cell_proportions_fc.pdf", width = 8,height = 6, device = cairo_pdf)

all <- gridExtra::grid.arrange(gplt, gplt1, ncol=2)
all

#ggsave(plot=all, filename="cell_proportions_fc_all.png", width = 16,height =7)
ggsave(plot=all, filename="cell_proportions_fc_all.pdf", width = 18,height = 7, device = cairo_pdf)


# stat.test$label <- factor(stat.test$label, levels=unique(stat.test$label))
# stat.test$label <- factor(stat.test$label, levels=unique(stat.test$label))
# 
# palette=c("#003049", "#669bbc")
# gplt = ggplot(df, aes(x=label, y=nr_cells)) + facet_grid(label , scales="free", space="free") + ylab("Mean percentage of total cells (%)")
# #gplt = gplt + facet_grid(major_label~, scales='free_x', space='free_x')
# gplt = gplt + geom_bar(stat="identity", position="dodge") #+ facet_wrap(vars(freeze))
# gplt = gplt + stat_pvalue_manual(stat.test,  label = "p.adj.signif", coord.flip = TRUE)
# gplt
# #gplt = gplt + scale_fill_manual(values=palette)#+ scale_y_reverse(breaks=c(-300, -200, -100, 0, 100,200,300), labels=c("300", "200","100", "0", "100", "200", "300"))
# gplt = gplt #+ geom_errorbar(aes(x=label, ymin=percentage, ymax=percentage+sd),  width=.2, position=position_dodge(.9)) 
# gplt = gplt + theme(axis.title.y = element_blank()) + theme(legend.title = element_blank()) + coord_flip()
# gplt = gplt 
# gplt


# + stat_pvalue_manual(stat.test,  label = "p.adj.signif") 
