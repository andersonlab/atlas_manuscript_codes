
#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/00055-study_overview/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/00055-study_overview/")


library(ggplot2)
f1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/cell_proportions_discovery.csv"
f2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/cell_proportions_replication.csv"
f3="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/cell_proportions_full.csv"


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


df=df[,c("sanger_sample_id", "disease_status", "nr_cells", "freeze", "label", "category")]
df=df[df$label %in% c("Monocytes", "Mac intermediate (2)"),]


df1=df[df$freeze=="Discovery",]
df2=df[df$freeze=="Replication",]


check <- df %>% group_by(freeze, disease_status, label) %>% summarise(mean_nr_cells = mean(nr_cells))

lala <- check %>% group_by(freeze, label) %>% dplyr::summarise(fc  = round(mean_nr_cells[disease_status == 'cd'] / mean_nr_cells[disease_status == 'healthy'],0) )

lala=lala[lala$freeze!="Full",]
lala
