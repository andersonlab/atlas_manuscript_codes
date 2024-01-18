
f="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/adata_dendrogram/cell_types_dendro.tsv"
dend=read.table(f, sep="\t", header=TRUE, stringsAsFactors = FALSE)


f3='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
ww=read.table(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)


p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Discovery/"
p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Replication/"
p2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/auto-annot-121/dge/MT_100-cohort_Full/"



f1=paste0(p1, "disease_status_merged-gsea_results.tsv")
f4=paste0(p4, "disease_status_merged-gsea_results.tsv")
f2=paste0(p2, "disease_status_merged-gsea_results.tsv")
df1=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df4=read.table(f4, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df2=read.table(f2, sep="\t", header=TRUE, stringsAsFactors = FALSE)

df1$freeze <- "Discovery"
df4$freeze <- "Replication"
df2$freeze <- "Full"
keepcol=c("label", "annot_id","signed_ranking","annot_coef", "qvalue_bh",  "ES", "freeze")
#keepcol=c("label", "annot_id","signed_ranking","annot_coef", "qvalue_bh",  "ES", "freeze", "leadingEdge")
df1=df1[,keepcol]
df4=df4[,keepcol]
df2=df2[,keepcol]
df = rbind(df1,df4, df2)
df=merge(df, ww, by = "label", all.x = TRUE)
df=df[grep("REACTOME", df$annot_id),]
df=df[df$signed_ranking=="True",]
df$minuslog10qval=-log10(df$qvalue_bh)
df$regulation <- ifelse(df$annot_coef >= 0, "up-regulated", "down-regulated")
df=df[!is.na(df$annot_id),]
df$annot_id=gsub("RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS_", "RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS", df$annot_id)
df$annot_id=gsub("REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS_", "REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR", df$annot_id)
df$annot_id=gsub("THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT", "CITRIC_ACID_CYCLE_RESPIRATORY_ELECTRON_TRANSPORT", df$annot_id)
df$annot_id=gsub("ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC", "FOLDING_ASSEMBLY_PEPTIDE_LOADING_OF_CLASS_I_MHC", df$annot_id)
df$annot_id=gsub("REACTOME_", "", df$annot_id)
df$significant <- df$qvalue_bh < 0.05

check=df[df$qvalue_bh < 0.05,]
check=check[check$freeze == "Discovery",]

keepme=c("label", "category", "annot_id", "significant",  "freeze", "qvalue_bh", "ES")
ll=df[,keepme]
ll <- reshape2::dcast(ll, annot_id + category + label ~ freeze , value.var = "significant")
ll$hit_type <- "Neither"
ll$hit_type[(ll$Discovery & ll$Replication)] <- "Discovery & Replication"
ll$Discovery <- NULL
ll$Replication <- NULL
ll$Full <- NULL
ll=ll[ll$hit_type=="Discovery & Replication",]

#####################
keepme=c("category", "label", "annot_id", "significant",  "freeze", "qvalue_bh", "ES")
ll1=df[,keepme]
ll1 <- reshape2::dcast(ll1, annot_id + category + label ~ freeze , value.var = "ES")
names(ll1)[names(ll1)=="Discovery"] <- "d_ES"
names(ll1)[names(ll1)=="Replication"] <- "r_ES"
names(ll1)[names(ll1)=="Full"] <- "f_ES"

ll2=df[,keepme]
ll2 <- reshape2::dcast(ll2, annot_id + category + label ~ freeze , value.var = "qvalue_bh")
names(ll2)[names(ll2)=="Discovery"] <- "d_qval"
names(ll2)[names(ll2)=="Replication"] <- "r_qval"
names(ll2)[names(ll2)=="Full"] <- "f_qval"

df_plt <- ll1 %>% right_join(ll, by=c("category", "label","annot_id"))
df_plt <- ll2 %>% right_join(df_plt, by=c("category", "label","annot_id"))

#############################

top.pathways=df_plt[df_plt$hit_type=="Discovery & Replication",] 
diff_sign=which(sign(top.pathways$d_ES) != sign(top.pathways$r_ES))

if((is.integer(diff_sign) && length(diff_sign) == 0L)==FALSE){
  top.pathways=top.pathways[-diff_sign,]
  print("Remove diff sign")
}

keep=c("annot_id", "label", "f_qval", "f_ES", "hit_type")
top.pathways=top.pathways[,keep]
length(unique(top.pathways$annot_id)) #number of unique pathways

##########

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

top.pathways$label <- factor(top.pathways$label, levels= oo$label)
top.pathways_label_f_ES=top.pathways[order(top.pathways$label, desc(top.pathways$f_ES)),]
top.pathways$annot_id=stringr::str_replace_all(top.pathways$annot_id, '_', ' ')

names(top.pathways_label_f_ES)[names(top.pathways_label_f_ES) %in% "annot_id"] <- "pathway"
names(top.pathways_label_f_ES)[names(top.pathways_label_f_ES) %in% "label"] <- "cell_type"


write.table(file="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/tables/TableS5_dge_gsea.tsv", top.pathways_label_f_ES,  sep="\t", row.names = FALSE)















