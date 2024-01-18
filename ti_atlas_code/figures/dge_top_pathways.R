#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/")

f="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/adata_dendrogram/cell_types_dendro.tsv"
dend=read.table(f, sep="\t", header=TRUE, stringsAsFactors = FALSE)


f3='/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-clean_annotation-full.csv'
ww=read.table(f3, sep=",", header=TRUE, stringsAsFactors = FALSE)

#p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy/analysis/dge/mast_random_effect_final/disease_status-run_8-mast-cp10k_greater_1-filter_mt_ribo_ig/"
#p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-old_samples_only/cell_prob_filter_onestd/"
#p5="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/differential_expression-new_samples_only/cell_prob_filter_onestd/"

#p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-old_samples_only/"
#p4="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-new_samples_only/"
#p2="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/data/gut-freeze003/ti-cd_healthy-auto_processed_v001-labels_freeze_v003/analysis/dge/cell_prob_filter_0pt5/differential_expression-all_samples/"

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


keepme=c("label", "category", "annot_id", "significant",  "freeze", "qvalue_bh", "ES")
ll=df[,keepme]
ll <- reshape2::dcast(ll, annot_id + category + label ~ freeze , value.var = "significant")
ll$hit_type <- "Neither"
ll$hit_type[(ll$Discovery & ll$Replication)] <- "Discovery & Replication"
ll$Discovery <- NULL
ll$Replication <- NULL
ll$Full <- NULL

length(unique(ll[ll$hit_type=="Discovery & Replication",]$annot_id))

# Select only discovery and replication
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

library(dplyr)
df_plt <- ll1 %>% right_join(ll, by=c("category", "label","annot_id"))
df_plt <- ll2 %>% right_join(df_plt, by=c("category", "label","annot_id"))


############################

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)


choose_celltypes=c("Stem cells", "Enterocyte", "Secretory",  "Myeloid")
#choose_celltypes=c("B Cell plasma", "B Cell")
#choose_celltypes=c("Stem cells", "Enterocyte", "Secretory") # "Myeloid", "T Cell", "Mesenchymal", "B Cell plasma", "B Cell"
#choose_celltypes=c("Myeloid" ) # "Myeloid",
#choose_celltypes=c("T Cell" )
df_plt_sub=df_plt[df_plt$category %in% choose_celltypes,]
oo=oo[oo$category %in% choose_celltypes,]
df_plt_sub$label <- factor(df_plt_sub$label, levels=oo$label)


length(unique(df_plt_sub[df_plt_sub$hit_type=="Discovery & Replication",]$annot_id))

df_plt_sub1=df_plt_sub[df_plt_sub$hit_type=="Discovery & Replication",] 
df_plt_sub1=df_plt_sub1[df_plt_sub1$d_ES > 0,]


##Take top pathways
###############
top_category <- df_plt_sub1 %>% group_by(category) %>% slice_min(f_qval, n =30) %>% arrange(f_qval) 
top_category1 <- df_plt_sub1 %>% group_by(category) %>% slice_max(f_ES, n =30) %>% arrange(f_ES) 
#top_category <- df_plt_sub1 %>% group_by(hit_type, category) %>% slice_min(r_qval, n =30) %>% arrange(r_qval) 
cap_annot = unique(top_category$annot_id)[1:10]


#cap_annot = c("CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION" ,"ANTIGEN_PROCESSING_CROSS_PRESENTATION","INNATE_IMMUNE_SYSTEM",                                
# "CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM","ER_PHAGOSOME_PATHWAY","ADAPTIVE_IMMUNE_SYSTEM",                              
#"INTERFERON_SIGNALING" , "IMMUNE_SYSTEM","INTERFERON_ALPHA_BETA_SIGNALING",                     
# "POST_TRANSLATIONAL_PROTEIN_MODIFICATION")          
###############

df_plt_sub1=df_plt_sub1[df_plt_sub1$annot_id %in% cap_annot,]
#df_plt_sub1=reshape2::melt(df_plt_sub1, id.vars=c("annot_id", "label", "category", 
#                                                  "hit_type", "d_qval","r_qval", "d_ES", "r_ES"),
#                           measure.vars=c( ))


df_plt_sub1 <- df_plt_sub1 %>% 
  dplyr::mutate(qvalue_bh = case_when(variable == 'd_ES' ~ d_qval,
                               variable == 'r_ES' ~ r_qval),
         plot_shape = case_when(variable == 'd_ES' ~ "◀",
                                variable == 'r_ES' ~ "▶"),
         plot_just = case_when(variable == 'd_ES' ~ 1,
                               variable == 'r_ES' ~ 0))

df_plt_sub1$category <- factor(df_plt_sub1$category , levels=c("Stem cells", "Enterocyte", "Secretory",  "Myeloid", "T Cell", "Mesenchymal", "B Cell plasma", "B Cell")) 
#df_plt_sub1$hit_type <- factor(df_plt_sub1$hit_type, levels=c("Discovery & Replication", "Discovery", "Replication", "Neither"))


#df_plt_sub1$variable <- factor(df_plt_sub1$variable, levels=c("d_ES", "r_ES"), labels=c("Discovery", "Replication"))
df_plt_sub1$annot_id=factor(df_plt_sub1$annot_id, levels=rev(cap_annot), labels=stringr::str_replace_all(rev(cap_annot), '_', ' '))
df_plt_sub1$label <- factor(df_plt_sub1$label, levels=oo$label)


######################

# plt <- ggplot(df_plt_sub1, aes(x=label, y=annot_id, 
#                            group=variable, color=value, shape=variable)) +
#   geom_point(alpha=0) +
#   geom_tile(alpha=0) +
#   geom_text(mapping=aes(size=-log10(qvalue_bh), label=plot_shape, hjust=plot_just),
#             position=position_dodge(width=0.3)) +
#   scale_color_gradient2(low="blue", mid="white",
#                         high="red",limits=c(-1,1)) +
#   scale_shape_manual(values=c("◀", "▶")) +
#  scale_x_discrete(drop=FALSE) + #  facet_grid(~category) 
#   labs(title = "", y = "", x = "", shape = 'Dataset',
#        color = "Enrichment score", size=expression(-log[10](P[FGSEA]))) +
#   theme_classic() +
#   theme( strip.text.x = element_text(size=12),
#          axis.text.x = element_text(size=12, angle=90, hjust=1, vjust=0.5),
#          axis.text.y = element_text(size=12),
#          #strip.text.x = element_text(angle = 0, size = 17),
#          strip.background = element_blank()) + 
#   guides(
#     size = guide_legend(override.aes = list(label = "▶"), order = 2),
#     shape = guide_legend(override.aes = list(alpha = 1, size = 5)), order = 1)
# 
# plt
# df_plt_sub2=df_plt_sub1[,c("annot_id", "label", "f_ES")]
# df_plt_sub3=dcast(df_plt_sub2, annot_id~label)
# annot_id=df_plt_sub3$annot_id
# df_plt_sub3$annot_id <- NULL
# df_plt_sub4=as.matrix(df_plt_sub3)
# rownames(df_plt_sub4) <- annot_id
# 
# 
# hc <- hclust(dist(df_plt_sub4), "ave")
# 
# library(ggdendro)
# p <- ggdendro::ggdendrogram(hc, rotate = FALSE, size = 2) 
# p <- p + coord_flip() + theme_dendro()
# p
# 
# #df_plt_sub4[is.na(df_plt_sub4)] <- 0
# intersect(as.character(unique(df_plt_sub2$label)),oo$label)
# 
# 
# sub_anno <- oo$label[colnames(df_plt_sub4) %in% oo$label]
# sub_anno1 <- c(                       "SCAVENGING OF HEME FROM PLASMA"  ,         "INTERFERON GAMMA SIGNALING"    ,             "IMMUNE SYSTEM"  ,                                     
#                         "ER PHAGOSOME PATHWAY"     ,                            "CYTOKINE SIGNALING IN IMMUNE SYSTEM"      ,            "INNATE IMMUNE SYSTEM"  ,                              
#                          "ANTIGEN PROCESSING CROSS PRESENTATION"  ,              "CLASS I MHC MEDIATED ANTIGEN PROCESSING PRESENTATION" ,"POST TRANSLATIONAL PROTEIN MODIFICATION" ,            
#                          "ADAPTIVE IMMUNE SYSTEM" )
# 
# library(heatmaply)
# p <- pheatmap::pheatmap(df_plt_sub4, cluster_rows=TRUE, cluster_cols=FALSE, na_col = "grey90", annotation_col = sub_anno, annotation_row = sub_anno1)#treeheight_row = 0, treeheight_col = 0, , cluster_rows = F
# p

##########
summary(df_plt_sub1$f_ES)

plt <- ggplot(df_plt_sub1, aes(x=label, y=annot_id, fill=f_ES)) + facet_grid(~category, scales="free", space="free") 
plt <- plt + geom_tile() + scale_fill_gradient2(low = "bisque",mid = "lightsalmon",  high = "red", midpoint = 0.7) #+ scale_fill_gradient(low = "cornsilk1", high = "red")  
plt <- plt + geom_tile(df_plt_sub1, aes(x=label, y=annot_id, fill=f_ES, colour = 'black')) #geom_point(shape=8, size=1)
plt <- plt + theme_classic()  +
  theme( strip.text.x = element_text(size=12),
         axis.text.x = element_text(size=12, angle=90, hjust=1, vjust=0.5),
         axis.text.y = element_text(size=12),
         axis.title.x = element_blank(), axis.title.y = element_blank(),
         #strip.text.x = element_text(angle = 0, size = 17),
         strip.background = element_blank()) + labs(fill="Full cohort ES")
plt


ggsave(plt, filename=paste("dge_top_pathways.pdf"), width=15, height =8, device = cairo_pdf)
ggsave(plt, filename=paste("dge_top_pathways_new1.png"), width=15, height =8)


cap_annot








