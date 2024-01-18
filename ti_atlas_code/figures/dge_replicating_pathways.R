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

#plot(sort(df1$ES), df1$gset_size)
#plot(df1$pvalue, df1$gset_size)

pp=data.frame("annot_id"= df1$annot_id, "gset_size"= df1$gset_size)
pp=unique(pp)
pp$annot_id=gsub("REACTOME_", "", pp$annot_id)


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


#choose_celltypes=c("Stem cells", "Enterocyte", "Secretory",  "Myeloid")
#choose_celltypes=c("T Cell", "Mesenchymal", "B Cell plasma", "B Cell") # "Myeloid", "T Cell", "Mesenchymal", "B Cell plasma", "B Cell"


df_plt_sub=df_plt[df_plt$category %in% choose_celltypes,]
oo=oo[oo$category %in% choose_celltypes,]
df_plt_sub$label <- factor(df_plt_sub$label, levels=oo$label)


length(unique(df_plt_sub[df_plt_sub$hit_type=="Discovery & Replication",]$annot_id))

df_plt_sub1=df_plt_sub[df_plt_sub$hit_type=="Discovery & Replication",] 
df_plt_sub1=df_plt_sub1[df_plt_sub1$d_ES > 0,]


# Average by categories
#check <- df_plt_sub1 %>% group_by(annot_id) %>% summarize(d_ES_mean=mean(d_ES), r_ES_mean=mean(r_ES), f_ES_mean=mean(f_ES))
check1 <- df_plt_sub1 %>% group_by(category, annot_id) %>% arrange(d_qval)
check2 <- df_plt_sub1 %>% group_by(category, annot_id) %>% arrange(r_qval)
check3 <- df_plt_sub1 %>% group_by(category, annot_id) %>% arrange(f_qval)


# high.level.pathways=c(
#   "Autophagy", "Cell Cycle", "Cell Cell communication", "Cellular responses to external stimuli", "Chromatin modifying enzymes", #"Cellular responses to stimuli", "Chromatin organization",
#   "Circadian Clock", "Developmental Biology", "Digestion and absorption", "Disease", "DNA Repair", 
#   "DNA Replication", "Drug ADME", "Extracellular matrix organization", "Gene expression (Transcription)",
#   "Hemostasis", "Immune System", "Metabolism", "Metabolism of proteins", "Metabolism of RNA", "Muscle contraction", 
#   "Neuronal System", "Organelle biogenesis and maintenance", "Programmed Cell Death", "Protein localization", 
#   "Reproduction", "Sensory Perception", "Signal Transduction", "Transport of small molecules", "Vesicle-mediated transport")
# 
# high.level.pathways=toupper(high.level.pathways)
# high.level.pathways=gsub(" ", "_", high.level.pathways)
# high.level.pathways

pa=c(
  "TCR_SIGNALING", 
  "COSTIMULATION_BY_THE_CD28_FAMILY",
  "SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR_",
  "CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
  "MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  "RAP1_SIGNALLING",
  "BUTYROPHILIN_BTN_FAMILY_INTERACTIONS")

pa1=c("INTERFERON_SIGNALING",
      "SIGNALING_BY_INTERLEUKINS",
      "GROWTH_HORMONE_RECEPTOR_SIGNALING",
      "PROLACTIN_RECEPTOR_SIGNALING",
      "TNFR2_NON_CANONICAL_NF_KB_PATHWAY",
      "FLT3_SIGNALING") #SIGNALLING BY CSF1 and CSF3

pa2=c("TOLL_LIKE_RECEPTOR_CASCADES",
      "COMPLEMENT_CASCADE",
      "NUCLEOTIDE_BINDING_DOMAIN_LEUCINE_RICH_REPEAT_CONTAINING_RECEPTOR_NLR_SIGNALING_PATHWAYS",
      "ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING",
      "DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA",
      "CYTOSOLIC_SENSORS_OF_PATHOGEN_ASSOCIATED_DNA_",
      "FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS",
      "DAP12_INTERACTIONS",
      "FC_EPSILON_RECEPTOR_FCERI_SIGNALING",
      "C_TYPE_LECTIN_RECEPTORS_CLRS_",
      "ANTIMICROBIAL_PEPTIDES",
      "NEUTROPHIL_DEGRANULATION",
      "ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES",
      "ALPHA_PROTEIN_KINASE_1_SIGNALING_PATHWAY")




high.level.pathways = c(pa, pa1, pa2)
#df_plt_sub1_cut=df_plt_sub1[df_plt_sub1$annot_id %in% high.level.pathways,]
df=df[df$freeze!="Full",]

df[df$annot_id %in% pa, "pathway_type"] <- "Adaptive Immune System"
df[df$annot_id %in% pa1, "pathway_type"] <- "Cytokine Signaling"
df[df$annot_id %in% pa2, "pathway_type"] <- "Innate Immune System"


#Add labels
file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
labels_used=oo[oo$label %in%df$label,]$label
df$label <- factor(df$label, levels=c(labels_used))


# Cap to high level pathways
df_plt_sub1_cut_nonsig=df[df$annot_id %in% high.level.pathways,]
df_plt_sub1_cut_nonsig=df_plt_sub1_cut_nonsig[df_plt_sub1_cut_nonsig$category %in% choose_celltypes,]
df_plt_sub1_cut_nonsig$freeze <- factor(df_plt_sub1_cut_nonsig$freeze, levels=c("Discovery", "Replication", "Full"))
df_plt_sub1_cut_nonsig$annot_id <- factor(df_plt_sub1_cut_nonsig$annot_id, levels=rev(high.level.pathways))
df_plt_sub1_cut_nonsig$annot_id <- factor(df_plt_sub1_cut_nonsig$annot_id, levels=rev(high.level.pathways))
df_plt_sub1_cut_sig=df_plt_sub1_cut_nonsig[df_plt_sub1_cut_nonsig$qvalue_bh <0.05,]
df_plt_sub1_cut_sig=df_plt_sub1_cut_sig[!is.na(df_plt_sub1_cut_sig$qvalue_bh),]
df_plt_sub1_cut_sig$annot_id <- factor(df_plt_sub1_cut_sig$annot_id, levels=rev(high.level.pathways))


#Add replicable
temp=df_plt_sub1_cut_sig[,c("label", "annot_id","ES", "freeze")]
annot_ids=as.character(unique(temp$annot_id))
cell_types=as.character(unique(temp$label))
for(i in 1:length(annot_ids)){#length(annot_ids)
  for(j in 1:length(cell_types)){
    #i=2
    #j=1
    #length(cell_types)
    print(paste0("i: ", i, " " , "j: ", j))
    temp1=temp[temp$annot_id==annot_ids[i] & temp$label==cell_types[j],]
    
    
    if(nrow(temp1)==2){
      #print(temp1)
      if(sign(temp1[temp1$freeze=="Discovery",]$ES)==sign(temp1[temp1$freeze=="Replication",]$ES))
        print(paste0("common: ", annot_ids[i]," ", cell_types[j])) 
      print(temp1)
      df_plt_sub1_cut_sig[df_plt_sub1_cut_sig$annot_id==annot_ids[i] & df_plt_sub1_cut_sig$label==cell_types[j],"replicable"] <- "yes"
      
    } else {
      #print(temp1)
      print(paste0("not common: ", annot_ids[i]," ", cell_types[j])) 
      df_plt_sub1_cut_sig[df_plt_sub1_cut_sig$annot_id==annot_ids[i] & df_plt_sub1_cut_sig$label==cell_types[j],"replicable"] <- "no"
    }
    
  }
}

#Make order
oderme=df_plt_sub1_cut_sig[,c("annot_id", "label")] 
myorder=rev(sort(table(oderme$annot_id))) 
myorder_new=c(rev(sort(myorder[pa])), rev(sort(myorder[pa1])), rev(sort(myorder[pa2])))
df_plt_sub1_cut_nonsig$annot_id <- gsub("_", " ", df_plt_sub1_cut_nonsig$annot_id)
df_plt_sub1_cut_sig$annot_id <- gsub("_", " ", df_plt_sub1_cut_sig$annot_id)
df_plt_sub1_cut_nonsig$annot_id <- factor(df_plt_sub1_cut_nonsig$annot_id, levels=gsub("_", " ", rev(names(myorder_new))))
df_plt_sub1_cut_sig$annot_id <- factor(df_plt_sub1_cut_sig$annot_id, levels=gsub("_", " ", rev(names(myorder_new))))

df_plt_sub1_cut_sig_replicable=df_plt_sub1_cut_sig[df_plt_sub1_cut_sig$replicable=="yes",]



length(unique(df_plt_sub1_cut_sig_replicable$annot_id))

###

epi=df_plt_sub1_cut_sig_replicable[df_plt_sub1_cut_sig_replicable$category %in% c("Enterocyte", "Secretory", "Stem cells", "Myeloid"),]
# epi=df_plt_sub1_cut_sig_replicable[df_plt_sub1_cut_sig_replicable$category %in% c("Enterocyte", "Secretory", "Stem cells"),]
# epi=df_plt_sub1_cut_sig_replicable[df_plt_sub1_cut_sig_replicable$category %in% c("Enterocyte", "Secretory", "Stem cells"),]
# myo=df_plt_sub1_cut_sig_replicable[df_plt_sub1_cut_sig_replicable$category %in% c("Myeloid"),]
# tc=df_plt_sub1_cut_sig_replicable[df_plt_sub1_cut_sig_replicable$category %in% c("T Cell"),]
# bc=df_plt_sub1_cut_sig_replicable[df_plt_sub1_cut_sig_replicable$category %in% c("B Cell"),]
file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)

pp= df_plt_sub1_cut_sig_replicable %>% group_by(label, .drop = FALSE) %>% dplyr::summarise(n_replicable=n()) %>% arrange(dplyr::desc(n_replicable))
oo=oo[c("label", "category")]
pp=merge(pp, oo, by="label")
#pp$category <- factor(pp$category , levels=rev(unique(oo$category)))   

#pp=data.frame(table(df_plt_sub1_cut_sig_replicable$label))
#pp$category <- factor(pp$category , levels=c("Stem cells", "Enterocyte", "Secretory", "Mesenchymal", "Myeloid", "T Cell", "B Cell", "B Cell plasma")) 

file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/palette.tsv"
pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
palette=pal$palette_label
names(palette) <- pal$label

pp[pp$category %in% c("Enterocyte","Stem cells", "Myeloid",  "Secretory"),"category_type"] <- "Epithelial and Myeloid"
pp[is.na(pp$category_type),"category_type"] <- "Immune and Mesenchymal"
pp$category_type <- factor(pp$category_type, levels=c("Epithelial and Myeloid", "Immune and Mesenchymal"), labels=c("Epithelial and Myeloid", "Immune and Mesenchymal"))

pp$category <- factor(pp$category , levels=rev(c("Enterocyte","Stem cells", "Myeloid",  "Secretory",  "B Cell","T Cell", "Mesenchymal",  "B Cell plasma")))

# file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/dge_power.tsv"
# power=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
# 
# pp2 <- pp %>% group_by(category) %>% dplyr::summarize(sum_rep=sum(n_replicable))
# pp1=merge(pp2, power, by="category")
file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/palette.tsv"
pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
palette=pal$palette_label
names(palette) <- pal$label
#palette

gplt = ggplot(pp, aes(x=category, y=n_replicable, fill=label)) #+ geom_jitter(cex=0.3)
gplt = gplt + facet_grid( vars(category_type), scales = "free", space = "free") # 
gplt = gplt + theme_bw() + ylab("Number of replicating pathways in both cohorts") + theme() + scale_fill_manual(values=palette)
gplt = gplt + geom_bar(stat="identity") + coord_flip() #+ scale_y_discrete(drop=F)
gplt = gplt + theme(legend.title = element_blank(), 
                  #strip.text = element_blank(angle=180),
                   axis.title.y=element_blank(), 
                   axis.title.x = element_text(size=13), 
                   strip.text.y = element_text(angle=270, size=12, hjust=0),
                   strip.text.x = element_text(size = 15), 
                   legend.text = element_text(size=15),
                   legend.position= "none",
                   axis.text.y=element_text(size=15),
                   axis.text.x=element_text(size=10),
                   strip.background = element_blank()) 
gplt

ggsave(plot=gplt, filename="dge_replicating_pathways.png", width = 7,height = 5) 
ggsave(plot=gplt, filename="dge_replicating_pathways.pdf", width = 7,height = 5, device = cairo_pdf, dpi=200)
# 
# file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/palette.tsv"
# pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
# palette=pal$palette_category
# names(palette) <- pal$category
# palette=palette[!duplicated(palette)]
# 
# pp$category <- factor(pp$category , levels=c("Enterocyte","Stem cells", "Myeloid",  "Secretory",  "B Cell","T Cell", "Mesenchymal",  "B Cell plasma"))
# 
# gplt = ggplot(pp, aes(x=label, y=n_replicable, fill=category)) #+ geom_jitter(cex=0.3)
# gplt = gplt + facet_grid( vars(category), scales = "free", space = "free") # 
# gplt = gplt + theme_bw() + ylab("Number of replicating pathways") + theme() #+ scale_fill_manual(values=palette)
# gplt = gplt + geom_bar(stat="identity") + coord_flip() #+ scale_y_discrete(drop=F)
# gplt = gplt + theme(legend.title = element_blank(), 
#                     #strip.text = element_blank(angle=180),
#                     axis.title.y=element_blank(), 
#                     axis.title.x = element_text(size=13), 
#                     strip.text.y = element_text(angle=270, size=15, hjust=0),
#                     strip.text.x = element_text(size = 15), 
#                     legend.text = element_text(size=15),
#                     legend.position= "none",
#                     axis.text.y=element_text(size=15),
#                     axis.text.x=element_text(size=10),
#                     strip.background = element_blank()) 
# gplt


# Number of replicable pathways
pp2 <- pp %>% group_by(category) %>% dplyr::summarize(sum_rep=sum(n_replicable))






file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/palette.tsv"
pal=read.table(file1, sep="\t", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
palette=pal$palette_category
names(palette) <- pal$category
palette=palette[!duplicated(palette)]

gplt = ggplot(pp1, aes(x=union_degs, y=sum_rep, colour=category)) #+ geom_jitter(cex=0.3)
gplt = gplt + facet_grid(~freeze, scales = "free", space = "free") # 
gplt = gplt + theme_bw() + ylab("Number of replicating pathways in both cohorts") + theme() + scale_colour_manual(values=palette)
gplt = gplt + geom_point() + coord_flip() #+ scale_y_discrete(drop=F)
#gplt
gplt = gplt + theme(legend.title = element_blank(), 
                    #strip.text = element_blank(angle=180),
                    axis.title.y=element_blank(), 
                    axis.title.x = element_text(size=13), 
                    strip.text.y = element_text(angle=270, size=15, hjust=0),
                    strip.text.x = element_text(size = 15), 
                    legend.text = element_text(size=15),
                   # legend.position= "none",
                    axis.text.y=element_text(size=15),
                    axis.text.x=element_text(size=10),
                    strip.background = element_blank()) 
gplt

#+ scale_fill_manual(values=palette) 
# gplt = gplt + theme(axis.title.x = element_blank(), axis.text.x = element_text(size=17, angle = 45, vjust = 1, hjust=1),
#                     strip.text.x = element_text(size = 17),
#                     #strip.text.x = element_text(size = 17), 
#                     axis.text.y = element_text(size = 15), 
#                     strip.background = element_blank(),
#                     #legend.position = "none",
#                     plot.title = element_text(size = 17)) + facet_wrap(~disease_status, scales = 'free_x') #+ theme(axis_text_x=element_text(angle=90))
# gplt = gplt + ggtitle("MHC-I classical and non-classical")+ stat_compare_means(aes(label = ..p.signif..), method = "t.test",comparisons = my_comparisons) #+ ggtitle("MHC-I classical (HLA-A,B,C)")
# gplt












length(unique(epi$annot_id))
# length(unique(myo$annot_id))
# length(unique(tc$annot_id))
# length(unique(bc$annot_id))


round(mean(epi$ES),2)
round(mean(myo$ES),2)


