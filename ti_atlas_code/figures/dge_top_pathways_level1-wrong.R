library(plyr)

#setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures_final_new/")
setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/protocol_paper/DGE/")

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

p1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/protocol_paper/DGE/"
f1=paste0(p1, "state_merged-gsea_results.tsv")
df1=read.table(f1, sep="\t", header=TRUE, stringsAsFactors = FALSE)

#plot(sort(df1$ES), df1$gset_size)
#plot(df1$pvalue, df1$gset_size)

pp=data.frame("annot_id"= df1$annot_id, "gset_size"= df1$gset_size)
pp=unique(pp)
pp$annot_id=gsub("REACTOME_", "", pp$annot_id)


df1$dataset <- "storage"
keepcol=c("category", "annot_id","signed_ranking","annot_coef", "qvalue_bh",  "ES", "dataset")
#keepcol=c("label", "annot_id","signed_ranking","annot_coef", "qvalue_bh",  "ES", "freeze", "leadingEdge")
df1=df1[,keepcol]
df = rbind(df1)
#df=merge(df, ww, by = "label", all.x = TRUE)
df=df[grep("REACTOME", df$annot_id),]
df=df[df$signed_ranking=="True",]
df$minuslog10qval=-log10(df$qvalue_bh)
df$regulation <- ifelse(df$annot_coef >= 0, "up-regulated", "down-regulated")
df=df[!is.na(df$annot_id),]
df$annot_id=gsub("RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS_", "RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS", df$annot_id)
df$annot_id=gsub("REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS_", "REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR", df$annot_id)
df$annot_id=gsub("THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT", "CITRIC_ACID_CYCLE_RESPIRATORY_ELECTRON_TRANSPORT", df$annot_id)
df$annot_id=gsub("ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC", "FOLDING_ASSEMBLY_PEPTIDE_LOADING_OF_CLASS_I_MHC", df$annot_id)
df$annot_id=gsub("IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL", "IMMUNOREGULATORY_INTERACTIONS", df$annot_id)
df$annot_id=gsub("NUCLEOTIDE_BINDING_DOMAIN_LEUCINE_RICH_REPEAT_CONTAINING_RECEPTOR_NLR_SIGNALING_PATHWAYS", "NUCLEOTIDE_BINDING_DOMAIN_NLR_SIGNALING_PATHWAYS", df$annot_id)


df$annot_id=gsub("REACTOME_", "", df$annot_id)
df$significant <- df$qvalue_bh < 0.05


keepme=c("category", "annot_id", "significant",  "dataset", "qvalue_bh", "ES")
ll=df[,keepme]
ll <- reshape2::dcast(ll, annot_id + category ~ dataset , value.var = "significant")
ll$hit_type <- "Neither"
ll$hit_type[(ll$storage==TRUE)] <- "Storage"

length(unique(ll[ll$hit_type=="Storage",]$annot_id))

# Select only discovery and replication
ll=ll[ll$hit_type=="Storage",]

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
choose_celltypes=c("Stem cells", "Enterocyte", "Secretory",  "Myeloid", "T Cell", "Mesenchymal", "B Cell plasma", "B Cell")
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
                      "IMMUNOREGULATORY_INTERACTIONS", #"IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
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
                       "NUCLEOTIDE_BINDING_DOMAIN_NLR_SIGNALING_PATHWAYS", #NUCLEOTIDE_BINDING_DOMAIN_LEUCINE_RICH_REPEAT_CONTAINING_RECEPTOR_NLR_SIGNALING_PATHWAYS",
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

df[df$category %in% c("Stem cells", "Enterocyte", "Secretory"), "compartment"] <- "Epithelial"
df[df$category %in% c("Myeloid"), "compartment"] <- "Myeloid"
df[df$category %in% c("T Cell", "Mesenchymal", "B Cell plasma", "B Cell"), "compartment"] <- "Rest"


df$compartment <- factor(df$compartment, levels=c("Epithelial", "Myeloid", "Rest"))
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
oderme=df_plt_sub1_cut_sig[,c("annot_id", "label", "replicable")] 
#myorder=rev(sort(table(oderme$annot_id))) 
oderme1=oderme[oderme$replicable=="yes",]
myorder1=ddply(oderme1,.(annot_id),nrow)
myorder1$nr_replicable=myorder1$V1/2
myorder1$V1 <- NULL
`%ni%` <- Negate(`%in%`)
all_annot_id=unique(as.character(df_plt_sub1_cut_nonsig$annot_id))
all_annot_id=all_annot_id[!is.na(all_annot_id)]
ll=data.frame("annot_id"=all_annot_id[all_annot_id %ni% unique(as.character(oderme1$annot_id))],
               "nr_replicable"=0)
myorder2=rbind(myorder1, ll)
myorder=myorder2$nr_replicable
names(myorder) <- myorder2$annot_id

#oderme2=oderme[oderme$replicable=="no",]
#myorder2=ddply(oderme2,.(annot_id),nrow)

myorder_new=c(rev(sort(myorder[pa])), rev(sort(myorder[pa1])), rev(sort(myorder[pa2])))
df_plt_sub1_cut_nonsig$annot_id <- gsub("_", " ", df_plt_sub1_cut_nonsig$annot_id)
df_plt_sub1_cut_sig$annot_id <- gsub("_", " ", df_plt_sub1_cut_sig$annot_id)
df_plt_sub1_cut_nonsig$annot_id <- factor(df_plt_sub1_cut_nonsig$annot_id, levels=gsub("_", " ", rev(names(myorder_new))))
df_plt_sub1_cut_sig$annot_id <- factor(df_plt_sub1_cut_sig$annot_id, levels=gsub("_", " ", rev(names(myorder_new))))

df_plt_sub1_cut_sig_replicable=df_plt_sub1_cut_sig[df_plt_sub1_cut_sig$replicable=="yes",]


mycols=data.frame("pathways"=rev(myorder_new), colours=c(rep("cornflowerblue", length(sort(myorder[pa2]))),
                                                    rep("navy", length(sort(myorder[pa1]))), 
                                                    rep("brown1", length(sort(myorder[pa])))))
ppalette <- as.character(mycols$colours)


#df_plt_sub1_cut_nonsig=df_plt_sub1_cut_nonsig[df_plt_sub1_cut_nonsig$annot_id %in% gsub("_", " ", myorder1$annot_id),]
#df_plt_sub1_cut_sig=df_plt_sub1_cut_nonsig[df_plt_sub1_cut_sig$annot_id %in% gsub("_", " ", myorder1$annot_id),]
#df_plt_sub1_cut_sig_replicable=df_plt_sub1_cut_sig_replicable[df_plt_sub1_cut_sig_replicable$annot_id %in% gsub("_", " ", myorder1$annot_id),]


library(ggplot2)
#+compartment
plt <- ggplot(df_plt_sub1_cut_nonsig, aes(x=label, y=annot_id, fill=ES)) + facet_grid(pathway_type~freeze, scales="free", space="free") 
plt <- plt + geom_tile() + scale_fill_gradient2(low = "steelblue", mid = "bisque", high = "red", midpoint = 0)#+ scale_fill_gradient2(low = "bisque",mid = "lightsalmon",  high = "red", midpoint = 0.7) #+ scale_fill_gradient(low = "cornsilk1", high = "red")  
plt <- plt + geom_tile(data = df_plt_sub1_cut_sig, aes(x=label, y=annot_id, fill=ES), colour = 'black', size=0.7)
plt <- plt + geom_point(data=df_plt_sub1_cut_sig_replicable, aes(x=label, y=annot_id), shape=8, size=1)
#plt <- plt + facet_grid(~freezcategory, scale="free")
plt <- plt + theme_classic()  +
  theme( strip.text.x = element_text(size=35),
         strip.text.y = element_text(size=30, angle=360, hjust=0),
         #strip.text.x = element_blank(),
         axis.text.x = element_text(size=18, angle=90, hjust=1, vjust=0.5),
         axis.text.y = element_text(size=18),
         legend.text=element_text(size=15),
         legend.title=element_text(size=15),
         axis.title.x = element_blank(), axis.title.y = element_blank(),
         #strip.text.x = element_text(angle = 0, size = 17),
         strip.background = element_blank()) + labs(fill="Enrichment score")
plt

ggsave(plt, filename=paste("dge_top_pathways_level1.png"), width=35, height =15)
ggsave(plt, filename=paste("dge_top_pathways_level1.pdf"), width=35, height =15, device = cairo_pdf, dpi=200)





###

