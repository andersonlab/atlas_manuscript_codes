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
#df$annot_id=gsub("SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR ", "SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR", df$annot_id)
# sub("\\_+$", "", df$annot_id)
df$annot_id <- gsub("\\_+$", "", df$annot_id)    # 
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
#choose_celltypes=c("Stem cells", "Enterocyte", "Secretory",  "Myeloid", "T Cell", "Mesenchymal", "B Cell plasma", "B Cell") # "Myeloid", "T Cell", "Mesenchymal", "B Cell plasma", "B Cell"


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


pa1=aa$annot_level1

# pa=c("FOLDING_ASSEMBLY_PEPTIDE_LOADING_OF_CLASS_I_MHC", 
#       "ANTIGEN_PROCESSING_CROSS_PRESENTATION"
#       , "ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION")
#  
# pa1=c("INTERFERON_ALPHA_BETA_SIGNALING",
#        "INTERFERON_GAMMA_SIGNALING",
#        "ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES") #SIGNALLING BY CSF1 and CSF3
# 
# pa2=c("CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
#       "INTERFERON_SIGNALING",) 


pa=c("CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
     "FOLDING_ASSEMBLY_PEPTIDE_LOADING_OF_CLASS_I_MHC", 
     "ANTIGEN_PROCESSING_CROSS_PRESENTATION", 
     "ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION",
     "INTERFERON_SIGNALING",
      "INTERFERON_ALPHA_BETA_SIGNALING",
      "INTERFERON_GAMMA_SIGNALING",
      "ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES") #SIGNALLING BY CSF1 and CSF3



w1=c("CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
     "FOLDING_ASSEMBLY_PEPTIDE_LOADING_OF_CLASS_I_MHC", 
     "ANTIGEN_PROCESSING_CROSS_PRESENTATION", 
     "ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION")
     
w2=c("INTERFERON_SIGNALING",
     "INTERFERON_ALPHA_BETA_SIGNALING",
     "INTERFERON_GAMMA_SIGNALING",
     "ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES") #SIGNALLING BY CSF1 and CSF3


high.level.pathways = pa
#df_plt_sub1_cut=df_plt_sub1[df_plt_sub1$annot_id %in% high.level.pathways,]
df=df[df$freeze!="Full",]
# 
df[df$annot_id %in% w1, "pathway_type"] <- "Level0"
df[df$annot_id %in% w2, "pathway_type"] <- "Level1"

#level0=c("CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION","INTERFERON_SIGNALING") 
#df[df$annot_id %in% level0, "pathway_type"] <- "Level0"
#df[is.na(df$pathway_type),"pathway_type"] <- "level1"

#df[df$annot_id %in% pa2, "pathway_type"] <- "Innate Immune System"


#Add labels
file1="/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/clean_annotation/data-order.csv"
oo=read.table(file1, sep=",", header=TRUE, comment.char = "",check.names = FALSE, stringsAsFactors = FALSE)
labels_used=oo[oo$label %in%df$label,]$label
df$label <- factor(df$label, levels=c(labels_used))

df[df$category %in% c("Stem cells", "Enterocyte", "Secretory"), "compartment"] <- "Epithelial"
df[df$category %in% c("Myeloid"), "compartment"] <- "Myeloid"

intersect(df$annot_id, high.level.pathways)
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

# Order
high.level.pathways

 df_plt_sub1_cut_nonsig$annot_id <- gsub("_", " ", df_plt_sub1_cut_nonsig$annot_id)
 df_plt_sub1_cut_sig$annot_id <- gsub("_", " ", df_plt_sub1_cut_sig$annot_id)
 df_plt_sub1_cut_nonsig$annot_id <- factor(df_plt_sub1_cut_nonsig$annot_id, levels=gsub("_", " ", rev(high.level.pathways)))
 df_plt_sub1_cut_sig$annot_id <- factor(df_plt_sub1_cut_sig$annot_id, levels=gsub("_", " ", rev(high.level.pathways)))

df_plt_sub1_cut_sig_replicable=df_plt_sub1_cut_sig[df_plt_sub1_cut_sig$replicable=="yes",]


# mycols=data.frame("pathways"=high.level.pathways, 
#                   colours=c("navy", "brown1", "brown1", "brown1", "navy", "brown1", "brown1", "brown1"))
# ppalette <- as.character(mycols$colours)

#ppalette <- as.character(rev(c("navy", "brown1", "brown1", "brown1", "navy", "brown1", "brown1", "brown1")))
ppalette <- as.character(rev(c("navy", "steelblue", "steelblue", "steelblue", "navy", "steelblue", "steelblue", "steelblue")))
library(paletteer)
paletter=paletteer_c("ggthemes::Red-Blue-White Diverging", 50)
paletter=as.vector(paletter)


# ll=c(3,2,1)
# ll1=c(3,2,1)
# pp=cbind(ll, ll1)
# pp=as.matrix(pp)
# model <- hclust(dist(pp), "ave")
# dhc <- as.dendrogram(model)
# ddata <- dendro_data(dhc, type = "rectangle")
# p <- ggplot(segment(ddata)) + 
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#   #coord_flip() + 
#   #scale_y_reverse(expand = c(0.2, 0))+ 
#   coord_flip() + 
#   theme_dendro()
# p
library(ggh4x)
??facet_nested

plt <- ggplot(df_plt_sub1_cut_nonsig, aes(x=label, y=annot_id, fill=ES)) + facet_nested(pathway_type~freeze+compartment, scales="free", space="free") 
plt <- plt + geom_tile() 
#plt <- plt + scale_fill_gradient2(low = "steelblue", mid = "bisque", high = "red", midpoint = 0)#+ scale_fill_gradient2(low = "bisque",mid = "lightsalmon",  high = "red", midpoint = 0.7) #+ scale_fill_gradient(low = "cornsilk1", high = "red")  
plt <- plt + scale_fill_gradientn(colours = rev(paletter),space = "Lab", na.value = "white", limit = c(-1,1))
plt <- plt + geom_tile(data = df_plt_sub1_cut_sig, aes(x=label, y=annot_id, fill=ES), colour = 'black', size=0.7)
plt <- plt + geom_point(data = df_plt_sub1_cut_sig_replicable, aes(x=label, y=annot_id), shape=8, size=1)
plt <- plt #+ geom_point(shape=8, size=1)
#plt
plt <- plt + theme_classic()  +
  theme( strip.text.x = element_text(size=25),
         strip.text.y = element_blank(),
         axis.text.x = element_text(size=15, angle=90, hjust=1, vjust=0.5),
         axis.text.y = element_text(size=19, colour = ppalette),
         legend.text=element_text(size=15),
         legend.title=element_text(size=15),
         axis.title.x = element_blank(), axis.title.y = element_blank(),
         #strip.text.x = element_text(angle = 0, size = 17),
         strip.background = element_blank()) + labs(fill="Enrichment\nscore")
#plt <- plt + geom_vline(xintercept = c(15) +1.5, linetype="dashed", colour="black", size=2)

plt

ggsave(plt, filename=paste("dge_top_pathways_level2.png"), width=23, height =9)
ggsave(plt, filename=paste("dge_top_pathways_level2.pdf"), width=23, height =9, device = cairo_pdf, dpi=200)
















dd1=df_plt_sub1_cut_sig[((df_plt_sub1_cut_sig$freeze=="Discovery") & 
                           (df_plt_sub1_cut_sig$annot_id == "CLASS I MHC MEDIATED ANTIGEN PROCESSING PRESENTATION")),]


dd2=df_plt_sub1_cut_sig[((df_plt_sub1_cut_sig$freeze=="Replication") & 
                           (df_plt_sub1_cut_sig$annot_id == "CLASS I MHC MEDIATED ANTIGEN PROCESSING PRESENTATION")),]


round(mean(dd1$ES),2)
round(mean(dd2$ES),2)




##Take top pathways
###############
top_category <- df_plt_sub1 %>% group_by(category) %>% slice_min(d_qval, n =30) %>% arrange(d_qval) 
top_category1 <- df_plt_sub1 %>% group_by(category) %>% slice_max(d_ES, n =30) %>% arrange(d_ES) 
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

##########
summary(df_plt_sub1$f_ES)

plt <- ggplot(df_plt_sub1, aes(x=label, y=annot_id, fill=f_ES)) + facet_grid(~category, scales="free", space="free") 
plt <- plt + geom_tile() + scale_fill_gradient2(low = "bisque",mid = "lightsalmon",  high = "red", midpoint = 0.7) #+ scale_fill_gradient(low = "cornsilk1", high = "red")  
plt <- plt + geom_point(shape=8, size=1)
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








