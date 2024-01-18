url="http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/Pathways/Human_Reactome_April_02_2023_symbol.gmt"
ll=download.file(url, destfile="lala.gmt")
library(GSA)
df=GSA.read.gmt("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/lala.gmt")
pathways=df$geneset.names
df=df[[1]]

names=lapply(pathways, function(x) strsplit(x, "%"))
names=lapply(names, "[[", 1)  
names=lapply(names, "[[", 1)  
names=unlist(names)

names(df) <- names

ww=data.frame(names)

temp=list()
for(i in 1:length(df)){
  #i=1
  temp[[i]]=data.frame("gene_symbol"=df[[i]], "annot_id"=names[i])
  #df_total=rbind(df_total, temp)
}
df_total = do.call(rbind, temp)
df_total$gene_symbol=as.character(df_total$gene_symbol)
df_total=df_total[which(df_total$gene_symbol!=""),]


df_total$annot_id <- gsub("\\(|\\)", "", df_total$annot_id)
df_total$annot_id <- gsub(" &", "", df_total$annot_id)
df_total$annot_id <- gsub("-", " ", df_total$annot_id)
df_total$annot_id <- gsub("BCR", "BCR ", df_total$annot_id)
df_total$annot_id <- gsub(" ", "_", df_total$annot_id)

df_total_dcast = reshape2::dcast(df_total, annot_id ~ gene_symbol)

# ee=read.table("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/second.level.tsv", sep="\t")
# names(ee) <- c("annot_category", "annot_id")
# length(unique(ee$annot_id))
# paths_level2=gsub("_", " ", unique(ee$annot_id))

#df_total_dcast$annot_id=gsub(" ", "_", df_total_dcast$annot_id)
#df_total_dcast$pathway1 <- toupper(df_total_dcast$pathway1)
#df_total_dcast$pathway2 <- toupper(df_total_dcast$pathway2)


cap_annot =c("IMMUNE_SYSTEM"    , "ADAPTIVE_IMMUNE_SYSTEM"        ,
             "CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION","MHC_CLASS_II_ANTIGEN_PROCESSING_PRESENTATION",
              "ANTIGEN_PROCESSING_CROSS_PRESENTATION" ,    "INTERFERON_SIGNALING" ,      "SIGNALING_BY_INTERLEUKINS" ,                             
              "CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM" , "INTERFERON_GAMMA_SIGNALING", "INTERFERON_ALPHA_BETA_SIGNALING")

X1=df_total_dcast[df_total_dcast$annot_id %in% cap_annot,]
paths=X1$annot_id
X1$annot_id <- NULL
X2=as.matrix(X1)
rownames(X2) <- paths

#ee=ee[ee$annot_id %in% df_plt_sub1$annot_id,]
#df_plt_sub1=merge(df_plt_sub1, ee, by="annot_id")

hc <- hclust(dist(X2), "complete")
plot(hc)
#plot(as.dendrogram(hc),horiz=T)


