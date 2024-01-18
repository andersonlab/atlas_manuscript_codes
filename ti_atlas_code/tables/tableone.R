setwd("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/figures-ti_cd_singlecell/figures/02000-deg/tables/")


library(UpSetR)
library(stringr)
library(xtable)

# Read in the data
dat_meta <- read.csv(
  paste0("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/ti_atlas_code/figures/data/GUT_scRNAseq-cleaned.csv"),
  
  #paste0( "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/GUT_scRNAseq-cleaned.csv"),
  #paste0( "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/repo/scrna_cellranger/sync_status/samples_metainfo.tsv"),
  
  sep = ",", stringsAsFactors = FALSE)

# Read in the TI sample ids
dat_smpls_ti <- read.csv(
  paste0("/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/sample_ids.tsv"),
  sep = "\t", stringsAsFactors = FALSE)

# Get the metadata on these samples
dat_ti <- subset(
  dat_meta,
  sanger_sample_id %in% dat_smpls_ti$experiment_id &
    #(biopsy_type == "ti" | biopsy_type == "neoti")
  biopsy_type == "TI" 
)

dat_meta$biopsy_type

dat_ti[dat_ti$sanger_sample_id %in% dat_smpls_ti[dat_smpls_ti$cohort == "Discovery","sanger_sample_id"],"cohort"] <- "Discovery"
dat_ti[dat_ti$sanger_sample_id %in% dat_smpls_ti[dat_smpls_ti$cohort == "Replication","sanger_sample_id"],"cohort"] <- "Replication"

# There are two samples that came on the same day with paired ti and rectum
# but for whom we do not know the original participant id
# (both women similar age). Fix that for sample counts
#dat_ti$patient_id[is.na(dat_ti$patient_id)] <- c("unknown1", "unknown2")

# Merge the data
dat <- dat_ti

# Clean up data labels
#dat$biopsy_type <- gsub("ti", "Terminal ileum", dat$biopsy_type)
dat$biopsy_type <- gsub("TI", "Terminal ileum", dat$biopsy_type)
dat$disease_status <- gsub("cd", "CD", dat$disease_status)
dat$disease_status <- gsub("healthy", "Healthy", dat$disease_status)


table(dat$disease_status, dat$cohort)

#dat_ti[is.na(dat$medication_class),]$medication_class <- "unknown"

#dat$medication_class=as.character(dat$medication_class)
# Read in the data
dat_medication <- read.csv(
  #paste0( "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/samples_metainfo.csv"),
  paste0( "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/manuscript/figures/data/GUT_scRNAseq-cleaned.csv"),
  #paste0( "/home/ubuntu/ubuntu/sc_ibd_project/ti_atlas/repo/scrna_cellranger/sync_status/samples_metainfo.tsv"),
  sep = ",", stringsAsFactors = FALSE)

dat_medication=data.frame("sanger_sample_id"=dat_medication$sanger_sample_id, 
                          "medication_class"=dat_medication$medication_class)

table(dat$medication_class)
table(dat$medication_class, dat$disease_status)

dat=merge(dat, dat_medication, by="sanger_sample_id")
table(dat$medication_class, dat$disease_status)

# dat$medication_class <- unlist(lapply(
#   dat$medication_class,
#   FUN = function(x) {
#     strsplit(x, "\\(")[[1]][1]
#   }
# ))
library(stringr)
dat$medication_class=str_replace(dat$medication_class, " \\s*\\([^\\)]+\\)", "")

#dat[dat$disease_status == "CD" & is.na(dat$medication_class), "sanger_sample_id"]
table(is.na(dat$medication_class))
#dat$medication_class[dat$disease_status == "CD" & is.na(dat$medication_class)] <- "Unknown"
#dat$medication_class[is.na(dat$medication_class)] <- "Healthy"

#dat$medication_class <- gsub("Anti-TNF+Immunosuppressant", "Anti-TNF + Immunosuppressants", dat$medication_class)
table(dat$medication_class, dat$disease_status)

dat[dat$medication_class=="Anti-TNF+Immunosuppressant", "medication_class"] <- "Anti-TNF + Immunosuppressants"
dat[dat$medication_class=="Anti-TNF ", "medication_class"] <- "Anti-TNF"
dat[dat$medication_class=="Immunosuppressants ", "medication_class"] <- "Immunosuppressants"
dat[dat$medication_class=="Corticosteriods ", "medication_class"] <- "Corticosteriods"
#dat[dat$medication_class=="No Meds", "medication_class"] <- "No medication"
#dat[dat$medication_class=="Anti-TNF (adalimumab-infliximab-golimumab)", "medication_class"] <- "Anti-TNF"
#dat[dat$medication_class=="Immunosuppressants (azathioprine-mercaptopurine-methotrexate)", "medication_class"] <- "Immunosuppressants"
#dat[dat$medication_class=="5-ASA drugs ", "medication_class"] <- "Other medication"
#dat[dat$medication_class=="Corticosteriods (prednisolone-budesonide)", "medication_class"] <- "Other medication"
#dat[dat$medication_class=="Corticosteriods ", "medication_class"] <- "Other medication"
#dat[dat$medication_class=="Ustekinumab", "medication_class"] <- "Other medication"
table(dat$medication_class, dat$disease_status)


# Make table 1
## traits to summarize
var_factor <- c(
  "disease_status",
  "sex",
  "medication_status",
  "medication_class"
  
  # "medications_details",
  # "bead_version",
)
var_continuous <- c(
  # "hours_to_chromium_processing",
  "age"
)
var <- c(var_factor, var_continuous)

dat$disease_status <-  factor(dat$disease_status, levels=c("Healthy", "CD"), labels=c("Healthy", "CD"))

table(dat$medication_class)

dat$medication_class <-  factor(dat$medication_class, levels=c(
  "No Meds",
  "Anti-TNF",
  "Immunosuppressants", 
  "Anti-TNF + Immunosuppressants", 
  "Corticosteriods",
 "5-ASA drugs",
 "Ustekinumab",
"Healthy")   )

dat[dat$medication_class=="No Meds","medication_status"] <- "No medication"
dat[is.na(dat$medication_status),"medication_status" ] <- "On medication"

dat[dat$medication_class=="No Meds","medication_class"] <- NA
                          
