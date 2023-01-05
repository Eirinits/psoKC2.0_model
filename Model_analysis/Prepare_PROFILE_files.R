library(biomaRt)

################# PROFILE pipeline preparations ########################## 
#PROFILE pipeline: https://github.com/sysbio-curie/PROFILE
#Prep files: 
#make "HUGO_Entrez.txt" file:
#Get ensembl-id:
psoRNAseq <- read_excel("Data_analysis/data/psoRNAseq 1.xlsx")
psoRNAseq_Ensemble_ID = as.character(row.names(psoRNAseq)) 
psoRNAseq_Ensemble = gsub("\\..*", "", psoRNAseq_Ensemble_ID)

#psoRNAseq_Symbol_Entrez = bitr(psoRNAseq_Ensemble, fromType = "ENSEMBL",
#                    toType = c("SYMBOL", "ENTREZID"), OrgDb = "org.Hs.eg.db") 


biomart <- useEnsembl(biomart ="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", mirror = "useast")
HUGO_Entrez = getBM(mart = biomart,
                    attributes=c("ensembl_gene_id","hgnc_id","hgnc_symbol", "entrezgene_id"), #PROFILE expects this
                    filter = "ensembl_gene_id", values = psoRNAseq_Ensemble, uniqueRows = T)
HUGO_Entrez = data.frame(
  psoRNAseq_Ensemble_ID[match(HUGO_Entrez$ensembl_gene_id, psoRNAseq_Ensemble)],
  HUGO_Entrez)
colnames(HUGO_Entrez) = c("ensemble_id", c("ensembl_gene_id", " HGNC ID", "Approved Symbol", "Entrez Gene ID"))


HUGO_Entrez_IDs = HUGO_Entrez 
HUGO_Entrez_IDs$` HGNC ID` = gsub(".*:","",HUGO_Entrez_IDs$` HGNC ID`) #hgnc-id

HUGO_Entrez_removed = HUGO_Entrez_IDs[!(is.na(HUGO_Entrez_IDs$`Entrez Gene ID`) & (HUGO_Entrez_IDs$` HGNC ID`=="") 
                                        & (HUGO_Entrez_IDs$`Approved Symbol`=="")), ]
#saving:
all_symbols = HUGO_Entrez_removed

###RNA-expression dataframe (patients):
#check for negatives
#code from:https://www.geeksforgeeks.org/select-rows-from-r-dataframe-that-contain-both-positive-and-negative-values/
negatives = subset(psoRNAseq_corrected,(rowSums(sign(psoRNAseq_corrected)<0)>0) & (rowSums(sign(psoRNAseq_corrected)>0)>0))
#code from: https://www.youtube.com/watch?v=zm9Iv6C4k60
psoRNAseq_corrected_1 = psoRNAseq_corrected
psoriasis_cohort = psoRNAseq_corrected_1[ , colnames(psoRNAseq_corrected_1) %in% rownames(meta_les_groups)]

#Keep rows in psoriasis_cohort that exist in HUGO_Entrez_removed:
psoriasis_cohort_in = psoriasis_cohort[rownames(psoriasis_cohort) %in% HUGO_Entrez_removed$ensemble_id, ]

#Approved symbols:
#code from: https://www.youtube.com/watch?v=6FjgGxm9sww 
psoriasis_cohort_symbol = psoriasis_cohort_in
psoriasis_cohort_symbol$Hugo_Symbol = HUGO_Entrez_removed$`Approved Symbol`[match(rownames(psoriasis_cohort_symbol), HUGO_Entrez_removed$ensemble_id)]

#Entrez-id:
psoriasis_cohort_symbol$Entrez_Gene_Id = HUGO_Entrez_removed$`Entrez Gene ID`[match(rownames(psoriasis_cohort_symbol), HUGO_Entrez_removed$ensemble_id)]

#Entrez_Gene_Id first column. Code from: https://statisticsglobe.com/move-column-to-first-position-of-data-frame-in-r
psoriasis_cohort_symbol_reordering = psoriasis_cohort_symbol %>% 
  dplyr::select("Entrez_Gene_Id", everything())
#HUGO_Symbol first column:
psoriasis_cohort_symbol_reordering_Hugo = psoriasis_cohort_symbol_reordering %>% 
  dplyr::select("Hugo_Symbol", everything())

#Pipeline input (patients):
psoriasis_RNA_Seq_expression = data.frame(psoriasis_cohort_symbol_reordering_Hugo, row.names = NULL)

#RNA-expression (controls):
#code line from: https://stackoverflow.com/questions/24104445/subsetting-a-large-txt-file-before-reading-it-into-the-variable-in-r?lq=1
meta_healthy = meta_RNAseq[meta_RNAseq$lesional == "healthy_control", ]
healthy_cohort = psoRNAseq_corrected_1[ , colnames(psoRNAseq_corrected_1) %in% rownames(meta_healthy)]
#keep rows in healthy_cohort that exist in HUGO_Entrez_removed:
healthy_cohort_in = healthy_cohort[rownames(healthy_cohort) %in% HUGO_Entrez_removed$ensemble_id, ]

#Approved symbols. code from: https://www.youtube.com/watch?v=6FjgGxm9sww 
healthy_cohort_symbol = healthy_cohort_in
healthy_cohort_symbol$Hugo_Symbol= HUGO_Entrez_removed$`Approved Symbol`[match(rownames(healthy_cohort_symbol), HUGO_Entrez_removed$ensemble_id)]
#Entrez-id:
healthy_cohort_symbol$Entrez_Gene_Id = HUGO_Entrez_removed$`Entrez Gene ID`[match(rownames(healthy_cohort_symbol), HUGO_Entrez_removed$ensemble_id)]

#Entrez_Gene_Id first column. Code from: https://statisticsglobe.com/move-column-to-first-position-of-data-frame-in-r
healthy_cohort_symbol_reordering = healthy_cohort_symbol %>%
  dplyr::select("Entrez_Gene_Id", everything())
#HUGO_Symbol first column:
healthy_cohort_symbol_reordering_Hugo = healthy_cohort_symbol_reordering %>%
  dplyr::select("Hugo_Symbol", everything())

#Pipeline input (controls):
healthy_RNA_Seq_expression = data.frame(healthy_cohort_symbol_reordering_Hugo, row.names = NULL)

#Remove "ensemble_Id" + "ensemble_gene_id":
HUGO_Entrez_removed$ensemble_id = NULL
HUGO_Entrez_removed$ensembl_gene_id = NULL

#Re-save model nodes:
setwd("~/Desktop/Laurence_approach")
library(readxl)
psoKC_namesToHugo_curated = read.delim("psoKC_v2.0/psoKC_namesToHugo_curated.txt")
setwd("~/Desktop/Laurence_approach/Data")
#code: https://www.rdocumentation.org/packages/pgirmess/versions/1.7.0/topics/write.delim
#write.delim(psoKC_namesToHugo_curated, file = "psoKC_namesToHugo_curated.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#Save HUGO_Entrez_removed: 
setwd("~/Desktop/Laurence_approach/Data/Try_negative")
#write.delim(HUGO_Entrez_removed, file = "HUGO_Entrez.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#Save healthy_RNA_Seq_expression:
setwd("~/Desktop/Laurence_approach/Data/Try_negative")
#write.delim(healthy_RNA_Seq_expression, file = "healthy_RNA_Seq_expression.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#Save psoriasis_RNA_Seq_expression:
setwd("~/Desktop/Laurence_approach/Data/Try_negative")
#write.delim(psoriasis_RNA_Seq_expression, file = "psoriasis_RNA_Seq_expression.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#Metadata (patient-id, group-nr):
data_clinical_patient = meta_les_groups
#code from: https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
data_clinical_patient = cbind(rownames(data_clinical_patient), 
                              data.frame(data_clinical_patient, row.names=NULL))
#Remove columns:
data_clinical_patient$GSE = NULL
data_clinical_patient$diagnosis = NULL
data_clinical_patient$lesional = NULL
#Rename columns:
colnames(data_clinical_patient) = c("PATIENT_ID", "GROUP_NR")
#Check the order:
all(data_clinical_patient$PATIENT_ID %in% colnames(psoriasis_cohort_in)) #TRUE
all(data_clinical_patient$PATIENT_ID == colnames(psoriasis_cohort_in)) #TRUE
#Save data_clinical_patient.txt:
setwd("~/Desktop/Laurence_approach/Data")
#write.delim(data_clinical_patient, file = "data_clinical_patient.txt", quote = FALSE, row.names = FALSE, sep = "\t" )

############Solving pipeline error
#Same length by Entrez_Gene_ID:
setwd("~/Desktop/Laurence_approach/Data/original_set")
Original_HUGO_Entrez = read.delim("HUGO_Entrez.txt", header = T, sep = "\t")
Original_psoriasis_RNA_Seq_expression = read.delim("psoriasis_RNA_Seq_expression.txt", header = T, sep = "\t")
#code from https://stackoverflow.com/questions/9126840/delete-rows-with-blank-values-in-one-particular-column
Original_HUGO_Entrez_notblank = Original_HUGO_Entrez[!(Original_HUGO_Entrez$Approved.Symbol == ""), ]
Original_HUGO_Entrez_notblank2 = Original_HUGO_Entrez_notblank[Original_HUGO_Entrez_notblank$Entrez.Gene.ID %in% Original_psoriasis_RNA_Seq_expression$Entrez_Gene_Id, ]
Original_psoriasis_RNA_Seq_expression2 = Original_psoriasis_RNA_Seq_expression[Original_psoriasis_RNA_Seq_expression$Entrez_Gene_Id %in% Original_HUGO_Entrez_notblank2$Entrez.Gene.ID, ]
#write the files: 
setwd("~/Desktop/Laurence_approach/Data/same_length")
#write.delim(Original_psoriasis_RNA_Seq_expression2, file = "Psoriasis_RNA_same.txt", quote = FALSE, row.names = FALSE, sep = "\t")
#write.delim(Original_HUGO_Entrez_notblank2, file = "HUGO_Entrez.txt", quote = FALSE, row.names = FALSE, sep = "\t")