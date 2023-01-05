library(tidyverse)
library(PCAtools)
library(limma) 
library(vsn) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(genefilter)
library(openxlsx)
library(factoextra)
library(cluster)
library(randomForest)
library(reshape2)

# This code was produced by Eir Aker.

# Helper functions --------------------------------------------------------

datamerging = function(pathname){
  name_files=list.files(path=pathname, full.names=TRUE) 
  list_of_dat = lapply(name_files, function(x){read.delim(file=x, header=T)}) 
  Reduce(function(x,y) {merge(x,y)}, list_of_dat)
}


# Load data and metadata --------------------------------------------------
annotation_RNASeq <- read.delim("Data_analysis/data/annotation_RNASeq.txt")

psoRNAseq = datamerging("Data_analysis/data/Federico_et_al_data/RNAseq_norm/")
psoRNAseq = column_to_rownames(psoRNAseq,"X") #contains samples for removal

#Make metadata available:

metadata_path = "Data_analysis/data/Federico_et_al_data/Study_metadata/"
meta_list = list.files(path=metadata_path, pattern="*.txt", full.names = T)
meta_list = lapply(meta_list, function(x) {read.delim(file=x, header=T, sep = ";")})
meta_RNAseq= do.call("rbind", lapply(meta_list, as.data.frame))
meta_RNAseq = data.frame(meta_RNAseq, row.names=2)
table(meta_RNAseq$lesional) 

#Make a new psoRNAseq frame which contains the metadata samples 

psoRNAseq = psoRNAseq[ ,colnames(psoRNAseq) %in% rownames(meta_RNAseq)]

matching = match(rownames(meta_RNAseq), colnames(psoRNAseq))
psoRNAseq = psoRNAseq[, matching]

# Check and correct for batch effects -------------------------------------

psoRNAseq_log = log2(psoRNAseq +1) #log2 transformation
pca_not_removed_batch = pca(psoRNAseq_log, metadata=meta_RNAseq, center = TRUE, #normalized + log-transformed + scaled
                            scale = TRUE)

biplot(pca_not_removed_batch, lab = paste0(pca_not_removed_batch$metadata$GSE), 
       vline = 0, colby = "GSE", legendPosition = "right", hline = 0) 

biplot(pca_not_removed_batch, lab = paste0(pca_not_removed_batch$metadata$lesional), 
       vline = 0, colby = "lesional", legendPosition = "right", hline =0) 

batch_corrected = removeBatchEffect(psoRNAseq_log, meta_RNAseq$GSE) #GSE study number as a batch
psoRNAseq_corrected = as.data.frame(batch_corrected) #normalized + log2 + batch corrected

norm_matrix = data.matrix(psoRNAseq) 
norm_log_batch = data.matrix(psoRNAseq_corrected)
meanSdPlot(norm_matrix)#normalized
meanSdPlot(norm_log_batch)#normalized + log2-transformed + batch-corrected

pcaTools_batch_corrected = pca(psoRNAseq_corrected, metadata=meta_RNAseq, center = TRUE, 
                               scale = TRUE)

biplot(pcaTools_batch_corrected, 
       lab = paste0(pcaTools_batch_corrected$metadata$GSE), 
       vline = 0, legendPosition = "right", hline = 0, colby = "GSE")

biplot(pcaTools_batch_corrected, 
       lab = paste0(pcaTools_batch_corrected$metadata$lesional),  
       vline = 0, legendPosition = "right", hline = 0, colby = "lesional")

#Hierarchical clustering (HC)

scaled_psoRNAseq = scale(psoRNAseq_corrected)
distance = dist(t(scaled_psoRNAseq))
h_clustering = hclust(distance, method = "complete")
plot(h_clustering, labels =F)


# Differential expression analysis ----------------------------------------

raw_RNAseq = datamerging("Data_analysis/data/Federico_et_al_data/RNAseq_raw/")
raw_RNAseq = data.frame(raw_RNAseq, row.names=1) 

#Remove treated samples: 
raw_RNAseq = raw_RNAseq[ ,colnames(raw_RNAseq) %in% rownames(meta_RNAseq)]
#Order the samples:
matching2 = match(rownames(meta_RNAseq), colnames(raw_RNAseq))
raw_RNAseq = raw_RNAseq[, matching2] 

library(edgeR)

#Remove "non_lesional" from lesional column in metadata:
meta_RNAseq_edgeR = meta_RNAseq[!(meta_RNAseq$lesional == "non_lesional"), ]
table(meta_RNAseq_edgeR$lesional)
#subset raw data to include same samples as meta_RNAseq_edgeR
raw_RNAseq_EdgeR = raw_RNAseq[ ,colnames(raw_RNAseq) %in% rownames(meta_RNAseq_edgeR)]
#Check if same ordering between raw data and meta data:
all(rownames(meta_RNAseq_edgeR) %in% colnames(raw_RNAseq_EdgeR)) #TRUE
all(rownames(meta_RNAseq_edgeR) == colnames(raw_RNAseq_EdgeR))#TRUE
disease_group = factor(meta_RNAseq_edgeR$lesional)#healhty_control, lesional
correction_batch = factor(meta_RNAseq_edgeR$GSE)
healthy_lesional_DGE = DGEList(counts = raw_RNAseq_EdgeR, group = disease_group)
healthy_lesional_DGE = calcNormFactors(healthy_lesional_DGE) #TMM normalization
healthy_lesional_DGE$samples$group = relevel(healthy_lesional_DGE$samples$group, ref = "healthy_control")
healthy_lesional_DGE$samples$group
EdgeR_healthy_design = model.matrix(~0 + correction_batch + healthy_lesional_DGE$samples$group)
healthy_lesional_DGE = estimateDisp(healthy_lesional_DGE, EdgeR_healthy_design)

#likelihood ratio tests:
healthy_fitting_DGE = glmFit(healthy_lesional_DGE, EdgeR_healthy_design)
likelihood_ratio_test = glmLRT(healthy_fitting_DGE)
#Results:
result_lesional_vs_healthy_edgeR= likelihood_ratio_test$table 
sign_results_edgeR = decideTestsDGE(likelihood_ratio_test, adjust.method = "BH", p.value = 0.05, lfc = 1)
DEGs_lesional_vs_healthy_edgeR = data.frame("Ensembl" = row.names(result_lesional_vs_healthy_edgeR)[which(sign_results_edgeR != 0)])
DEGs_lesional_vs_healthy_edgeR$symbol = annotation_RNASeq$Gene.Symbol[match(DEGs_lesional_vs_healthy_edgeR$Ensembl, annotation_RNASeq$Transcript.ID)]
#DEGs psoriasis vs healthy edgeR 
DEGs_psoriasis_vs_healthy = result_lesional_vs_healthy_edgeR[row.names(result_lesional_vs_healthy_edgeR) %in% DEGs_lesional_vs_healthy_edgeR$Ensembl, ]
DEGs_psoriasis_vs_healthy$symbol = annotation_RNASeq$Gene.Symbol[match(rownames(DEGs_psoriasis_vs_healthy), annotation_RNASeq$Transcript.ID)]

underexpressed_DEGs <- DEGs_psoriasis_vs_healthy %>% 
  filter(logFC < 0, PValue < 0.05)

overexpressed_DEGs <- DEGs_psoriasis_vs_healthy %>% 
  filter(logFC > 0, PValue < 0.05)

# Pathway enrichment ------------------------------------------------------

keytypes(org.Hs.eg.db)

genes_ensemble_under = gsub("\\..*", "", rownames(underexpressed_DEGs))

genes_entrez_under = bitr(genes_ensemble_under, fromType = "ENSEMBL",
                    toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

Path_enrich_under = enrichPathway(gene = genes_entrez_under$ENTREZID, pvalueCutoff =  0.05, 
                            readable = T)
head(Path_enrich) 

#Pathway enrichment (over-expressed genes). Reactome:
genes_ensemble_over = gsub("\\..*", "", rownames(overexpressed_DEGs))
genes_entrez_over = bitr(genes_ensemble_over, fromType = "ENSEMBL",
                         toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
Path_enrich_over = enrichPathway(gene = genes_entrez_over$ENTREZID, pvalueCutoff =  0.05, 
                            readable = T)
head(Path_enrich_over)
#enriched_over = data.frame(Path_enrich_over@result[["Description"]], Path_enrich_over@result[["p.adjust"]])

barplot(Path_enrich_over, showCategory = 25) 
barplot(Path_enrich_under, showCategory = 25) 

meta_les = meta_RNAseq[meta_RNAseq$lesional=="lesional" & is.na(meta_RNAseq$treatment) ,]
psoRNAseq_lesional_notscaled = psoRNAseq_corrected[, colnames(psoRNAseq_corrected) %in% rownames(meta_les)]
psoRNAseq_lesional=scale(psoRNAseq_lesional_notscaled) #batch-corrected + scaled + log2 + normalized


# Patient clustering ------------------------------------------------------


# Select the most 200 variable genes: 
variation_200_lesionals= head(order(rowVars(psoRNAseq_lesional), decreasing = T), 200) 
tops_lesional = psoRNAseq_lesional[variation_200_lesionals, ]
tops_lesional = tops_lesional - rowMeans(tops_lesional)
variable_200 = as.data.frame(rownames(tops_lesional)) #kun 200 gener

#Gene symbols: 
genes_ensemble_var200 = gsub("\\..*", "", variable_200$`rownames(tops_lesional)`)
genes_entrez_var200 = bitr(genes_ensemble_var200, fromType = "ENSEMBL",
                          toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
#Removing sex-related genes:
sex_specific = data.frame(c("ENSG00000129824.16", "ENSG00000229807.12",
                                            "ENSG00000067048.17", "ENSG00000012817.15", "ENSG00000114374.13", "ENSG00000183878.15",
                                            "ENSG00000131002.12","ENSG00000067646.12","ENSG00000099725.14","ENSG00000198692.10",
                                            "ENSG00000165246.14", "ENSG00000176728.9"))
variable_lesionals = filter(variable_200, !(rownames(tops_lesional) %in% sex_specific)) 

#Get symbols:
variable_lesionals_symbol = variable_lesionals
variable_lesionals_symbol$symbol = annotation_RNASeq$Gene.Symbol[match(variable_lesionals_symbol$`rownames(tops_lesional)`, annotation_RNASeq$Transcript.ID)]

#Gender genes:

colnames(sex_specific) = c("ENSEMBL")
sex_specific$symbol = annotation_RNASeq$Gene.Symbol[match(sex_specific$ENSEMBL, annotation_RNASeq$Transcript.ID)]

#variable_lesionals er kun 188 gener
###188 variable genes: 
var188_lesional = psoRNAseq_lesional[rownames(psoRNAseq_lesional) %in% variable_lesionals$`rownames(tops_lesional)`, ]

#Silhouette + Elbow for grouping-process: 
set.seed(234)
fviz_nbclust(t(var188_lesional), kmeans, method = "wss") #4 groups (Elbow)
fviz_nbclust(t(var188_lesional), kmeans, method = "silhouette") #2 groups (Silhouette)
 
#HC for subgrouping:
#Code from: https://compgenomr.github.io/book/clustering-grouping-samples-based-on-their-similarity.html 
#and from: https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html#Lesson_Objectives 
distance_lesional_var188 = dist(t(var188_lesional))
h_clustering_lesional_var188 = hclust(distance_lesional_var188 , method = "ward.D2")
plot(h_clustering_lesional_var188, labels =F)

set.seed(234)
RF_var188_les = randomForest(x = t(var188_lesional), ntree = 2000, proximity = T)
RF_var188_les
#PAM subgrouping:
get_proximity_188var = RF_var188_les$proximity
pam.RF_var188_les = pam(get_proximity_188var, k=2) #2 groups
#Subgroups:
group_RF_var188 = as.data.frame(pam.RF_var188_les[["clustering"]])
table(group_RF_var188$`pam.RF_var188_les[["clustering"]]`) 

#Save PAM subgroups:
#setwd("~/Desktop/Files_Eirs_supplementary")
#RF_PAM_subgroups = group_RF_var188
#code from: https://stackoverflow.com/questions/29511215/convert-row-names-into-first-column
#RF_PAM_subgroups = cbind(rownames(RF_PAM_subgroups),
#                         data.frame(RF_PAM_subgroups, row.names=NULL))
#colnames(RF_PAM_subgroups) = c("PATIENT","SUBGROUP")
#write.csv(RF_PAM_subgroups,"RF_PAM_subgroups.csv")

#Save matrix for clustering:
#setwd("~/Desktop/Files_Eirs_supplementary")
#psoriasis_188variable = data.frame(var188_lesional)
#write.csv(psoriasis_188variable, "psoriasis_188variable.csv")

#extract subgroups:
groupings1_RF_var188= as.data.frame(group_RF_var188[group_RF_var188$`pam.RF_var188_les[["clustering"]]`== 1 , , drop=FALSE])
groupings2_RF_var188= as.data.frame(group_RF_var188[group_RF_var188$`pam.RF_var188_les[["clustering"]]`== 2 , , drop=FALSE])
#group 1:
G1_les_var188 = var188_lesional[ ,colnames(var188_lesional) %in% rownames(groupings1_RF_var188)] #83 samples
#group 2:
G2_les_var188 = var188_lesional[, colnames(var188_lesional) %in% rownames(groupings2_RF_var188)] #102 samples

######################    Boxplot analysis  #########################
#Look at DEGs:
##melt and boxplot code from: https://github.com/hamidghaedi/RNA-seq-differential-expression
#and from: https://www.biostars.org/p/469559/
var188_lesional

melt_lesional_var188 = data.frame(melt(var188_lesional))
colnames(melt_lesional_var188) = c("gene_name", "patient_sample", "count_values")
melt_lesional_var188$group_nr = ifelse(melt_lesional_var188$patient_sample %in% rownames(groupings1_RF_var188), "1", "2")

check= melt_lesional_var188 %>%   #checking correct grouping
  group_by(patient_sample, group_nr) %>%
  summarize(numbers=n())

all(check$patient_sample %in% rownames(group_RF_var188)) #TRUE
all(check$patient_sample == rownames(group_RF_var188)) #TRUE
all(check$group_nr == group_RF_var188$`pam.RF_var188_les[["clustering"]]`) #TRUE
table(check$group_nr) #correct
table(group_RF_var188$`pam.RF_var188_les[["clustering"]]`) #correct

        #over-expressed
#60 over-expressed genes:
over60_exp=head(DEG_result_0.05padj[order(DEG_result_0.05padj$log2FoldChange, decreasing = TRUE), ], 60)
over60_DEGs = data.frame(over60_exp@rownames) 
#Matrix of over-expressed:
var188groups_over60DEGs = psoRNAseq_lesional_notscaled[rownames(psoRNAseq_lesional_notscaled) %in% over60_DEGs$over60_exp.rownames, ]

var188groups_over60DEGs$symbol = annotation_RNA$Gene.Symbol[match(rownames(var188groups_over60DEGs), annotation_RNA$Transcript.ID)]
var188groups_over60DEGs = data.frame(var188groups_over60DEGs, row.names ="symbol")
var188groups_over60DEGs = as.matrix(var188groups_over60DEGs) #matrix for melt

melt_var188groups_over60DEGs = data.frame(melt(var188groups_over60DEGs))
colnames(melt_var188groups_over60DEGs) = c("gene_name", "patient_sample", "count_value")
#adding group-nr (1 or 2):
melt_var188groups_over60DEGs$Subgroup = ifelse(melt_var188groups_over60DEGs$patient_sample %in% rownames(groupings1_RF_var188), "1", "2")
#code from: https://www.biostars.org/p/469559/
ggplot(melt_var188groups_over60DEGs, aes(x = as.factor(gene_name), y = count_value)) +
  geom_boxplot(aes(fill = Subgroup), position = position_dodge(0.9), notch = T) +
  scale_fill_manual(values = c("red", "green"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(y = "Expression value", x = "Gene symbol")

    #under-expressed
#60 under-expressed genes:
under60_exp= head(DEG_result_0.05padj[order(DEG_result_0.05padj$log2FoldChange), ], 60)
under_60DEGs = data.frame(under60_exp@rownames)
#Matrix 60 under-expressed:
var188groups_under_60DEGs=psoRNAseq_lesional_notscaled[rownames(psoRNAseq_lesional_notscaled) %in% under_60DEGs$under60_exp.rownames, ]  
#Gene symbols:
var188groups_under_60DEGs$symbol = annotation_RNA$Gene.Symbol[match(rownames(var188groups_under_60DEGs), annotation_RNA$Transcript.ID)]
#Symbol to rowname
var188groups_under_60DEGs = data.frame(var188groups_under_60DEGs, row.names = "symbol")
var188groups_under_60DEGs = as.matrix(var188groups_under_60DEGs)#matrix for melt

melt_var188groups_under_60DEGs= data.frame(melt(var188groups_under_60DEGs)) 
colnames(melt_var188groups_under_60DEGs) = c("gene_name", "patient_sample", "count_value")
#add group-nr (1 or 2):
melt_var188groups_under_60DEGs$Subgroup = ifelse(melt_var188groups_under_60DEGs$patient_sample %in% rownames(groupings1_RF_var188), "1", "2")
#boxplot:
ggplot(melt_var188groups_under_60DEGs, aes(x = as.factor(gene_name), y = count_value)) +
  geom_boxplot(aes(fill = Subgroup), position = position_dodge(0.9), notch = T) +
  scale_fill_manual(values = c("red", "green"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(y = "Expression value", x = "Gene symbol")
# DEA between the two clusters --------------------------------------------

#DIFFERENTIAL EXPRESSION ANALYSIS(G1 & G2) 
#DESeq: 
#Matrix of lesionals (raw):
raw_lesionals = raw_RNAseq[, colnames(raw_RNAseq) %in% colnames(psoRNAseq_lesional)]
#Metadata (study+group):
meta_les_groups= meta_les
meta_les_groups$group_nr = ifelse(rownames(meta_les_groups) %in% rownames(groupings1_RF_var188), "G1", "G2")

all(rownames(meta_les_groups) %in% colnames(raw_lesionals)) #TRUE
all(rownames(meta_les_groups) == colnames(raw_lesionals)) #TRUE

all(rownames(meta_les_groups) == colnames(raw_lesionals)) #TRUE
psoriasis_group = factor(meta_les_groups$group_nr) #G1 and G2
correct_batch = factor(meta_les_groups$GSE) 
psoriasis_DGE = DGEList(counts = raw_lesionals, group = psoriasis_group) 
psoriasis_DGE = calcNormFactors(psoriasis_DGE) #TMM normalization
psoriasis_DGE$samples$group = relevel(psoriasis_DGE$samples$group, ref = "G1")
psoriasis_DGE$samples$group #baseline = G1
EdgeR_design = model.matrix(~0 + correct_batch + psoriasis_DGE$samples$group)
psoriasis_DGE = estimateDisp(psoriasis_DGE, EdgeR_design)
#Likelihood ratio tests:
fitting_DGE = glmFit(psoriasis_DGE, EdgeR_design)
L_ratio_test = glmLRT(fitting_DGE)
#Results:
result_G2vsG1_edgeR = L_ratio_test$table
significant_results = decideTestsDGE(L_ratio_test, adjust.method = "BH", p.value = 0.05, lfc = 1)
DEGs_G2vsG1_edgeR = data.frame(row.names(result_G2vsG1_edgeR)[which(significant_results != 0)])
DEGs_G2vsG1_edgeR$symbol = annotation_RNASeq$Gene.Symbol[match(DEGs_G2vsG1_edgeR$row.names.result_G2vsG1_edgeR..which.significant_results....0.., annotation_RNASeq$Transcript.ID)]

###
ncol(EdgeR_design)
qr(EdgeR_design)$rank
###

DEGs_subgroup2_1 = result_G2vsG1_edgeR[row.names(result_G2vsG1_edgeR) %in% DEGs_G2vsG1_edgeR$row.names.result_G2vsG1_edgeR..which.significant_results....0.., ]
DEGs_subgroup2_1$symbol = annotation_RNASeq$Gene.Symbol[match(rownames(DEGs_subgroup2_1), annotation_RNASeq$Transcript.ID)]

#over-expressed (G2):
edgeR_G2vsG1_up = data.frame(row.names(result_G2vsG1_edgeR)[which(significant_results == 1)])
edgeR_G2vsG1_up$symbol = annotation_RNASeq$Gene.Symbol[match(edgeR_G2vsG1_up$row.names.result_G2vsG1_edgeR..which.significant_results....1.., annotation_RNASeq$Transcript.ID)]
#under-expressed (G2):
edgeR_G2vsG1_down = data.frame(row.names(result_G2vsG1_edgeR)[which(significant_results == -1)])
edgeR_G2vsG1_down$symbol = annotation_RNASeq$Gene.Symbol[match(edgeR_G2vsG1_down$row.names.result_G2vsG1_edgeR..which.significant_results.....1.., annotation_RNASeq$Transcript.ID)]


entrez_up_edgeR_G2vsG1 = bitr(edgeR_G2vsG1_up$symbol, fromType = "SYMBOL", 
                           toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
up_pathways_edgeR_G2vsG1 = enrichPathway(gene = entrez_up_edgeR_G2vsG1$ENTREZID, pvalueCutoff =  0.05,
                                         readable = T)
head(up_pathways_edgeR_G2vsG1)
barplot(up_pathways_edgeR_G2vsG1, showCategory = 10)
path_genes_up_reactome= data.frame(up_pathways_edgeR_G2vsG1@result[["Description"]], up_pathways_edgeR_G2vsG1@result[["geneID"]], up_pathways_edgeR_G2vsG1@result[["p.adjust"]], up_pathways_edgeR_G2vsG1@result[["GeneRatio"]])
colnames(path_genes_up_reactome) = c("pathway","genes", "p.adjust", "gene ratio")
#save for supplementary:
#setwd("~/Desktop/Files_Eirs_supplementary")
#write.csv(path_genes_up_reactome, "upregulated_pathways.csv")

#Reactome down-regulated pathways:
entrez_down_edgeR_G2vsG1  = bitr(edgeR_G2vsG1_down$symbol, fromType = "SYMBOL", 
                                 toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
down_pathways_edgeR_G2vsG1 = enrichPathway(gene = entrez_down_edgeR_G2vsG1$ENTREZID, pvalueCutoff =  0.05,
                                           readable = T)
head(down_pathways_edgeR_G2vsG1)
barplot(down_pathways_edgeR_G2vsG1)
path_genes_down_reactome = data.frame(down_pathways_edgeR_G2vsG1@result[["Description"]], down_pathways_edgeR_G2vsG1@result[["geneID"]], down_pathways_edgeR_G2vsG1@result[["p.adjust"]], down_pathways_edgeR_G2vsG1@result[["GeneRatio"]])
colnames(path_genes_down_reactome) = c("pathways", "genes", "p.adjust", "gene ratio")

#Reactome pathway enrichment:
entrez_allDEGs_edgeR_G2vsG1 = bitr(DEGs_G2vsG1_edgeR$symbol, fromType = "SYMBOL",
                                   toType = "ENTREZID", OrgDb= "org.Hs.eg.db")
allDEGs_pathways_edgeR_G2vsG1 = enrichPathway(gene = entrez_allDEGs_edgeR_G2vsG1$ENTREZID, pvalueCutoff = 0.05,
                                              readable = T)
head(allDEGs_pathways_edgeR_G2vsG1)
barplot(allDEGs_pathways_edgeR_G2vsG1)
