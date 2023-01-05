library(tidyverse)
library(grid)
library(tidyverse)

all_states <- read.csv("Model_analysis/PROFILE/Results/Simulations/all_final_states.csv", sep="")
all_states <- split(all_states, ~Perturbation)

all_states <- lapply(all_states, function(x) x %>%
                       select(-Perturbation) %>% 
                       remove_rownames() %>% 
                       column_to_rownames("Patient")) 

all_states_noWT <- lapply(all_states, function(x) x %>% filter(rownames(.) != "WT"))

PSORIASIS_clin <- read.csv("Model_analysis/Patient_annotation.csv") 

genenames <- read.table("Model_analysis/psoKC_namesToHugo_curated.txt",header=T,sep="\t")

geneindex <- strsplit(as.character(genenames[,2]), split = ",") %>% sapply(function(l){gsub(" ","",l)})
geneindex <- data.frame(V1 = rep(genenames[,1], sapply(geneindex, length)), V2 = unlist(geneindex))
model_nodes_HUGO <- unique(geneindex[,2]) %>% sub("^\\s+", "", .)


WT <- read.csv("Model_analysis/PROFILE/Results/Simulations/untreated_allPatients_phenotypes.csv", row.names = 1)

common_patients <- intersect(substr(rownames(WT),1,12), PSORIASIS_clin$PATIENT)
data_plot_extended <- inner_join(PSORIASIS_clin, rownames_to_column(WT,var="PATIENT") %>%
                                   mutate(PATIENT=substr(PATIENT,1,12)),
                                 by="PATIENT") 
data_plot_extended[is.na(data_plot_extended)] <- 0

data_plot <- gather(data_plot_extended, key = "Phenotype", value = "Score", -c(PATIENT, SUBGROUP,PASI,Sex,Severity)) %>% 
  mutate(Score = as.numeric(Score)) %>% 
  mutate(SUBGROUP = factor(SUBGROUP)) %>% 
  mutate(Phenotype = gsub("_WT","",Phenotype)) %>% 
  mutate(Phenotype = factor(Phenotype, levels = c("Apoptosis","Proliferation","Differentiation","Inflammation","Th1","Th17","Neutrophils","Immune_cells"))) %>% 
  filter(!Phenotype %in% c("Differentiation","Immune_cells"))

GROUPcolours <- c("1" = "lightyellow", "2" = "darkred")

p <-  data_plot %>%
  filter(Phenotype %in% c("Apoptosis","Proliferation","Differentiation","Inflammation")) %>% 
  ggplot() + geom_density(aes(x = Score, y= ..scaled.., fill=SUBGROUP),) + 
  scale_fill_manual(values=GROUPcolours, name="Group", guide = "none") +
  scale_y_continuous(breaks=c( 0.5 ,1)) +
  facet_grid(Phenotype~SUBGROUP, scales = "free") +
  xlim(-0.001,1.05)+
 # ylab("")+
  theme_minimal(base_size = 15) +
  theme(panel.spacing.x = unit(2,"lines"))

print(p)

pdf("~/Documents/Git/Personalized_psoKC_Eir/final_maboss_clusters.pdf")
for (i in 1:length(all_states)) {
  
  simulations <- all_states[[i]]
  sim_name <- names(all_states)[i]
  #Prepare data
  common_patients <- intersect(substr(rownames(simulations),1,12), PSORIASIS_clin$PATIENT)
  data_plot_extended <- inner_join(PSORIASIS_clin, rownames_to_column(simulations,var="PATIENT") %>%
                                     mutate(PATIENT=substr(PATIENT,1,12)),
                                   by="PATIENT") 
  
  data_plot <- gather(data_plot_extended, key = "Phenotype", value = "Score", -c(PATIENT, SUBGROUP, PASI, Sex, Severity)) %>% 
    mutate(SUBGROUP = factor(SUBGROUP)) %>% 
    mutate(Score = as.numeric(Score))
  
  SUBGROUPcolours <- c("1" = "lightyellow", "2" = "darkred")
  
  p <-  data_plot %>%
    ggplot() + geom_density(aes(x = Score, y=..scaled.., fill=SUBGROUP)) + 
    scale_fill_manual(values=SUBGROUPcolours, name="Patient cluster") +
    facet_grid(Phenotype~SUBGROUP, scales = "free") +
    xlim(-0.001,1.05)+
    labs(title = sim_name) +
    theme_minimal() 
  
  print(p)
  ggsave(
    paste(sim_name,"density.png"),
    plot = last_plot(),
    scale = 1,
    height = 25,
    width = 15,
    units = "cm",
  #  units = c("in", "cm", "mm", "px"),
    dpi = 300
  )
}

dev.off()
