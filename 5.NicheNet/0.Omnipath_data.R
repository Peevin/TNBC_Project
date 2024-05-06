library(OmnipathR)
library(nichenetr)
library(tidyverse)
library(mlrMBO)
library(parallelMap)

setwd('/Users/liupeiwen/NicheNet/')

# ----------------------Ligands to Receptors---------------------
# We first define a function to transform the format of interactions in the Omnipath to the one used in the NicheNet method:
interactionFormatTransf <- function(InputDf, InteractionType){
  
  OutputInt <- tibble(from = character(), to = character(), 
                      source = character(), database = character())  
  
  n <- nrow(InputDf)
  sources <- dplyr::pull(InputDf, sources)
  sourceNodes <- dplyr::pull(InputDf, from)
  targetNodes <- dplyr::pull(InputDf, to)
  
  for (i in seq(n)){
    currentSources <- unlist(strsplit(sources[i],";"))
    for (j in seq(length(currentSources))){
      OutputInt <- add_row(OutputInt, 
                           from = sourceNodes[i] , 
                           to = targetNodes[i],  
                           # source = paste(currentSources[j], InteractionType, sep="_"),
                           source = currentSources[j],
                           database = currentSources[j]) 
    }
  }
  
  return(OutputInt)
}

# Generating the Omnipath ligand-receptor network
# Omnipath possess a dedicated dataset storing these type of interactions ( LigrecExtra ). Therefore, we get these interactions:
## We remove self-interacitons and duplicated records
lr_Interactions_Omnipath <- import_ligrecextra_interactions() %>%
  dplyr::select(source_genesymbol,target_genesymbol,sources) %>%
  dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>% 
  dplyr::filter(from != to) %>% 
  dplyr::distinct()


## We import Omnipath Inter cellular annotations
InterCell_Annotations <- import_omnipath_intercell() 

## We filter those proteins which are mainly annotated as receptor or ligand
Ligands_Receptors <- InterCell_Annotations %>%
  dplyr::filter(category %in% c("receptor","ligand"))

## There are also some complexes. We are going to deal with them by including
## each of its individual proteins in our list
Ligand_Receptors_class <- character()
Ligand_Receptors_name <- character()
for (i in seq(nrow(Ligands_Receptors))){
  if (Ligands_Receptors$entity_type[i] == "complex"){
    Genescomplex <-unlist(strsplit(gsub("COMPLEX:", "", 
                                        Ligands_Receptors$genesymbol[i]),"_"))
    class <- rep(Ligands_Receptors$category[i],length(Genescomplex))
    Ligand_Receptors_name <- c(Ligand_Receptors_name,Genescomplex)
    Ligand_Receptors_class <- c(Ligand_Receptors_class,class)
    
  } else {
    Ligand_Receptors_name <- 
      c(Ligand_Receptors_name, Ligands_Receptors$genesymbol[i]) 
    Ligand_Receptors_class <- 
      c(Ligand_Receptors_class, Ligands_Receptors$category[i]) 
  }
}

## We create a vector with all the ligands and another with all the receptors.
Ligand_Receptors_df <- data.frame(GeneSymbol = Ligand_Receptors_name, 
                                  Class = Ligand_Receptors_class, stringsAsFactors = FALSE) %>%
  dplyr::distinct()
AllLigands_vec <- 
  dplyr::filter(Ligand_Receptors_df, Class == "ligand") %>%
  dplyr::pull(GeneSymbol)
AllReceptors_vec <- 
  dplyr::filter(Ligand_Receptors_df, Class == "receptor") %>%
  dplyr::pull(GeneSymbol)

## We next get protein-protein interactions from the different datasets availabe
## in Omnipath
AllInteractions <- 
  import_post_translational_interactions(exclude = "ligrecextra") %>% 
  dplyr::select(source_genesymbol, target_genesymbol, sources) %>% 
  dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>% 
  dplyr::filter(from != to) %>% 
  dplyr::distinct() 

## I finally match interactions and annotations.
Matching_Interactions_Annotations <- AllInteractions %>%
  dplyr::filter(from %in% AllLigands_vec) %>%
  dplyr::filter(to %in% AllReceptors_vec) %>%
  dplyr::distinct()


## We access to the Omnipath webservice usign the OmnipathR package and we 
## set the number of sources reporting an interactions as its weight 
lr_Network_Omnipath <- 
  bind_rows(lr_Interactions_Omnipath, Matching_Interactions_Annotations) %>%
  dplyr::distinct() %>%
  interactionFormatTransf(InteractionType="LigrecExtra") %>%
  dplyr::distinct() 

## I have to remove self-interactions
lr_Network_Omnipath <- lr_Network_Omnipath %>% 
  dplyr::filter(from != to) 

## I also have to remove interactions going to ligands. See Methods Nichenet 
## paper
ligands <- unique(dplyr::pull(lr_Network_Omnipath, from))
lr_Network_Omnipath <- lr_Network_Omnipath %>% 
  dplyr::filter(!(to %in% ligands))

## There are in addition some records containing not input gene, we remove them
## since they are giving problems with running the model.
lr_Network_Omnipath <- lr_Network_Omnipath %>% 
  dplyr::filter(from != "") %>% 
  dplyr::filter(to != "")
nrow(lr_Network_Omnipath)

saveRDS(lr_Network_Omnipath, 
        "OmniNetworks_NNformat/lr_Network_Omnipath.rds")




# --------------------Receptore to TF-----------------------
# Generating the Omnipath signalling network
# We generate a signaling network using Omnipath resources. In Omnipath , 
# we can find different datasets describing protein interactions ( https://github.com/saezlab/pypath ).
# We will merge these datasts to generate a signaling network.

## Original Omnipath interactions
sig_Network_Omnipath <- 
  interactionFormatTransf(AllInteractions, InteractionType="Signalling") %>%
  dplyr::distinct() 

## I have to remove self-interactions in the signaling network
sig_Network_Omnipath <- sig_Network_Omnipath %>% 
  dplyr::filter(from != to)

## I also have to remove interactions going to ligands. See Methods Nichenet 
## paper
sig_Network_Omnipath <- sig_Network_Omnipath %>% 
  dplyr::filter(!(to %in% ligands))

## There are in addition some records containing not input gene, we remove them
## since they are giving problems with running the model.
sig_Network_Omnipath <- sig_Network_Omnipath %>% 
  dplyr::filter(from != "") %>% 
  dplyr::filter(to != "")


## We also remove signaling interactions that are already in the lig-receptor 
## network. 
sig_Network_Omnipath <- dplyr::anti_join(
  sig_Network_Omnipath, 
  lr_Network_Omnipath, 
  by = c("from" = "from", "to" = "to"))

nrow(sig_Network_Omnipath)

saveRDS(sig_Network_Omnipath, 
        "OmniNetworks_NNformat/sig_Network_Omnipath.rds")

# Generating the Omnipath gene regulatory network
# In this section, we generate a GRN network using the DoRothEA regulons which are available through Omnipath :
gr_Interactions_Omnipath <- 
  import_dorothea_interactions(dorothea_levels = c("A","B","C")) %>%  
  dplyr::select(source_genesymbol, target_genesymbol, sources) %>%
  dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>% 
  dplyr::filter(from != to) %>%
  dplyr::distinct()  

gr_Network_Omnipath <- 
  interactionFormatTransf(
    gr_Interactions_Omnipath, 
    InteractionType="Dorothea") %>%
  dplyr::distinct() 
nrow(gr_Network_Omnipath)

saveRDS(gr_Network_Omnipath,
        "OmniNetworks_NNformat/gr_Network_Omnipath.rds")

# Parameter optimization via mlrMBO
# We are going to optimice Nichenet method, 
# i.e. PageRank parameters and source weights, 
# basedon a collection of experiments where the effect of a ligand on gene expression was measured. 
# We have to remove the experiments involving the ligand IFNA1 because it is not present in our network.
expression_settings_validation <- 
  readRDS('OmniNetworks_NNformat/expression_settings.rds')

index <- which(!unlist(lapply(expression_settings_validation, 
                              function(x) any(x$from != "IFNA1"))))

expression_settings_validation <- expression_settings_validation[-index]

lr_Network_Omnipath <- readRDS('OmniNetworks_NNformat/lr_Network_Omnipath.rds')
sig_Network_Omnipath <- readRDS('OmniNetworks_NNformat/sig_Network_Omnipath.rds')
gr_Network_Omnipath <- readRDS('OmniNetworks_NNformat/gr_Network_Omnipath.rds')

All_sources <- unique(c(lr_Network_Omnipath$source,
                        sig_Network_Omnipath$source, gr_Network_Omnipath$source))
my_source_weights_df <- 
  tibble(source = All_sources, weight = rep(1,length(All_sources)))

additional_arguments_topology_correction <- 
  list(source_names = my_source_weights_df$source %>% unique(), 
       algorithm = "PPR", 
       correct_topology = FALSE,
       lr_network = lr_Network_Omnipath, 
       sig_network = sig_Network_Omnipath, 
       gr_network = gr_Network_Omnipath, 
       settings = lapply(expression_settings_validation, 
                         convert_expression_settings_evaluation), 
       secondary_targets = FALSE, 
       remove_direct_links = "no", 
       cutoff_method = "quantile")

nr_datasources <- additional_arguments_topology_correction$source_names %>% 
  length()

obj_fun_multi_topology_correction = makeMultiObjectiveFunction(name = "nichenet_optimization",
                                                               description = "data source weight and hyperparameter optimization: expensive black-box function", 
                                                               fn = model_evaluation_optimization, 
                                                               par.set = makeParamSet(
                                                                 makeNumericVectorParam("source_weights", len = nr_datasources, 
                                                                                        lower = 0, upper = 1, tunable = FALSE), 
                                                                 makeNumericVectorParam("lr_hub", len = 1, lower = 0, upper = 1, 
                                                                                        tunable = TRUE),  
                                                                 makeNumericVectorParam("sig_hub", len = 1, lower = 0, upper = 1, 
                                                                                        tunable = TRUE), 
                                                                 makeNumericVectorParam("gr_hub", len = 1, lower = 0, upper = 1, 
                                                                                        tunable = TRUE),  
                                                                 makeNumericVectorParam("ltf_cutoff", len = 1, lower = 0.9, 
                                                                                        upper = 0.999, tunable = TRUE),  
                                                                 makeNumericVectorParam("damping_factor", len = 1, lower = 0.01, 
                                                                                        upper = 0.99, tunable =TRUE)), 
                                                               has.simple.signature = FALSE,
                                                               n.objectives = 4, 
                                                               noisy = FALSE,
                                                               minimize = c(FALSE,FALSE,FALSE,FALSE))

optimization_results = 
  lapply(1,mlrmbo_optimization, obj_fun = obj_fun_multi_topology_correction, 
         niter = 8, ncores = 8, nstart = 160, 
         additional_arguments = additional_arguments_topology_correction)

saveRDS(optimization_results, "OmniNetworks_NNformat/Optimization_results.rds")


setwd('/Users/liupeiwen/BC/New_analysis/5.NicheNet_rec_target/data/')

