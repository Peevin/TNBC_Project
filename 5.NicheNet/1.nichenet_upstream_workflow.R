devtools::install_github('Peevin/Receptor-TF-nichenetr')


library(nichenetr)
library(tidyverse)

setwd('/Users/liupeiwen/NicheNet/')

## We load the networks generated in the ParameterOptimization script
lr_network_Omnipath = readRDS("OmniNetworks_NNformat/lr_Network_Omnipath.rds")
sig_network_Omnipath = readRDS("OmniNetworks_NNformat/sig_Network_Omnipath.rds")
gr_network_Omnipath = readRDS("OmniNetworks_NNformat/gr_Network_Omnipath.rds")

## The interactions are weighted only based in the number of data sources that
## report them
All_sources <- unique(c(lr_network_Omnipath$source,
                        sig_network_Omnipath$source, gr_network_Omnipath$source))

my_source_weights_df <- 
  tibble(source = All_sources, weight = rep(1,length(All_sources)))

weighted_networks <- construct_weighted_networks(
  lr_network = lr_network_Omnipath, 
  sig_network = sig_network_Omnipath, 
  gr_network = gr_network_Omnipath, 
  source_weights_df = my_source_weights_df,
  )

## We read the results of the optimization 
resultsOptimization <- readRDS("OmniNetworks_NNformat/Optimization_results.rds")

optimized_parameters = resultsOptimization %>% 
  process_mlrmbo_nichenet_optimization(
    source_names = my_source_weights_df$source %>% unique())

weighted_networks <- apply_hub_corrections(
  weighted_networks = weighted_networks, 
  lr_sig_hub = optimized_parameters$lr_sig_hub, 
  gr_hub = optimized_parameters$gr_hub)

saveRDS(weighted_networks,"OmniNetworks_NNformat/weighted_networksNonSourceWeights.rds")


# We now generate the matrix containing the ligand-target regulatory potential scores based on the weighted integrated networks.
ligand_signaling_network = weighted_networks$lr_sig
regulatory_network = weighted_networks$gr

allgenes = c(ligand_signaling_network$from, ligand_signaling_network$to, regulatory_network$from, regulatory_network$to) %>% unique() %>% sort()
allgenes_integer = allgenes %>% factor() %>% as.numeric()
allgenes_id_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
id2allgenes = mapper(allgenes_id_tbl,"allgenes_integer","allgenes")

receptor <- unique(lr_network_Omnipath$to)
receptor <- intersect(receptor, names(id2allgenes))
receptor <- as.list(receptor)

ligand_target_matrix <- construct_ligand_target_matrix(
  weighted_networks = weighted_networks, 
  ligands = receptor, 
  algorithm = "PPR", 
  damping_factor = optimized_parameters$damping_factor, 
  ltf_cutoff = optimized_parameters$ltf_cutoff)
saveRDS(ligand_target_matrix,"OmniNetworks_NNformat/Receptor_TF_matrixNoweights.rds")

# Now, we create an alternative model using the optimized weights for the different sources of data.
## Here we also take into account the optimized source weights
weighted_networks <- construct_weighted_networks(
  lr_network = lr_network_Omnipath, 
  sig_network = sig_network_Omnipath, 
  gr_network = gr_network_Omnipath,
  source_weights_df = optimized_parameters$source_weight_df)

weighted_networks <- apply_hub_corrections(
  weighted_networks = weighted_networks, 
  lr_sig_hub = optimized_parameters$lr_sig_hub, 
  gr_hub = optimized_parameters$gr_hub)

ligand_target_matrix = construct_ligand_target_matrix(
  weighted_networks = weighted_networks, 
  ligands = receptor, 
  algorithm = "PPR", 
  damping_factor = optimized_parameters$damping_factor, 
  ltf_cutoff = optimized_parameters$ltf_cutoff)

saveRDS(ligand_target_matrix,"OmniNetworks_NNformat/receptor_TF_matrixWithweights.rds")
saveRDS(weighted_networks,"OmniNetworks_NNformat/weighted_networksWithSourceWeights.rds")

predict_ligand_activities









