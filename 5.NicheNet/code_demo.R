# construct_ligand_tf_matrix = function(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, ligands_as_cols = FALSE) {
#   
#   # input check
#   if (!is.list(weighted_networks))
#     stop("weighted_networks must be a list object")
#   if (!is.data.frame(weighted_networks$lr_sig))
#     stop("lr_sig must be a data frame or tibble object")
#   if (!is.data.frame(weighted_networks$gr))
#     stop("gr must be a data frame or tibble object")
#   
#   if (!is.numeric(weighted_networks$lr_sig$weight))
#     stop("lr_sig must contain a column named data source weights")
#   if (!is.numeric(weighted_networks$gr$weight))
#     stop("gr must contain a column named data source weights")
#   
#   if (!is.list(ligands))
#     stop("ligands must be a list object")
#   if ( sum((unique(unlist(ligands)) %in% unique(c(lr_network$from,lr_network$to))) == FALSE) > 0)
#     warning("One or more ligands of interest not present in the ligand-receptor network 'lr_network'. You can possibly ignore this warning if you provided your own ligand_receptor network to the weighted networks." )
#   
#   if(is.null(ltf_cutoff)){
#     if( algorithm == "PPR" | algorithm == "SPL" )
#       warning("Did you not forget to give a value to ltf_cutoff?")
#   } else {
#     if (ltf_cutoff < 0 | ltf_cutoff > 1)
#       stop("ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
#   }
#   
#   if (algorithm != "PPR" & algorithm != "SPL" & algorithm != "direct")
#     stop("algorithm must be 'PPR' or 'SPL' or 'direct'")
#   if(algorithm == "PPR"){
#     if (damping_factor < 0 | damping_factor >= 1)
#       stop("damping_factor must be a number between 0 and 1 (0 included, 1 not)")
#   }
#   if (!is.logical(ligands_as_cols) | length(ligands_as_cols) != 1)
#     stop("ligands_as_cols must be a logical vector of length 1")
#   
#   requireNamespace("dplyr")
#   
#   # load in weighted networks
#   ligand_signaling_network = weighted_networks$lr_sig
#   regulatory_network = weighted_networks$gr
#   
#   # convert ids to numeric for making Matrix::sparseMatrix later on
#   allgenes = c(ligand_signaling_network$from, ligand_signaling_network$to, regulatory_network$from, regulatory_network$to) %>% unique() %>% sort()
#   allgenes_integer = allgenes %>% factor() %>% as.numeric()
#   allgenes_id_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
#   mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
#   id2allgenes = mapper(allgenes_id_tbl,"allgenes_integer","allgenes")
#   
#   ligand_signaling_network = ligand_signaling_network %>% mutate(from_allgenes = id2allgenes[from], to_allgenes = id2allgenes[to]) %>% arrange(from_allgenes) %>% dplyr::select(from_allgenes,to_allgenes,weight)
#   
#   if (algorithm == "PPR"){
#     # Make Matrix::sparse signaling weighted matrix and graph to apply personalized pagerank
#     ligand_signaling_network_matrix = Matrix::sparseMatrix(ligand_signaling_network$from_allgenes %>% as.integer, ligand_signaling_network$to_allgenes %>% as.integer, x=ligand_signaling_network$weight %>% as.numeric, dims = c(length(allgenes), length(allgenes)))
#     signaling_igraph = igraph::graph_from_adjacency_matrix(ligand_signaling_network_matrix, weighted=TRUE, mode="directed")
#     # personalized pagerank
#     # prepare preference vector E
#     E = rep(0,times = length(igraph::V(signaling_igraph)))
#     # ppr for every ligand individual
#     complete_matrix = lapply(ligands,PPR_wrapper,E,signaling_igraph,damping_factor,id2allgenes,ltf_cutoff)
#   } else if (algorithm == "SPL"){
#     # Make Matrix::sparse signaling weighted matrix and graph to apply shortest path length
#     ligand_signaling_network = ligand_signaling_network %>% mutate(weight = 1/weight) # reverse the weight because spl: find shortest path, with lowest weight, so original weights needed to be reversed
#     ligand_signaling_network_matrix = Matrix::sparseMatrix(ligand_signaling_network$from_allgenes %>% as.integer, ligand_signaling_network$to_allgenes %>% as.integer, x=ligand_signaling_network$weight %>% as.numeric, dims = c(length(allgenes), length(allgenes)))
#     signaling_igraph = igraph::graph_from_adjacency_matrix(ligand_signaling_network_matrix, weighted=TRUE, mode="directed")
#     # SPL
#     complete_matrix = lapply(ligands,SPL_wrapper,signaling_igraph,id2allgenes,ltf_cutoff)
#   } else if (algorithm == "direct"){
#     ligand_ligand_network = tibble::tibble(from_allgenes = id2allgenes[unlist(ligands)], to_allgenes = id2allgenes[unlist(ligands)]) %>% inner_join(ligand_signaling_network %>% filter(from_allgenes %in% id2allgenes[unlist(ligands)]) %>% group_by(from_allgenes) %>% top_n(1,weight) %>% ungroup() %>% distinct(from_allgenes,weight), by = "from_allgenes")
#     ligand_signaling_network = ligand_signaling_network %>% bind_rows(ligand_ligand_network)
#     ligand_signaling_network_matrix = Matrix::sparseMatrix(ligand_signaling_network$from_allgenes %>% as.integer, ligand_signaling_network$to_allgenes %>% as.integer, x=ligand_signaling_network$weight %>% as.numeric, dims = c(length(allgenes), length(allgenes)))
#     signaling_igraph = igraph::graph_from_adjacency_matrix(ligand_signaling_network_matrix, weighted=TRUE, mode="directed")
#     
#     complete_matrix = lapply(ligands,direct_wrapper,ligand_signaling_network_matrix,id2allgenes,ltf_cutoff)
#   }
#   ltf_matrix = matrix(unlist(complete_matrix), ncol = length(igraph::V(signaling_igraph)), byrow = TRUE)
#   
#   rownames(ltf_matrix) = sapply(ligands,function(x){paste0(x,collapse = "-")})
#   colnames(ltf_matrix) = allgenes
#   
#   if (ligands_as_cols == TRUE){
#     tf_matrix = ltf_matrix %>% as.matrix()
#     ltf_matrix = ltf_matrix %>% t()
#   }
#   
#   return(ltf_matrix)
# }

library(dplyr)
ligand_signaling_network = weighted_networks$lr_sig
regulatory_network = weighted_networks$gr
receptor <- unique(lr_network_Omnipath$to)
receptor <- intersect(receptor, names(id2allgenes))
receptor <- as.list(receptor)

# convert ids to numeric for making Matrix::sparseMatrix later on
allgenes = c(ligand_signaling_network$from, ligand_signaling_network$to, regulatory_network$from, regulatory_network$to) %>% unique() %>% sort()
allgenes_integer = allgenes %>% factor() %>% as.numeric()
allgenes_id_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
id2allgenes = mapper(allgenes_id_tbl,"allgenes_integer","allgenes")
names(id2allgenes)

ligand_signaling_network = ligand_signaling_network %>% mutate(from_allgenes = id2allgenes[from], to_allgenes = id2allgenes[to]) %>% arrange(from_allgenes) %>% dplyr::select(from_allgenes,to_allgenes,weight)

ligand_signaling_network_matrix = Matrix::sparseMatrix(ligand_signaling_network$from_allgenes %>% as.integer, ligand_signaling_network$to_allgenes %>% as.integer, x=ligand_signaling_network$weight %>% as.numeric, dims = c(length(allgenes), length(allgenes)))
signaling_igraph = igraph::graph_from_adjacency_matrix(ligand_signaling_network_matrix, weighted=TRUE, mode="directed")
# personalized pagerank
# prepare preference vector E
E = rep(0,times = length(igraph::V(signaling_igraph)))
# ppr for every ligand individual



damping_factor = optimized_parameters$damping_factor
ltf_cutoff = optimized_parameters$ltf_cutoff

PPR_wrapper = function(ligand,E,G,delta,id2allgenes,ppr_cutoff) {
  partial_matrix = lapply(ligand,single_ligand_ppr_wrapper,E,G,delta,id2allgenes)
  ppr_matrix = matrix(unlist(partial_matrix), ncol = length(E), byrow = TRUE)
  if (delta == 0 & is.null(ppr_cutoff)){
    ppr_cutoff = 0
  }
  # put cutoff
  if (ppr_cutoff > 0){
    ppr_matrix_TRUE = apply(ppr_matrix,1,function(x){x <= quantile(x,ppr_cutoff)}) %>% t()
    ppr_matrix[ppr_matrix_TRUE] = 0
  }
  # if multiple ligands: total ligand-tf score = maximum score a particular ligand-tf interaction; if only one ligand: just the normal score
  return(apply(ppr_matrix, 2, mean))
}

single_ligand_ppr_wrapper = function(l,E,G,delta,id2allgenes){
  # set ligand as seed in preference vector E
  E[id2allgenes[l]] = 1
  return(igraph::page_rank(G, algo = c("prpack"), vids = igraph::V(G),directed = TRUE, damping = delta, personalized = E) %>% .$vector)
}
igraph::V(signaling_igraph)


E[id2allgenes[receptor[[1]]]] = 1
G <- signaling_igraph
delta <- damping_factor
ppr_cutoff <- optimized_parameters$ltf_cutoff
g <- igraph::page_rank(G, algo = c("prpack"), vids = igraph::V(G),directed = TRUE, damping = delta, personalized = E) %>% .$vector
ppr_matrix = matrix(unlist(g), ncol = length(E), byrow = TRUE)
ppr_matrix_TRUE = apply(ppr_matrix,1,function(x){x <= quantile(x,ppr_cutoff)}) %>% t()
ppr_matrix[ppr_matrix_TRUE] = 0

complete_matrix = apply(ppr_matrix, 2, mean)
ltf_matrix = matrix(unlist(complete_matrix), ncol = length(igraph::V(signaling_igraph)), byrow = TRUE)
rownames(ltf_matrix) = sapply(rec,function(x){paste0(x,collapse = "-")})




complete_matrix = lapply(receptor,PPR_wrapper,E,signaling_igraph,damping_factor,id2allgenes,ltf_cutoff)
