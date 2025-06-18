library(shiny)
library(Seurat)

#Load data from hpc2n mystore/cartdata/data/Car_pos_TCells.rds
#carPosT <- LoadSeuratRds('Car_pos_TCells.rds')

T_cells_noCC <- LoadSeuratRds('T_cells_noCC.rds')
#Load deseq data for volcanoes
#cd4CARPseudoBulk <- readRDS("./cd4_bulk_noCC_Car+.rds")

#Format comparisons for volcanoes
# comparisons <- resultsNames(cd4CARPseudoBulk)[grepl('hypoxia|CAR',resultsNames(cd4CARPseudoBulk))] 
# comparisons <- unlist(lapply(comparisons, FUN = function(x){
#   gsub('_', ' ', x)
# }))

#Function for Comapre Two Groups tab to get pseudobulk object for DEG analysis
getPseudoBulkObject <- function(dat, designVars, intercept = "include", return.Seurat = FALSE){
  bulk <- AggregateExpression(
    dat,
    return.seurat = T,
    assays = "RNA",
    group.by = designVars)
  
  # Add in number of cells by sample and celltype
  n_cells <- dat[[]] %>% dplyr::select(all_of(designVars)) %>% 
    dplyr::group_by(across(all_of(designVars))) %>% 
    dplyr::count() 
  
  bulk[[]] <- n_cells
  
  if(return.Seurat){
    return(bulk)
  }
  # Get count matrix
  cluster_counts <- FetchData(bulk, layer="counts", vars=rownames(bulk))
  
  # Create DESeq2 object
  if(intercept == "include"){
    designForm <- reformulate(termlabels = designVars)
    dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                  colData = bulk[[]],
                                  design = designForm)
  }else{
    designForm <- reformulate(termlabels = c(0,designVars))
    dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                  colData = bulk[[]],
                                  design = designForm)
  }
  
  dds
}
