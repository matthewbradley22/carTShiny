#Load libraries
library(DT)
library(stringr)
library(dplyr)
library(xlsx)

shinyServer(function(input, output, session) {
  updateSelectizeInput(session,
                       'feature',
                       choices = rownames(T_cells_noCC[['RNA']]),
                       server = TRUE)
  
  #Subset dataset when user presses submit on sidebar
  selected <- eventReactive(input$Submit, {
    (
      subset(
        T_cells_noCC,
        CD_pred %in% input$CD_groups & CAR %in% input$CAR_groups &
          hypoxia %in% input$Hypoxia_groups &
          day %in% input$day_Groups &
          carExpression %in% input$CAR_presence
      )
    )
  })
  
  #Create pseudobulk comparison when user clicks run
  compareGroupsDat <- eventReactive(input$runPseudo, {
    comparison <- input$comparisonChoice
    columnOfImportance <- case_when(startsWith(comparison, 'hypoxia')~'hypoxia',
                                    startsWith(comparison, 'day')~'day',
                                    startsWith(comparison, 'CD') ~ 'CD_pred',
                                    startsWith(comparison, 'carExpression')~'carExpression',
                                    startsWith(comparison, 'CAR_M')~'CAR',
                                    startsWith(comparison, 'CAR_u')~'CAR')
    
    referenceVar <- str_extract(comparison, "_(?!.*_).*")
    referenceVar <- sub('.', '', referenceVar)
    dat <- selected()
    dat[[]][[columnOfImportance]] <- relevel(factor(dat[[]][[columnOfImportance]]), ref = referenceVar)
    varList <- dat[[]] %>% select('hypoxia', 'day', 'CAR', 'CD_pred', 'carExpression', 'donor_id') %>% 
      sapply(unique) %>% lapply(length)
    designVars <- names(which(varList>1))
    dat_bulk <- getPseudoBulkObject(dat, c(designVars))
    dat_bulk <- DESeq(dat_bulk)
    results(dat_bulk, name=input$comparisonChoice)
  })
  
  #Create data for UMAP
  UMAPDat <- eventReactive(input$SubmitUMAP, {
    dat <- selected()
    dat <- NormalizeData(dat)
    dat <-FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
    dat <- RunPCA(dat, features = VariableFeatures(object = dat))
    dat <- FindNeighbors(dat, dims = 1:input$numNeighbors)
    dat <- FindClusters(dat, resolution = 0.5)
    dat <- RunUMAP(dat, dims = 1:input$numNeighbors)
    dat
  })
  
  # selectedVolc <- eventReactive(input$Submit, { #This won't work with the resultsNames inputs for volcData below
  #   (
  #     cd4CARPseudoBulk[cd4CARPseudoBulk$CD_pred %in% input$CD_groups & 
  #                        cd4CARPseudoBulk$CAR == input$CAR_groups,
  #                      cd4CARPseudoBulk$hypoxia %in% input$Hypoxia_groups & 
  #                        cd4CARPseudoBulk$day %in% input$day_Groups]
  #   )
  # })
  # 
  # volcData <- eventReactive( input$RunVolcano, {
  #   coefChoice <- gsub(' ', '_', input$volcanoComparison)
  #   lfcShrink(cd4CARPseudoBulk, coef = coefChoice, type = "apeglm")
  # })
  
  
  output$dimPlot1 <- renderPlot({ 
    dat <- UMAPDat()
    DimPlot(dat)
  })
  
  observeEvent(input$SubmitUMAP,{
    dat = UMAPDat()
    updateSelectInput(session, "findMarkersIdent", choices = sort(unique((dat$seurat_clusters))))
  })
  

  output$SeuratClusterDEGs <- renderDataTable({
      dat <- UMAPDat()
      markers <- FindMarkers(dat, group.by = 'seurat_clusters', ident.1 = as.numeric(input$findMarkersIdent))
      markers[markers$p_val_adj < 0.01,]
  })
  
  output$violinPlot <- renderPlot({
    VlnPlot(selected(),
            features = input$feature,
            group.by = input$vlnGroup)
  })
  
  
  # output$VolcanoPlot <- renderPlot({
  #   EnhancedVolcano(
  #     volcData(),
  #     lab = rownames(volcData()),
  #     x = 'log2FoldChange',
  #     y = 'pvalue',
  #     title = 'filler',
  #     max.overlaps =  15,
  #     drawConnectors = TRUE
  #   )
  # })

  output$GeneList <- renderDataTable({
    comparison <- input$variableComparison
    compDirectory <- paste0(gsub(' ','', comparison), 'Comps/')
    geneFiles = list.files(compDirectory)
    geneFiles = paste0(compDirectory, geneFiles)
    geneTables = lapply(geneFiles, read.csv)
    
    geneListNames = list.files(compDirectory)
    names(geneTables) <- geneListNames
    
    geneTables <- lapply(geneTables, FUN = function(x){
      head(x, n = 20)
      })
    
    for(i in 1:length(geneTables)){
      geneTables[[i]]$comp = gsub('.csv','',geneListNames[[i]])
    }
    
    geneTables <- do.call(rbind, geneTables)
    tableDat = as.data.frame(sort(table(geneTables$X), decreasing = TRUE))
    colnames(tableDat) = c('Gene', 'Count')
    tableDat
  })
  
  output$comparisonsList <- renderDataTable({
    comparison <- input$variableComparison
    compDirectory <- paste0(gsub(' ','', comparison), 'Comps/')
    geneFiles = list.files(compDirectory)
    geneFiles = paste0(compDirectory, geneFiles)
    geneTables = lapply(geneFiles, read.csv)
    
    geneListNames = list.files(compDirectory)
    names(geneTables) <- geneListNames
    
    geneTables <- lapply(geneTables, FUN = function(x){
      head(x, n = 20)
    })
    
    present <- lapply(geneTables, FUN = function(x){
      input$geneForListSearch %in% x$X
    })
    
    unname(unlist(present))
    presenceDF <- data.frame(names = names(present), presence =  unname(unlist(present)))
    presenceDF <- presenceDF[presenceDF$presence == TRUE,]
    geneListDF <- as.data.frame(str_remove(presenceDF$names, '.csv'))
    colnames(geneListDF) <- 'Gene Tables'
    geneListDF
  })
  
  output$DifferentialGeneList <- renderDataTable({
    dat <- compareGroupsDat()
    if(input$removeTrivial){
      dat <- dat[!grepl('^MT', rownames(dat)),] #Remove mitochondrial genes
    }
    dat[order(dat$padj),]  %>% as.data.frame() %>% head(n = 1000)
  })

  output$curSubset <- renderText({
    dat <- selected()
    hypoxiaVals <- paste(as.character(unique(dat$hypoxia)), collapse = ' ')
    CARVals <-  paste(as.character(unique(dat$CAR)), collapse = ' ')
    dayVals <- paste(as.character(unique(dat$day)), collapse = ' ')
    CDVals <- paste(as.character(unique(dat$CD_pred)), collapse = ' ')
    CARPredVals <- paste(as.character(unique(dat$carExpression)), collapse = ' ')
    
    paste('Current subset:<br> <B>hypoxia values:</B> ', hypoxiaVals, 
          ' <br><B>CAR values:</B> ', CARVals,
          ' <br><B>Day values:</B> ', dayVals,
          ' <br><B>CD values:</B> ', CDVals,
          ' <br><B>CAR presence values:</B> ', CARPredVals)
  })
  
  output$downloadExcel <- downloadHandler(
       filename = function() {
        paste(input$comparisonChoice, Sys.Date(), '.xlsx', sep='')
      },
     content = function(con) {
       dat <- compareGroupsDat()
       dat <- dat[order(dat$padj),]  %>% as.data.frame() %>% head(n = 1000)
       write.xlsx(dat, con, sheetName = "Sheet1",
                  col.names = TRUE, row.names = TRUE, append = FALSE)
     }
    )
    
})