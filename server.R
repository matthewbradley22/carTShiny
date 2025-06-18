#library(EnhancedVolcano)
#library(apeglm)
library(DT)
library(stringr)
library(dplyr)
library(xlsx)
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  updateSelectizeInput(session,
                       'feature',
                       choices = rownames(T_cells_noCC[['RNA']]),
                       server = TRUE)
  
  selected <- eventReactive(input$Submit, {
    (
      subset(
        T_cells_noCC,
        CD_pred %in% input$CD_groups & CAR %in% input$CAR_groups &
          hypoxia %in% input$Hypoxia_groups &
          day %in% input$day_Groups &
          CAR_pred %in% as.numeric(input$CAR_presence)
      )
    )
  })
  
  compareGroupsDat <- eventReactive(input$runPseudo, {
    comparison <- input$comparisonChoice
    columnOfImportance <- case_when(startsWith(comparison, 'hypoxia')~'hypoxia',
                                    startsWith(comparison, 'day')~'day',
                                    startsWith(comparison, 'CD') ~ 'CD_pred',
                                    startsWith(comparison, 'CAR_pred')~'CAR_pred',
                                    startsWith(comparison, 'CAR_M')~'CAR',
                                    startsWith(comparison, 'CAR_u')~'CAR')
    
    referenceVar <- str_extract(comparison, "_(?!.*_).*")
    referenceVar <- sub('.', '', referenceVar)
    dat <- selected()
    dat[[]][[columnOfImportance]] <- relevel(factor(dat[[]][[columnOfImportance]]), ref = referenceVar)
    varList <- dat[[]] %>% select('hypoxia', 'day', 'CAR', 'CD_pred', 'CAR_pred', 'donor_id') %>% 
      sapply(unique) %>% lapply(length)
    designVars <- names(which(varList>1))
    dat_bulk <- getPseudoBulkObject(dat, c(designVars))
    dat_bulk <- DESeq(dat_bulk)
    results(dat_bulk, name=input$comparisonChoice)
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
    #selected <- subset(carPosT, CD_pred == input$CD_groups & CAR %in% input$CAR_groups &
    # hypoxia %in% input$Hypoxia_groups)
    FeaturePlot(selected(), features = input$feature, slot = 'data')
  })
  
  output$violinPlot <- renderPlot({
    VlnPlot(selected(),
            features = input$feature,
            group.by = input$vlnGroup)
  })
  
  output$dotPlot <- renderPlot({
    DotPlot(selected(),
            features = input$feature,
            group.by = input$vlnGroup,
            scale = FALSE)
  })
  
  
  
  output$VolcanoPlot <- renderPlot({
    EnhancedVolcano(
      volcData(),
      lab = rownames(volcData()),
      x = 'log2FoldChange',
      y = 'pvalue',
      title = 'filler',
      max.overlaps =  15,
      drawConnectors = TRUE
    )
  })

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
    dat[order(dat$padj),]  %>% as.data.frame() %>% head(n = 100)
  })

  output$curSubset <- renderText({
    dat <- selected()
    hypoxiaVals <- paste(as.character(unique(dat$hypoxia)), collapse = ' ')
    CARVals <-  paste(as.character(unique(dat$CAR)), collapse = ' ')
    dayVals <- paste(as.character(unique(dat$day)), collapse = ' ')
    CDVals <- paste(as.character(unique(dat$CD_pred)), collapse = ' ')
    CARPredVals <- paste(as.character(unique(dat$CAR_pred)), collapse = ' ')
    
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
       dat <- dat[order(dat$padj),]  %>% as.data.frame() %>% head(n = 100)
       write.xlsx(dat, con, sheetName = "Sheet1",
                  col.names = TRUE, row.names = TRUE, append = FALSE)
     }
    )
    
})