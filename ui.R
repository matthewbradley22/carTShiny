library(shiny)
library(Seurat)
library(DESeq2)
library(shinycssloaders)

# Define UI for application that plots random distributions 
shinyUI(fluidPage(
  
  # Application title
  headerPanel("Car-T Shiny App"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    selectizeInput("feature", "gene to plot", choices = NULL),
    checkboxGroupInput('CD_groups', 'CD Group', choices = c('CD4', 'CD8'), selected = c('CD4','CD8')),
    checkboxGroupInput('CAR_groups', 'CAR Groups', choices = c('MBBz', 'M28z', 'M1XX', 'untransduced'),
                       selected = c('MBBz', 'M28z', 'M1XX', 'untransduced')),
    checkboxGroupInput('Hypoxia_groups', 'Hypoxia Groups', choices = c('HH', 'NH', 'NN'),
                       selected = c('HH', 'NH', 'NN', 'untransduced')),
    checkboxGroupInput('day_Groups', 'Day', choices = c('D7', 'D13'),
                       selected = c('D7', 'D13')),
    checkboxGroupInput('CAR_presence', 'CAR presence', choices = c('carNeg', 'carPos'),
                       selected = c('carNeg', 'carPos')),
    radioButtons('vlnGroup', 'Group Violin-Plot', c('CAR', 'CD_pred', 'hypoxia', 'day')),
    actionButton("Submit", "Submit", class = "btn-lg btn-success")
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel('Summary Plots',
                 plotOutput("violinPlot", height=400, width = 600)
               ),
      tabPanel('Gene Lists',
               h2('This is now outdated, need to update with new dataset'),
               p(em("For each gene, these tables show the number of sub comparisons that have this gene
                 as a top DEG (in the top 20). For example, in the HH vs NN comparison, SBF2 shows up in 13 of the
                 30 subcomparisons as a top 20 DEG. You can view which subcomparisons by typing SBF2 in the box below
                 labeled 'gene' and see that it is a top DEG in D13_CD4_M1XX_0_HHNN, D13_CD4_M1XX_1_HHNN ...")),
               selectInput('variableComparison', 'comparison', choices = c('HH vs NN', 'Car Neg Vs Pos')),
               DT::dataTableOutput("GeneList"),
               tableOutput("OverlappingGenes"),
               textInput('geneForListSearch', 'Gene'),
               DT::dataTableOutput("comparisonsList")),
      tabPanel('Compare Two Groups',
               p("Make sure to click 'Submit' on the sidebar before clicking 'Run', whether you subset the data or not."),
               selectInput('comparisonChoice', 'Choose comparison', choices = c("hypoxia_NH_vs_HH",
                                                                               "hypoxia_NN_vs_HH",
                                                                               "hypoxia_NH_vs_NN",
                                                                               "CAR_M28z_vs_M1XX",
                                                                               "CAR_MBBz_vs_M1XX",
                                                                               "CAR_untransduced_vs_M1XX",
                                                                               "CAR_untransduced_vs_M28z",
                                                                               "CAR_untransduced_vs_MBBz",
                                                                               "CAR_M28z_vs_MBBz",
                                                                               "day_D13_vs_D7",
                                                                               "CD_pred_CD8_vs_CD4",
                                                                               "carExpression_carPos_vs_carNeg")),
               checkboxInput('removeTrivial', 'Remove Trivial Genes', value = FALSE),
               actionButton("runPseudo","Run"),
               withSpinner(DT::dataTableOutput("DifferentialGeneList"), color="#0dc5c1"),
               htmlOutput('curSubset'),
               downloadButton('downloadExcel', label = "Download", class = NULL)),
      tabPanel('UMAP', 
               sliderInput('numNeighbors', 'Choose a number of neighbors', min = 1, max = 40, value = 20),
               actionButton("SubmitUMAP", "Submit", class = "btn-lg btn-success"),
               withSpinner(plotOutput("dimPlot1", height=400, width = 400), color="#0dc5c1"),
               selectInput('findMarkersIdent', 'Select Group', choices = NULL),
               DT::dataTableOutput("SeuratClusterDEGs"),
               actionButton("GenerateDEGs", "Generate DEGs", class = "btn-lg btn-success")),
    ),
 
  )
))
