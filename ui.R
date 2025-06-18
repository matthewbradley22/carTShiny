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
    checkboxGroupInput('CAR_presence', 'CAR presence', choices = c('0', '1'),
                       selected = c('0', '1')),
    radioButtons('vlnGroup', 'Group Violin-Plot', c('CAR', 'CD_pred', 'hypoxia', 'day')),
    actionButton("Submit", "Submit", class = "btn-lg btn-success")
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel('Summary Plots',
                 plotOutput("violinPlot", height=400, width = 600),
                 plotOutput("dotPlot", height=400, width = 400)
               ),
      tabPanel('UMAP', plotOutput("dimPlot1", height=400, width = 400)),
      tabPanel('Gene Lists',
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
                                                                               "CAR_pred_1_vs_0")),
               actionButton("runPseudo","Run"),
               p("Make sure to click 'submit' on the sidebar before clicking 'Run', whether you subset the data or not"),
               withSpinner(DT::dataTableOutput("DifferentialGeneList"), color="#0dc5c1"))
    ),
 
  )
))
