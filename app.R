gc()

library(shiny)
library(shinythemes)
library(rio)
library(here)
library(plotly)
library(dplyr)
library(biomaRt)
library(DT)

##for biomaRt API usage##
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# snpmart <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
##for biomaRt API usage##

variants_data<- import(here("phewas-data-karla.xlsx"))

names(variants_data)[names(variants_data) == "chromosome"] <- "chr"
names(variants_data)[names(variants_data) == "cases"] <- "Cases"
names(variants_data)[names(variants_data) == "phewas phenotype"] <- "Phenotype"
names(variants_data)[names(variants_data) == "p-value"] <- "P_value"
names(variants_data)[names(variants_data) == "gene_name"] <- "Gene_name"
names(variants_data)[names(variants_data) == "gwas-associations"] <- "Gwas_Association"

ui <- 
  navbarPage("App", collapsible = TRUE, inverse = TRUE, theme = shinytheme("united"),
             tabPanel("All-data",
                      sidebarLayout(
                        sidebarPanel(
                          fileInput(
                            'variants_file',
                            'File:',
                            multiple=FALSE,
                            accept=".csv",
                            buttonLabel = "Browse...",
                            placeholder = "No file selected"
                          )),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Manhattan Plot",
                                     br(),br(),p("This graph represents a Manhattan Plot of some of the phenotypes associated to 
                                     variants (rsID) across the whole human genome reported in PheWAS Catalog.  For detailed phenotype information please look at the tables"),
                                     plotlyOutput(outputId = "all_data"),br(),br(),br(),p("The color of the marker is different for each chromosome."),
                                     p("The size of the marker is determined by the number of cases reported of a phenotype in that variant."),
                                     uiOutput("link")),
                            tabPanel("All-data table",dataTableOutput(outputId = "all_data_table")),
                            type="tab"
                            )
                          ))),
             tabPanel("PheWAS-browser",
                                titlePanel("PheWAS-browser"),
                                # Sidebar with a slider input for number of bins 
                                sidebarLayout(
                                  sidebarPanel(
                                    selectizeInput(
                                      'gene_name',
                                      'Gene name:',
                                      choices=NULL,
                                      selected=NULL,
                                      multiple=FALSE
                                    ),
                                    uiOutput("secondSelection")
                                  ),
                                  mainPanel(
                                    plotlyOutput(outputId = "phenotypes" )
                                  ))),
             tabPanel("Gene-browser",
                      fluidPage(
                        tabsetPanel(
                          tabPanel("Gene data",
                                   sidebarLayout(
                                     sidebarPanel(
                                       selectizeInput(
                                         'gene_name_data',
                                         'Gene name:',
                                         choices=NULL,
                                         selected=NULL,
                                         multiple=FALSE
                                       ),
                                     ),
                                     mainPanel(
                                       dataTableOutput(outputId = "gene_data" )
                                     )
                                   )),
                          tabPanel("Transcript data",
                                   sidebarLayout(
                                     sidebarPanel(
                                       uiOutput("transcriptChoices")
                                     ),
                                     mainPanel(
                                       dataTableOutput(outputId = "transcript_seq" )
                                     )
                                  ))
                          # tabPanel("Variant data",
                          #          sidebarLayout(
                          #            sidebarPanel(
                          #              uiOutput("variantChoices")
                          #            ),
                          #            mainPanel(
                          #              dataTableOutput(outputId = "variant_datatable")
                          #            )
                          #         ))
                        )))
  )

server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'gene_name', choices=c("",variants_data$Gene_name), server=TRUE)
  updateSelectizeInput(session, 'gene_name_data', choices=c("",variants_data$Gene_name), server=TRUE)
  
  url <- a("PheWAS Catalog", href="https://urldefense.com/v3/__https://phewascatalog.org__;!!Nmw4Hv0!0RPhOdwdvpSVrpMMxprAMi83uk3zmFk9GGVzmvbBWioR29Oxfa2xKvQN38Ero-yrgilgnf-wymOJfw0yYlLru-c8PcIUIGnrhyid$ ")
  
  output$link <- renderUI({
    tagList("Link to ",url,".")
  })
  
  output$all_data <- renderPlotly ({
    plot_ly() %>%
      add_markers(data=variants_data, x = ~chr, y = ~logPvalue,
                  customdata = ~Gene_name, color = ~chr,
                  marker = list(opacity=0.25,size = ~Cases/200), text = ~snp,
                  hoverinfo=FALSE,
                  hovertemplate = paste("<b>Gene: </b>%{customdata}",
                                        "<br><b>snp: </b>%{text}",
                                        "<extra></extra>")) %>%
      hide_colorbar() %>%
      layout(xaxis = list(title='Chromosome',range=c(0,23),autotick=F,dtick=1),
             yaxis = list(title='- Log10(P-value)'))
  })
  
  output$all_data_table <- renderDataTable({
    DT::datatable(variants_data[c("chr","pos","snp","Gene_name","Phenotype",
                                  "Cases","P_value", "Gwas_Association")],
                  filter='top', rownames=FALSE, options = list(pageLength = 10),
                  selection=list(mode='multiple',target='row'))
  })
  
  gene.choices = reactive({
    if(is.null(variants_data)){return()}
    else{
      variants_data %>%
        filter(Gene_name==input$gene_name)}
  })
  
  dataFiltGenes <- reactive({
    rowsIN <- which(variants_data$Gene_name == input$gene_name)
    variants_data[rowsIN,]
  })
  
  dataFiltered <- reactive({
    rowsIN <- which(variants_data$snp == input$snp)
    variants_data[rowsIN,]
  })
  
  ## Reactive input snp (depends on the first input selection)
  output$secondSelection <- renderUI({selectizeInput(
    'snp',
    'SNP:',
    choices=c("",gene.choices()$snp),
    selected=NULL,
    multiple=FALSE
  )})
  
  output$phenotypes <-renderPlotly({
    if (nrow(dataFiltGenes())==0) {
      return(NULL)
    }
    else if ((nrow(dataFiltGenes())) != 0 && (nrow(dataFiltered()) == 0))  {
      plot_ly(dataFiltGenes()) %>%
        add_trace(x = ~Cases, y = ~Phenotype,
                  type='bar', orientation='h',
                  marker = ~list(color = c('rgb(158,202,225)','rgb(165,139,196)','rgb(139,196,79)','rgb(247,236,13)','rgb(242,205,247)',
                                           'rgb(247,104,65)','rgb(103,213,183)','rgb(158,202,225)','rgb(165,139,196)','rgb(139,196,79)',
                                           'rgb(247,236,13)','rgb(242,205,247)','rgb(247,104,65)','rgb(103,213,183)'), 
                                 line = list(color = 'rgb(0,0,0)', width = 1.5)))%>%
        layout(xaxis = list(title='Number of cases'),
               yaxis = list(title='Phenotype'))
    }
    else if ((nrow(dataFiltGenes())) != 0 && (nrow(dataFiltered()) != 0))   {
      plot_ly(dataFiltered()) %>%
        add_trace(x = ~Cases, y = ~Phenotype,
                  type='bar', orientation='h',
                  marker = ~list(color = c('rgb(158,202,225)','rgb(165,139,196)','rgb(139,196,79)','rgb(247,236,13)','rgb(242,205,247)',
                                           'rgb(247,104,65)','rgb(103,213,183)','rgb(158,202,225)','rgb(165,139,196)','rgb(139,196,79)',
                                           'rgb(247,236,13)','rgb(242,205,247)','rgb(247,104,65)','rgb(103,213,183)'), 
                                line = list(color = 'rgb(0,0,0)', width = 1.5)))%>%
        
        layout(xaxis = list(title='Number of cases'),
               yaxis = list(title='Phenotype'))
    }
    
    
  })
  
  ###################################################################################################################################################
  
  ##### USING biomaRt API #####
  
  gene_data <- reactive ({
    getBM(attributes=c('chromosome_name','start_position','end_position',
                       'ensembl_gene_id','ensembl_transcript_id'), 
          filters ='external_gene_name',
          values =input$gene_name_data, mart = ensembl)
  })
  
  output$gene_data <- DT::renderDataTable({
    DT::datatable(gene_data(), rownames=FALSE,
                  options = list(lengthChange=FALSE, info=FALSE, scrollY=F, dom= 't',
                  fillContainer=TRUE))
  })
  
  ## Reactive input transcript_data
  
  output$transcriptChoices <- renderUI({selectizeInput(
    'trans',
    'Transcript:',
    choices=c("",gene_data()$ensembl_transcript_id),
    selected="ENST00000425614",
    multiple=FALSE
  )})
  
  current_selection <- reactiveVal(NULL)
  
  observeEvent(input$trans, {
    current_selection(input$trans)
  })
  
  seq_aa <- reactive ({
    getSequence(id = input$trans,
    type = "ensembl_transcript_id",
    seqType = "peptide",mart = ensembl)
  })
  
  output$transcript_seq <- DT::renderDataTable({
    DT::datatable(seq_aa(),rownames=FALSE,
                  options=list(
                    scrollX = TRUE,
                    autoWidth=TRUE,
                    columnDefs=list(list(width='10px', targets="_all"))
                  ))
  })
  
  gene_data.choices = reactive({
     variants_data %>%
      filter(Gene_name==input$gene_name_data)
    })
  
  # ## Reactive input variant_data
  # output$variantChoices <- renderUI({selectInput(
  #   'snp_data',
  #   'SNP:',
  #   choices=c("",gene_data.choices()$snp),
  #   selected=NULL,
  #   multiple=FALSE
  #   )})
  # 
  # current_selection2 <- reactiveVal(NULL)
  # 
  # observeEvent(input$snp_data, {
  #   current_selection2(input$snp_data)
  #   })
  # 
  # variant_data <- reactive ({
  #   getBM(attributes=c('refsnp_id','chr_name','chrom_start','chrom_end',
  #                       'chrom_strand','allele','allele_1','minor_allele'),
  #          filters ='snp_filter', values =input$snp_data, mart = snpmart)
  #  })
  # 
  #  output$variant_datatable <- DT::renderDataTable({
  #    DT::datatable(variant_data(), rownames=FALSE,
  #                  options=list(
  #                    scrollX = TRUE,
  #                    autoWidth=TRUE,
  #                    columnDefs=list(list(width='10px', targets="_all"))
  #                  ))
  #  })
  
}

shinyApp(ui = ui, server = server)
