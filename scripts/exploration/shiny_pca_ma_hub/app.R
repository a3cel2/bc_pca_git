library(igraph)
library(devtools)
library(shinyapps)
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)
devtools::load_all('../../../packages/bcPcaAnalysis')
devtools::document('../../../packages/bcPcaAnalysis')

pca_universe = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.6.txt'
expression_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.14.txt'
protein_abundance_file = "/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/data/paxdb_abundance.tsv"
output_path <- 'www'#/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/scripts/exploration/shiny_pca_ma/www'

my_predictions <- pca_ma_prediction(pca_universe
                                    ,protein_abundance_file,
                                    expression_file,'ethanol',
                                    expression_condition_regexp='Ethanol.4h')

ui <- fluidPage(
  sidebarPanel( 
    sliderInput('node_size', 'Node Size', 1,50,30,step = 1),
    sliderInput('edge_size', 'Edge Size', 1,50,30,step = 1),
    textInput('hub', 'Hub', 'HXT1')
  ),
  mainPanel(
    htmlOutput('file1'))
)





server <- function(input, output) {
  
  generate_file <- reactive({
    reactive_element <- input$xsize;
    outfile <- paste(c(output_path,'plot1.pdf'),collapse='/')
    print(outfile)
    #tempfile(fileext = ".png")
    CairoPDF(outfile, width=12, height=6)
    hub_comparison_graph(
      my_predictions,
      hub_name=input$hub,
      node_size=input$node_size,
      edge_width=input$edge_size,
      output_path = '',filename='',draw = T
    )
      
    dev.off()
    return('plot1.pdf')
  }
  )
  
  output$file1 <- renderUI({tags$iframe(style="height:900px; width:100%; scrolling=yes", 
                          src=generate_file())})
  
  
  
}

shinyApp(ui = ui, server = server)
