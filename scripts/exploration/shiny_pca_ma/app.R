library(devtools)
library(shinyapps)
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)
devtools::load_all('../../../packages/bcPcaAnalysis')
devtools::document('../../../packages/bcPcaAnalysis')

pca_universe = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.6.txt'
expression_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.14.txt'
protein_abundance_file = "/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/data/paxdb_abundance.tsv"
output_path <- '/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/scripts/exploration/shiny_pca_ma/www'

my_predictions <- pca_ma_prediction(pca_universe
                                    ,protein_abundance_file,
                                    expression_file,'ethanol',
                                    expression_condition_regexp='Ethanol.4h')

ui <- fluidPage(
  sidebarPanel( 
    sliderInput('label_size', 'Label Size', 0,5,1.23,step = 1/100),
    sliderInput('xsize', 'Point Size', 0,10,0.8,step = 0.1),
    sliderInput('point_red', 'Point colour red', 0,1,0.15,step = 1/255),
    sliderInput('point_green', 'Point colour green', 0,1,0.22,step = 1/255),
    sliderInput('point_blue', 'Point colour blue', 0,1,0.33,step = 1/255),
    sliderInput('point_transparent', 'Point colour opacity', 0,1,0.37,step = 1/255),
    sliderInput('outline_transparent', 'Outline opacity', 0,1,0,step = 1/255)
    
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
    CairoPDF(outfile, width=5, height=5)
    pca_ma_prediction_plot(
      my_predictions,output_path = '',draw = T,
      point_size = input$xsize,
      point_colours = rgb(
        input$point_red,
        input$point_green,
        input$point_blue,
        input$point_transparent
      ),
      outline_colours = rgb(
        0,0,0,input$outline_transparent
      ),
      label_size = input$label_size
    )
    dev.off()
    return('plot1.pdf')
  }
  )
  
  output$file1 <- renderUI({tags$iframe(style="height:900px; width:100%; scrolling=yes", 
                          src=generate_file())})
  
  
  
}

shinyApp(ui = ui, server = server)
