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
    sliderInput('label_size', 'Label Size', 0,5,1.23,step = 1/100),
    sliderInput('xsize', 'Point Size', 0,10,0.8,step = 0.1),
    
    sliderInput('enhanced_red', 'enhanced_red', 0,1,0.15,step = 1/255),
    sliderInput('enhanced_green', 'enhanced_green', 0,1,0.22,step = 1/255),
    sliderInput('enhanced_blue', 'enhanced_blue', 0,1,0.33,step = 1/255),
    
    sliderInput('depleted_red', 'depleted_red', 0,1,0.15,step = 1/255),
    sliderInput('depleted_green', 'depleted_green', 0,1,0.22,step = 1/255),
    sliderInput('depleted_blue', 'depleted_blue', 0,1,0.33,step = 1/255),
    
    sliderInput('polygon_transparent', 'Polygon transparency', 0,1,0.37,step = 1/255),
    sliderInput('legend_x', 'Legend X position', -4,4,-0.9,step = 1/255),
    sliderInput('legend_y', 'Legend Y position', 0,1,0.85,step = 1/255)
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
    CairoPDF(outfile, width=6, height=6)
    pca_ma_precision_plot(
      my_predictions,output_path = '',draw = T,
      bootstrap_iters=10,
      enhanced_colour=rgb(input$enhanced_red,input$enhanced_green,input$enhanced_blue),
      depleted_colour=rgb(input$depleted_red,input$depleted_green,input$depleted_blue),
      polygon_transparency=input$polygon_transparent,
      legend_x_position=input$legend_x,
      legend_y_position=input$legend_y,
      label_size=1.3)
    dev.off()
    return('plot1.pdf')
  }
  )
  
  output$file1 <- renderUI({tags$iframe(style="height:900px; width:100%; scrolling=yes", 
                          src=generate_file())})
  
  
  
}

shinyApp(ui = ui, server = server)
