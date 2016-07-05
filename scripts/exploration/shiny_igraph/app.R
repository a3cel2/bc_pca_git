
ui <- fluidPage(
  
  sidebarPanel(
    actionButton('edit', 'Grab layout from editor')
  ),
  mainPanel(
    plotOutput('file1')
    )
  
  
)

server <- function(input, output) {
  my_graph <- erdos.renyi.game(50,0.1)
  tkid <- tkplot(my_graph)

  
  my_layout <- reactive({
    test <- input$edit
    l <- tkplot.getcoords(tkid)
    return(l)
  }
  )
  
  
  output$file1 <- renderPlot({
    plot(my_graph,layout=my_layout())
  })
}

shinyApp(ui = ui, server = server)
