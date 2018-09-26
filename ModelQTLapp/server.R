#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
   X <- reactive({ inFile <- input$X
   
   if (is.null(inFile))
     return(NULL)
   X <-  read.csv(inFile$datapath, header = input$header, sep = input$sep)
   })
   
   Y <- reactive({
     inFile <- input$Y
     
     if (is.null(inFile))
       return(NULL)
     
     Y <- read.csv(inFile$datapath, header = input$header, sep = input$sep)
     Y
   })
   
   
  output$Regressors <- renderDataTable({
        X()
  })
 
  
  output$Responses <- renderDataTable({
    Y()
  })

  model <- reactive({
    input$go
    isolate({
      model <- QTLmod_group_univ$new(X(), Y())
      model$estime()
      model
    })
  })
  
  output$Model1 <- renderPlot({
       if(!input$plot)  return(NULL)
         
         model()$plot_coef(input$tresh)
  })
  
  
  output$CV1 <- renderPlot({
    if(!input$CV)  return(NULL)
    
    model()$plot_cv(s =input$t_cv)
  })
  

  
})
