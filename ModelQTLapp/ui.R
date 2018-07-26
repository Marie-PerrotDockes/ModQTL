#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Old Faithful Geyser Data"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      downloadButton("data", "Download data :", class =".Rdata"),
      
      selectInput("type", "Type of Criterion:", choices = c(Lasso = "lasso",
                                                         Group =" group")),
      selectInput("uni_mult", "Type of Variate:", choices = c(Multivariate =" multi",
                                                              Univariate = "uni"
                                                            )),
      selectInput("grp","Type of group", choices = c(Both =" both", 
                                                     Trait =" trait",
                                                     Marker ="marker"
                                                     )),
      sliderInput("tresh",
                  "Index of the threshold in the stability path",
                  min = 1,
                  max = 100,
                  value = 30),
      selectInput("color", "Histogramme color", 
                  choices= c(Red ="red", Orange = "orange", Blue = "blue")),
      textInput("main" ,"Histogramme title",value = "Histo"),
      radioButtons("col",label = "Lapin:",choices = colnames(faithful))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot"),
      textOutput("numberofbins"),
      verbatimTextOutput("summary"),
      dataTableOutput("table")
    )
  )
))
