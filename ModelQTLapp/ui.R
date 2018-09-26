#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(ModQTL)
library(glmnet)
require(gganimate)
require(RColorBrewer)
# Define UI for application that draws a histogram
shinyUI(
  navbarPage("Model QTL",
             tabPanel("Upload Data",
                     column(width=3 ,
                            wellPanel(
                      fluidPage(
                        fileInput("X", "Choose CSV File for regressors"
                        ),
                        fileInput("Y", "Choose CSV File for responses"
                        ),
                        tags$hr(),
                        checkboxInput("header", "Header", TRUE),
                        selectInput("sep", "Separator:",
                                    c("," = ",",
                                      ";" = ";",
                                      "tabulation" = "\t"))
                        
                        
                      ))),
                      column(width=9,
                      dataTableOutput("Regressors"),
                      dataTableOutput("Responses")
                      )
             ),
             tabPanel("Univariate Lasso",
                      column(width=3 ,
                             wellPanel(
                               fluidPage(
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
                                             "Index of the threshold in the regression path",
                                             min = 1,
                                             max = 100,
                                             value = 30),
                                 selectInput("t_cv","Choice for Cross-validation", choices = c(LambdaMin =" lambda.min", 
                                                                                Lambda1se =" Lambda.1se"
                                 )),
                                 checkboxInput("plot","Plot", F),
                                 checkboxInput("CV","Perform Cross Validation", F),
                                 checkboxInput("anim","Animate", F),
                                 actionButton("go", "Run Model")
                                 
                               ))),
                      column(width=9,
                             plotOutput("Model1"),
                             plotOutput("CV1")
                      )
             )
  ))