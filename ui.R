
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
options(shiny.maxRequestSize = 50 * 1024 ^ 2)
library(cttools)
library(shiny)
library(shinyBS)
library(shinyjs)

shinyUI(fluidPage(
  shinyjs::useShinyjs(),  
  
  # Application title
  titlePanel("ICH Segmentation"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      fileInput('img_fname', 'Choose NIfTI image to upload',
                accept = c(
                  '.nii',
                  '.nii.gz'
                )
      ),
      bsButton("ss", "Skull Strip Image", disabled = TRUE),
      br(),
      br(),
      bsButton("reg", "Register Skull-Stripped Image", disabled = TRUE),
      br(),
      br(),
      bsButton("make_pred", "Create Predictors", disabled = TRUE),
      br(),
      br(),
      bsButton("predict", "Create Prediction Images", disabled = TRUE),
      br(),
      br(),
      downloadButton('download', 'Download')
      # actionButton("process", "Process Skull-stripped Image")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("origPlot"),
      plotOutput("ssPlot"),
      plotOutput("regPlot"),
      plotOutput("candPlot"),
      plotOutput("predPlot")
    )
  )
))
