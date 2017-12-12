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
shinyUI(
  fluidPage(
    
    # Application title
    titlePanel("Large Simulation Results"),
    
    selectInput("physicalSystem", 
                "Physical System",  
                physicalSystem.choices,
                physicalSystem.choices[1],
                width="400%"),
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        selectInput("nobs", 
                    "Number of Observations",  
                    nobs.candidates,
                    51),
        hr(),
        selectInput("maternDf", 
                    "Degree of Freedom for Matern kernel",  
                    c(2.01, 2.5),
                    2.01),
        hr(),
        selectInput("noise", 
                    "simulated noise level",  
                    noise.candidates,
                    0.1),
        hr(),
        selectInput("phaseType", 
                    "what mu to plug-in",
                    c("phase1", "phase2", "trueMu", "priorTempered", "priorTemperedPhase2"),
                    "priorTemperedPhase2"),
        hr(),
        checkboxGroupInput("variablePrint", "Variables to print:",
                           c("avg 2.5% quantile" = "2.5%",
                             "avg median" = "50%",
                             "avg 97.5% quantile" = "97.5%",
                             "avg mean" = "mean",
                             "avg coverage" = "coverage"),
                           "coverage"),
        hr(),
        textInput("pdfBaseDir", "directory for pdf results (must be absolute)", 
                  "/Users/shihaoyang/GoogleDrive/results")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        tabPanel("Summary Statistics", 
                 h4("Summary Statistics"),
                 DT::dataTableOutput('performTable')),
        hr(),
        tabPanel("repetition size", 
                 h4('repetition size'),
                 DT::dataTableOutput('repSizeTable'))
      )
      
    )
  )
)
