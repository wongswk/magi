#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  # total performance
  performTable <- reactive({
    
    outTab <- organizeOutput(maternDf = input$maternDf, 
                             noise = input$noise, 
                             nobs = input$nobs)
    
    outTabCoverage <- outTab$coverage[,,input$phaseType,]
    outTabCoverage <- round(outTabCoverage, 4)
    
    outTabCoveragePrint <- lapply(1:dim(outTabCoverage)[3], function(itab){
      outEachTab <- format(outTabCoverage[,,itab][input$variablePrint,,drop=FALSE], nsmall=4)
      rownames(outEachTab) <- input$variablePrint
      rownames(outEachTab) <- paste0(dimnames(outTabCoverage)[[3]][itab],
                                     " || ", rownames(outEachTab))
      outEachTab
    })
    outTabCoveragePrint <- do.call(rbind, outTabCoveragePrint)
    data.frame(outTabCoveragePrint)
  })
  output$performTable <- DT::renderDataTable({
    DT::datatable(performTable(), 
                  options = list(destroy = TRUE,paging=FALSE,searching=FALSE, bInfo=FALSE))
    })
  
  repSizeTable <- reactive({
    outTab <- organizeOutput(maternDf = input$maternDf, 
                             noise = input$noise, 
                             nobs = input$nobs)
    outTab$repetitionSize[,input$phaseType]
    data.frame(outTab$repetitionSize)
  })
  output$repSizeTable <- DT::renderDataTable({
    DT::datatable(repSizeTable(), 
                  options = list(destroy = TRUE,paging=FALSE,searching=FALSE, bInfo=FALSE))
  })
})
