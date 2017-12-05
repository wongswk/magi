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
    baseInfoTab <- data.frame(outTab$repetitionSize)
    
    if(input$maternDf == 2.01){
      kernelType = "generalMatern"
    }else if(input$maternDf == 2.5){
      kernelType = "matern"
    }else{
      stop("invalid kernel df")
    }
    pdfBaseLink <- input$pdfBaseDir
    subDir <- paste0("withmeanBand-", kernelType, "-nobs", input$nobs, "-noise", input$noise, "-ndis")
    subDir <- list.dirs(pdfBaseLink)[grep(subDir, list.dirs(pdfBaseLink))]
    ndis <- as.numeric(gsub(".*-ndis([0-9]+)", "\\1", subDir))
    subDir <- subDir[order(ndis)]
    ndis <- sort(ndis)
    urls <- paste0("file://", subDir)
    refs <- paste0("<a href='",  urls, "' target='_blank'>pdf visualization ndis-",ndis,"</a>")
    baseInfoTab$visualizations <- refs
    baseInfoTab
    
  })
  output$repSizeTable <- DT::renderDataTable(
     repSizeTable(), escape = FALSE, options = list(destroy = TRUE,paging=FALSE,searching=FALSE, bInfo=FALSE))
})
