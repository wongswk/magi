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
shinyServer(function(input, output, session) {
  
  observe({
    physicalSystem <- input$physicalSystem
    
    if(physicalSystem == "FitzHugh-Nagumo (FN) Model"){
      load("data/FN-largeExperimentSummary.rda", envir = .GlobalEnv)
      assign("temperPrior.candidates", c(TRUE, FALSE), envir = .GlobalEnv)
      assign("phaseType.candidates", c("phase1", "phase2"), envir = .GlobalEnv)
    }else {
      stop("only FN system")
    }
    
    updateSelectInput(session, "nobs", choices=nobs.candidates, selected=nobs.candidates[4])
    updateSelectInput(session, "temperPrior", choices=temperPrior.candidates, selected=temperPrior.candidates[1])
    updateSelectInput(session, "noise", choices=noise.candidates, selected=noise.candidates[2])
    updateSelectInput(session, "phaseType", choices=phaseType.candidates, selected=phaseType.candidates[1])
  })
  # total performance
  performTable <- reactive({
    outTab <- organizeOutput(temperPrior = input$temperPrior, 
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
    print("data.frame(outTabCoveragePrint)")
    data.frame(outTabCoveragePrint)
  })
  output$performTable <- DT::renderDataTable({
    DT::datatable(performTable(), 
                  options = list(destroy = TRUE,paging=FALSE,searching=FALSE, bInfo=FALSE))
    })
  
  
  repSizeTable <- reactive({
    outTab <- organizeOutput(temperPrior = input$temperPrior, 
                             noise = input$noise, 
                             nobs = input$nobs)
    baseInfoTab <- data.frame(outTab$repetitionSize)
    
    
    kernelType = "generalMatern"
    
    pdfBaseLink <- "35.237.17.250/results/"
    if(input$physicalSystem=="FitzHugh-Nagumo (FN) Model"){
      folderPrefix <- "FN-"
    }else if(input$physicalSystem=="Oscillatory expression of the Hes1"){
      folderPrefix <- "hes1-"
    }else if(input$physicalSystem=="HIV model"){
      folderPrefix <- "HIV-"
    }
    subDir <- paste0(folderPrefix, 
                     "withmeanBand-", kernelType, 
                     "-nobs", input$nobs, 
                     "-noise", input$noise, "_", input$noise)
    ndis <- substr(rownames(baseInfoTab), 12, nchar(rownames(baseInfoTab)))
    subDir <- paste0(subDir, "-ndis", ndis)
    if(input$temperPrior){
      subDir <- paste0(subDir, "-temperPrior")
    }else{
      subDir <- paste0(subDir, "-unitHeatPrior")
    }
    subDir <- file.path(pdfBaseLink, subDir)
    urls <- paste0("http://", subDir)
    refs <- paste0("<a href='",  urls, "' target='_blank'>pdf visualization ndis-",ndis,"</a>")
    baseInfoTab$visualizations <- refs
    baseInfoTab
    
  })
  output$repSizeTable <- DT::renderDataTable(
     repSizeTable(), escape = FALSE, options = list(destroy = TRUE,paging=FALSE,searching=FALSE, bInfo=FALSE))
})
