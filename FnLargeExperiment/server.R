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
      load("data/largeExperimentSummary.rda", envir = .GlobalEnv)
      assign("maternDf.candidates", c(2.01, 2.5), envir = .GlobalEnv)
      assign("phaseType.candidates", c("phase1", "phase2", "trueMu", "priorTempered", "priorTemperedPhase2"), envir = .GlobalEnv)
    }else if(physicalSystem == "Oscillatory expression of the Hes1"){
      load("data/hes1-largeExperimentSummary.rda", envir = .GlobalEnv)
      assign("maternDf.candidates", 2.01, envir = .GlobalEnv)
      assign("kernel.candidates", "generalMatern", envir = .GlobalEnv)
      assign("phaseType.candidates", c("trueMu", "priorTempered", "priorTemperedPhase2"), envir = .GlobalEnv)
      noise.candidates <- sapply(noise.candidates, function(noi) paste(round(noi * c(4, 1, 8), 3), collapse = "_"))
      assign("noise.candidates", noise.candidates, envir = .GlobalEnv)
    }else if(physicalSystem == "HIV model"){
      load("data/HIV-largeExperimentSummary.rda", envir = .GlobalEnv)
      assign("maternDf.candidates", 2.01, envir = .GlobalEnv)
      assign("kernel.candidates", "generalMatern", envir = .GlobalEnv)
      assign("phaseType.candidates", c("trueMu", "priorTempered", "priorTemperedPhase2"), envir = .GlobalEnv)
    }
    
    updateSelectInput(session, "nobs", choices=nobs.candidates, selected=nobs.candidates[4])
    updateSelectInput(session, "maternDf", choices=maternDf.candidates, selected=maternDf.candidates[1])
    updateSelectInput(session, "noise", choices=noise.candidates, selected=noise.candidates[2])
    updateSelectInput(session, "phaseType", choices=phaseType.candidates, selected="priorTemperedPhase2")
  })
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
    pdfBaseLink <- "35.227.57.42/results"
    if(input$physicalSystem=="FitzHugh-Nagumo (FN) Model"){
      folderPrefix <- ""
    }else if(input$physicalSystem=="Oscillatory expression of the Hes1"){
      folderPrefix <- "hes1-"
    }else if(input$physicalSystem=="HIV model"){
      folderPrefix <- "HIV-"
    }
    subDir <- paste0(folderPrefix, "withmeanBand-", kernelType, "-nobs", input$nobs, "-noise", input$noise)
    ndis <- substr(rownames(baseInfoTab), 12, nchar(rownames(baseInfoTab)))
    subDir <- paste0(subDir, "-ndis", ndis)
    subDir <- file.path(pdfBaseLink, subDir)
    urls <- paste0("http://", subDir)
    refs <- paste0("<a href='",  urls, "' target='_blank'>pdf visualization ndis-",ndis,"</a>")
    baseInfoTab$visualizations <- refs
    baseInfoTab
    
  })
  output$repSizeTable <- DT::renderDataTable(
     repSizeTable(), escape = FALSE, options = list(destroy = TRUE,paging=FALSE,searching=FALSE, bInfo=FALSE))
})
