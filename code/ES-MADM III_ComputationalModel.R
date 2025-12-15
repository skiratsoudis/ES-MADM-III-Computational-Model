###############################################################################
# Filename: ESMADMIII_FinalAlphaVariation.R
# ES-MADM III: Fuzzy Preference-Enhanced MADM with robust code & per-criterion α
###############################################################################

library(shiny)
library(shinythemes)
library(shinyjs)
library(readxl)
library(writexl)
library(ggplot2)
library(dplyr)
library(DT)
library(gridExtra)
library(reshape2)
library(ggrepel)

########################################
# 1) Helper Functions
########################################

round3 <- function(x) round(x, 3)
format3 <- function(x) sprintf("%.3f", x)
log2safe <- function(x, eps=1e-15) ifelse(x > eps, log2(x), 0)

preferenceFunction <- function(d, funcType, q=0, p=0, s=0) {
  ft <- tolower(funcType)
  eps <- 1e-15
  if(ft=="usual"){
    return(ifelse(d>0,1,0))
  } else if(ft=="u-shape"){
    return(ifelse(d<=q,0,1))
  } else if(ft=="v-shape"){
    if(p<=eps) return(ifelse(d>0,1,0))
    if(d<=0) return(0)
    if(d<p) return(d/p)
    return(1)
  } else if(ft=="level"){
    if(d<=q) return(0)
    if(d<=p) return(0.5)
    return(1)
  } else if(ft=="linear"){
    if(d<=q) return(0)
    if(abs(p-q)<=eps) return(ifelse(d>q,1,0))
    if(d<=p) return((d-q)/(p-q))
    return(1)
  } else if(ft=="gaussian"){
    if(d<=0) return(0)
    if(abs(s)<=eps) return(ifelse(d>0,1,0))
    return(1 - exp(-(d^2)/(2*s^2)))
  } else {
    if(p<=eps) return(ifelse(d>0,1,0))
    if(d<=0) return(0)
    if(d<p) return(d/p)
    return(1)
  }
}

computePreferenceEnhancedRow <- function(row, funcType, q=0, p=0, s=0){
  N <- length(row)
  F <- numeric(N)
  for(j in 1:N){
    diffs <- row[j]-row[-j]
    prefs <- sapply(diffs, function(d) preferenceFunction(d, funcType, q, p, s))
    F[j] <- mean(prefs)
  }
  return(F)
}

fuzzyNormalizationVector <- function(dataMat, deltaMat, typeVec, alphaVec, eps=1e-15){
  M <- nrow(dataMat)
  N <- ncol(dataMat)
  normData <- matrix(0, nrow=M, ncol=N)
  for(i in 1:M){
    alpha <- alphaVec[i]
    lower <- dataMat[i,] - (1-alpha)*deltaMat[i,]
    upper <- dataMat[i,] + alpha*deltaMat[i,]
    if(tolower(typeVec[i])=="benefit"){
      maxUpper <- max(upper)+eps
      normData[i,] <- lower/maxUpper
    } else {
      maxVal <- max(dataMat[i,]+deltaMat[i,])+eps
      inv <- (maxVal - lower)
      normData[i,] <- inv/(sum(inv)+eps)
    }
  }
  return(normData)
}

fuzzySubjectiveWeights <- function(sbjDF, alphaVec){
  M <- nrow(sbjDF)
  lower <- numeric(M)
  upper <- numeric(M)
  for(i in 1:M){
    cVal <- as.numeric(sbjDF$Central[i])
    dVal <- as.numeric(sbjDF$Delta[i])
    aVal <- alphaVec[i]
    lower[i] <- cVal - (1-aVal)*dVal
    upper[i] <- cVal + aVal*dVal
  }
  names(lower) <- sbjDF$Criterion
  names(upper) <- sbjDF$Criterion
  list(lower=lower, upper=upper)
}

computeESMADMIII <- function(dataMat, deltaMat, sbjDF, prefParams, benefitCost, alphaVec){
  M <- nrow(dataMat)
  N <- ncol(dataMat)
  
  critTypes <- sapply(rownames(dataMat), function(x){
    idx <- which(benefitCost$Criterion==x)
    if(length(idx)==0) return("benefit") else return(tolower(benefitCost$Type[idx]))
  })
  
  normData <- fuzzyNormalizationVector(dataMat, deltaMat, critTypes, alphaVec)
  
  prefMatrix <- matrix(0, nrow=M, ncol=N)
  rownames(prefMatrix) <- rownames(dataMat)
  colnames(prefMatrix) <- colnames(dataMat)
  for(i in 1:M){
    row_i <- normData[i,]
    critName <- rownames(dataMat)[i]
    rowParam <- prefParams[prefParams$Criterion==critName, ]
    if(nrow(rowParam)==0){
      funcType <- "v-shape"
      q_val <- 0; p_val <- 0.1; s_val <- 0
    } else {
      funcType <- as.character(rowParam$PreferenceFunction)
      q_val <- as.numeric(rowParam$Threshold_q)
      p_val <- as.numeric(rowParam$Threshold_p)
      s_val <- as.numeric(rowParam$Threshold_s)
      if(is.na(q_val)) q_val <- 0
      if(is.na(p_val)) p_val <- 0
      if(is.na(s_val)) s_val <- 0
    }
    F_i <- computePreferenceEnhancedRow(row_i, funcType, q_val, p_val, s_val)
    product <- row_i*F_i
    prefMatrix[i,] <- product/(sum(product)+1e-15)
  }
  condProb <- prefMatrix
  
  h <- numeric(M)
  for(i in 1:M){
    row_i <- condProb[i,]
    nonz <- row_i>0
    val <- -sum(row_i[nonz]*log2safe(row_i[nonz]))/log2(N)
    h[i] <- min(max(val,0),1)
  }
  d <- 1-h
  objW <- d/sum(d+1e-15)
  
  sbjFz <- fuzzySubjectiveWeights(sbjDF, alphaVec)
  prod_lower <- objW * sbjFz$lower
  prod_upper <- objW * sbjFz$upper
  sum_upper <- sum(objW*sbjFz$upper)+1e-15
  sum_lower <- sum(objW*sbjFz$lower)+1e-15
  intW_lower <- prod_lower/sum_upper
  intW_upper <- prod_upper/sum_lower
  intW <- (intW_lower+intW_upper)/2
  
  condEntropy <- numeric(M)
  for(i in 1:M){
    row_i <- condProb[i,]
    nonz <- row_i>0
    cval <- -sum(row_i[nonz]*log2safe(row_i[nonz]))
    condEntropy[i] <- min(max(cval,0),log2(N))
  }
  
  altScores <- numeric(N)
  for(j in 1:N){
    altScores[j] <- sum(condProb[,j]*intW)
  }
  altNonzero <- altScores[altScores > 0]
  valY <- -sum(altNonzero * log2safe(altNonzero))
  overallEntropyY <- min(max(valY,0), log2(N))
  
  totalCondEntropy <- sum(intW*condEntropy)
  normTotalCondEntropy <- totalCondEntropy/(overallEntropyY+1e-15)
  
  jointP <- matrix(0, nrow=M, ncol=N)
  for(i in 1:M){
    for(j in 1:N){
      jointP[i,j] <- condProb[i,j]*intW[i]
    }
  }
  jNonz <- jointP[jointP>0]
  valXY <- -sum(jNonz*log2safe(jNonz))
  SXY <- min(max(valXY,0),log2(M*N))
  
  Sx <- -sum(intW[intW>0]*log2safe(intW[intW>0]))
  Sx <- max(0,Sx)
  
  Jxy <- Sx + overallEntropyY - SXY
  Jxy <- min(max(Jxy,0), min(Sx, overallEntropyY))
  
  IXY <- SXY / (Sx + overallEntropyY + 1e-15)
  IXY <- pmin(1, pmax(0, IXY))
  
  partialMI <- numeric(M)
  for(i in 1:M){
    partialMI[i] <- max(0, overallEntropyY - condEntropy[i])
  }
  CES <- (1/M)*sum(partialMI/(overallEntropyY+1e-15))
  CSF <- 1-(totalCondEntropy/(overallEntropyY+1e-15))
  ADI <- 1-(overallEntropyY/(log2(N)+1e-15))
  NMI <- Jxy/(overallEntropyY+1e-15)
  
  sumAll <- NMI + CES + ADI
  NMGI <- 0
  if(sumAll>1e-15){
    p_nmi <- NMI/sumAll
    p_ces <- CES/sumAll
    p_adi <- ADI/sumAll
    K <- 3
    k <- 1/log2(K)
    e_nmi <- ifelse(p_nmi>0, -k*p_nmi*log2(p_nmi),0)
    e_ces <- ifelse(p_ces>0, -k*p_ces*log2(p_ces),0)
    e_adi <- ifelse(p_adi>0, -k*p_adi*log2(p_adi),0)
    d_nmi <- 1-e_nmi
    d_ces <- 1-e_ces
    d_adi <- 1-e_adi
    d_sum <- d_nmi+d_ces+d_adi
    if(d_sum>1e-15){
      w_nmi <- d_nmi/d_sum
      w_ces <- d_ces/d_sum
      w_adi <- d_adi/d_sum
      valNMGI <- w_nmi*NMI + w_ces*CES + w_adi*ADI
      NMGI <- pmin(1, pmax(0, valNMGI))
    }
  }
  
  condEntropyNorm <- condEntropy/(log2(N)+1e-15)
  critImport <- intW * (1 - condEntropyNorm)
  
  list(
    M=M, N=N,
    xOBJ=objW,
    xINT=intW,
    critImport=critImport,
    altScores=altScores,
    overallEntropyY=overallEntropyY,
    totalCondEntropy=totalCondEntropy,
    normTotalCondEntropy=normTotalCondEntropy,
    SXY=SXY,
    Sx=Sx,
    Jxy=Jxy,
    IXY=IXY,
    NMI=NMI,
    CES=CES,
    CSF=CSF,
    ADI=ADI,
    NMGI=NMGI
  )
}

generateDetailedEntropyInsights <- function(SXY, SX, SY, IXY, SYX, IYX, NMI, ADI, CES, CSF, NMGI){
  insights <- c()
  if(SXY>1.5){
    insights <- c(insights, "<li>A high joint entropy S(X,Y) indicates complex interactions with significant diversity.</li>")
  } else {
    insights <- c(insights, "<li>A low joint entropy S(X,Y) suggests limited informational diversity among criteria & alternatives.</li>")
  }
  if(SX>1.0){
    insights <- c(insights, "<li>Criteria entropy S(X) is well-distributed, supporting robust decision-making.</li>")
  } else {
    insights <- c(insights, "<li>A lower S(X) indicates that a few criteria dominate.</li>")
  }
  if(SY>1.0){
    insights <- c(insights, "<li>The alternatives entropy S(Y) highlights a diverse set of options.</li>")
  } else {
    insights <- c(insights, "<li>A lower S(Y) indicates limited distinction among alternatives.</li>")
  }
  if(IYX>0.5){
    insights <- c(insights, "<li>Mutual information I(X;Y) shows criteria effectively explain the alternatives.</li>")
  } else {
    insights <- c(insights, "<li>A lower I(X;Y) suggests criteria do not fully capture alternatives&#39; variability.</li>")
  }
  if(SYX<0.5){
    insights <- c(insights, "<li>The conditional entropy I(Y|X) demonstrates strong uncertainty reduction given the criteria.</li>")
  } else {
    insights <- c(insights, "<li>A higher I(Y|X) implies residual ambiguity after considering the criteria.</li>")
  }
  if(IXY>0.7){
    insights <- c(insights, "<li>A high normalized joint entropy I(X,Y) indicates well-balanced information distribution.</li>")
  } else {
    insights <- c(insights, "<li>A lower I(X,Y) suggests redundancy or imbalance in the criteria&#39;s discrimination.</li>")
  }
  if(NMI>0.7){
    insights <- c(insights, "<li>High NMI indicates strong correlation between criteria & alternatives, reflecting robust explanatory power.</li>")
  } else {
    insights <- c(insights, "<li>A lower NMI indicates weaker correlation, suggesting refinements needed.</li>")
  }
  if(ADI>0.8){
    insights <- c(insights, "<li>A high ADI reflects strong distinctions among alternatives, aiding clarity.</li>")
  } else {
    insights <- c(insights, "<li>A lower ADI indicates overlapping alternatives, complicating the final choice.</li>")
  }
  if(CES>0.7){
    insights <- c(insights, "<li>A high CES shows criteria effectively reduce uncertainty in the alternatives.</li>")
  } else {
    insights <- c(insights, "<li>A lower CES suggests the criteria have limited impact on uncertainty reduction.</li>")
  }
  if(CSF>0.75){
    insights <- c(insights, "<li>A high CSF signifies stable decisions with minimal residual uncertainty.</li>")
  } else {
    insights <- c(insights, "<li>A lower CSF implies potential instability, requiring further refinements.</li>")
  }
  if(NMGI>0.8){
    insights <- c(insights, "<li>A high NMGI indicates a robust, well-integrated, and effective decision process.</li>")
  } else {
    insights <- c(insights, "<li>A lower NMGI highlights areas for synergy & stability improvements.</li>")
  }
  paste("<ul>", paste(insights, collapse=""), "</ul>")
}

########################################
# 2) UI
########################################

ui <- fluidPage(
  theme=shinytheme("slate"),
  useShinyjs(),
  tags$head(
    tags$style(HTML("
      table.dataTable tbody td { color: #FFFFFF !important; }
      table.dataTable thead th { color: #FFFF00 !important; }
      .dataTables_wrapper .dataTables_length label,
      .dataTables_wrapper .dataTables_filter label,
      .dataTables_wrapper .dataTables_info {
        color: #FFFF00 !important;
      }
      .dataTables_wrapper .dataTables_paginate .paginate_button {
        color: #000000 !important;
        background-color: #CCCCCC !important;
      }
      .dataTables_wrapper .dataTables_length select,
      .dataTables_wrapper .dataTables_filter input {
        background-color: #FFFFFF !important;
        color: #000000 !important;
      }
    "))
  ),
  
  titlePanel("ES-MADM III: Fuzzy Preference-Entropy Synergy MADM Model (Version III)"),
  
  navbarPage(
    title="ES-MADM III",
    id="tabs",
    
    tabPanel("1. Instructions",
             fluidRow(
               column(12,
                      h3("ES-MADM III with Per-Criterion α Variation"),
                      tags$ul(
                        tags$li("Load 6 sheets in Tab 2 (DataMatrix, FuzzyDeviations, FuzzySubjectiveWeights, PreferenceParams, BenefitCost, FuzzyAlpha_Criteria)."),
                        tags$li("If FuzzyAlpha_Criteria is missing, we default to α=0.5 for all criteria."),
                        tags$li("Then see baseline results in Tabs 3–5."),
                        tags$li("Sensitivity Analysis in Tab 6 includes data & weight changes, plus an 'Alpha Variation' sub-tab for custom α scenarios.")
                      ),
                      hr(),
                      h4("Key Notations"),
                      DTOutput("tableNotation")
               )
             )
    ),
    
    tabPanel("2. Data Import",
             sidebarLayout(
               sidebarPanel(
                 fileInput("fileExcel","Upload Excel File:", accept=c(".xls",".xlsx")),
                 textInput("sheetMatrix","Sheet for DataMatrix:",value="DataMatrix"),
                 textInput("sheetDelta","Sheet for FuzzyDeviations:",value="FuzzyDeviations"),
                 textInput("sheetWeights","Sheet for FuzzySubjectiveWeights:",value="FuzzySubjectiveWeights"),
                 textInput("sheetPrefParams","Sheet for PreferenceParams:",value="PreferenceParams"),
                 textInput("sheetBenefitCost","Sheet for BenefitCost:",value="BenefitCost"),
                 textInput("sheetAlpha","Sheet for FuzzyAlpha_Criteria:",value="FuzzyAlpha_Criteria"),
                 actionButton("loadData","Load Data",class="btn-primary")
               ),
               mainPanel(
                 h4("DataMatrix Preview"),
                 DTOutput("tableMatrixPreview"),
                 h4("FuzzyDeviations Preview"),
                 DTOutput("tableDeltaPreview"),
                 h4("FuzzySubjectiveWeights Preview"),
                 DTOutput("tableWeightsPreview"),
                 h4("PreferenceParams Preview"),
                 DTOutput("tablePrefParamsPreview"),
                 h4("BenefitCost Preview"),
                 DTOutput("tableBCPreview"),
                 h4("FuzzyAlpha_Criteria Preview"),
                 DTOutput("tableAlphaPreview")
               )
             )
    ),
    
    tabPanel("3. Alternative Results",
             sidebarLayout(
               sidebarPanel(),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Alternatives Table", DTOutput("tableAltScores")),
                   tabPanel("Alternatives Plot", plotOutput("plotAltScores",height="600px"))
                 )
               )
             )
    ),
    
    tabPanel("4. Criteria Results",
             sidebarLayout(
               sidebarPanel(),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Criteria Table", DTOutput("tableCriteriaResults")),
                   tabPanel("Criteria Plots",
                            tabsetPanel(
                              tabPanel("a) xSBJ Weights", plotOutput("plotxSBJWeights", height="600px")),
                              tabPanel("b) xOBJ Weights", plotOutput("plotxOBJWeights", height="600px")),
                              tabPanel("c) xINT Weights", plotOutput("plotxINTWeights", height="600px")),
                              tabPanel("d) Comparative Weights",
                                       checkboxGroupInput("selectWeights","Select Weights to Display:",
                                                          choices=list("xSBJ"="xSBJ","xOBJ"="xOBJ","xINT"="xINT"),
                                                          selected=c("xSBJ","xOBJ","xINT")),
                                       plotOutput("plotComparativeWeights", height="800px")
                              ),
                              tabPanel("e) Integrated Criteria Importance", plotOutput("plotIntegratedCI",height="600px"))
                            )
                   )
                 )
               )
             )
    ),
    
    tabPanel("5. Entropy Measures & Indices",
             sidebarLayout(
               sidebarPanel(helpText("View baseline global entropies & key indices.")),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Entropy Measures", DTOutput("tableGlobalEntropies")),
                   tabPanel("Entropy Indices", DTOutput("tableIndices")),
                   tabPanel("NMGI", DTOutput("tableNMGI")),
                   tabPanel("Entropy Measures Plot", plotOutput("plotEntropyMeasures",height="600px")),
                   tabPanel("Entropy Indices Plot", plotOutput("plotIndices",height="600px"))
                 )
               )
             )
    ),
    
    tabPanel("6. Sensitivity Analysis",
             sidebarLayout(
               sidebarPanel(
                 h4("Sensitivity Controls (Data/Weights/Preferences)"),
                 tabsetPanel(
                   tabPanel("Criteria Weight",
                            selectInput("sens_weight_criterion","Select Criterion:", choices=NULL),
                            numericInput("sens_weight_value","New Fuzzy SBJ Weight (Central):",1,step=0.1),
                            numericInput("sens_weight_delta","New Fuzzy SBJ Weight (Delta):",0.05,step=0.01),
                            actionButton("apply_weight","Apply Weight Change")
                   ),
                   tabPanel("Data Matrix Element",
                            selectInput("sens_data_criterion","Select Criterion:",choices=NULL),
                            selectInput("sens_data_alternative","Select Alternative:",choices=NULL),
                            numericInput("sens_data_value","New Data Value (ξ):",1,step=0.1),
                            actionButton("apply_data","Apply Data Change")
                   ),
                   tabPanel("Criteria Type",
                            selectInput("sens_type_criterion","Select Criterion:",choices=NULL),
                            radioButtons("sens_type_value","Select Type:",choices=c("benefit","cost"),inline=TRUE),
                            actionButton("apply_type","Apply Type Change")
                   ),
                   tabPanel("Preference Function",
                            selectInput("sens_pf_criterion","Select Criterion:",choices=NULL),
                            selectInput("sens_pf_type","Select Preference Function:",
                                        choices=c("Usual","U-Shape","V-Shape","Level","Linear","Gaussian")),
                            actionButton("apply_pf","Apply Preference Function Change")
                   ),
                   tabPanel("Threshold",
                            selectInput("sens_threshold_criterion","Select Criterion:",choices=NULL),
                            numericInput("sens_threshold_q","New q:",0,step=0.01),
                            numericInput("sens_threshold_p","New p:",0,step=0.01),
                            numericInput("sens_threshold_s","New s:",0,step=0.01),
                            actionButton("apply_threshold","Apply Threshold Change")
                   )
                 ),
                 hr(),
                 actionButton("reset_sensitivity","Reset to Baseline", class="btn-danger")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Modified Results",
                            DTOutput("tableModAltScores"),
                            DTOutput("tableModCriteriaResults"),
                            plotOutput("plotModAltScores",height="600px"),
                            plotOutput("plotModEntropy",height="600px")
                   ),
                   tabPanel("Comparison",
                            DTOutput("tableCompIndices"),
                            plotOutput("plotCompIndices",height="600px")
                   ),
                   tabPanel("Alpha Variation",
                            fluidRow(
                              column(12,
                                     h4("Adjust α-cuts per Criterion"),
                                     helpText("Edit the α column (0..1). Then 'Apply α Changes' to re-run."),
                                     DTOutput("tableAlpha"),
                                     numericInput("universalAlpha","Universal α:",0.5,min=0,max=1,step=0.05),
                                     actionButton("applyUniversalAlpha","Set All Rows to Universal α"),
                                     hr(),
                                     actionButton("applyAlphaChanges","Apply α Changes", class="btn-info"),
                                     hr(),
                                     h4("Alpha Variation Results Comparison"),
                                     DTOutput("tableAlphaComp"),
                                     plotOutput("plotAlphaNMGI", height="500px")
                              )
                            )
                   )
                 )
               )
             )
    ),
    
    tabPanel("7. Summary",
             fluidRow(
               column(3,
                      div(style="text-align:center;",
                          img(src="https://cdn.thecollector.com/wp-content/uploads/2024/07/thinker-auguste-rodin-what-so-special.jpg?width=1400&quality=70",
                              width="160px", height="160px", style="border-radius:50%; border:2px solid black;"),
                          downloadButton("downloadAllSummary","Download Results (Excel)", class="btn-primary", style="margin-top:15px;")
                      )
               ),
               column(9,
                      h3("Summary of Baseline Results"),
                      uiOutput("summaryResults")
               )
             )
    ),
    
    tabPanel("8. Settings",
             fluidRow(
               column(12,
                      h3("Settings"),
                      fluidRow(
                        column(4,
                               h4("Font Size Adjustment"),
                               sliderInput("fontSize","Select Font Size:",min=12,max=24,value=16,step=1)
                        )
                      )
               )
             )
    )
  )
)

########################################
# 4) SERVER
########################################

server <- function(input, output, session){
  
  storedDataMat <- reactiveVal(NULL)
  storedDeltaMat <- reactiveVal(NULL)
  storedSBJw <- reactiveVal(NULL)
  storedPrefParams <- reactiveVal(NULL)
  storedBenefitCost <- reactiveVal(NULL)
  storedAlphaCrit <- reactiveVal(NULL)
  
  scenario <- reactiveValues(
    dataMat=NULL,
    SBJw=NULL,
    PrefParams=NULL,
    BenefitCost=NULL
  )
  
  alphaDF <- reactiveVal(NULL)
  
  baselineResults <- reactive({
    req(storedDataMat(), storedDeltaMat(), storedSBJw(), storedPrefParams(), storedBenefitCost())
    mat <- storedDataMat()
    M <- nrow(mat)
    if(!is.null(storedAlphaCrit())){
      dfA <- storedAlphaCrit()
      rownames(dfA) <- dfA$Criterion
      alphaVec <- numeric(M)
      for(i in 1:M){
        cr <- rownames(mat)[i]
        if(cr %in% dfA$Criterion) alphaVec[i] <- as.numeric(dfA[cr,"Alpha"]) else alphaVec[i] <- 0.5
      }
    } else {
      alphaVec <- rep(0.5, M)
    }
    computeESMADMIII(
      dataMat=mat,
      deltaMat=storedDeltaMat(),
      sbjDF=storedSBJw(),
      prefParams=storedPrefParams(),
      benefitCost=storedBenefitCost(),
      alphaVec=alphaVec
    )
  })
  
  observeEvent(input$loadData, {
    req(input$fileExcel)
    safeReadSheet <- function(sheetRef, desc="unknown"){
      out <- NULL
      tryCatch({
        out <- read_excel(input$fileExcel$datapath, sheet=sheetRef, col_names=TRUE)
      }, error=function(e){
        showNotification(paste("Could not read sheet for",desc,"->",sheetRef,"Error:",e$message), type="error")
      })
      out
    }
    dfMatrix <- safeReadSheet(input$sheetMatrix,"DataMatrix")
    if(is.null(dfMatrix)) return(NULL)
    critNames <- dfMatrix[[1]]
    dataMatrix <- as.matrix(dfMatrix[,-1])
    rownames(dataMatrix) <- critNames
    colnames(dataMatrix) <- colnames(dfMatrix)[-1]
    
    dfDelta <- safeReadSheet(input$sheetDelta,"FuzzyDeviations")
    if(is.null(dfDelta)) return(NULL)
    deltaMatrix <- as.matrix(dfDelta[,-1])
    rownames(deltaMatrix) <- dfDelta[[1]]
    colnames(deltaMatrix) <- colnames(dfDelta)[-1]
    
    dfSbj <- safeReadSheet(input$sheetWeights,"FuzzySubjectiveWeights")
    if(is.null(dfSbj)) return(NULL)
    
    dfPref <- safeReadSheet(input$sheetPrefParams,"PreferenceParams")
    if(is.null(dfPref)) return(NULL)
    
    dfBC <- safeReadSheet(input$sheetBenefitCost,"BenefitCost")
    if(is.null(dfBC)) return(NULL)
    
    dfAlpha <- safeReadSheet(input$sheetAlpha,"FuzzyAlpha_Criteria")
    if(is.null(dfAlpha)){
      dfAlpha <- data.frame(Criterion=critNames,Alpha=rep(0.5,length(critNames)),stringsAsFactors=FALSE)
    }
    
    storedDataMat(dataMatrix)
    storedDeltaMat(deltaMatrix)
    storedSBJw(dfSbj)
    storedPrefParams(dfPref)
    storedBenefitCost(dfBC)
    storedAlphaCrit(dfAlpha)
    
    scenario$dataMat <- dataMatrix
    scenario$SBJw <- dfSbj
    scenario$PrefParams <- dfPref
    scenario$BenefitCost <- dfBC
    
    alphaDF(dfAlpha)
    
    updateSelectInput(session,"sens_weight_criterion", choices=critNames)
    updateSelectInput(session,"sens_data_criterion", choices=critNames)
    updateSelectInput(session,"sens_type_criterion", choices=critNames)
    updateSelectInput(session,"sens_pf_criterion", choices=critNames)
    updateSelectInput(session,"sens_threshold_criterion", choices=critNames)
    updateSelectInput(session,"sens_data_alternative", choices=colnames(dataMatrix))
    
    showNotification("Data loaded successfully!",type="message")
  })
  
  output$tableMatrixPreview <- renderDT({
    req(storedDataMat())
    df <- as.data.frame(storedDataMat())
    datatable(df, options=list(pageLength=5), rownames=TRUE) %>%
      formatRound(columns=1:ncol(df),digits=3)
  })
  output$tableDeltaPreview <- renderDT({
    req(storedDeltaMat())
    df <- as.data.frame(storedDeltaMat())
    datatable(df,options=list(pageLength=5),rownames=TRUE) %>%
      formatRound(columns=1:ncol(df),digits=3)
  })
  output$tableWeightsPreview <- renderDT({
    req(storedSBJw())
    df <- storedSBJw()
    datatable(df,options=list(pageLength=5),rownames=FALSE) %>%
      formatRound(columns=c("Central","Delta"),digits=3)
  })
  output$tablePrefParamsPreview <- renderDT({
    req(storedPrefParams())
    df <- storedPrefParams()
    datatable(df,options=list(pageLength=5),rownames=FALSE) %>%
      formatRound(columns=c("Threshold_q","Threshold_p","Threshold_s"),digits=3)
  })
  output$tableBCPreview <- renderDT({
    req(storedBenefitCost())
    df <- storedBenefitCost()
    datatable(df,options=list(pageLength=5),rownames=FALSE)
  })
  output$tableAlphaPreview <- renderDT({
    req(storedAlphaCrit())
    datatable(storedAlphaCrit(),options=list(pageLength=5),rownames=FALSE) %>%
      formatRound(columns="Alpha",digits=3)
  })
  
  output$tableNotation <- renderDT({
    df <- data.frame(
      Symbol=c("ξ(μν)","x(μ)^SBJ","x(μ)^OBJ","x(μ)^INT",
               "S(X,Y)","S(X)","S(Y)","I(X;Y)","I(Y|X)","I(X,Y)",
               "NMI","ADI","CES","CSF","NMGI"),
      Meaning=c(
        "Fuzzy performance (central) of criterion μ on alt ν",
        "Fuzzy subjective weight (Central±Delta)",
        "Objective weight (entropy-based)",
        "Integrated weight (fuzzy SBJ + OBJ)",
        "Joint entropy of X & Y",
        "Entropy of criteria distribution",
        "Entropy of alternatives distribution",
        "Mutual information between X & Y",
        "Normalized conditional entropy of Y given X",
        "Normalized joint entropy of X & Y",
        "Normalized Mutual Information",
        "Alternatives Distinction Index",
        "Criteria Effectiveness Score",
        "Conditional Stability Factor",
        "Net Mutual Growth Index"
      )
    )
    datatable(df, options=list(pageLength=20))
  })
  
  output$tableAltScores <- renderDT({
    r <- baselineResults(); req(r)
    altNames <- colnames(storedDataMat())
    df <- data.frame(Alternative=altNames, Score=r$altScores)
    df <- df[order(-df$Score),]
    datatable(df,options=list(pageLength=10),rownames=FALSE) %>%
      formatRound(columns="Score",digits=3)
  })
  output$plotAltScores <- renderPlot({
    r <- baselineResults(); req(r)
    altNames <- colnames(storedDataMat())
    df_alt <- data.frame(Alternative=altNames, Score=r$altScores)
    df_alt <- df_alt[order(-df_alt$Score),]
    ggplot(df_alt, aes(x=reorder(Alternative,-Score), y=Score))+
      geom_bar(stat="identity", fill="coral")+
      geom_text(aes(label=format3(Score)), vjust=-0.3, color="white", size=5)+
      labs(title="Baseline Alternative Scores", x="Alternative", y="Score")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  
  output$tableCriteriaResults <- renderDT({
    r <- baselineResults(); req(r)
    critNames <- rownames(storedDataMat())
    df <- data.frame(
      Criterion=critNames,
      xOBJ_Weight=r$xOBJ,
      xINT_Weight=r$xINT,
      Integrated_Criteria_Importance=r$critImport
    )
    datatable(df,options=list(pageLength=10),rownames=FALSE) %>%
      formatRound(columns=c("xOBJ_Weight","xINT_Weight","Integrated_Criteria_Importance"),digits=3)
  })
  
  output$plotxSBJWeights <- renderPlot({
    req(storedSBJw())
    sbjDF <- storedSBJw()
    df_subj <- data.frame(Criterion=sbjDF$Criterion, Weight=sbjDF$Central)
    df_subj <- df_subj[order(-df_subj$Weight),]
    ggplot(df_subj, aes(x=reorder(Criterion,-Weight), y=Weight))+
      geom_bar(stat="identity", fill="darkgreen")+
      geom_text(aes(label=format3(Weight)), vjust=-0.3, color="white", size=5)+
      labs(title="a) xSBJ Criteria Weights (Central Values)", x="Criterion", y="Weight")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  output$plotxOBJWeights <- renderPlot({
    r <- baselineResults(); req(r)
    critNames <- rownames(storedDataMat())
    df_obj <- data.frame(Criterion=critNames, Weight=r$xOBJ)
    df_obj <- df_obj[order(-df_obj$Weight),]
    ggplot(df_obj, aes(x=reorder(Criterion,-Weight), y=Weight))+
      geom_bar(stat="identity", fill="dodgerblue")+
      geom_text(aes(label=format3(Weight)), vjust=-0.3, color="white", size=5)+
      labs(title="b) xOBJ Criteria Weights (Entropy-Based)", x="Criterion", y="Weight")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  output$plotxINTWeights <- renderPlot({
    r <- baselineResults(); req(r)
    critNames <- rownames(storedDataMat())
    df_int <- data.frame(Criterion=critNames, Weight=r$xINT)
    df_int <- df_int[order(-df_int$Weight),]
    ggplot(df_int, aes(x=reorder(Criterion,-Weight), y=Weight))+
      geom_bar(stat="identity", fill="orange")+
      geom_text(aes(label=format3(Weight)), vjust=-0.3, color="white", size=5)+
      labs(title="c) xINT Criteria Weights (Integrated)", x="Criterion", y="Weight")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  output$plotComparativeWeights <- renderPlot({
    r <- baselineResults(); req(r)
    sbjDF <- storedSBJw()
    critNames <- rownames(storedDataMat())
    xSBJ <- sapply(critNames, function(cn){
      rowIdx <- which(sbjDF$Criterion==cn)
      if(length(rowIdx)==0) return(0) else return(sbjDF$Central[rowIdx])
    })
    df_comp <- data.frame(
      Criterion=critNames,
      xSBJ=xSBJ,
      xOBJ=r$xOBJ,
      xINT=r$xINT
    )
    df_melt <- melt(df_comp, id.vars="Criterion", variable.name="Weight_Type", value.name="Weight")
    selected_weights <- input$selectWeights
    if(is.null(selected_weights)) selected_weights<-c("xSBJ","xOBJ","xINT")
    df_melt <- df_melt %>% filter(Weight_Type %in% selected_weights)
    ggplot(df_melt, aes(x=Criterion, y=Weight, color=Weight_Type, group=Weight_Type))+
      geom_line(size=1.2)+
      geom_point(size=3)+
      geom_text_repel(aes(label=format3(Weight)),size=3,color="white",
                      box.padding=0.1,point.padding=0.1,segment.color="transparent",
                      max.overlaps=Inf,show.legend=FALSE)+
      labs(title="d) Comparative Criteria Weights", x="Criterion", y="Weight")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        legend.title=element_text(color="#FFFF00",size=16),
        legend.text=element_text(color="#FFFFFF",size=14),
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )+
      scale_color_manual(values=c("xSBJ"="darkgreen","xOBJ"="dodgerblue","xINT"="orange"))+
      coord_cartesian(clip='off')
  })
  output$plotIntegratedCI <- renderPlot({
    r <- baselineResults(); req(r)
    critNames <- rownames(storedDataMat())
    df_ci <- data.frame(Criterion=critNames, ICI=r$critImport)
    df_ci <- df_ci[order(-df_ci$ICI),]
    ggplot(df_ci, aes(x=reorder(Criterion,-ICI), y=ICI))+
      geom_bar(stat="identity", fill="cyan")+
      geom_text(aes(label=format3(ICI)),vjust=-0.3,color="white",size=5)+
      labs(title="e) Integrated Criteria Importance", x="Criterion", y="ICI")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  
  output$tableGlobalEntropies <- renderDT({
    r <- baselineResults(); req(r)
    df <- data.frame(
      "S(Y)"=r$overallEntropyY,
      "S(Y|X)"=r$totalCondEntropy,
      "I(Y|X)"=r$normTotalCondEntropy,
      "S(X,Y)"=r$SXY,
      "S(X)"=r$Sx,
      "I(X,Y)"=r$IXY
    )
    df_t <- as.data.frame(t(df))
    colnames(df_t) <- "Value"
    datatable(df_t, options=list(pageLength=10), rownames=TRUE) %>%
      formatRound(columns="Value",digits=3)
  })
  output$tableIndices <- renderDT({
    r <- baselineResults(); req(r)
    df <- data.frame(
      NMI=r$NMI,
      CES=r$CES,
      CSF=r$CSF,
      ADI=r$ADI
    )
    datatable(df, options=list(pageLength=10), rownames=FALSE) %>%
      formatRound(columns=c("NMI","CES","CSF","ADI"),digits=3)
  })
  output$tableNMGI <- renderDT({
    r <- baselineResults(); req(r)
    df <- data.frame(NMGI=r$NMGI)
    datatable(df, options=list(pageLength=10), rownames=FALSE) %>%
      formatRound(columns="NMGI",digits=3)
  })
  output$plotEntropyMeasures <- renderPlot({
    r <- baselineResults(); req(r)
    df_entropy <- data.frame(
      Measure=c("S(Y)","S(Y|X)","I(Y|X)","S(X,Y)","S(X)","I(X,Y)"),
      Value=c(r$overallEntropyY,r$totalCondEntropy,r$normTotalCondEntropy,r$SXY,r$Sx,r$IXY)
    )
    df_entropy$Measure <- factor(df_entropy$Measure, levels=df_entropy$Measure)
    ggplot(df_entropy, aes(x=Measure, y=Value, fill=Measure))+
      geom_bar(stat="identity")+
      geom_text(aes(label=format3(Value)), vjust=-0.3,color="white",size=5)+
      labs(title="Baseline Entropy Measures", x="Measure", y="Value")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        legend.position="none",
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  output$plotIndices <- renderPlot({
    r <- baselineResults(); req(r)
    df_idx <- data.frame(
      Index=c("NMI","CES","CSF","ADI","NMGI"),
      Value=c(r$NMI,r$CES,r$CSF,r$ADI,r$NMGI)
    )
    df_idx <- df_idx[order(-df_idx$Value),]
    ggplot(df_idx, aes(x=reorder(Index, -Value), y=Value))+
      geom_bar(stat="identity", fill="purple")+
      geom_text(aes(label=format3(Value)),vjust=-0.3,color="white",size=5)+
      labs(title="Baseline Entropy Indices (Descending)", x="Index", y="Value")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  
  modifiedResults <- reactive({
    req(scenario$dataMat, scenario$SBJw, scenario$PrefParams, scenario$BenefitCost)
    M <- nrow(scenario$dataMat)
    alphaVec <- rep(0.5, M)
    computeESMADMIII(
      dataMat=scenario$dataMat,
      deltaMat=storedDeltaMat(),
      sbjDF=scenario$SBJw,
      prefParams=scenario$PrefParams,
      benefitCost=scenario$BenefitCost,
      alphaVec=alphaVec
    )
  })
  
  output$tableModAltScores <- renderDT({
    r <- modifiedResults(); req(r)
    altNames <- colnames(scenario$dataMat)
    df <- data.frame(Alternative=altNames, Score=r$altScores)
    df <- df[order(-df$Score),]
    datatable(df,options=list(pageLength=10),rownames=FALSE) %>%
      formatRound(columns="Score",digits=3)
  })
  output$tableModCriteriaResults <- renderDT({
    r <- modifiedResults(); req(r)
    critNames <- rownames(scenario$dataMat)
    df <- data.frame(
      Criterion=critNames,
      xOBJ_Weight=r$xOBJ,
      xINT_Weight=r$xINT,
      Integrated_Criteria_Importance=r$critImport
    )
    datatable(df,options=list(pageLength=10),rownames=FALSE) %>%
      formatRound(columns=c("xOBJ_Weight","xINT_Weight","Integrated_Criteria_Importance"),digits=3)
  })
  output$plotModAltScores <- renderPlot({
    r <- modifiedResults(); req(r)
    altNames <- colnames(scenario$dataMat)
    dfA <- data.frame(Alternative=altNames, Score=r$altScores)
    dfA <- dfA[order(-dfA$Score),]
    ggplot(dfA, aes(x=reorder(Alternative,-Score), y=Score))+
      geom_bar(stat="identity", fill="coral")+
      geom_text(aes(label=format3(Score)),vjust=-0.3,color="white",size=5)+
      labs(title="Modified Alternative Scores", x="Alternative", y="Score")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  output$plotModEntropy <- renderPlot({
    r <- modifiedResults(); req(r)
    df_entropy <- data.frame(
      Metric=c("S(Y)","S(Y|X)","I(Y|X)","S(X,Y)","S(X)","NMI","CES","CSF","ADI","NMGI"),
      Value=c(r$overallEntropyY,r$totalCondEntropy,r$normTotalCondEntropy,
              r$SXY,r$Sx,r$NMI,r$CES,r$CSF,r$ADI,r$NMGI)
    )
    ggplot(df_entropy, aes(x=reorder(Metric,-Value), y=Value, fill=Metric))+
      geom_bar(stat="identity")+
      geom_text(aes(label=format3(Value)), vjust=-0.3,color="white",size=4)+
      labs(title="Modified Entropy Measures & Indices", x="Metric", y="Value")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        legend.position="none",
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )+
      coord_flip()
  })
  output$tableCompIndices <- renderDT({
    base <- baselineResults()
    mod <- modifiedResults()
    req(base,mod)
    df <- data.frame(
      Metric=c("S(Y)","S(Y|X)","I(Y|X)","S(X,Y)","S(X)",
               "J(X;Y)","I(X,Y)","NMI","CES","CSF","ADI","NMGI"),
      Baseline=c(base$overallEntropyY,base$totalCondEntropy,base$normTotalCondEntropy,
                 base$SXY,base$Sx,base$Jxy,base$IXY,base$NMI,base$CES,base$CSF,base$ADI,base$NMGI),
      Modified=c(mod$overallEntropyY,mod$totalCondEntropy,mod$normTotalCondEntropy,
                 mod$SXY,mod$Sx,mod$Jxy,mod$IXY,mod$NMI,mod$CES,mod$CSF,mod$ADI,mod$NMGI)
    )
    datatable(df,options=list(pageLength=12),rownames=FALSE) %>%
      formatRound(columns=c("Baseline","Modified"),digits=3)
  })
  output$plotCompIndices <- renderPlot({
    base <- baselineResults()
    mod <- modifiedResults()
    metrics <- c("S(Y)","S(Y|X)","I(Y|X)","S(X,Y)","S(X)","J(X;Y)","I(X,Y)","NMI","CES","CSF","ADI","NMGI")
    df <- data.frame(
      Metric=rep(metrics,2),
      Value=c(
        base$overallEntropyY, base$totalCondEntropy, base$normTotalCondEntropy, base$SXY, base$Sx,
        base$Jxy, base$IXY, base$NMI, base$CES, base$CSF, base$ADI, base$NMGI,
        mod$overallEntropyY, mod$totalCondEntropy, mod$normTotalCondEntropy, mod$SXY, mod$Sx,
        mod$Jxy, mod$IXY, mod$NMI, mod$CES, mod$CSF, mod$ADI, mod$NMGI
      ),
      Scenario=rep(c("Baseline","Modified"),each=length(metrics))
    )
    df$Metric <- factor(df$Metric, levels=metrics)
    ggplot(df,aes(x=Metric,y=Value,fill=Scenario))+
      geom_bar(stat="identity",position=position_dodge(width=0.7))+
      labs(title="Comparison of Baseline vs. Modified Indices",x="Metric",y="Value")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  
  output$tableAlpha <- renderDT({
    req(alphaDF())
    datatable(alphaDF(), editable=TRUE, rownames=FALSE, options=list(pageLength=10)) %>%
      formatRound(columns="Alpha",digits=3)
  })
  proxyAlpha <- dataTableProxy("tableAlpha")
  observeEvent(input$tableAlpha_cell_edit, {
    info <- input$tableAlpha_cell_edit
    i <- info$row
    j <- info$col
    v <- info$value
    dfA <- alphaDF()
    if(j==2){
      valNum <- suppressWarnings(as.numeric(v))
      if(!is.na(valNum) && valNum>=0 && valNum<=1){
        dfA[i,j] <- valNum
      } else {
        showNotification("Alpha must be numeric in [0,1]. Reverting.", type="error")
      }
    }
    alphaDF(dfA)
    replaceData(proxyAlpha, alphaDF(), resetPaging=FALSE)
  })
  observeEvent(input$applyUniversalAlpha, {
    dfA <- alphaDF()
    dfA$Alpha <- input$universalAlpha
    alphaDF(dfA)
    replaceData(proxyAlpha, alphaDF(), resetPaging=FALSE)
  })
  
  alphaVarResults <- reactiveVal(NULL)
  observeEvent(input$applyAlphaChanges, {
    req(scenario$dataMat, storedSBJw(), storedPrefParams(), storedBenefitCost(), alphaDF())
    dfA <- alphaDF()
    M <- nrow(scenario$dataMat)
    rownames(dfA) <- dfA$Criterion
    alphaVec <- numeric(M)
    for(i in 1:M){
      cr <- rownames(scenario$dataMat)[i]
      if(cr %in% dfA$Criterion) alphaVec[i] <- as.numeric(dfA[cr,"Alpha"]) else alphaVec[i] <- 0.5
    }
    res <- computeESMADMIII(
      dataMat=scenario$dataMat,
      deltaMat=storedDeltaMat(),
      sbjDF=scenario$SBJw,
      prefParams=scenario$PrefParams,
      benefitCost=scenario$BenefitCost,
      alphaVec=alphaVec
    )
    alphaVarResults(res)
    showNotification("Alpha Variation scenario computed!", type="message")
  })
  output$tableAlphaComp <- renderDT({
    base <- baselineResults()
    mod <- alphaVarResults()
    req(mod)
    df <- data.frame(
      Metric=c("S(Y)","S(Y|X)","I(Y|X)","S(X,Y)","S(X)","NMI","CES","CSF","ADI","NMGI"),
      Baseline=c(base$overallEntropyY, base$totalCondEntropy, base$normTotalCondEntropy,
                 base$SXY, base$Sx, base$NMI, base$CES, base$CSF, base$ADI, base$NMGI),
      Modified=c(mod$overallEntropyY, mod$totalCondEntropy, mod$normTotalCondEntropy,
                 mod$SXY, mod$Sx, mod$NMI, mod$CES, mod$CSF, mod$ADI, mod$NMGI)
    )
    datatable(df,options=list(pageLength=10),rownames=FALSE) %>%
      formatRound(columns=c("Baseline","Modified"),digits=3)
  })
  output$plotAlphaNMGI <- renderPlot({
    req(alphaVarResults())
    base <- baselineResults()
    mod <- alphaVarResults()
    df <- data.frame(
      Scenario=c("Baseline","Alpha Variation"),
      NMGI=c(base$NMGI, mod$NMGI)
    )
    ggplot(df, aes(x=Scenario, y=NMGI, fill=Scenario))+
      geom_bar(stat="identity", width=0.4)+
      geom_text(aes(label=format3(NMGI)), vjust=-0.3, color="white", size=5)+
      labs(title="NMGI Comparison: Baseline vs. Alpha Variation", x="", y="NMGI")+
      theme_minimal()+
      theme(
        plot.title=element_text(color="#FFFF00",size=18,face="bold"),
        axis.title=element_text(color="#FFFFFF",size=16),
        axis.text=element_text(color="#FFFFFF",size=14),
        legend.position="none",
        panel.background=element_rect(fill="#222222"),
        plot.background=element_rect(fill="#222222")
      )
  })
  
  output$summaryResults <- renderUI({
    r <- baselineResults(); req(r)
    altNames <- colnames(storedDataMat())
    dfAlt <- data.frame(Alternative=altNames, Score=r$altScores)
    dfAlt <- dfAlt[order(-dfAlt$Score),]
    topAlt <- dfAlt$Alternative[1]
    topScore <- dfAlt$Score[1]
    critNames <- rownames(storedDataMat())
    dfICI <- data.frame(Criterion=critNames, ICI=r$critImport)
    dfICI <- dfICI[order(-dfICI$ICI),]
    topCrit <- dfICI$Criterion[1]
    topVal <- dfICI$ICI[1]
    
    insights <- generateDetailedEntropyInsights(
      SXY=r$SXY,
      SX=r$Sx,
      SY=r$overallEntropyY,
      IXY=r$IXY,
      SYX=r$normTotalCondEntropy,
      IYX=r$Jxy,
      NMI=r$NMI,
      ADI=r$ADI,
      CES=r$CES,
      CSF=r$CSF,
      NMGI=r$NMGI
    )
    
    HTML(paste0(
      "<table style='width:100%; border-collapse:collapse; border:1px solid black;'>",
      "<tr style='background-color:#ffff00;'>",
      "<th colspan='2' style='text-align:left; padding:8px; border:1px solid black; color:black; font-weight:bold;'>Metric</th></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'><strong>Optimal Alternative</strong></td>",
      "<td style='padding:8px; border:1px solid black;'>", topAlt, " (Score: <strong>", format3(topScore), "</strong>)</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'><strong>Most Significant Criterion</strong></td>",
      "<td style='padding:8px; border:1px solid black;'>", topCrit, " (ICI: <strong>", format3(topVal), "</strong>)</td></tr>",
      "<tr style='background-color:#ffff00;'>",
      "<th colspan='2' style='text-align:left; padding:8px; border:1px solid black; color:black; font-weight:bold;'>Key Indices</th></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>S(X,Y) - Joint Entropy</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$SXY), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>S(X) - Criteria Entropy</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$Sx), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>S(Y) - Alternatives Entropy</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$overallEntropyY), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>I(X,Y) - Normalized Joint Entropy</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$IXY), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>I(Y|X) - Conditional Entropy</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$normTotalCondEntropy), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>I(X;Y) - Mutual Information</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$Jxy), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>Normalized Mutual Information (NMI)</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$NMI), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>Alternatives Distinction Index (ADI)</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$ADI), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>Criteria Effectiveness Score (CES)</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$CES), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>Conditional Stability Factor (CSF)</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$CSF), "</td></tr>",
      "<tr><td style='padding:8px; border:1px solid black;'>Net Mutual Growth Index (NMGI)</td>",
      "<td style='padding:8px; border:1px solid black;'>", format3(r$NMGI), "</td></tr>",
      "<tr style='background-color:#ffff00;'>",
      "<th colspan='2' style='text-align:left; padding:8px; border:1px solid black; color:black; font-weight:bold;'>Detailed Insights</th></tr>",
      "<tr><td colspan='2' style='padding:8px; border:1px solid black;'>", insights, "</td></tr>",
      "</table>"
    ))
  })
  
  observe({
    fs <- input$fontSize
    runjs(sprintf("$('body').css('font-size','%spx');", fs))
  })
  
  output$downloadAllSummary <- downloadHandler(
    filename=function(){"ESMADMIII_Summary_Results.xlsx"},
    content=function(file){
      r <- baselineResults(); req(r)
      critNames <- rownames(storedDataMat())
      df_obj <- data.frame(
        Criterion=critNames,
        Objective_Weight=r$xOBJ
      )
      sbj <- storedSBJw()
      subj_central <- sapply(critNames, function(cn){
        idx <- which(sbj$Criterion==cn)
        if(length(idx)==0) NA_real_ else sbj$Central[idx[1]]
      })
      subj_delta <- sapply(critNames, function(cn){
        idx <- which(sbj$Criterion==cn)
        if(length(idx)==0) NA_real_ else sbj$Delta[idx[1]]
      })
      df_int <- data.frame(
        Criterion=critNames,
        Subjective_Weight_Central=subj_central,
        Subjective_Weight_Delta=subj_delta,
        Objective_Weight=r$xOBJ,
        Integrated_Weight=r$xINT,
        Integrated_Criteria_Importance=r$critImport
      )
      altNames <- colnames(storedDataMat())
      df_alt <- data.frame(Alternative=altNames, Score=r$altScores)
      df_alt <- df_alt[order(-df_alt$Score),]
      df_indices <- data.frame(
        S_Y=r$overallEntropyY,
        Total_Cond_Entropy=r$totalCondEntropy,
        Norm_Cond_Entropy=r$normTotalCondEntropy,
        S_X=r$Sx,
        S_XY=r$SXY,
        J_XY=r$Jxy,
        I_XY=r$IXY,
        NMI=r$NMI,
        CES=r$CES,
        CSF=r$CSF,
        ADI=r$ADI,
        NMGI=r$NMGI
      )
      writexl::write_xlsx(list(
        Criteria_Objective=df_obj,
        Criteria_Integrated=df_int,
        Alternatives=df_alt,
        Indices_Entropy=df_indices
      ), path=file)
    }
  )
}

########################################
# 5) Run the App
########################################

shinyApp(ui=ui, server=server)
