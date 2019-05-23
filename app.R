library(shiny)
library(shiny.semantic)
library(shinydashboard)

# Define UI for miles per gallon app ----
#ui <- semanticPage(
# ui <- fluidPage(
col1="grey"
col2="teal"

main_title=h1("Prior explorer for coloc", class="ui header",color=col2)

  condpanel_marg=conditionalPanel(
                       condition = "input.priors == 'Marginal/Conditional'",
                       title = "Marginal/Conditional priors",
                       numericInput("nsnp", label="Number of SNPs in region",
                                    value = 1000),br(),
                       numericInput("q12", 
                                    "q1|2",
                                    value = 1e-2),br(),
                       numericInput("q1", 
                                    "q1", 
                                    value = 1e-4),br(),
                       numericInput("q2", 
                                    "q2", 
                                    value = 1e-4),br(),
                       helpText("nsnps > 0",br(),"0 < q1|2 < 1",br(),"0 < q1 < 1/nsnps",br(),"0 < q2 < 1/nsnps"))
  
  condpanel_raw= conditionalPanel(
                       condition = "input.priors == 'Raw'",
                       title = "Raw priors",
                       numericInput("nsnp", 
                                    "Number of SNPs in region", 
                                    value = 1000),br(),
                       numericInput("p12", 
                                    "p12", 
                                    value = 1e-6),br(),
                       numericInput("p1", 
                                    "p1", 
                                    value = 1e-4),br(),
                       numericInput("p2", 
                                    "p2", 
                                    value = 1e-4),br(),
                       helpText("nsnps > 0",br(),"p1 * p2 < p12 < min(p1,p2)",br(),"0 < p1 < 1/nsnps",br(),"0 < p2 < 1/nsnps") )
  
ui <- dashboardPage(
  
  # App title ----
  title = "Prior explorer for coloc",
  # suppressDependencies("bootstrap"),
  dashboardHeader(color = col2,disable=TRUE,
                  h4("Prior explorer for coloc",class="ui center aligned basic segment")),
  
  
  dashboardSidebar(side = "left", size = "thin", color = col1, #inverted=TRUE,
                   sidebarMenu(
                     menuItem(tabName = "output", "App"),
                     menuItem(tabName = "docs", "Documentation"))
  ),
  
  # Main panel for displaying outputs ----
  # mainPanel(
  dashboardBody(tabItems(
    tabItem(tabName = "output", 
            
            div(class="ui basic segment",
                
                main_title,
                
                div(class="ui two column doubling stackable grid container",
                    div(class="column",
                        # Inputs ----
                        # Type of input ----  
                        box(title="Input: Type of Per SNP priors",collapsible=FALSE,color=col1,title_side = "top left",
                            selectInput(inputId="priors", choices=c("Raw", "Marginal/Conditional"), label="select one",selected="Raw",selectize=FALSE,size=2)),
                        box(title="Input: Parameter values",collapsible=FALSE,color=col1,title_side = "top left",
                            sidebarPanel(
                              condpanel_raw, 
                              condpanel_marg
                    ))),
                    div(class="column",
                        # div(class="ui basic segment",
                        box(title="Result",collapsible = FALSE,color = col2,
                            # Output: Tables ----
                            h4("Per-SNP priors",  class="ui horizontal divider header"),
                            tableOutput("priors"),
                            
                            h4("Per-hypothesis priors",  class="ui horizontal divider header"),
                            tableOutput("table"),
                            
                            # Output: Barplot ----
                            plotOutput(outputId = "plot")))))),
            
    tabItem(tabName = "docs", 
            
                main_title,
            
            # titlePanel("Prior explorer for coloc"),
            div(p(a("Coloc",href="http://chr1swallace.github.io/coloc"),
              "can be used to decide whether two traits share a common causal genetic variant in an LD-defined genomic region.",
              "As a Bayesian method, it requires the user think carefully about their prior beliefs before running the analysis.",
              "This app is designed to convert the per-SNP raw priors described in", a("Giambartolomei et al 2014",href="https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383"), 
              "or the per SNP marginal and conditional priors described in [manuscript in progress], to per-hypothesis priors."),
            
            p("This should help investigators better understand the values they supply to coloc, and ensure they are consistent with prior beliefs across a range of levels.",
              "The different parameters which can be input (raw or marginal/conditiona per SNP priors) and their interpretation in terms of per hypothesis priors are described below."),
           
            div(class="ui basic segment", 
              tags$h5("Per SNP priors - raw", class="ui horizontal divider header"),
              tags$ul(
                tags$li("p12 = Prior probability a random SNP in the region is jointly causal for both traits"),
                tags$li("p1 = Prior probability a random SNP in the region is causally associated with trait 1 and not trait 2"),
                tags$li("p2 = Prior probability a random SNP in the region is causally associated with trait 2 and not trait 1")),
            
            tags$h5("Per SNP priors - marginal and conditional", class="ui horizontal divider header"),
            tags$ul(
              tags$li( "q1|2 = Prior probability a random SNP in the region is jointly causal for both traits given it is causal for trait 2"),
              tags$li("q1 = Prior probability a random SNP in the region is causally associated with trait 1 (whether or not causal for trait 2)"),
              tags$li("q2 = Prior probability a random SNP in the region is causally associated with trait 2 (whether or not causal for trait 1)")),
           
            tags$h5("Per-hypothesis priors",class="ui horizontal divider header"),
            tags$ul(
              tags$li("H0 = no association with either trait in region"),
              tags$li("H1 = association with trait 1 in region, but not trait 2"),
              tags$li("H2 = association with trait 2 in region, but not trait 1"),
              tags$li("H3 = association with both traits in region, but seperate causal variants"),
              tags$li("H4 = association with both traits in region, same causal variants"))),
              
            div(class="ui right aligned basic segment",
                tags$i("View the source code on",
                a("github",href="http://github.com/chr1swallace/coloc-priors"))  ))
    ))))
    
# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  
  r2c <- function(nsnp,p12,p1,p2) {
    list(nsnp=nsnp,
         q12=p12/(p2+p12),
         q1=p12+p1,
         q2=p12+p2)
  }
  c2r <- function(nsnp,q12,q1,q2) {
    p12=q12 * q2
    list(nsnp=nsnp,
         p12=p12,
         p1=q1-p12,
         p2=q2-p12)
  }
  prior.raw <- function(nsnp,p12,p1,p2) {
    if(p12<p1*p2 || p12 > p1 || p12 > p2)
      return(NULL)
    tmp <- c(nsnp * p1,
             nsnp * p2,
             nsnp * (nsnp-1) * p1 * p2,
             nsnp * p12)
    tmp <- c(1-sum(tmp),tmp)
    names(tmp) <- paste0("H",0:4)
    tmp
  }
  prior.convert <- function(input) {
    # if(input$priors=="Raw")
    prior.raw(input$nsnp,input$p12,input$p1,input$p2)
    # else {
    #   newinput=c2r(input$nsnp,input$q12,input$q1,input$q2)
    #   prior.raw(newinput$nsnp,newinput$p12,newinput$p1,newinput$p2)
    # }
  }
  
  isValid <- reactive({
    validate(
      need(vals$p1 > 0 & vals$p1 < 1/vals$nsnp, "single trait prob p1 must be between 0 and 1/nsnp"),
      need(vals$q1 > 0 & vals$q1 < 1/vals$nsnp, "single trait prob q1 must be between 0 and 1/nsnp"),
      need(vals$q2 > 0 & vals$q2 < 1/vals$nsnp, "single trait prob q2 must be between 0 and 1/nsnp"),
      need(vals$p2 > 0 & vals$p2 < 1/vals$nsnp, "single trait prob p2 must be between 0 and 1/nsnp"),
      need(vals$nsnp > 1, "nsnp must be > 1"),
      need(vals$p12>=vals$p1 * vals$p2 &&  vals$p12 <= min(vals$p1,vals$p2), "p12 must be between p1*p2 and min(p1,p2)"),
      need((vals$p1+vals$p2+vals$p12)*vals$nsnp + vals$p1*vals$p2*vals$nsnp*(vals$nsnp-1) < 1, "sum of P(H_i), i=1..4, is greater than 1"))
  })
  
  
  isValidNomsg <- reactive({
    req(vals$p1 > 0 & vals$p1 < 1/vals$nsnp)
    req(vals$q1 > 0 & vals$q1 < 1/vals$nsnp)
    req(vals$q2 > 0 & vals$q2 < 1/vals$nsnp)
    req(vals$p2 > 0 & vals$p2 < 1/vals$nsnp)
    req(vals$p12>=vals$p1 * vals$p2 &&  vals$p12 <= min(vals$p1,vals$p2))
    req((vals$p1+vals$p2+vals$p12)*vals$nsnp + vals$p1*vals$p2*vals$nsnp*(vals$nsnp-1) < 1)
  })
  
  
  vals <- reactiveValues()
  observe({
    if(input$priors=="Raw") {
      newinput=r2c(input$nsnp,input$p12,input$p1,input$p2)
      vals$nsnp=input$nsnp
      vals$p12=input$p12
      vals$p1=input$p1
      vals$p2=input$p2
      vals$q12=newinput$q12
      vals$q1=newinput$q1
      vals$q2=newinput$q2
    } else {
      newinput=c2r(input$nsnp,input$q12,input$q1,input$q2)
      vals$nsnp=input$nsnp
      vals$p12=newinput$p12
      vals$p1=newinput$p1
      vals$p2=newinput$p2
      vals$q12=input$q12
      vals$q1=input$q1
      vals$q2=input$q2
    }
    # output$valid = isValid()
  })
  
  # vals <- reactiveValues(getvals())
  
  output$plot = renderPlot({
    isValidNomsg()
    hyp.prior <- prior.convert(vals)
    barplot(hyp.prior,names.arg=names(hyp.prior))
  })
  
  output$priors <- renderTable({
    isValid()
    cp=c("q1|2"=vals$q12,q1=vals$q1,q2=vals$q2)
    gp=c(p12=vals$p12,p1=vals$p1,p2=vals$p2)
    t(c(gp,cp))
  },digits=-1)
  
  # output$priors <- renderTable({
  #   if(input$priors=="Raw") {
  #     newinput=r2c(input$nsnp,input$p12,input$p1,input$p2)
  #     cp=c("q1|2"=newinput$q12,q1=newinput$q1,q2=newinput$q2)
  #     gp=c(p12=input$p12,p1=input$p1,p2=input$p2)
  #   } else {
  #     cp=c("q1|2"=input$q12,q1=input$q1,q2=input$q2)
  #     newinput=c2r(input$nsnp,input$q12,input$q1,input$q2)
  #     gp=c(p12=newinput$p12,p1=newinput$p1,p2=newinput$p2)
  #   }
  #   t(c(gp,cp))
  # },digits=-1)
  
  output$table <- renderTable({
    isValidNomsg()
    hyp.prior <- prior.convert(vals)
    t(hyp.prior)
  },digits=3)
}

shinyApp(ui, server)

