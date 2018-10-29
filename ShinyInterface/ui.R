# Libraries ####
libs =c('shiny','parallel','rstan','FitOCTLib',
        'inlmisc','shinycssloaders')
for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}

navbarPage("FitOCT",
           # Config ####
           theme = shinythemes::shinytheme(c("cosmo","cerulean","spacelab","yeti")[2]),
           tabPanel("Noise",
                    sidebarLayout(
                      sidebarPanel(
                        fileInput(
                          inputId = 'dataFile',
                          label   = 'Select data file',
                          multiple= FALSE,
                          accept  = c('.txt','.csv')
                        ),
                        hr( style="border-color: #666;"),
                        sliderInput(inputId = 'smooth_df',
                                    label   = 'Smoothing level',
                                    min     = 2,
                                    max     = 30,
                                    value   = 15),
                        hr( style="border-color: #666;"),
                        verbatimTextOutput("resMonoExp")
                      ),
                      mainPanel(
                        wellPanel(
                          h4("Noise estimation"),
                          plotOutput("plotNoise",height='300px')
                        ),
                        wellPanel(
                          h4("Monoexponential fit"),
                          plotOutput("plotMonoExp",height='300px')
                        )
                      )
                    )
           ),
           tabPanel("Modulation",
                    sidebarLayout(
                      sidebarPanel(
                        fluidRow(
                          column(7,
                                 radioButtons(
                                   inputId = 'method',
                                   label   = 'Bayesian Inference',
                                   choices = c(M.A.P. = 'optim',
                                               Sample = 'sample'),
                                   selected = 'optim',
                                   inline = FALSE
                                 )
                          ),
                          column(3,
                                 actionButton("runExpGP",
                                              strong("Run"),
                                              icon=icon('gear')
                                 ),
                                 tags$style(type='text/css',
                                            "#runALS { width:100%; margin-top: 25px;}")
                          )
                        ),
                        hr( style="border-color: #666;"),
                        verbatimTextOutput("resExpGP")
                      ),
                      mainPanel(
                        wellPanel(
                          h4("Modulated exponential fit"),
                          tabsetPanel(
                            tabPanel('Fit',
                                     withSpinner(
                                       plotOutput("plotExpGP",height='600px')
                                       ,type = 4
                                     )
                            ),
                            tabPanel(
                              "Statistics",
                              DT::dataTableOutput("summaryOut")
                            ),
                            tabPanel('Traces',
                                    withSpinner(
                                      plotOutput("tracesExpGP",height='600px')
                                      ,type = 4
                                    )
                            )
                          )
                        )
                      )
                    )
           ),
           # About ####
           tabPanel(
             title="About",
             sidebarPanel(
               h5("Author      : P. Pernot"),
               h5("Affiliation : CNRS"),
               h5("Version     : 0.1"),
               h5("Date        : 2018/10/28"),
               hr( style="border-color: #666;"),
               a(href="https://github.com/ppernot/FitOCT","How to cite..."),
               br(),
               a(href="https://github.com/ppernot/FitOCT","code@github"),
               br(),
               a(href="https://github.com/ppernot/FitOCT/issues","Bugs report, Features request")
             )
           )

)

