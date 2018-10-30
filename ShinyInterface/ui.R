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
           theme = shinythemes::shinytheme(c("cosmo","cerulean",
                                             "spacelab","yeti")[2]),
           tabPanel("Noise + Exp.",
                    sidebarLayout(
                      sidebarPanel(
                        fileInput(
                          inputId = 'dataFile',
                          label   = 'Select data file',
                          multiple= FALSE,
                          accept  = c('.txt','.csv')
                        ),
                        hr( style="border-color: #666;"),
                        sliderInput("depthSel",
                                    "Select depth range",
                                    min = 0,
                                    max = 1,
                                    value = c(0,1),
                                    sep=""
                        ),
                        hr( style="border-color: #666;"),
                        sliderInput(inputId = 'smooth_df',
                                    label   = 'Smoothing level',
                                    min     = 2,
                                    max     = 30,
                                    value   = 15),
                        hr( style="border-color: #666;"),
                        h5("Fit results"),
                        verbatimTextOutput("resMonoExp")
                      ),
                      mainPanel(
                        wellPanel(
                          h4("Noise estimation: uy(x) = a_1 * exp( -depth / a_2 )"),
                          plotOutput("plotNoise",height='300px')
                        ),
                        wellPanel(
                          h4("Monoexponential fit: y(x) = b_1 + b_2  *exp( -depth / b_3 )"),
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
                        conditionalPanel(
                          condition = "input.method == 'sample'",
                          fluidRow(
                            column(5,
                                   numericInput("nb_warmup",
                                                label='nb_warmup',
                                                value = 100,
                                                min = 100,
                                                max = 1000,
                                                step = 100,
                                                width = "100px")
                            ),
                            column(5,
                                   numericInput("nb_sample",
                                                label='nb_sample',
                                                value = 100,
                                                min = 100,
                                                max = 10000,
                                                step = 100,
                                                width = "100px")
                            )
                          )
                        ),
                        hr( style="border-color: #666;"),
                        h4("Constraints"),
                        fluidRow(
                          column(5,
                                 numericInput("Nn",
                                              label='Nb. ctrl',
                                              value = 10,
                                              min = 5,
                                              max = 20,
                                              step = 1)
                          ),
                          column(5,
                                 numericInput("ru_theta",
                                              label='ru_theta',
                                              value = 0.1,
                                              min = 0.01,
                                              max = 0.5,
                                              step = 0.01)
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

