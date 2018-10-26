library(shiny)
library(markdown)

navbarPage("FitOCT",
           tabPanel("Noise",
                    sidebarLayout(
                      sidebarPanel(
                        sliderInput(inputId = 'set_df',
                                    label   = 'df',
                                    min     = 2,
                                    max     = 30,
                                    value   = 15)
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
           tabPanel("Summary",
                    verbatimTextOutput("summary")
           )
)

