# Libraries ####
libs =c('shiny','parallel','rstan',
        'inlmisc','shinycssloaders')
for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}

navbarPage(
  title = 'FitOCT',
  # Config ####
  theme = shinythemes::shinytheme(c("cosmo","cerulean",
                                    "spacelab","yeti")[2]),
  tabPanel(
    # Noise + Exp ####
    title = 'Noise + Exp.',
    sidebarLayout(
      sidebarPanel(
        fileInput(
          inputId = 'dataFile',
          label   = 'Select data file',
          multiple= FALSE,
          accept  = c('.txt','.csv')
        ),
        hr(style="border-color: #666;"),
        sliderInput(
          inputId = 'depthSel',
          label   = 'Select depth range',
          min = 0,
          max = 1,
          value = c(0,1),
          sep = ''
        ),
        sliderInput(
          inputId = 'subSample',
          label   = 'Sub-sample factor',
          min   = 1,
          max   = 10,
          value = 1,
          step  = 1,
          sep = ''
        ),
        hr(style="border-color: #666;"),
        sliderInput(
          inputId = 'smooth_df',
          label   = 'Smoothing level',
          min     = 2,
          max     = 30,
          value   = 15
        ),
        hr(style='border-color: #666;'),
        h5('Fit results'),
        verbatimTextOutput(
          'resMonoExp'
        )
      ),
      mainPanel(
        wellPanel(
          h4("Noise estimation: uy(x) = a_1 * exp( -2*depth / a_2 )"),
          plotOutput(
            outputId = 'plotNoise',
            height   = '300px'
          )
        ),
        wellPanel(
          h4("Monoexponential fit: y(x) = b_1 + b_2  *exp( -2*depth / b_3 )"),
          plotOutput(
            outputId = 'plotMonoExp',
            height = '300px'
          )
        )
      )
    )
  ),
  tabPanel(
    # Modulation ####
    title = 'Modulation',
    sidebarLayout(
      sidebarPanel(
        h4('Bayesian Inference'),
        wellPanel(
          fluidRow(
            column(
              width = 7,
              radioButtons(
                inputId = 'method',
                label   = 'Method',
                choices = c(M.A.P. = 'optim',
                            Sample = 'sample'),
                selected = 'optim',
                inline = FALSE
              )
            ),
            column(
              width = 3,
              actionButton(
                inputId = 'runExpGP',
                strong("Run"),
                icon=icon('gear')
              ),
              tags$style(type='text/css',
                         "#runExpGP { width:200%; margin-top: 0px;}")
            )
          ),
          conditionalPanel(
            condition = "input.method == 'sample'",
            fluidRow(
              column(
                width = 5,
                numericInput(
                  "nb_warmup",
                  label='nb_warmup',
                  value = 500,
                  min = 100,
                  max = 5000,
                  step = 100,
                  width = "100px"
                )
              ),
              column(
                width = 5,
                numericInput(
                  "nb_sample",
                  label='nb_sample',
                  value = 500,
                  min = 100,
                  max = 10000,
                  step = 100,
                  width = "100px"
                )
              )
            )
          )
        ),
        hr( style="border-color: #666;"),
        h4('Controls'),
        tabsetPanel(
          tabPanel(
            title = 'Prior',
            wellPanel(
              fluidRow(
                column(
                  width = 6,
                  numericInput(
                    "ru_theta",
                    label='ru_theta',
                    value = 0.05,
                    min = 0.01,
                    max = 0.5,
                    step = 0.01,
                    width = "100px"
                  )
                ),
                column(
                  width = 6,
                  actionButton(
                    inputId = 'runPriExpGP',
                    strong("Simul."),
                    icon=icon('gear')
                  ),
                  tags$style(type='text/css',
                             "#runPriExpGP { width:100%; margin-top: 25px;}"
                  )
                )
              ),
              fluidRow(
                column(
                  width = 10,
                  sliderInput(
                    'lambda_rate',
                    label='lambda_rate',
                    value = 0.1,
                    min   = 0.0,
                    max   = 0.5,
                    step  = 0.01
                  )
                )
              )
            )
          ),
          tabPanel(
            title = 'GP',
            wellPanel(
              fluidRow(
                column(
                  width = 6,
                  numericInput(
                    "Nn",
                    label='Nb. ctrl',
                    value = 10,
                    min = 5,
                    max = 20,
                    step = 1
                  )
                ),
                column(
                  width = 6,
                  radioButtons(
                    inputId = 'gridType',
                    label   = 'Grid Type',
                    choices = c(Extremal = 'extremal',
                                Internal = 'internal'),
                    selected = 'internal',
                    inline = FALSE
                  )
                )
              ),
              fluidRow(
                column(
                  width = 10,
                  sliderInput(
                    'rho_scale',
                    label='rho',
                    value = 0.0,
                    min   = 0.0,
                    max   = 1.0,
                    step  = 0.01
                  )
                ),
                column(
                  width = 6
                )
              )
            )
          ),
          tabPanel(
            title = 'Graphics',
            wellPanel(
              fluidRow(
                column(
                  width = 12,
                  sliderInput(
                    "modRange",
                    "Scale Modulation range",
                    min   = 0.1,
                    max   = 5,
                    value = 0.3,
                    step  = 0.1
                  )
                )
              )
            )
          )
        ),
        hr( style="border-color: #666;"),
        h4("Fit results"),
        verbatimTextOutput(
          'resExpGP'
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            title = 'Posterior',
            wellPanel(
              tabsetPanel(
                tabPanel(
                  title = 'Predictions',
                  wellPanel(
                    withSpinner(
                      plotOutput(
                        'plotExpGP',
                        height = '600px'
                      ),
                      type = 4
                    )
                  )
                ),
                tabPanel(
                  title = 'Statistics',
                  wellPanel(
                    DT::dataTableOutput(
                      'summaryOut'
                    )
                  )
                ),
                tabPanel(
                  title = 'Traces',
                  wellPanel(
                    withSpinner(
                      plotOutput(
                        'tracesExpGP',
                        height = '600px'
                      ),
                      type = 4
                    )
                  )
                ),
                tabPanel(
                  title = 'Output',
                  withSpinner(
                    verbatimTextOutput(
                      'outExpGP'
                    ),
                    type = 4
                  )
                )
              )
            )
          ),
          tabPanel(
            title = 'Prior',
            wellPanel(
              tabsetPanel(
                tabPanel(
                  title = 'Predictions',
                  wellPanel(
                    withSpinner(
                      plotOutput(
                        'plotPriExpGP',
                        height = '300px'
                      ),
                      type = 4
                    )
                  )
                ),
                tabPanel(
                  title = 'Statistics',
                  wellPanel(
                    DT::dataTableOutput(
                      'summaryPriOut'
                    )
                  )
                ),
                tabPanel(
                  title = 'Traces',
                  wellPanel(
                    withSpinner(
                      plotOutput(
                        'tracesPriExpGP',
                        height = '600px'
                      ),
                      type = 4
                    )
                  )
                )
              )
            )
          ),
          tabPanel(
            title = 'Identification',
            wellPanel(
              tabsetPanel(
                tabPanel(
                  title = 'Pri/Post',
                  wellPanel(
                    withSpinner(
                      plotOutput(
                        'priPostExpGP',
                        height = '600px'
                      ),
                      type = 4
                    )
                  )
                )
              )
            )
          ),
          type='pills'
        )
      )
    )
  ),
  tabPanel(
    # Test GP ####
    title = "Test GP",
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(
            width = 5,
            numericInput(
              inputId = 'NnTest',
              label='Nb. ctrl',
              value = 10,
              min = 5,
              max = 20,
              step = 1
            ),
            radioButtons(
              inputId = 'gridTypeTest',
              label   = 'Grid Type',
              choices = c(Extremal = 'extremal',
                          Internal = 'internal'),
              selected = 'internal',
              inline = FALSE
            )
          ),
          column(
            width = 5,
            numericInput(
              'alpha_scaleTest',
              label='alpha_scale',
              value = 0.1,
              min = 0.01,
              max = 1.0,
              step = 0.05
            ),
            sliderInput(
              'rho_scaleTest',
              label='rho_scale',
              value = 0.0,
              min   = 0.0,
              max   = 1.0,
              step  = 0.01
            )
          )
        )
      ),
      mainPanel(
        wellPanel(
          h4('Simulate GP for test modulation curve'),
          plotOutput(
            'plotGP',
            height = '600px'
          )
        )
      )
    )
  ),
  # About ####
  tabPanel(
    title="Save",
    sidebarPanel(
      wellPanel(
        h4('Generate report'),
        fluidRow(
          column(
            width = 6,
            downloadButton(
              outputId = 'report',
              label    = 'Download  (Ctrl+Click)'
            )
          )
        ),
        h4('Save parameters'),
        fluidRow(
          column(
            width = 6,
            downloadButton(
              outputId = 'params',
              label    = 'Download  (Ctrl+Click)'
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

