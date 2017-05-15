library(shiny)

shinyUI(pageWithSidebar(
  headerPanel('Estimation Modell by Shiny'),
  
  sidebarPanel(
    fileInput('file', 'Choose A CSV File', 
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    numericInput('obs', 'Number of observations to view', 10),
    conditionalPanel(
      condition = 'output.element || output.plot1',
      textAreaInput('exp_input', 'Enter the expression of system', resize='vertical'),
      helpText('Parallelsystem is with P() evaluated.',
               'Seriensystem is with S() evaluated.',
               'For example: P(1, 2, 3, S(4, 5))'),
    
      radioButtons('distribution', 'Choose a distribution function',
                   c(weibull='w',
                     exp='e',
                     lognormal='l'),
                   'w'),
      numericInput('par1', 'Set the value of the first parameter', NULL),
      numericInput('par2', 'Set the value of the second parameter', NULL),
      numericInput('sys_t', 'Set a timepoint of observation', 0),
      radioButtons('selectedGraph', 'Choose a graph',
                   c(static='s',
                     dynamic='d'),
                   's')
    ),
    checkboxInput('hd', 'header', TRUE),
    radioButtons('sep', 'Separator',
                c(Comma=',',
                  Semicolon=';',
                  Tab='\t')),
    uiOutput('selectableData'),
    submitButton('Confirm')
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel('Table', tableOutput('table')),
      tabPanel('Summary', verbatimTextOutput('summary')),
      tabPanel('Parameters Estimation',
               h4('Weibull-distribution with two parameters'),
               tableOutput('weib2.pm.table'),
               
               h4('Weibull-distribution with three parameters'),
               tableOutput('weib3.pm.table'),
               
               h4('Mixed weibull-distribution'),
               tableOutput('mixedWeib.pm.table'),
               
               h4('Exponential distribution'),
               tableOutput('exp.pm.table'),
               
               h4('Lognormal distribution'),
               tableOutput('logNormal.pm.table'),
               
               h4('Gumbel distribution'),
               tableOutput('gumbel.pm.table'),
               
               h4('Gamma distribution'),
               tableOutput('gamma.pm.table')
               ),
      
      
      #tabPanel('Dynamic Graph',
      #         uiOutput('html')),

      tabPanel('Plot', plotOutput('plot1'), 
               plotOutput('plot2'), 
               plotOutput('plot3'),
               plotOutput('plot4')
               ),
      tabPanel('Avalibility estimation and graph of systems', 
               tableOutput('element'),
               
               h4('Expression:'),
               textOutput('exp'),
               
               h4('Density Curve:'),
               plotOutput('dplot'),
               
               h4('Result with entered t:'),
               tableOutput('expTb'),
               
               h4('The diagram of system'),
               uiOutput('d_graph'),
               plotOutput('s_graph')
      )
    )
  )
))