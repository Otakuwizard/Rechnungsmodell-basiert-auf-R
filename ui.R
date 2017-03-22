library(shiny)

shinyUI(pageWithSidebar(
  headerPanel('ML Schaetzung der Parameter unter verschiedenen Verteilungsmodelen'),
  
  sidebarPanel(
    fileInput('file', 'Choose CSV File', 
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    conditionalPanel(
      condition = 'output.element || output.plot1',
      textAreaInput('exp_input', 'Enter the expression of system', resize='vertical'),
      helpText('Parallelsystem is with P() evaluated.',
               'Seriensystem is with S() evaluated.',
               'For example: P(1, 2, 3, S(4, 5))')
    ),
    numericInput('sys_t', 'Set a timepoint of observation', 0),
    numericInput('obs', 'Number of observations to view', 10),
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
      tabPanel('Parameter Estimate',
               h4('Zwei parameterische Weibull-Verteilung'),
               tableOutput('weib2.pm.table'),
               conditionalPanel(
                 condition = 'input.sys_t != 0',
                 h4('Unter t zensierte Daten:'),
                 tableOutput('z_weib2.pm.table')
               ),
               h4('Drei parameterische Weibull-Verteilung'),
               tableOutput('weib3.pm.table'),
               conditionalPanel(
                 condition = 'input.sys_t != 0',
                 h4('Unter t zensierte Daten:'),
                 tableOutput('z_weib3.pm.table')
               ),
               
               h4('Gemischt Weibull-Verteilung'),
               tableOutput('mixedWeib.pm.table'),
               conditionalPanel(
                 condition = 'input.sys_t != 0',
                 h4('Unter t zensierte Daten:'),
                 tableOutput('z_mixedweib.pm.table')
               ),
               
               h4('Exponential Verteilung'),
               tableOutput('exp.pm.table'),
               
               h4('Log-normal Verteilung'),
               tableOutput('logNormal.pm.table'),
               
               h4('Gumbel Verteilung'),
               tableOutput('gumbel.pm.table'),
               
               h4('Gamma Verteilung'),
               tableOutput('gamma.pm.table')
               ),
      tabPanel('Schaetzung der Verfuegbatkeit von Systemen', 
               tableOutput('element'),
               
               h4('Expression:'),
               textOutput('exp'),
               
               h4('Density Expression:'),
               textOutput('dexp'),
               
               h4('Result with entered t:'),
               tableOutput('expTb')),
      
      tabPanel('Summary', verbatimTextOutput('summary')),
      tabPanel('Plot', plotOutput('plot1'), 
               plotOutput('plot2'), 
               conditionalPanel(
                condition = 'input.sys_t != 0',          
                plotOutput('plot3')
                )
               )
    )
  )
))