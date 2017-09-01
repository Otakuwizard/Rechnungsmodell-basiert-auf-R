library(shiny)
#library(Ryacas)
library(Deriv)
library(datasets)
library(car)
library(lattice)
library(FAdist)
library(Renext)
source('parameters_ml.R')
source('system.R')


shinyServer(function(input, output){
  datasetInput = reactive({
    inFile = input$file
    if(is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header=input$hd, sep=input$sep)
  })
  
  get.numeric.colnames = reactive({
    if(is.null(datasetInput()))
      return(NULL)
    data = datasetInput()
    colnames = NULL
    for(coln in names(data)){
      i = which.names(coln, names(data))
      if(is.numeric(data[,i]))
        colnames = append(colnames, coln)
    }
    return(colnames)
  })
  
  get.selected.data = reactive({
    data = datasetInput()
    coln = input$seletedData
    i = which.names(coln, names(data))
    r = data[,i]
    return(r)
  })
  
  get.selected.status = reactive({
    coln = input$seletedStatus
    data = datasetInput()
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return(NULL)
    r = data[,i]
    return(r)
  })

  
  get.sys.exp = reactive({
    i_str = input$exp_input
    if(i_str == '')
      return('Please enter a expression.')
    return(c(get_sys_expression(i_str), get_sys_d_expression(i_str)))
  })
  
  get.parameter.fn = reactive({
    fn = switch(input$distribution,
                'w'=get.weib2.result,
                'e'=get.exp.result,
                'l'=get.logNormal.result,
                get.weib2.result)
    
    return(fn)
  })
  
  get.parameter.table = function(censoring=FALSE, ...){
    if(is.null(datasetInput()))
      return(NULL)
    data = datasetInput()
    estimation.list = list('distribution'=NULL, 'par1'=NULL, 'par2'=NULL, 'par3'=NULL, 'par4'=NULL, 'par5'=NULL, 'N'=NULL, 'censored'=NULL, 'value'=NULL)
    fns = c(get.exp.result, get.gamma.result, get.gumbel.result, get.logNormal.result, get.mixedWeib.result, get.weib2.result, get.weib3.result)
    coln = input$selectedData
    i = which.names(coln, names(data))
    if(is.numeric(data[,i])){
      coln2 = input$selectedStatus
      j = which.names(coln2, names(data))
      for(fn in fns){
        if(censoring && is.numeric(data[,j])){
          result = fn(data[,i], q_data=data[,j])
        }else{
          result = fn(data[,i])
        }
        parameters = round(result$par, digits=3)
        for(n in c(2:6)){
          if(is.na(parameters[n-1])){
            estimation.list[[n]] = append(estimation.list[[n]], '')
          }else{
            estimation.list[[n]] = append(estimation.list[[n]], parameters[n-1])
          }
        }
        estimation.list$N = append(estimation.list$N, result$N)
        estimation.list$censored = append(estimation.list$censored, result$Zensierung)
        estimation.list$value = append(estimation.list$value, round(-result$value, digits=3))
        estimation.list$distribution = append(estimation.list$distribution, result$distribution)
        
      }
    }
    
    pm.dataset = data.frame(estimation.list)
    pm.dataset = pm.dataset[order(-pm.dataset$value),]
    return(pm.dataset)
  }
  
  output$table = renderTable({
    head(datasetInput(), n=input$obs)
  })
  
  output$summary = renderPrint({
    summary(datasetInput())
  })
  
  output$pm.table = renderTable({
    if(is.null(input$selectedStatus)){
      return()
    }else if(input$selectedStatus == 'None'){
      get.parameter.table()
    }else{
      get.parameter.table(censoring=TRUE)
    }
  })
  
  output$selectableData = renderUI({
    if(is.null(get.numeric.colnames()))
      return(NULL)
    tagList(
      radioButtons('selectedData', 'Choose a time data',
                   get.numeric.colnames()),
      radioButtons('selectedStatus', 'Choose a status data',
                   append('None', get.numeric.colnames()))
    )
  })

  output$plot1 = renderPlot({
    if(is.null(datasetInput()))
      return()
    data = datasetInput()
    coln = input$selectedData
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return()
    r = get.weib2.result(data[,i])
    par(mfrow=c(1,2))
    hist(data[,i], col='darkgray', border='white', main=paste('Histogram of', coln), prob=T, xlab=coln)
    rug(jitter(data[,i]))
    weibplot(data[,i], shape=c(1, r$par[1]), scale=c(mean(data[,i]), r$par[2]), labels=c('before', 'after'), mono=F, main=paste('Weib Plot of ', coln))
  })
  
  output$plot2 = renderPlot({
    if(is.null(datasetInput()))
      return()
    data = datasetInput()
    coln = input$selectedData
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return()
    r = get.weib2.result(data[,i])
    beta = r$par[1]
    tau = r$par[2]
    par(mfrow=c(1,2))
    boxplot(data[,i], ylab=coln, main=paste('boxplot of ', coln))
    rug(jitter(data[,i]), side=2)
    abline(h=mean(data[,i], na.rm=T), lty=2)
    plot(density(data[,i]), type='l', col='red', lty=2, ylim=c(0,0.02), main='')
    lines(data[,i], dweibull(data[,i], beta, tau), col='green', lty=1, type='p', pch=4)
    title('Weibul distribution and densitiy line')
    legend(x=min(data[,i]), legend=c('density', 'dweibull'), col=c('red', 'green'), lty=c(2, 1), pch=c(-1, 4), merge=TRUE, bg='gray90')
    
  })
  
  output$plot3 = renderPlot({
    if(is.null(datasetInput()))
      return()
    data = datasetInput()
    coln = input$selectedData
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return()
    r = get.weib2.result(data[,i])
    beta = r$par[1]
    tau = r$par[2]
    plot(x=data[,i], y=pweibull(data[,i], beta, tau), xlab=coln, ylab='Probability', main=paste('Weibull distribution curve of', coln), col='blue')
  })
  
  output$plot4 = renderPlot({
    if(is.null(datasetInput()))
      return()
    data = datasetInput()
    coln1 = input$selectedData
    coln2 = input$selectedStatus
    i = which.names(coln1, names(data))
    j = which.names(coln2, names(data))
    if(length(i) == 0 || length(j) == 0)
      return()
    #bwplot(data[,1]~data[,i], data=data, ylab=names(data)[1], xlab=coln)
    require(quantreg)
    require(survival)
    y = data[,i]
    d = data[,j]
    if(length(d[d == TRUE]) == 0){
      return('Invalid status data')
    }
    plot(survfit(Surv(y, d)~1), main='Kaplan-Meier Plot')
    f = crq(Surv(y, d)~1, method='Portnoy', grid='pivot')
    x = f$sol[2,]
    p = 1 - f$sol[1,]
    p = c(p, p[length(p)])
    par(col='red')
    fs = plot(stepfun(x, p), do.points=F, add=T)
  })

  output$element = renderTable({
    if(is.null(datasetInput()))
      return()
    data = datasetInput()
    coln = input$selectedData
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return('Choose a data')
    exp_str = input$exp_input
    dis_fn = input$distribution
    fn = get.parameter.fn()
    r = fn(data[,i])
    par1 = r$par[1]
    par2 = r$par[2]
    t = data[,i]
    if(exp_str == ''){
      value = rep(NA, times=length(t))
      values = list(value=value, density=value)
    }else if(is.na(input$par1) && is.na(input$par2)){  
      values = get.sys.value(exp_str, t, par1, par2, dis_fn)
    }else{
      values = get.sys.value(exp_str, t, input$par1, input$par2, dis_fn)
    }
    distribution = switch(input$distribution,
                                'w'=pweibull,
                                'e'=pexp,
                                'l'=plnorm,
                                pweibull)
    if(input$distribution == 'e'){
      distribution.value = distribution(t, par1)
    }else{
      distribution.value = distribution(t, par1, par2)
    }
    id = seq(1, length(t))
    tb = data.frame('id'=id, 'data'=data[,i], 'distribution'=distribution.value, 'Survive Probability'=1-distribution.value, 'System Survive'=1-values$value, 'Density'=values$density)
    names(tb)[2] = coln
    head(tb, n=input$obs)
  },
  digits = 3)
  
  output$exp = renderText({
    i_str = input$exp_input
    dis_fn = input$distribution
    if(i_str == '')
      return('Please enter a expression.')
    re = get.expression(i_str, dis_fn)
    exp_sys = re$expression
    simp_exp = Simplify(exp_sys)
    return(as.character(simp_exp))
  })
  
  output$dplot = renderPlot({
    if(is.null(datasetInput()))
      return()
    data = datasetInput()
    coln = input$selectedData
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return('Choose a data')
    exp_str = input$exp_input
    if(exp_str == ''){
      return()
    }
    dis_fn = input$distribution
    fn = get.parameter.fn()
    r = fn(data[,i])
    par1 = r$par[1]
    par2 = r$par[2]
    end_t = max(data[,i])
    t = seq(0, end_t, 0.5)
    if(is.na(input$par1) && is.na(input$par2)){  
      values = get.sys.value(exp_str, t, par1, par2, dis_fn)
    }else{
      values = get.sys.value(exp_str, t, input$par1, input$par2, dis_fn)
    }
    plot(x=t, y=values$density, type='l', col="red", main=paste('Density Curve of ', coln), xlab=coln, ylab='density')
  })
  
  output$expTb = renderTable({
    if(input$exp_input == '')
      return()
    if(is.null(datasetInput()))
      return()
    exp_str = input$exp_input
    t = input$sys_t
    dis_fn = input$distribution
    data = datasetInput()
    coln = input$selectedData
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return('Choose a data')
    exp_str = input$exp_input
    fn = get.parameter.fn()
    r = fn(data[,i])
    if(is.na(input$par1) && is.na(input$par2)){
      values = get.sys.value(exp_str, t, r$par[1], r$par[2], dis_fn)
    }else{
      values = get.sys.value(exp_str, t, input$par1, input$par2, dis_fn)
    }
    tb = data.frame('System Survive'=1-values$value, 'Density'=values$density)
    return(tb)
  },
  digits = 3)
  
  output$d_graph = renderText({
    i_str = input$exp_input
    dis_fn = input$distribution
    if(i_str == '')
      return()
    if(input$selectedGraph != 'd')
      return()
    require(d3Network)
    re = get.expression(i_str, dis_fn)
    links = re$links
    names(links) = c('source', 'target')
    shiny.path = path.package('shiny')
    dir.www = paste0(shiny.path, '/www-dir/elements_links.html')
    d3SimpleNetwork(Data=links,Source = "source", Target = "target",
                    height=600, width=900, fontsize = 13,
                    linkDistance=150, file=dir.www)
    return('<iframe src="elements_links.html" height=1308 width=1890 scrolling="no" frameborder="0"></iframe>')
  })
  
  output$s_graph =renderPlot({
    i_str = input$exp_input
    dis_fn = input$distribution
    if(i_str == '')
      return()
    if(input$selectedGraph != 's')
      return()
    re = get.expression(i_str, dis_fn)
    require(igraph)
    g = graph_from_data_frame(re$links)
    #E(g)$curved = 0.2 
    plot(g, layout=layout_as_tree, vertex.size=6,
         vertex.color='#3182bd', vertex.label.cex=1.5,
         vertex.label.dist=1, vertex.label.color='red')
  })
})