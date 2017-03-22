library(shiny)
library(Ryacas)
library(Deriv)
#library(grDevices)
library(datasets)
library(car)
library(lattice)
library(FAdist)
library(Renext)
source('parameters_ml.R')


shinyServer(function(input, output){
  datasetInput = reactive({
    inFile = input$file
    if(is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header=input$hd, sep=input$sep)
  })
  
  getColnames = reactive({
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
  
  #get_exp_str = reactive({
  #  return(eval(parse(text=input$exp_input)))
  #})
  
  get.sys.exp = reactive({
    i_str = input$exp_input
    if(i_str == '')
      return('Please enter a expression.')
    return(c(get_sys_expression(i_str), get_sys_d_expression(i_str)))
  })
  
  get.parameter.table = function(fn, distribution, ...){
    if(is.null(datasetInput()))
      return(NULL)
    data = datasetInput()
    dataId = NULL
    par1 = NULL
    par2 = NULL
    par3 = NULL
    par4 = NULL
    par5 = NULL
    value = NULL
    N = NULL
    zens = NULL
    #for(coln in names(data)){
    coln = input$selectedData
      i = which.names(coln, names(data))
      if(is.numeric(data[,i])){
        result = (fn(data[,i], ...))
        parameters = result$par
        if(!is.na(parameters[1]))
          par1 = append(par1, parameters[1])
        if(!is.na(parameters[2]))
          par2 = append(par2, parameters[2])
        if(!is.na(parameters[3]))
          par3 = append(par3, parameters[3])
        if(!is.na(parameters[4]))
          par4 = append(par4, parameters[4])
        if(!is.na(parameters[5]))
          par5 = append(par5, parameters[5])
        N = append(N, result$N)
        zens = append(zens, result$Zensierung)
        value = append(value, -result$value)
        dataId = append(dataId, coln)
      }
    #}
    if(distribution == 'weib2')
      pm.dataset = data.frame('data'=dataId, 'beta'=par1, 'tau'=par2, 'value'=value)
    if(distribution == 'z_weib2')
      pm.dataset = data.frame('data'=dataId, 'beta'=par1, 'tau'=par2, 'N'=N, 'Zensierung'=zens, 'value'=value)
    if(distribution == 'weib3')
      pm.dataset = data.frame('data'=dataId, 'beta'=par1, 'tau'=par2, 'sigma'=par3, 'value'=value)
    if(distribution == 'z_weib3')
      pm.dataset = data.frame('data'=dataId, 'beta'=par1, 'tau'=par2, 'sigma'=par3, 'N'=N, 'Zensierung'=zens, 'value'=value)
    if(distribution == 'mixedWeib')
      pm.dataset = data.frame('data'=dataId, 'beta1'=par1, 'tau1'=par2, 'beta2'=par3, 'tau2'=par4, 'p'=par5, 'value'=value)
    if(distribution == 'z_mixedWeib')
      pm.dataset = data.frame('data'=dataId, 'beta1'=par1, 'tau1'=par2, 'beta2'=par3, 'tau2'=par4, 'p'=par5, 'N'=N, 'Zensierung'=zens, 'value'=value)
    if(distribution == 'exp')
      pm.dataset = data.frame('data'=dataId, 'lambda'=par1, 'value'=value)
    if(distribution == 'logNormal')
      pm.dataset = data.frame('data'=dataId, 'miu'=par1, 'sigma'=par2, 'value'=value)
    if(distribution == 'gumbel')
      pm.dataset = data.frame('data'=dataId, 'miu'=par1, 'beta'=par2, 'value'=value)
    if(distribution == 'gamma')
      pm.dataset = data.frame('data'=dataId, 'alpha'=par1, 'beta'=par2, 'value'=value)
    return(pm.dataset)
  }
  
  output$table = renderTable({
    head(datasetInput(), n=input$obs)
  })
  
  output$summary = renderPrint({
    summary(datasetInput())
  })
  
  output$weib2.pm.table = renderTable({
    data = datasetInput()
    coln = input$selectedStatus
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return()
    get.parameter.table(get.weib2.zensiert.result, 'z_weib2', q_data=data[,i])
  })
  
  output$z_weib2.pm.table = renderTable({
    get.parameter.table(get.weib2.zensiert.result, 'z_weib2', input$sys_t)
  })
  
  output$weib3.pm.table = renderTable({
    get.parameter.table(get.weib3.result, 'weib3')
  })
  
  output$z_weib3.pm.table = renderTable({
    get.parameter.table(get.weib3.zensiert.result, 'z_weib3', input$sys_t)
  })
  
  output$mixedWeib.pm.table = renderTable({
    get.parameter.table(get.mixedWeib.result, 'mixedWeib')
  })
  
  output$z_mixedweib.pm.table = renderTable({
    get.parameter.table(get.mixedweib.zensiert.result, 'z_mixedWeib', input$sys_t)
  })
  
  output$exp.pm.table = renderTable({
    get.parameter.table(get.exp.result, 'exp')
  })
  
  output$logNormal.pm.table = renderTable({
    get.parameter.table(get.logNormal.result, 'logNormal')
  })
  
  output$gumbel.pm.table = renderTable({
    get.parameter.table(get.gumbel.result, 'gumbel')
  })
  
  output$gamma.pm.table = renderTable({
    get.parameter.table(get.gamma.result, 'gamma')
  })
  
  output$selectableData = renderUI({
    if(is.null(getColnames()))
      return(NULL)
    tagList(
      radioButtons('selectedData', 'Choose a time data',
                   getColnames()),
      radioButtons('selectedStatus', 'Choose a status data',
                   getColnames())
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
    #lines(density(data[,i], n=length(data[,i])), col='red')
    rug(jitter(data[,i]))
    #qq.plot(data[,i], ylab=coln, main=paste('Normal QQ plot of', coln))
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
    r = get.weib3.result(data[,i])
    beta = r$par[1]
    tau = r$par[2]
    par(mfrow=c(1,2))
    boxplot(data[,i], ylab=coln, main=paste('boxplot of ', coln))
    rug(jitter(data[,i]), side=2)
    abline(h=mean(data[,i], na.rm=T), lty=2)
    #line_data = matrix(c(density(data[,i]), pweibull(data[,i], beta, tau)), length(data[,i]), 2)
    #line_data = data.frame('density'=density(data[,i]), 'weibull'=pweibull(data[,i], beta, tau))
    #matplot(line_data, type='l', lty=1:2, col=2:3, bg=4:5, main=paste('Density and weibull distribution line of ', coln))
    plot(density(data[,i]), type='l', col='red', ylim=0:1, lty=2, main='')
    lines(data[,i], dweibull(data[,i], beta, tau), col='green', lty=1, type='b', pch=4)
    lines(data[,i], pweibull(data[,i], beta, tau), col='blue', lty=3, type='b', pch=3)
    title('Weibul distribution and densitiy line')
    legend(x=min(data[,i]), y=0.5, legend=c('density', 'dweibull', 'pweibull'), col=c('red', 'green', 'blue'), lty=c(2, 1, 3), pch=c(-1, 4, 3), merge=TRUE, bg='gray90')
    
  })
  
  output$plot3 = renderPlot({
    if(is.null(datasetInput()))
      return()
    data = datasetInput()
    coln = input$selectedData
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return()
    #bwplot(data[,1]~data[,i], data=data, ylab=names(data)[1], xlab=coln)
    require(quantreg)
    require(survival)
    y = pmin(data[,i], input$sys_t)
    d = (y < input$sys_t)
    if(length(d[d == TRUE]) == 0){
      return('Invalid observation time')
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
    seriensystem = NULL
    parallelsystem = NULL
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return('Choose a data')
    exp_str = input$exp_input
    r = get.weib3.result(data[,i])
    beta = r$par[1]
    tau = r$par[2]
    #t = input$sys_t
    t = data[,i]
    if(exp_str == ''){
      value = rep(NA, times=length(t))
      values = rep(value, times=2)
    }else{  
      values = get_sys_value(exp_str, t, beta, tau)
    }
    value_weibull = pweibull(t, beta, tau)
    id = seq(1, length(t))
    tb = data.frame('id'=id, 'data'=data[,i], 'Weibull'=as.character(value_weibull), 'Survive Probability'=as.character(1-value_weibull), 'System Survive'=as.character(values[1]), 'Density'=as.character(values[2]))
    names(tb)[2] = coln
    head(tb, n=input$obs)
  })
  
  #get_exp_value = reactive({
  #  exp_sys = parse(text=eval(parse(text=input$exp_input)))
  #  t = input$sys_t
  #  return(eval(exp_sys))
  #})
  
  output$exp = renderText({
    i_str = input$exp_input
    if(i_str == '')
      return('Please enter a expression.')
    exp_sys = get_sys_expression(i_str)
    simp_exp = Simplify(exp_sys)
    return(as.character(simp_exp))
  })
  
  output$dexp = renderText({
    i_str = input$exp_input
    if(i_str == '')
      return('Please enter a expression.')
    dexp = get_sys_d_expression(i_str)
    simp_exp = Simplify(dexp)
    return(as.character(simp_exp))
  })
  
  output$expTb = renderTable({
    exp_str = input$exp_input
    t = input$sys_t
    if(exp_str == '')
      return()
    if(is.null(datasetInput()))
      return()
    data = datasetInput()
    coln = input$selectedData
    seriensystem = NULL
    parallelsystem = NULL
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return('Choose a data')
    exp_str = input$exp_input
    r = get.weib3.result(data[,i])
    beta = r$par[1]
    tau = r$par[2]
    values = get_sys_value(exp_str, t, beta, tau)
    tb = data.frame('System Survive'=as.character(values[1]), 'Density'=as.character(values[2]))
    return(tb)
  })
  #output$exp2 = renderText({
  #  return(get_exp_value())
  #})
})