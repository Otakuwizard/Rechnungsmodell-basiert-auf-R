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
  
  get.parameter.table = function(fn, distribution, censoring=FALSE, ...){
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
    coln = input$selectedData
    i = which.names(coln, names(data))
    if(is.numeric(data[,i])){
      coln2 = input$selectedStatus
      j = which.names(coln2, names(data))
      if(censoring && is.numeric(data[,j])){
        result = fn(data[,i], q_data=data[,j])
      }else{
        result = fn(data[,i])
      }
      parameters = round(result$par, digits=3)
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
      value = append(value, round(-result$value, digits=3))
      dataId = append(dataId, coln)
    }
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
    if(distribution == 'z_exp')
      pm.dataset = data.frame('data'=dataId, 'lambda'=par1, 'N'=N, 'Zensierung'=zens, 'value'=value)
    if(distribution == 'logNormal')
      pm.dataset = data.frame('data'=dataId, 'miu'=par1, 'sigma'=par2, 'value'=value)
    if(distribution == 'z_logNormal')
      pm.dataset = data.frame('data'=dataId, 'miu'=par1, 'sigma'=par2, 'N'=N, 'Zensierung'=zens, 'value'=value)
    if(distribution == 'gumbel')
      pm.dataset = data.frame('data'=dataId, 'miu'=par1, 'beta'=par2, 'value'=value)
    if(distribution == 'z_gumbel')
      pm.dataset = data.frame('data'=dataId, 'miu'=par1, 'beta'=par2, 'N'=N, 'Zensierung'=zens, 'value'=value)
    if(distribution == 'gamma')
      pm.dataset = data.frame('data'=dataId, 'alpha'=par1, 'beta'=par2, 'value'=value)
    if(distribution == 'z_gamma')
      pm.dataset = data.frame('data'=dataId, 'alpha'=par1, 'beta'=par2, 'N'=N, 'Zensierung'=zens, 'value'=value)
    return(pm.dataset)
  }
  
  output$table = renderTable({
    head(datasetInput(), n=input$obs)
  })
  
  output$summary = renderPrint({
    summary(datasetInput())
  })
  
  output$weib2.pm.table = renderTable({
    if(is.null(input$selectedStatus)){
      return()
    }else if(input$selectedStatus == 'None'){
      get.parameter.table(get.weib2.result, 'weib2')
    }else{
      get.parameter.table(get.weib2.result, 'z_weib2', censoring=TRUE)
    }
  })
  
  output$weib3.pm.table = renderTable({
    if(is.null(input$selectedStatus)){
      return()
    }else if(input$selectedStatus == 'None'){
      get.parameter.table(get.weib3.result, 'weib3')
    }else{
      get.parameter.table(get.weib3.result, 'z_weib3', censoring=TRUE)
    }
  })
  
  output$mixedWeib.pm.table = renderTable({
    if(is.null(input$selectedStatus)){
      return()
    }else if(input$selectedStatus == 'None'){
      get.parameter.table(get.mixedWeib.result, 'mixedWeib')
    }else{
      get.parameter.table(get.mixedWeib.result, 'z_mixedWeib', TRUE)
    }
  })
  
  output$exp.pm.table = renderTable({
    if(is.null(input$selectedStatus)){
      return()
    }else if(input$selectedStatus == 'None'){
      get.parameter.table(get.exp.result, 'exp')
    }else{
      get.parameter.table(get.exp.result, 'z_exp', TRUE)
    }
  })
  
  output$logNormal.pm.table = renderTable({
    if(is.null(input$selectedStatus)){
      return()
    }else if(input$selectedStatus == 'None'){
      get.parameter.table(get.logNormal.result, 'logNormal')
    }else{
      get.parameter.table(get.logNormal.result, 'z_logNormal', TRUE)
    }
  })
  
  output$gumbel.pm.table = renderTable({
    if(is.null(input$selectedStatus)){
      return()
    }else if(input$selectedStatus == 'None'){
      get.parameter.table(get.gumbel.result, 'gumbel')
    }else{
      get.parameter.table(get.gumbel.result, 'z_gumbel', TRUE)
    }
  })
  
  output$gamma.pm.table = renderTable({
    if(is.null(input$selectedStatus)){
      return()
    }else if(input$selectedStatus == 'None'){
      get.parameter.table(get.gamma.result, 'gamma')
    }else{
      get.parameter.table(get.gamma.result, 'z_gamma', TRUE)
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