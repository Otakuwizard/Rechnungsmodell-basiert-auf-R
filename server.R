library(shiny)
library(Ryacas)
library(Deriv)
#library(grDevices)
library(datasets)
library(car)
library(lattice)
library(FAdist)
library(Renext)
library(d3Network)
source('parameters_ml.R')


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
    #if(length(i) == 0)
    #  return(NULL)
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
  
  #get_exp_str = reactive({
  #  return(eval(parse(text=input$exp_input)))
  #})
  
  get.sys.exp = reactive({
    i_str = input$exp_input
    if(i_str == '')
      return('Please enter a expression.')
    return(c(get_sys_expression(i_str), get_sys_d_expression(i_str)))
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
    #for(coln in names(data)){
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
    get.parameter.table(get.weib2.result, 'weib2')
  })
  
  output$z_weib2.pm.table = renderTable({
    get.parameter.table(get.weib2.result, 'z_weib2', censoring=TRUE)
  })
  
  output$weib3.pm.table = renderTable({
    get.parameter.table(get.weib3.result, 'weib3')
  })
  
  output$z_weib3.pm.table = renderTable({
    get.parameter.table(get.weib3.result, 'z_weib3', censoring=TRUE)
  })
  
  output$mixedWeib.pm.table = renderTable({
    get.parameter.table(get.mixedWeib.result, 'mixedWeib')
  })
  
  output$z_mixedweib.pm.table = renderTable({
    get.parameter.table(get.mixedWeib.result, 'z_mixedWeib', TRUE)
  })
  
  output$exp.pm.table = renderTable({
    get.parameter.table(get.exp.result, 'exp')
  })
  
  output$z_exp.pm.table = renderTable({
    get.parameter.table(get.exp.result, 'z_exp', TRUE)
  })
  
  output$logNormal.pm.table = renderTable({
    get.parameter.table(get.logNormal.result, 'logNormal')
  })
  
  output$z_logNormal.pm.table = renderTable({
    get.parameter.table(get.logNormal.result, 'z_logNormal', TRUE)
  })
  
  output$gumbel.pm.table = renderTable({
    get.parameter.table(get.gumbel.result, 'gumbel')
  })
  
  output$z_gumbel.pm.table = renderTable({
    get.parameter.table(get.gumbel.result, 'z_gumbel', TRUE)
  })
  
  output$gamma.pm.table = renderTable({
    get.parameter.table(get.gamma.result, 'gamma')
  })
  
  output$z_gamma.pm.table = renderTable({
    get.parameter.table(get.gamma.result, 'z_gamma', TRUE)
  })
  
  output$selectableData = renderUI({
    if(is.null(get.numeric.colnames()))
      return(NULL)
    tagList(
      radioButtons('selectedData', 'Choose a time data',
                   get.numeric.colnames()),
      radioButtons('selectedStatus', 'Choose a status data',
                   get.numeric.colnames())
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
    i = which.names(coln, names(data))
    if(length(i) == 0)
      return('Choose a data')
    exp_str = input$exp_input
    dis_fn = input$distribution
    r = get.weib3.result(data[,i])
    beta = r$par[1]
    tau = r$par[2]
    t = data[,i]
    if(exp_str == ''){
      value = rep(NA, times=length(t))
      values = list(value=value, density=value)
    }else if(is.na(input$par1) && is.na(input$par2)){  
      values = get_sys_value(exp_str, t, beta, tau)
    }else{
      values = get_sys_value(exp_str, t, input$par1, input$par2, dis_fn)
    }
    value_weibull = pweibull(t, beta, tau)
    id = seq(1, length(t))
    tb = data.frame('id'=id, 'data'=data[,i], 'Weibull'=value_weibull, 'Survive Probability'=1-value_weibull, 'System Survive'=values$value, 'Density'=values$density)
    names(tb)[2] = coln
    head(tb, n=input$obs)
  },
  digits = 3)
  
  #get_exp_value = reactive({
  #  exp_sys = parse(text=eval(parse(text=input$exp_input)))
  #  t = input$sys_t
  #  return(eval(exp_sys))
  #})
  
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
  
  output$dexp = renderText({
    i_str = input$exp_input
    dis_fn = input$distribution
    if(i_str == '')
      return('Please enter a expression.')
    dexp = get_sys_d_expression(i_str, dis_fn)
    simp_exp = Simplify(dexp)
    return(as.character(simp_exp))
  })
  
  output$expTb = renderTable({
    exp_str = input$exp_input
    t = input$sys_t
    dis_fn = input$distribution
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
    if(is.na(input$par1) && is.na(input$par2)){
      r = get.weib3.result(data[,i])
      values = get_sys_value(exp_str, t, r$par[1], r$par[2])
    }else{
      values = get_sys_value(exp_str, t, input$par1, input$par2, dis_fn)
    }
    tb = data.frame('System Survive'=values$value, 'Density'=values$density)
    return(tb)
  },
  digits = 3)
  
  output$links = renderPlot({
    i_str = input$exp_input
    dis_fn = input$distribution
    if(i_str == '')
      return('Please enter a expression.')
    #require(d3Network)
    re = get.expression(i_str, dis_fn)
    links = re$links
    #names(links) = c('source', 'target', 'value')
    nodes = data.frame(name=unique(c(links$N1, links$N2)), stringsAsFactors=FALSE)
    nodes$seq = 0:(nrow(nodes)-1)
    links = merge(links, nodes, by.x='N1', by.y='name')
    names(links)[4] = 'source'
    links = merge(links, nodes, by.x='N2', by.y='name')
    names(links)[5] = 'target'
    names(links)[3] = 'value'
    
    links = subset(links, select=c('source', 'target', 'value'))
    nodes = subset(nodes, select=c('name'))
    
    d3ForceNetwork(Links=links, Nodes=nodes, Source = "source",  
                    Target = "target", Value='value', NodeID='name',
                    height=1200, width=1800, fontsize = 13,
                    linkDistance=300, file='./elements_links.html')
    require(igraph)
    g = graph_from_data_frame(re$links)
    #E(g)$curved = 0.2 
    plot(g, layout=layout_as_tree, vertex.size=6,
         vertex.color='#3182bd', vertex.label.cex=1.5,
         vertex.label.dist=0.5, vertex.label.color='black')
  })
  
})