ent_nan = function(x){
  if(is.nan(x)){
    reutrn(0)}else{
      return(x)}
}

#get.gamma.result = function(data){
#  gammaFn = function(theta, x){
#    alpha = theta[1]
#    beta = theta[2]
#    n = length(x)
#    lnL = sum(dgamma(x, alpha, beta, log=TRUE)) #sum(log((x**(alpha-1))/(gamma(alpha)*(beta**alpha)))) - sum(x/beta)
#    return(-lnL)
#  }
  
#  r = optim(c(1, 1), gammaFn, x=data)
#  return(r)
#}

get.gamma.result = function(x, q_data=NULL){
  if(is.null(q_data)){
    gammaFn = function(theta, x){
      alpha = theta[1]
      beta = theta[2]
      n = length(x)
      lnL = sum(dgamma(x, alpha, beta, log=TRUE)) #sum(log((x**(alpha-1))/(gamma(alpha)*(beta**alpha)))) - sum(x/beta)
      return(-lnL)
    }
    r = optim(c(1, 1), gammaFn, x=x)
  }else{
    gammaFn = function(theta, x, q){
      alpha = theta[1]
      beta = theta[2]
      lnL = sum(dgamma(x, alpha, beta, log=TRUE)*q + log(1-pgamma(x, alpha, beta))*(1-q))
      return(-lnL)
    }
    r = optim(c(5, 1), gammaFn, x=x, q=q_data)
  }
  zensiert = q_data[q_data == 0]
  r['N'] = length(x)
  r['Zensierung'] = length(zensiert)
  return(r)
}

#get.gumbel.result = function(data){
#  gumbelFn = function(theta, x){
#    miu = theta[1]
#    beta = theta[2]
#    n = length(x)
#    lnL = n*log(1/beta) - sum((x-miu)/beta + exp(-((x-miu)/beta)))
#    return(-lnL)
#  }
#  r = optim(c(0, 1), gumbelFn, x=data)
#  return(r)
#}

get.gumbel.result = function(x, q_data=NULL){
  if(is.null(q_data)){
    gumbelFn = function(theta, x){
      miu = theta[1]
      beta = theta[2]
      n = length(x)
      lnL = n*log(1/beta) - sum((x-miu)/beta + exp(-((x-miu)/beta)))
      return(-lnL)
    }
    r = optim(c(0, 1), gumbelFn, x=x)
  }else{
    gumbelFn = function(theta, x, q){
      miu = theta[1]
      beta = theta[2]
      n = length(x)
      lnL = sum((log(1/beta) - (x-miu)/beta + exp(-((x-miu)/beta)))*q + log(1-exp(-exp(-(x-miu)/beta)))*(1-q))
      return(-lnL)
    }
    r = optim(c(3, 5), fn=gumbelFn, x=x, q=q_data)
  }
  
  zensiert = q_data[q_data == 0]
  r['N'] = length(x)
  r['Zensierung'] = length(zensiert)
  return(r)
}

#get.logNormal.result = function(x){
#  n = length(x)
#  miu = sum(log(x))/n
#  sigm = sqrt(sum((log(x) - miu)**2)/n)
#  value = -(sum(log(1/(x*sigm*sqrt(2*pi))*exp(-((log(x)-miu)**2/(2*(sigm**2)))))))
#  r = list('par'=c(miu, sigm), 'value'=value)
#  return(r)
#}

get.logNormal.result = function(x, q_data=NULL){
  if(is.null(q_data)){
    n = length(x)
    miu = sum(log(x))/n
    sigm = sqrt(sum((log(x) - miu)**2)/n)
    value = -(sum(log(1/(x*sigm*sqrt(2*pi))*exp(-((log(x)-miu)**2/(2*(sigm**2)))))))
    r = list('par'=c(miu, sigm), 'value'=value)
  }else{
    require(stats)
    lognormFn = function(theta, x, q){
      miu = theta[1]
      sigm = theta[2]
      lnL = sum(dlnorm(x, miu, sigm, log=TRUE)*q + log(1-plnorm(x, miu, sigm))*(1-q))
      return(-lnL)
    }
    r = optim(c(0, 1), fn=lognormFn, x=x, q=q_data)
  }
  zensiert = q_data[q_data == 0]
  r['N'] = length(x)
  r['Zensierung'] = length(zensiert)
  return(r)
}

#get.exp.result = function(x){
#  adv = mean(x)
#  lambda = 1/adv
#  value = -(sum(log(lambda*(exp(1)**-(lambda*x)))))
#  r = list('par'=lambda, 'value'=value)
#  return(r)
#}

get.exp.result = function(x, q_data=NULL){
  if(is.null(q_data)){
    adv = mean(x)
    lambda = 1/adv
    value = -(sum(log(lambda*(exp(1)**-(lambda*x)))))
    r = list('par'=lambda, 'value'=value)
  }else{
    expFn = function(theta, x, q){
      lamb = theta
      lnL = sum(log(lamb*exp(-lamb*x))*q + log(exp(-lamb*x))*(1-q))
      return(-lnL)
    }
    init_v = 1/mean(x)
    r = optim(init_v, fn=expFn, x=x, q=q_data)
  }
  zensiert = q_data[q_data == 0]
  r['N'] = length(x)
  r['Zensierung'] = length(zensiert)
  return(r)
}

#get.weib3.result = function(data){
#  adv = mean(data)
#  weib3Fn = function(theta, x){
#    n = length(x)
#    beta = theta[1]
#    tau = theta[2]
#    sigm = theta[3]
#    lnL = n*log(beta/(tau**beta)) + sum(log((x-sigm)**(beta-1))) - sum(((x-sigm)/tau)**beta)
#    return(-lnL)
#  }
  
#  r = optim(c(1, adv, 0), weib3Fn, x=data)
#  return(r)
#}

get.weib3.result = function(x, q_data=NULL){
  if(is.null(q_data)){
    q_data = rep(1, times=length(x))
  }
  adv = mean(x)
  weib3Fn = function(theta, x, q){
    beta = theta[1]
    tau = theta[2]
    sigm = theta[3]
    lnL = sum((log(beta/(tau**beta)) + log((x-sigm)**(beta-1)) - ((x-sigm)/tau)**beta)*q + log(1-pweibull3(x, beta, tau, sigm))*(1-q))
    return(-lnL)
  }
  
  r = optim(c(1, adv, 0), weib3Fn, x=x, q=q_data)
  zensiert = q_data[q_data == 0]
  r['N'] = length(x)
  r['Zensierung'] = length(zensiert)
  return(r)
}

#get.mixedWeib.result = function(data){
#  adv = mean(data)
#  mixedWeib = function(theta, x){
#    n = length(x)
#    beta1 = theta[1]
#    tau1 = theta[2]
#    beta2 = theta[3]
#    tau2 = theta[4]
#    p = theta[5]
#    if(p < 0 | p > 1){
#      return(NA)
#    }
#    lnL = sum(log((p/(tau1**beta1))*beta1*(x**(beta1-1))*exp(-((x/tau1)**beta1)) + ((1-p)/(tau2**beta2))*beta2*(x**(beta2-1))*exp(-((x/tau2)**beta2))))
#    return(-lnL)
#  }
#  r = optim(c(1, adv, 1, adv, 0), mixedWeib, x=data)
#  return(r)
#}

get.mixedWeib.result = function(x, q_data=NULL){
  if(is.null(q_data)){
    q_data = rep(1, times=length(x))
  }
  adv = mean(x)
  mixedWeib = function(theta, x, q){
    n = length(x)
    beta1 = theta[1]
    tau1 = theta[2]
    beta2 = theta[3]
    tau2 = theta[4]
    p = theta[5]
    if(p < 0 | p > 1){
      return(NA)
    }
    lnL = sum((log((p/(tau1**beta1))*beta1*(x**(beta1-1))*exp(-((x/tau1)**beta1)) + ((1-p)/(tau2**beta2))*beta2*(x**(beta2-1))*exp(-((x/tau2)**beta2))))*q + log(1-(p*pweibull(x, beta1, tau1) + (1-p)*pweibull(x, beta2, tau2)))*(1-q))
    return(-lnL)
  }
  r = optim(c(1, adv, 1, adv, 0), mixedWeib, x=x, q=q_data)
  zensiert = q_data[q_data == 0]
  r['N'] = length(x)
  r['Zensierung'] = length(zensiert)
  return(r)
}

#get.weib2.result = function(data){
#  adv = mean(data)
#  weibull.ml = function(theta, x){
#    beta = theta[1]
#    tau = theta[2]
#    n = length(x)
#    lnL = n*log(beta/(tau**beta)) + sum(log(x**(beta-1))) - sum((x/tau)**beta)
#    return(-lnL)
#  }

#  r = optim(c(1, adv), weibull.ml, x=data)
  
#  return(r)
#}

get.weib2.result = function(x, q_data=NULL){
  if(is.null(q_data)){
    q_data = rep(1, times=length(x))
  }
  
  adv = mean(x)
  weibull.ml = function(theta, x, q){
    beta = theta[1]
    tau = theta[2]
    n = length(x)
    lnL = sum(log((beta/(tau**beta)) * (x**(beta-1)) * exp(- (x/tau)**beta)) * q + log(1-pweibull(x, beta, tau))*(1-q))
    return(-lnL)
  }
  r = optim(c(1, adv), weibull.ml, x=x, q=q_data)
  zensiert = q_data[q_data == 0]
  r['N'] = length(x)
  r['Zensierung'] = length(zensiert)
  return(r)
}

serienSys = function(e, fn){
  r = fn(e)
  beta = r$par[1]
  tau = r$par[2]
  p = 1 - pweibull(e, beta, tau)
  Pn = prod(p)
  return(as.character(Pn))
}

parallelSys = function(e, fn){
  r = fn(e)
  beta = r$par[1]
  tau = r$par[2]
  f = pweibull(e, beta, tau)
  Pn = 1 - prod(f)
  return(as.character(Pn))
}

get.expression = function(fn_str, dis_fn='w'){
  distributions = list(w='(1-exp(-(t/tau)**beta))',
                       e='(1-exp(-lambda*t))',
                       l='(1/(t*sigm*(2*pi)^(1/2))*exp(-(log(t)-miu)^2/(2*sigm^2)))'
                       )
  exp.distribution = distributions[dis_fn]
  node.id = 1
  P = function(...){
    els = list(...)
    #exp.distribution = '(1-exp(-(t/tau)**beta))'
    exp1 = '1'
    from = NULL
    to = NULL
    for(ele in els){
      if(is.numeric(ele)){
        if(length(ele) > 1){
          for(i in ele){
            exp1 = paste(exp1, '*', exp.distribution)
            from = append(from, c('prev', paste0('P', i)))
            to = append(to, c(paste0('P', i), 'next'))
          }
        }else{
          exp1 = paste(exp1, '*', exp.distribution)
          from = append(from, c('prev', paste0('P', ele)))
          to = append(to, c(paste0('P', ele), 'next'))
        }
      }else if(is.list(ele)){
        exp1 = paste(exp1, '* (1-', ele$expression, ')')
        from = append(from, ele$from)
        to = append(to, ele$to)
      }
    }
    exp2 = paste('(1 -', exp1, ')')
    result = list(expression=exp2, from=from, to=to, type='P')
    return(result)
    #exp3 = Simplify(exp2)
  }
  
  S = function(...){
    els = list(...)
    #exp_survive = 'exp(-(t/tau)**beta)'
    exp.survive = paste('(1-', exp.distribution, ')')
    exp1 = '1'
    from = NULL
    to = NULL
    for(ele in els){
      if(is.numeric(ele)){
        if(length(ele) > 1){
          for(i in ele){
            exp1 = paste(exp1, '*', exp.survive)
            if(is.null(from) && is.null(to)){
              from = append(from, c('prev', paste0('S', i)))
              to = append(to, c(paste0('S', i), 'next'))
            }else{
              node = paste0('S', i)
              replace.index = which(to == 'next')
              to = replace(to, replace.index, node)
              from = append(from, node)
              to = append(to, 'next')
            }
          }
        }else{
          exp1 = paste(exp1, '*', exp.survive)
          if(is.null(from) && is.null(to)){
            from = append(from, c('prev', paste0('S', ele)))
            to = append(to, c(paste0('S', ele), 'next'))
          }else{
            node = paste0('S', ele)
            replace.index = which(to == 'next')
            to = replace(to, replace.index, node)
            from = append(from, node)
            to = append(to, 'next')
          }
        }
      }else if(is.list(ele)){
        exp1 = paste(exp1, '*', ele$expression)
        if(is.null(from) && is.null(to)){
          from = append(from, ele$from)
          to = append(to, ele$to)
        }else{
          if(ele$type == 'S'){
            node = ele$to[1]
            replace.index = which(to == 'next')
            replace(to, replace.index, node)
            ele$to = ele$to[2:length(ele$to)]
            ele$from = ele$from[2:length(ele$from)]
            from = append(from, ele$from)
            to = append(to, ele$to)
          }else{
            if(startsWith(from[length(from)], 'S')){
              node = from[length(from)]
              replace.index = which(ele$from == 'prev')
              ele$from = replace(ele$from, replace.index, node)
              from = from[1:length(from)-1]
              to = to[1:length(to)-1]
              from = append(from, ele$from)
              to = append(to, ele$to)
            }else{
              node = paste0('node', node.id)
              node.id = node.id + 1
              replace.index1 = which(to == 'next')
              replace.index2 = which(ele$from == 'prev')
              to = replace(to, replace.index1, node)
              ele$from = replace(ele$from, replace.index2, node)
              from = append(from, ele$from)
              to = append(to, ele$to)
            }
          }
        }
      }
    }
    result = list(expression=paste('(', exp1, ')'), from=from, to=to, type='S')
    return(result)
  }
  eval.result = eval(parse(text=fn_str))
  exp.sys = parse(text=eval.result$expression)
  replace1 = which(eval.result$from == 'prev')
  replace2 = which(eval.result$to == 'next')
  from = replace(eval.result$from, replace1, 'start')
  to = replace(eval.result$to, replace2, 'end')
  value = rnorm(length(from), 5, 2)
  links = data.frame(N1=from, N2=to, Value=value, stringsAsFactors=FALSE)
  r = list(expression=exp.sys, links=links)
  return(r)
}

#get_sys_expression = function(fn_str){
#  #exp_sys = parse(text=eval(parse(text=fn_str)))
#  exp_sys = get.expression(parse(text=fn_str))
#  #return(yacas(exp_sys)$text)
#  return(exp_sys)
#}

get_sys_d_expression = function(fn_str, dis_fn){
  #return(D(parse(text=eval(parse(text=fn_str))), 't'))
  re = get.expression(fn_str, dis_fn)
  return(D(re$expression, 't'))
}

get_sys_value = function(fn_str, lebensdauer, par1, par2, dis_fn='w'){
  if(dis_fn == 'w'){
    beta = par1
    tau = par2
  }else if(dis_fn == 'e'){
    lambda = par1
  }else if(dis_fn == 'l'){
    miu = par1
    sigm = par2
  }
  t = lebensdauer
  #exp_sys = parse(text=eval(parse(text=fn_str)))
  re = get.expression(fn_str, dis_fn)
  exp_sys = re$expression
  sys_value = eval(exp_sys)
  gradient_value = eval(D(exp_sys, 't'))
  return(list(value=sys_value, density=gradient_value))
}
