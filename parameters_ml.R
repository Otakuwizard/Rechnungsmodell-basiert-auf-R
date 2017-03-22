ent_nan = function(x){
  if(is.nan(x)){
    reutrn(0)}else{
      return(x)}
}

get.gamma.result = function(data){
  gammaFn = function(theta, x){
    alpha = theta[1]
    beta = theta[2]
    n = length(x)
    lnL = sum(log((x**(alpha-1))/(gamma(alpha)*(beta**alpha)))) - sum(x/beta)
    return(-lnL)
  }
  
  r = optim(c(1, 1), gammaFn, x=data)
  return(r)
}

get.gamma.zensiert.result = function(data, ob_t){
  z_data = pmin(data, ob_t)
  q_data = NULL
  for(ob in z_data){
    if(ob < ob_t){
      q_data = append(q_data, 1)
    }else{
      q_data = append(q_data, 0)
    }
  }
  gammaFn = function(theta, x, q){
    alpha = theta[1]
    beta = theta[2]
    lnL = sum((log((x**(alpha-1))/(gamma(alpha)*(beta**alpha))) - x/beta)*q + log(1-pgamma(x, alpha, beta))*(1-q))
    return(-lnL)
  }
  r = optim(c(1, 1), gammaFn, x=z_data, q=q_data)
  return(r)
}

get.gumbel.result = function(data){
  gumbelFn = function(theta, x){
    miu = theta[1]
    beta = theta[2]
    n = length(x)
    lnL = n*log(1/beta) - sum((x-miu)/beta + exp(-((x-miu)/beta)))
    return(-lnL)
  }
  r = optim(c(0, 1), gumbelFn, x=data)
  return(r)
}

get.gumbel.zensiert.result = function(data, ob_t){
  z_data = pmin(data, ob_t)
  q_data = NULL
  for(ob in z_data){
    if(ob < ob_t){
      q_data = append(q_data, 1)
    }else{
      q_data = append(q_data, 0)
    }
  }
  gumbelFn = function(theta, x, q){
    miu = theta[1]
    beta = theta[2]
    n = length(x)
    lnL = sum((log(1/beta) - (x-miu)/beta + exp(-((x-miu)/beta)))*q + log(1-exp(-exp(-(x-miu)/beta)))*(1-q))
    return(-lnL)
  }
  r = optim(c(0, 1), fn=gumbelFn, x=data, q=q_data)
  return(r)
}

get.logNormal.result = function(x){
  n = length(x)
  miu = sum(log(x))/n
  sigm = sqrt(sum((log(x) - miu)**2)/n)
  value = -(sum(log(1/(x*sigm*sqrt(2*pi))*exp(-((log(x)-miu)**2/(2*(sigm**2)))))))
  r = list('par'=c(miu, sigm), 'value'=value)
  return(r)
}

get.exp.result = function(x){
  n = length(x)
  adv = mean(x)
  lambda = 1/adv
  value = -(sum(log(lambda*(exp(1)**-(lambda*x)))))
  r = list('par'=lambda, 'value'=value)
  return(r)
}

get.weib3.result = function(data){
  adv = mean(data)
  weib3Fn = function(theta, x){
    n = length(x)
    beta = theta[1]
    tau = theta[2]
    sigm = theta[3]
    lnL = n*log(beta/(tau**beta)) + sum(log((x-sigm)**(beta-1))) - sum(((x-sigm)/tau)**beta)
    return(-lnL)
  }
  
  r = optim(c(1, adv, 0), weib3Fn, x=data)
  return(r)
}

get.weib3.zensiert.result = function(x, ob_t=0, q_data=NULL){
  if(is.null(q_data)){
    z_data = pmin(x, ob_t)
    for(ob in z_data){
      if(ob < ob_t){
        q_data = append(q_data, 1)
      }else{
        q_data = append(q_data, 0)
      }
    }
  }else{
    z_data = x
  }
  adv = mean(z_data)
  weib3Fn = function(theta, x, q){
    beta = theta[1]
    tau = theta[2]
    sigm = theta[3]
    lnL = sum((log(beta/(tau**beta)) + log((x-sigm)**(beta-1)) - ((x-sigm)/tau)**beta)*q + log(1-pweibull3(x, beta, tau, sigm))*(1-q))
    return(-lnL)
  }
  
  r = optim(c(1, adv, 0), weib3Fn, x=z_data, q=q_data)
  zensiert = q_data[q_data == 0]
  r['N'] = length(z_data)
  r['Zensierung'] = length(zensiert)
  return(r)
}

get.mixedWeib.result = function(data){
  adv = mean(data)
  mixedWeib = function(theta, x){
    n = length(x)
    beta1 = theta[1]
    tau1 = theta[2]
    beta2 = theta[3]
    tau2 = theta[4]
    p = theta[5]
    if(p < 0 | p > 1){
      return(NA)
    }
    lnL = sum(log((p/(tau1**beta1))*beta1*(x**(beta1-1))*exp(-((x/tau1)**beta1)) + ((1-p)/(tau2**beta2))*beta2*(x**(beta2-1))*exp(-((x/tau2)**beta2))))
    return(-lnL)
  }
  r = optim(c(1, adv, 1, adv, 0), mixedWeib, x=data)
  return(r)
}

get.mixedweib.zensiert.result = function(x, ob_t=0, q_data=NULL){
  if(is.null(q_data)){
    z_data = pmin(x, ob_t)
    for(ob in z_data){
      if(ob < ob_t){
        q_data = append(q_data, 1)
      }else{
        q_data = append(q_data, 0)
      }
    }
  }else{
    z_data = x
  }
  adv = mean(z_data)
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
  r = optim(c(1, adv, 1, adv, 0), mixedWeib, x=z_data, q=q_data)
  zensiert = q_data[q_data == 0]
  r['N'] = length(z_data)
  r['Zensierung'] = length(zensiert)
  return(r)
}

get.weib2.result = function(data){
  adv = mean(data)
  weibull.ml = function(theta, x){
    beta = theta[1]
    tau = theta[2]
    n = length(x)
    lnL = n*log(beta/(tau**beta)) + sum(log(x**(beta-1))) - sum((x/tau)**beta)
    return(-lnL)
  }
  
  r = optim(c(1, adv), weibull.ml, x=data)
  
  return(r)
}

get.weib2.zensiert.result = function(x, ob_t=0, q_data=NULL){
  if(is.null(q_data)){
    z_data = pmin(x, ob_t)
    for(ob in z_data){
      if(ob < ob_t){
        q_data = append(q_data, 1)
      }else{
        q_data = append(q_data, 0)
      }
    }
  }else{
    z_data = x
  }
  adv = mean(z_data)
  weibull.ml = function(theta, x, q){
    beta = theta[1]
    tau = theta[2]
    n = length(x)
    lnL = sum(log((beta/(tau**beta)) * (x**(beta-1)) * exp(- (x/tau)**beta)) * q + log(1-pweibull(x, beta, tau))*(1-q))
    return(-lnL)
  }
  r = optim(c(1, adv), weibull.ml, x=z_data, q=q_data)
  zensiert = q_data[q_data == 0]
  r['N'] = length(z_data)
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

P = function(...){
  els = c(...)
  exp.distribution = '(1-exp(-(t/tau)**beta))'
  exp1 = '1'
  for(ele in els){
    
    if(is.numeric(ele) || !startsWith(ele, '(')){
      exp1 = paste(exp1, '*', exp.distribution)
    }else if(is.character(ele) && startsWith(ele, '(')){
      exp1 = paste(exp1, '* (1-', ele, ')')
    }
  }
  exp2 = paste('(1 -', exp1, ')')
  return(exp2)
  #exp3 = Simplify(exp2)
}

S = function(...){
  els = c(...)
  exp_survive = 'exp(-(t/tau)**beta)'
  exp1 = '1'
  for(ele in els){
    if(is.numeric(ele) || !startsWith(ele, '(')){
      exp1 = paste(exp1, '*', exp_survive)
    }else if(is.character(ele) && startsWith(ele, '(')){
      exp1 = paste(exp1, '*', ele)
    }
  }
  return(paste('(', exp1, ')'))
}

get_sys_expression = function(fn_str){
  exp_sys = parse(text=eval(parse(text=fn_str)))
  #return(yacas(exp_sys)$text)
  return(exp_sys)
}

get_sys_d_expression = function(fn_str){
  return(D(parse(text=eval(parse(text=fn_str))), 't'))
}

get_sys_value = function(fn_str, lebensdauer, shape, scale){
  beta = shape
  tau = scale
  t = lebensdauer
  exp_sys = parse(text=eval(parse(text=fn_str)))
  sys_value = eval(exp_sys)
  gradient_value = eval(D(exp_sys, 't'))
  return(c(sys_value, gradient_value))
}