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
    if(p < 0 || p > 1){
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

get.weib2.result = function(x, q_data=NULL){
  if(is.null(q_data)){
    q_data = rep(1, times=length(x))
  }
  
  adv = mean(x)
  weibull.ml = function(theta, x, q){
    beta = theta[1]
    tau = theta[2]
    lnL = sum(log((beta/(tau**beta)) * (x**(beta-1)) * exp(- (x/tau)**beta)) * q + log(1-pweibull(x, beta, tau))*(1-q))
    return(-lnL)
  }
  r = optim(c(1, adv), weibull.ml, x=x, q=q_data)
  zensiert = q_data[q_data == 0]
  r['N'] = length(x)
  r['Zensierung'] = length(zensiert)
  return(r)
}
