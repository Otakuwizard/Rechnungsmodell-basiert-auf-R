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
        node.prev = paste0('node', node.id)
        node.next = paste0('node', node.id+1)
        node.id = node.id + 2
        if(is.null(from) && is.null(to)){
          from = append(from, 'prev')
          to = append(to, 'next')
        }
        
        replace.index1 = which(to == 'next')
        replace.index2 = which(ele$from == 'prev')
        to = replace(to, replace.index1, node.prev)
        ele$from = replace(ele$from, replace.index2, node.prev)
        replace.index3 = which(ele$to == 'next')
        ele$to = replace(ele$to, replace.index3, node.next)
        from = append(from, ele$from)
        to = append(to, ele$to)
        
        from = append(from, node.next)
        to = append(to, 'next')
      }
    }
    result = list(expression=paste('(', exp1, ')'), from=from, to=to, type='S')
    return(result)
  }
  eval.result = eval(parse(text=fn_str))
  exp.sys = parse(text=paste('(1-', eval.result$expression, ')'))
  replace1 = which(eval.result$from == 'prev')
  replace2 = which(eval.result$to == 'next')
  from = replace(eval.result$from, replace1, 'start')
  to = replace(eval.result$to, replace2, 'end')
  links = data.frame(N1=from, N2=to, stringsAsFactors=FALSE)
  r = list(expression=exp.sys, links=links)
  return(r)
}

get.sys.d.expression = function(fn_str, dis_fn='w'){
  #return(D(parse(text=eval(parse(text=fn_str))), 't'))
  re = get.expression(fn_str, dis_fn)
  return(D(re$expression, 't'))
}

get.sys.value = function(fn_str, lebensdauer, par1, par2, dis_fn='w'){
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