# non-symmetric correspondence analysis (nsca)
# regularized nsca with linear restictions on the row points (reg.nsca)
library(MASS) 
library(multichull) # used in path.nsca
library(ggplot2) # used in path.nsca and boot.nsca
library(egg) # used in path.nsca
library(nnet) # used in boot.nsca

nsca = function(P, X = NULL, S = 2, maxiter = 65536, dcrit = 1e-6){
  # if frequency data is given
  if(sum(P) != 1){P = P/sum(P)}
  
  R = nrow(P)
  C = ncol(P)
  Dr = diag(rowSums(P))
  Dc = diag(colSums(P))
  ones.r = matrix(1,R,1)
  ones.c = matrix(1,C,1)
  
  PI = solve(Dr) %*% P - ones.r %*% t(ones.c) %*% Dc
  
  # GSVD of A
  if(is.null(X)){
    udv = svd(sqrt(Dr) %*% PI, nu = S, nv = S)
    U = solve(sqrt(Dr)) %*% udv$u
    V = udv$v
    d = udv$d
  }
  else{
    # not implemented yet
    udv = svd(sqrt(Dr) %*% PI, nu = S, nv = S)
    U = solve(sqrt(Dr)) %*% udv$u
    V = udv$v
    d = udv$d
  }
  
  #D = ifelse(S == 1, matrix(udv$d[1], 1, 1), diag(d[1:S], nrow = S, ncol = S))
  D = diag(d[1:S], nrow = S, ncol = S)
  
  # loss
  loss = sum( diag( t(PI - U %*% D %*% t(V)) %*% Dr %*% (PI - U %*% D %*% t(V)) ))
  
  # 
  output = list(
    P = P,
    Dr = Dr,
    Dc = Dc,
    S = S,
    PI = PI,
    U = U,
    V = V,
    d = d,
    UD = U %*% D,
    loss = loss
  )
  return(output)  
}

nsca.it = function(P, X = NULL, S = 2, maxiter = 65536, trace = TRUE, dcrit = 1e-6){
  # if frequency data is given
  if(sum(P) != 1){P = P/sum(P)}
  
  R = nrow(P)
  C = ncol(P)
  Dr = diag(rowSums(P))
  Dc = diag(colSums(P))
  ones.r = matrix(1,R,1)
  ones.c = matrix(1,C,1)
  
  PI = solve(Dr) %*% P - ones.r %*% t(ones.c) %*% Dc
  
  # random starting value for V
  V = svd(matrix(rnorm(C*S), C, S), nu = S, nv = S)$u
  
  if(is.null(X)){
    U = PI %*% V
  }
  else{
    B = solve(t(X) %*% Dr %*% X) %*% t(X) %*% Dr %*% PI %*% V
    U = X %*% B
  }

  loss.old = sum( diag( t(PI - U %*% t(V)) %*% Dr %*% (PI - U %*% t(V)) ))
  
  for (iter in 1:maxiter){
    # update U: U = A %*% V %*% solve(t(V) %*% V) = A %*% V - because of orthonormality
    # U = A %*% V
    if(is.null(X)){
      U = PI %*% V
      B = NULL
    }
    else{
      B = solve(t(X) %*% Dr %*% X) %*% t(X) %*% Dr %*% PI %*% V
      U = X %*% B
    }
    
    # update V
    udv = svd( t(PI) %*% Dr %*% U, nu = S, nv = S)
    V = udv$u %*% t(udv$v)    
    
    # loss
    loss.new = sum( diag( t(PI - U %*% t(V)) %*% Dr %*% (PI - U %*% t(V)) ))
    
    if (trace) {cat(iter, loss.old, loss.new, (2*(loss.old - loss.new)/(loss.old + loss.new)), "\n")}
    if ( (2*(loss.old - loss.new)/(loss.old + loss.new)) < dcrit ) break
    if (iter == maxiter) warning("Maximum number of iterations reached - not converged (yet)")
    loss.old = loss.new
    
  }
  # due to random starts solutions differ
  if(is.null(X)){
    udv = svd(U %*% t(V))
    U = udv$u %*% diag(udv$d)
    V = udv$v
  }
  else{
    udv = svd(B %*% t(V))
    B = udv$u %*% diag(udv$d)
    U = X %*% B
    V = udv$v
  }

  output = list(
    P = P,
    Dr = Dr,
    Dc = Dc,
    S = S,
    PI = PI,
    U = U,
    B = B,
    V = V,
    loss = loss.new,
    iter = iter
  )
  return(output)  
}

reg.nsca = function(P, X, S = 2, penalties = c(0,0,0), start = NULL, maxiter = 65536, trace = TRUE, dcrit = 1e-6){
  # Minimizes L(B,V) = || PI - XBV'||_{Dr,I} ^2, s.t., V'V = I 
  # and penalties on B (lasso, ridge, group.lasso)
  # ------------------------------------------------------------------
  
  # if frequency data is given
  if(sum(P) != 1){N = sum(P); P = P/sum(P)}
  else{N = NULL}
  
  # some constant
  eps = 1e-10
  
  R = nrow(P)
  C = ncol(P)
  Dr = diag(rowSums(P))
  Dc = diag(colSums(P))
  ones.r = matrix(1,R,1)
  ones.c = matrix(1,C,1)
  
  # response matrix
  PI = solve(Dr) %*% P - ones.r %*% t(ones.c) %*% Dc
  vpi = Vec(PI)
  
  # AA = kronecker(diag(S), t(X) %*% X)
  H = kronecker(diag(S), t(X) %*% Dr %*% X)
  IDr = kronecker(diag(C), Dr)
  

  if(is.null(start)){
    # smart starting value for V - using nsca on PI
    V = svd(sqrt(Dr) %*% PI, nu = S, nv = S)$v
    B = ginv( t(X) %*% Dr %*% X ) %*% t(X) %*% Dr %*% PI %*% V
  }
  else{
    B = start$B
    V = start$V
  }
  
  U = X %*% B
  theta = X %*% B %*% t(V)
  
  # starting loss value
  pen1 = penalties[1] * sum(abs(B))
  pen2 = penalties[2] * sum(B^2)
  pen3 = penalties[3] * sum(sqrt(rowSums(B^2))) 
  penalty = pen1 + pen2 + pen3
  r = PI - theta # residuals
  loss.old = sum( diag( t(r) %*% Dr %*% r )) + penalty
  
  for (iter in 1:maxiter){
    
    # ----------------------------------------------------
    # update B
    # ----------------------------------------------------
    bold = Vec(B)
    absb = abs(bold) # necessary for lasso
    sqrtb2 = rep( sqrt(rowSums(B^2)), S) 
    # D1 = diag(drop((penalties[1]/2) / absb)) # lasso
    D1 = penalties[1]/2 * diag(1 / drop(absb)) # lasso
    D2 = penalties[2] * diag(ncol(X)*S) # ridge 
    # D3 = diag(drop((penalties[3]/2) / sqrtb2)) # group lasso
    D3 = penalties[3]/2 * diag(1 / drop(sqrtb2)) # group lasso
    Dtot = D1 + D2 + D3 # combined
    
    # update 
    if(any(Dtot > (1/eps))){
      idx = which(Dtot > (1/eps))
      Dtot[idx] = 1/eps # upper bound on elements of D
      warning = TRUE
    }
    
    b = solve(H + Dtot, t(kronecker(V, X)) %*% IDr %*% vpi)
    if(any(Dtot > (1/eps))) b[idx] = bold[idx]
    B = matrix(b, ncol = S)
    
    U = X %*% B
    
    # ----------------------------------------------------
    # update V
    # ----------------------------------------------------

    # update V
    udv = svd( t(PI) %*% Dr %*% U, nu = S, nv = S)
    V = udv$u %*% t(udv$v)    

    theta = X %*% B %*% t(V)
    
    # compute loss
    pen1 = penalties[1] * sum(abs(B))
    pen2 = penalties[2] * sum(B^2)
    pen3 = penalties[3] * sum(sqrt(rowSums(B^2))) 
    penalty = pen1 + pen2 + pen3
    
    r = PI - theta # residuals

    loss.new = sum( diag( t(r) %*% Dr %*% r )) + penalty
    
    if (trace) {cat(iter, loss.old, loss.new, (2*(loss.old - loss.new)/(loss.old + loss.new)), "\n")}
    if ( (2*(loss.old - loss.new)/(loss.old + loss.new)) < dcrit ) break
    if (iter == maxiter) warning("Maximum number of iterations reached - not converged (yet)")
    loss.old = loss.new
  }
  
  SS = sum( diag( t(r) %*% Dr %*% r ))
  
  output = list(
    N = N,
    P = P,
    X = X,
    penalties = penalties,
    Dr = Dr,
    S = S,
    PI = PI,
    B = B,
    V = V,
    U = X %*% B,
    iter = iter,
    pen.loss = loss.new,
    SS = SS,
    penalty = penalty
  )
  return(output)
}

Vec = function(X){
  vecx = matrix(X, ncol = 1)
  return(vecx)
}

boot.pnsca = function(N, X, B = 200, S = 2, lambda.seq, ptype = 1, trace = FALSE){
  # N: contingency table
  # X: external info
  # B: number of bootstraps
  # S: dimensionality
  # lambda.seq: sequence of penalty parameters
  # ptype: 1 = lasso; 2 = ridge; 3 = group.lasso
  
  # make data frame from contingency table
  n = sum(N)
  df <-as.data.frame(matrix( nrow = n, ncol = 2 ))
  id = 0
  for(i in 1:nrow(N)){
    for(j in 1:ncol(N)){
      if(N[i,j] > 0){
        df[(id+1):(id + N[i,j]), ] = outer(rep(1, N[i,j] ), c(i, j))
        id = id + N[i, j]
      }
    }
  }
  
  # vecN = matrix(N, ncol = 1)
  # vecP = vecN/n
  
  if(ptype == 1){
    penalties = cbind(lambda.seq, 0, 0)
  }
  else if(ptype ==2){
    penalties = cbind(0, lambda.seq, 0)
  }
  else if(ptype ==3){
    penalties = cbind(0, 0, lambda.seq)
  }
  
  PE = SS = COM = matrix(NA, B, length(lambda.seq))
  
  for(b in 1:B){
    if(trace){cat("This is bootstrap sample", b, "from", B, "\n")}
    # bootstrap sample from df
    id = sample(1:n, n, replace = TRUE)
    
    # make contingency tables for training and test set
    train = t(class.ind(df[id, 1])) %*% class.ind(df[id, 2])
    test = t(class.ind(df[-id, 1])) %*% class.ind(df[-id, 2])
    
    while(!all(dim(train) == dim(test))){
      # new bootstrap sample
      id = sample(1:n, n, replace = TRUE)
      # make contingency tables for training and test set
      train = t(class.ind(df[id, 1])) %*% class.ind(df[id, 2])
      test = t(class.ind(df[-id, 1])) %*% class.ind(df[-id, 2])
    }

    
    # train = matrix(rmultinom(n = 1, size = n, prob = vecP), ncol = ncol(N))
    # test = pmax(N - train, 0)
    # 
    
    for(l in 1:length(lambda.seq)){
      # fit on training
      if(l ==1){
        out = reg.nsca(train, X = X, S = S, penalties = penalties[l, ], trace = FALSE)
        B0 = out$B
      }
      else{
        ini = list(B = out$B, V = out$V)
        out = reg.nsca(train, X = X, S = S, penalties = penalties[l, ], start = ini, trace = FALSE)
      }
      
      # HERE we can do something with convex hull
      # SS[b,l] = out$pen.loss
      SS[b,l] = out$SS
      # 
      if(ptype == 1) COM[b,l] = sum(abs(out$B)) / sum(abs(B0))
      if(ptype == 2) COM[b,l] = sum(out$B^2) / sum(B0^2)
      if(ptype == 3) COM[b,l] = sum(sqrt(rowSums(out$B^2))) / sum(sqrt(rowSums(B0^2))) 

      # test set
      P = test/sum(test)
      Dr = diag(rowSums(P))
      Dc = diag(colSums(P))
      PI = solve(Dr) %*% P - matrix(1,nrow(P),ncol(P)) %*% Dc
      
      # estimated values
      theta = out$X %*% out$B %*% t(out$V)
      r = PI - theta 
      # out-of-sample predictions
      PE[b, l] = sum( diag( t(r) %*% Dr %*% r ))
    }
  }
  
  Complexity = colMeans(COM)
  SSE = colMeans(SS)
  complexity.fit = cbind(Complexity, SSE)
  rownames(complexity.fit) = paste0("m", 1:length(lambda.seq))
  chull.sse = CHull(complexity.fit, bound = "lower")
  # plot(chull.sse)
  
  PSE = colMeans(PE)
  complexity.fit = cbind(Complexity, PSE)
  rownames(complexity.fit) = paste0("m", 1:length(lambda.seq))
  chull.pse = CHull(complexity.fit, bound = "lower")
  # plot(chull.pse)
  
  PEdf = data.frame(penalty = lambda.seq, 
                    m.error = colMeans(PE),
                    se.error = apply(PE, 2, sd))
  
  lambda.min = PEdf[which.min(PEdf$m.error), 1]
  lambda.1se = PEdf[max(which((PEdf$m.error - PEdf$se.error) < min(PEdf$m.error))), 1]  
  cat("The minimum cross validated error is attained at", lambda.min, "\n")
  cat("The 1 SE rule is attained at", lambda.1se, "\n")
  
  plt = ggplot(PEdf, aes(x = penalty, y = m.error)) + 
    geom_point() + 
    geom_errorbar(aes(ymax = m.error + se.error, ymin = m.error - se.error)) +
    geom_vline(xintercept = lambda.min, colour = "darkblue", linetype = "dotted") + 
    geom_vline(xintercept = lambda.1se, colour = "darkblue") + 
    # geom_hline(yintercept = min(PEdf$m.error), colour = "red") +
    xlab("Penalty") + 
    ylab("Out of Bag Prediction Error") 
  
  print(plt)
  
  
  # complexity = colMeans(COM/outer(COM[, 1], rep(1, length(lambda.seq)))) # norm(B) / norm(B:OLS)
  # complexity = lambda.seq
  # complexity.fit = cbind(complexity, colMeans(PE))
  # chull.out = CHull(complexity.fit, bound = "upper")
  # chull.out = CHull(complexity.fit, bound = "lower")
  # plot(chull.out, col = "green", pch = ".")
  
  output = list(
    PE = PE,
    PEdf = PEdf,
    lamda.seq = lambda.seq,
    lambda.min = lambda.min,
    lambda.1se = lambda.1se,
    plot = plt,
    SS = SS,
    COM = COM,
    CHss = chull.sse,
    CHps = chull.pse
  )
  return(output)
}

path.pnsca = function(N, X, S = 2, ptype = 1){
  # makes a path of the regression weights B for increasing values of penalty
  # also determines max value of penalty
  # !! strategy does not work well for ridge penalty
  PP = N/sum(N)
  
  lambda.seq = seq(0, 100, by = 0.001)
  
  if(ptype == 1){
    penalties = cbind(lambda.seq, 0, 0)
  }
  else if(ptype == 2){
    penalties = cbind(0, lambda.seq, 0)
  }
  else if(ptype == 3){
    penalties = cbind(0, 0, lambda.seq)
  }
  
  Blist = list()
  Vlist = list()
  SSE = rep(NA, length(lambda.seq))
  PSE = rep(NA, length(lambda.seq))
  Complexity = rep(NA, length(lambda.seq))
  
  for(l in 1:length(lambda.seq)){
    # fit on training
    if(l == 1){
      out = reg.nsca(PP, X = X, S = S, penalties = penalties[l, ], trace = FALSE)
      B0 = out$B
    }
    else{
      ini = list(B = out$B, V = out$V)
      out = reg.nsca(PP, X = X, S = S, penalties = penalties[l, ], start = ini, trace = FALSE)
    }
    Blist[[l]] = out$B
    Vlist[[l]] = out$V
    
    SSE[l] = out$SS
    PSE[l] = out$pen.loss
    if(ptype == 1) Complexity[l] = sum(abs(out$B)) / sum(abs(B0))
    if(ptype == 2) Complexity[l] = sum(out$B^2) / sum(B0^2)
    if(ptype == 3) Complexity[l] = sum(sqrt(rowSums(out$B^2))) / sum(sqrt(rowSums(B0^2))) 

   if( all(abs(out$B) < 1e-6) ) break 
  }
  
  # plot
  lambda.seq = lambda.seq[1:l]
  P = ncol(X)
  L = l
  
  SSE = SSE[1:L]
  Complexity = Complexity[1:L]
  # Complexity = max(lambda.seq) - lambda.seq
  complexity.fit = cbind(Complexity, SSE)
  
  rownames(complexity.fit) = paste0("m", 1:L)
  chull.out = CHull(complexity.fit, bound = "lower")
  # plot(chull.out)
  
  
  Bslong = matrix(NA, (P*L), S + 2); 
  Bslong[ , 1] = rep(lambda.seq, each = P)
  Bslong[ , 2] = rep(1:P, L)
  # if(is.null(colnames(X))){Bslong[ , 2] = rep(1:P, L)}
  # else{Bslong[ , 2] = rep(colnames(X), L) }
  for(l in 1:L){
    Bslong[((l-1)*P + 1):(l*P), 3:(2 + S)] = Blist[[l]]
  }
  Bslong = as.data.frame(Bslong)
  colnames(Bslong) = c("Penalty", "Variable", paste("dim", 1:S, sep = ""))
  Bslong$Variable = factor(Bslong$Variable, 
                           levels = 1:P, 
                           labels = colnames(X))
  
  plts = list()
  for(s in 1:S){
    if(s == 1){
      plt = ggplot(Bslong, aes(x = Penalty, 
                               y = dim1, 
                               colour = Variable)) + geom_line() 
      
      plt = plt + geom_label_repel(data = Bslong[1:P, ], 
                            aes(x = 0, y = dim1, label = Variable, hjust = 0), size = 1.5, fontface = "bold") 
                            # colour = "black", check_overlap = TRUE)

      plt = plt  + 
        labs(x = "Penalty",
             y = "Estimates",
             title = paste("Dimension", s))
      
      plt = plt + theme(legend.position = "none")
#      plt = plt + theme_set(theme_bw()) + theme_update(legend.position = "none")
#      plt = plt + theme_minimal()
      #if(s < S){plt = plt + theme(legend.position = "none")}

    }
    else if(s == 2){
      plt = ggplot(Bslong, aes(x = Penalty, 
                               y = dim2, 
                               colour = Variable)) + geom_line()      
      
      plt = plt + geom_label_repel(data = Bslong[1:P, ], 
                                   aes(x = 0, y = dim2, label = Variable, hjust = 0), size = 1.5, fontface = "bold") 

      plt = plt  + 
        labs(x = "Penalty",
             y = "Estimates",
             title = paste("Dimension", s))
      
      plt = plt + theme(legend.position = "none")
      # if(s < S){plt = plt + theme(legend.position = "none")}
    }
    else if(s == 3){
      plt = ggplot(Bslong, aes(x = Penalty, 
                               y = dim3, 
                               colour = Variable)) + geom_line()      
      
      plt = plt + geom_label_repel(data = Bslong[1:P, ], 
                                   aes(x = 0, y = dim3, label = Variable, hjust = 0), size = 1.5, fontface = "bold") 
      
      plt = plt  + 
        labs(x = "Penalty",
             y = "Estimates",
             title = paste("Dimension", s))
      plt = plt + theme(legend.position = "none")
    }
    
    if(s > 3){cat("Plots for Dimensionality > 3, not implemented yet")}
    plts[[s]] = plt
  }
  
  myplot = egg::ggarrange(plots = plts, ncol = S)

  # create output
  output = list(
    Bslong = Bslong,
    Bs = Blist,
    Vs = Vlist,
    lambda.seq = lambda.seq,
    figures = plts,
    plot = myplot,
    SSE = SSE,
    PSE = PSE,
    Complexity = Complexity,
    complexity.fit = complexity.fit,
    CH = chull.out
  )
  return(output)
}

boot.regnsca = function(object, B = 200, trace = TRUE){
  # get inpt from object
  N = round(object$N * object$P)
  X = object$X
  S = object$S
  penalties = object$penalties
  start = list(B = object$B, V = object$V)

  # make data frame from contingency table
  n = sum(N)
  df <-as.data.frame(matrix( nrow = object$N, ncol = 2 ))
  id = 0
  for(i in 1:nrow(N)){
    for(j in 1:ncol(N)){
      if(N[i,j] > 0){
        # print(c(i,j, id, N[i,j]))
        # if((i == 6) && (j = 1) ){
        #   df[(id+1):(id + N[i,j]), ] = outer(rep(1, 15 ), c(i, j))
        # }
        # else{
          df[(id+1):(id + N[i,j]), ] = outer(rep(1, N[i,j] ), c(i, j))
        # }
        id = id + N[i, j]
      }
    }
  }
  
  # PE = SS = COM = matrix(NA, B, length(lambda.seq))
  B.boot = array(NA, c(ncol(X), S, B))
    
  for(b in 1:B){
    if(trace){cat("This is bootstrap sample", b, "from", B, "\n")}
    # bootstrap sample from df
    id = sample(1:n, n, replace = TRUE)
    
    # make contingency tables for training and test set
    train = t(class.ind(df[id, 1])) %*% class.ind(df[id, 2])
    test = t(class.ind(df[-id, 1])) %*% class.ind(df[-id, 2])
    
    while(!all(dim(train) == dim(test))){
      # new bootstrap sample
      id = sample(1:n, n, replace = TRUE)
      # make contingency tables for training and test set
      train = t(class.ind(df[id, 1])) %*% class.ind(df[id, 2])
      test = t(class.ind(df[-id, 1])) %*% class.ind(df[-id, 2])
    }
    
    out = reg.nsca(train, X = X, S = S, penalties = penalties, start = start, trace = FALSE)
    B.boot[ , , b] = out$B
  }
  
  B.ave = apply(B.boot, c(1,2), mean)
  
  st.num = 0
  st.denom = 0
  for(b in 1:B){
    st.num = st.num + ssq(B.boot[, , b] - B.ave)
    st.denom = st.denom + ssq(B.boot)
  }
  
  ST = 1 - st.num/st.denom

  output = list(
    B.boot  = B.boot,
    B.ave = B.ave,
    ST = ST
  )
  return(output)
}

ssq = function(X){
  ssq = sum(X^2)
  return(ssq)
}