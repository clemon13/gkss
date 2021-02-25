library(igraph)
library(graphkernels)
library(MASS)

#######################
## Utility functions ##
#######################

sigmoid = function(x)1/(1+exp(-x))

logit = function(x)log(x/(1-x))

#The difference function for perturbing edge in 2-star term
del_t2s = function(x){
  deg = colSums(x)
  Deg_sum = outer(deg,deg,function(x,y)x+y)
  #Del = Deg_sum - x
  diag(Deg_sum) <- 0
  Deg_sum
}


#The difference function for perturbing edge in triangle term
del_tri = function(x){
  X2 = x%*%t(x)
  diag(X2) <- 0
  X2
}


#The difference function for perturbing edge
Delta.fun = function(x, p=0){
  # for E2ST graph
  # Del = (-2 + 0.01*del_tri(x))
  
  #for ER graph
  Del = x*p
}


################################
### Compute kernel statistic ###
################################

compute.transition.list = function(X){
  P=list()
  for (w in 1:length(X)){
    x = X
    x[w] = abs(1-X[w])
    G = graph_from_adjacency_matrix(x)
    P[[w]] = G
  }
  P[[length(X)+1]] = graph_from_adjacency_matrix(X)
  P
}


compute.sampled.list = function(X, sample.index){
  P=list()
  l = length(sample.index)
  for (w in 1:l){
    x = X
    x[sample.index[w]] = abs(1 - X[sample.index[w]])
    G = graph_from_adjacency_matrix(x)
    P[[w]] = G
  }
  P[[l+1]] = graph_from_adjacency_matrix(X)
  P
}


compute.normalized = function(K){
  V = diag(K)
  D.mat = diag(1.0/sqrt(V))
  D.mat%*%K%*%D.mat
}


## Weisfeiler_Lehman Graph Kernel
compute.wl.kernel=function(X, level=3, diag=1, normalize=TRUE){
  n = length(X)
  P = compute.transition.list(X)
  kernel.matrix = CalculateWLKernel(P, level)
  rm(P)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}


## Geometric Random Walk Kernel
compute.grw.kernel=function(X, level=3, diag=1, normalize=TRUE){
  n = length(X)
  P = compute.transition.list(X)
  kernel.matrix = CalculateGeometricRandomWalkKernel(P, level)
  rm(P)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}


## Shortest Path Kernel
compute.sp.kernel=function(X, level=1, diag=1, normalize=TRUE){
  n = length(X)
  P = compute.transition.list(X)
  kernel.matrix = CalculateShortestPathKernel(P)
  rm(P)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}


## Vertex-Edge Histogram Kernel
compute.veh.kernel=function(X, level=0.1, diag=1, normalize=TRUE){
  n = length(X)
  P = compute.transition.list(X)
  kernel.matrix = CalculateVertexEdgeHistGaussKernel(P,level)
  rm(P)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}


# Various special kernels

compute.linear.kernel=function(X,normalize = TRUE,diag=1){
  n = length(X)
  K = matrix(1,ncol=n,nrow = n) 
  if(diag==0)diag(K)<-0
  K
}

compute.constant.kernel=function(X,normalize = TRUE,diag=1){
  n = length(X)
  K = matrix(0,ncol=n,nrow = n) + 1
  if(diag==0)diag(K)<-0
  K
}

compute.flipsign.kernel=function(X,normalize = TRUE,diag=1){
  n = length(X)
  S = matrix(X, byrow =TRUE)*2 - 1
  K = S%*%t(S)
  if(diag==0)diag(K)<-0
  K
}




##################################################################
##### When single graphs are observed ##########################
##################################################################

# GKSS methods


# Using vvRKHS
generate.one.GKSS.condition=function(t_fun, X, kernel=compute.wl.kernel,diagonal=1)
{
  S=matrix(X,byrow=TRUE)
  n.sim=length(S)
  S.t =  t_fun(X)
  S.t.vec=abs(S - matrix(S.t, byrow=TRUE))
  S.mat = S.t.vec %*% t(S.t.vec)
  Ky.mat = (S*2-1)%*%t(S*2-1) 
  
  
  K = kernel(X, diag=diagonal)
  J.kernel = S.mat * Ky.mat * K
  
  W=rep(1/n.sim,n.sim)
  J.kernel.out=J.kernel
  if(diagonal==0)diag(J.kernel)=rep(0,n.sim)
  stats.value=n.sim*t(W)%*%J.kernel%*%W * sqrt(v.kernel)
  #Return:
  #stats.value: n times KSD
  #J.kernel: J kernel matrix for wild bootstrapt 
  list(stats.value=stats.value,J.kernel=J.kernel.out, K=K, S=S.mat)
}


generate.one.GKSS.linear=function(t_fun, X, kernel=compute.wl.kernel,diagonal=1)
{ S=matrix(X,byrow=TRUE)
  n.sim=length(S)
  S.t =  t_fun(X)
  S.t.vec=abs(S - matrix(S.t, byrow=TRUE))
  S.mat = S.t.vec %*% t(S.t.vec)
  Ky.mat = (S*2-1)%*%t(S*2-1) 
  
  J.kernel = S.mat * Ky.mat
    
  W=rep(1/n.sim,n.sim)
  J.kernel.out=J.kernel
  if(diagonal==0)diag(J.kernel)=rep(0,n.sim)
  stats.value=n.sim*t(W)%*%J.kernel%*%W
  #Return:
  #stats.value: n times KSD
  #J.kernel: J kernel matrix for wild bootstrapt 
  list(stats.value=stats.value,J.kernel=J.kernel.out, S=S.mat)
}


# Re-Sample algorithms 
generate.one.GKSS.sampled=function(t_fun, X, sample.index, g.kernel=CalculateWLKernel,level=3, diagonal=1, normalize=TRUE, v.scale=TRUE)
{
  S=matrix(X,byrow=TRUE)[sample.index]
  n.sim=length(S)
  S.t =  t_fun(X)[sample.index]
  S.t.vec=abs(S - S.t)
  S.mat = S.t.vec %*% t(S.t.vec)
  Ky.mat = (S*2-1)%*%t(S*2-1) 
  
  
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  kernel.matrix = g.kernel(P, level)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  
  J.kernel = S.mat * Ky.mat * K
  
  W=rep(1/n,n)
  J.kernel.out=J.kernel
  if(diagonal==0)diag(J.kernel)=rep(0,n.sim)
  v.kernel = var(S)
  stats.value=n*t(W)%*%J.kernel%*%W 
  if(v.scale)nKSD = stats.value* sqrt(v.kernel)
  #Return:
  #stats.value: n times KSD
  #J.kernel: J kernel matrix for wild bootstrapt 
  list(stats.value=stats.value,J.kernel=J.kernel.out, K=K, S=S.mat)
}



# Competing approaches


generate.degree.stats=function(t_fun, X, sample.index, mean=0, kernel=CalculateWLKernel,level=3, diagonal=1, normalize=TRUE, v.scale=TRUE)
{
  Deg = colSums(X)
  n = dim(X)[1]
  stats.value = mean((Deg-mean)^2)
  list(stats.value=stats.value)
}


generate.graphical.stats.degree=function(t_coef, X, sample.index, kernel=CalculateWLKernel,level=3, diagonal=1, normalize=TRUE, v.scale=TRUE)
{
  #X here is an ergm
  gest = ergm(X~ edges + kstar(2) + triangles, estimate = "MPLE")
  # gnull <- gof(gest, GOF=~distance + espartners + degree + dspartners+
  #                triadcensus, coef=t_coef)
  gnull <- gof(gest, GOF=~degree, coef=t_coef)
  obs = gnull$obs.deg
  obs = obs/sum(obs)
  sims = colMeans(gnull$sim.deg)
  sims = sims/sum(sims)
  stats.value = 0.5*sum(abs(obs-sims))
  list(stats.value=stats.value)
}

generate.graphical.stats.distance=function(t_coef, X, sample.index, kernel=CalculateWLKernel,level=3, diagonal=1, normalize=TRUE, v.scale=TRUE)
{
  #X here is an ergm
  gest = ergm(X~ edges + kstar(2) + triangles, estimate = "MPLE")
  # gnull <- gof(gest, GOF=~distance + espartners + degree + dspartners+
  #                triadcensus, coef=t_coef)
  gnull <- gof(gest, GOF=~distance, coef=t_coef)
  obs = gnull$obs.dist
  obs = obs/sum(obs)
  sims = colMeans(gnull$sim.dist)
  sims = sims/sum(sims)
  stats.value = 0.5*sum(abs(obs-sims))
  list(stats.value=stats.value)
}

generate.graphical.stats.espart=function(t_coef, X, sample.index, kernel=CalculateWLKernel,level=3, diagonal=1, normalize=TRUE, v.scale=TRUE)
{
  #X here is an ergm
  gest = ergm(X~ edges + kstar(2) + triangles, estimate = "CD")
  # gnull <- gof(gest, GOF=~distance + espartners + degree + dspartners+
  #                triadcensus, coef=t_coef)
  gnull <- gof(gest, GOF=~espartners, coef=t_coef)
  obs = gnull$obs.espart
  obs = obs/sum(obs)
  sims = colMeans(gnull$sim.espart)
  sims = sims/sum(sims)
  stats.value = 0.5*sum(abs(obs-sims))
  list(stats.value=stats.value)
}

generate.MD.degree=function(t_coef, X)
{
  gest = ergm(X~ edges + kstar(2) + triangles, estimate = "MPLE")
  gnull <- gof(gest, GOF=~degree, coef=t_coef)
  Ax = gnull$obs.deg
  sim.null = gnull$sim.deg
  mu = colMeans(sim.null)
  Sig = cov(sim.null)
  stats.value = t(Ax-mu) %*% ginv(Sig) %*%(Ax-mu) 
  list(stats.value=stats.value)
}

###perform MC-based test####

perform.test = function(generate.method, model.h1, coef.h1, n, sim.h0){
  g.sim.data <- simulate(model.h1, nsim=n,
                         coef=coef.h1,control=control.simulate(MCMC.burnin=1000, MCMC.interval=500))
  N = dim(sim.h0)[1]
  pval.list = matrix(0,n)
  test.stats.list = matrix(0,n)
  for (i in 1:n){
    X = g.sim.data[[i]][,]
    
    ##method
    test.out=generate.method(t_fun, X)
    test.stats = test.out[['stats.value']]; 
    # print(paste("Test stats:",test.stats))
    pval = mean(rep(test.stats, N)<sim.h0)
    # print(pval)
    pval.list[i] = pval
    test.stats.list[i] = test.stats
    if(i%%100==0)print(paste("Iter",i, (pval)))
  }
  list(pvalue=pval.list, stats=test.stats.list)
}

perform.test.sample = function(generate.method, model.h1, coef.h1, n, B, sim.h0){
  g.sim.data <- simulate(model.h1, nsim=n,
                         coef=coef.h1,control=control.simulate(MCMC.burnin=1000, MCMC.interval=50))
  N = dim(sim.h0)[1]
  pval.list = matrix(0,n)
  test.stats.list = matrix(0,n)
  d = dim(g.sim.data[[1]][,])[1]
  idx = sample.int(d^2, size = B, replace = TRUE)
  for (i in 1:n){
    X = g.sim.data[[i]][,]
    ##method
    test.out=generate.method(t_fun, X, idx)
    test.stats = test.out[['stats.value']]; 
    # print(paste("Test stats:",test.stats))
    pval = mean(rep(test.stats, N)<sim.h0)
    # print(pval)
    pval.list[i] = pval
    test.stats.list[i] = test.stats
    if(i%%100==0)print(paste("Iter",i, (pval)))
  }
  list(pvalue=pval.list, stats=test.stats.list)
}


perform.test.graphical = function(generate.method, model.h1, coef.h1, n, B, sim.h0, coef.h0=c(-2, -0.0, 0.01)){
  g.sim.data <- simulate(model.h1, nsim=n,
                         coef=coef.h1,control=control.simulate(MCMC.burnin=1000, MCMC.interval=500))
  N = dim(sim.h0)[1]
  pval.list = matrix(0,n)
  test.stats.list = matrix(0,n)
  d = dim(g.sim.data[[1]][,])[1]
  idx = sample.int(d^2, size = B, replace = TRUE)
  for (i in 1:n){
    X = g.sim.data[[i]]
    
    ##method
    test.out=generate.method(coef.h0, X, idx)
    test.stats = test.out[['stats.value']]; 
    # print(paste("Test stats:",test.stats))
    pval = mean(rep(test.stats, N)<sim.h0)
    # print(pval)
    pval.list[i] = pval
    test.stats.list[i] = test.stats
    if(i%%100==0)print(paste("Iter",i, (pval)))
  }
  list(pvalue=pval.list, stats=test.stats.list)
}

perform.test.size = function(generate.method, model.h1, coef.h1, n, sim.h0, B.list){
  g.sim.data <- simulate(model.h1, nsim=n,
                         coef=coef.h1,control=control.simulate(MCMC.burnin=1000, MCMC.interval=500))
  N = dim(sim.h0)[2]
  l = length(B.list)
  pval.list = matrix(0,l,n)
  test.stats.list = matrix(0,l,n)
  for (i in 1:n){
    X = g.sim.data[[i]][,]
    
    ##method
    test.out=generate.method(t_fun, X)
    J = test.out[['J.kernel']]
    for (w in 1:l){
      B = B.list[w]
      idx = sample.int(d^2, size = B, replace = TRUE)
      Mat = J[idx, idx]
      test.stats = mean(Mat)
      pval = mean(rep(test.stats, N)<sim.h0[w,])
      pval.list[w,i] = pval
      test.stats.list[w,i] = test.stats
    if(i%%100==0)print(paste("Iter",i, (pval)))
    }
  }
  list(pvalue=pval.list, stats=test.stats.list)
}

##################################################################
##### When multiple graphs are observed ##########################
##################################################################


#### Kernel Discrete Stein Discrepancy #######


compute.transition.by.dim = function(G.list, dim){
  P=list()
  n = length(G.list)
  for (w in 1:n){
    x = G.list[[w]][,]
    G = graph_from_adjacency_matrix(x)
    P[[w]] = G
    x[dim] = abs(1-x[dim])
    G = graph_from_adjacency_matrix(x)
    P[[w+n]] = G
  }
  P
}

compute.kxx = function(G.ist, kernel=CalculateWLKernel, par=1){
  d= dim(G.list[[1]][,])[1]; n=length(G.list) 
  K.ten = array(0,c(2*n, 2*n, d^2))
  for (i in 1:d^2){
    P = compute.transition.by.dim(G.list, i)
    K = kernel(P, par)
    K = compute.normalized(K)
    K.ten[,,i] = K
    gc()
  }
  K.ten
}


generate.nKDSD=function(Delta.fun, G.list, kernel=CalculateWLKernel, diagonal=0)
{
  d= dim(G.list[[1]][,])[1]; n=length(G.list) 
  Score.mat = matrix(0, nrow= d^2, ncol = n)
  for(i in 1:n){
    X = G.list[[i]][,]
    flip = X*2 - 1
    Score.mat[,i] = matrix(exp(- Delta.fun(X) * flip),nrow = d^2)[,1]
  }
  
  k.ten = compute.kxx(G.list, kernel = kernel) # 2n x 2n x d^2
  
  J.kernel = matrix(0, ncol = n, nrow = n)
  for (i in 1:d^2){
    K.use = k.ten[,,i]
    S.vec = Score.mat[i,]
    D.s = diag(S.vec, nrow=n)
    S.mat = outer(S.vec, S.vec, function(x,y)x*y)
    term1 = K.use[1:n,1:n]* S.mat
    term2 = K.use[1:n,(n+1):(2*n)] %*% D.s
    term3 = t(term2)
    term4 = K.use[(n+1):(2*n),(n+1):(2*n)] 
    H.mat = term1 - term2 - term3 +term4
    H.mat = H.mat - mean(diag(H.mat))
    J.kernel = H.mat
  }
  
  if (diagonal==0)diag(J.kernel)<-0
  stats.value=n * mean(H.mat)
  #Return:
  #stats.value: n times KSD
  #J.kernel: J kernel matrix for wild bootstrapt 
  list(stats.value=stats.value,J.kernel=J.kernel, k.ten = k.ten, Score.mat = Score.mat)
}

generate.one.GKSD=function(t_fun, X, kernel=compute.wl.kernel, diagonal=1)
{
  S=matrix(X,byrow=TRUE)
  n.sim=length(S)
  S1t0 = t_fun(X)
  
  P.vec = runif(n.sim)
  S1t0.vec=matrix(S1t0, byrow=TRUE)
  S.small = 1 * (P.vec < S1t0.vec)
  Flip = S.small*2 - 1
  
  S.vec <- abs(S - S1t0.vec) * Flip
  # J = outer(S.vec[,1], S.vec[,1], function(x,y)x*y)
  K = kernel(X,  normalize = TRUE, level = 1)
  
  J.W = S.vec %*% t(S.vec) 
  J.kernel = J.W * K
  
  W=rep(1/n.sim,n.sim)
  J.kernel.out=J.kernel
  if(diagonal==0)diag(J.kernel)=rep(0,n.sim)
  stats.value=n.sim*t(W)%*%J.kernel%*%W
  #Return:
  #stats.value: n times KSD
  #J.kernel: J kernel matrix for wild bootstrapt 
  list(stats.value=stats.value,J.kernel=J.kernel.out, K=K,W=J.W)
}

###################
##Wild Bootstrapt##
###################
Wild.Bootstrap=function(N.boot,J.kernel,diagonal=0,Opt=1,p=0.5)
{
  #Opt 1: Normal weights
  #Opt 2: Multinomial weights
  #Opt 3: Rademacher weights
  
  n.sim=dim(J.kernel)[1]
  
  W=rep(1/n.sim,n.sim)
  if(diagonal==0)diag(J.kernel)=rep(0,n.sim)
  nMMD.boot=c()
  
  for(i in 1:N.boot)
  {
    if(Opt==1)xi=rnorm(n.sim,0,1) #for gaussian samples
    else if(Opt==2)xi=rmultinom(1:n.sim,n.sim,rep(1/n.sim,n.sim))-1 #multinomial samples
    else if(Opt==3)
    {
      xi=sample(c(-1,1),n.sim,replace=TRUE,prob=c(p, 1-p)) #rademacher samples with 0.5 prob
    }
    
    W2=W*xi
    nMMD.boot=c(nMMD.boot,n.sim*t(W2)%*%J.kernel%*%W2)
    # nMMD.boot=c(nMMD.boot, t(W2)%*%J.kernel%*%W2)
  }
  nMMD.boot
}

