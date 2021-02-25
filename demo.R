library(ergm)
library(network)
library("xtable")
library("ggplot2")
library(igraph)
library(graphkernels)

setwd("~/Workspace/gksd/mc_gksd/code") #for ubuntu working directory
#setwd("C:/Users/user/Documents/GitHub/gksd/gksd_r/") #for windows working directory
source("Base_fun.R")


path.fig = "~/Workspace/gksd/mc_gksd/fig"
path.res = "~/Workspace/gksd/mc_gksd/res"


# The null model
# E2ST with edge  +  2-star + triangles 
coef.h0 = c(-2, -0.0, 0.01) # the setting in Yang et.al 2018 for KDSD

# coef.h1 is constructed by varying 2-star parameter
# list of parameters of 2-star
coef.2s = seq(-0.5, 0.5, 0.1)


t_fun = function(x){
  Del = coef.h0[1] + coef.h0[2]*del_t2s(x) + coef.h0[3] * del_tri(x)
  S = sigmoid(Del)
}


# construct the network  
d = 20
un=network(d, directed = FALSE)
model0<- un ~ edges + kstar(2) + triangles


###################
# gKSS basedtests #
###################


generate.method = generate.one.GKSS.sampled

#simluate the null distribution
N = 500
B.list = c(50, 100, 200)
sim <- matrix(1, length(B.list), N)
g.sim <- simulate(model0, nsim=N,
                  coef=coef.h0,control=control.simulate(MCMC.burnin=1000+10*d, MCMC.interval=d))

for (j in 1:length(B.list)){
  B = B.list[j]
  idx = sample.int(d^2, size = B, replace = TRUE)
  for (ii in 1:N){
    X.sim = g.sim[[ii]][,]
    method.out=generate.method(t_fun, X.sim,idx)
    sim[j, ii] = method.out[['stats.value']]
    if(ii%%100==0)print(paste("Iter",ii, (sim[j,ii])))
  }
}

sim.null = sim


n=200
alpha = c(0.01, 0.05, 0.1)
l = length(B.list)
l.2s = length(coef.2s)
test.stats.total = array(0, c(l,l.2s,n))
power = array(0, c(l, l.2s,length(alpha)))

for (j in 1:length(B.list)){
  B = B.list[j]
  sim.h0 = matrix(sim[j,])
  for (i in 1:l.2s){
  coef.h1 = c(-2, coef.2s[i], 0.01)
  res = perform.test.sample(generate.method, model0, coef.h1, n, B,sim.h0)
  pvalue = res[['pvalue']]
  test.stats.total[j,i,] = res[['stats']]
  for (c in 1:length(alpha)){
    power[j,i,c] = mean(pvalue<alpha[c])
    print(power[j,i,c])
  }
  }
  print(j)
}

##############################
# competing approaches tests #
##############################


generate.method = generate.degree.stats

#simluate the null distribution

N=1000
sim.degree <- matrix(1, N)


idx = sample.int(d^2, size = B, replace = TRUE)
for (ii in 1:N){
  X.sim = g.sim[[ii]][,]
  sim.degree[ii] = mean(colSums(X.sim))
  if(ii%%100==0)print(paste("Iter",ii, (sim.degree[ii])))
}

mean.degree = mean(sim.degree)

for (ii in 1:N){
  X.sim = g.sim[[ii]][,]
  method.out=generate.degree.stats(t_fun, X.sim,idx, mean = mean.degree)
  sim.degree[ii] = method.out[['stats.value']]
  if(ii%%100==0)print(paste("Iter",ii, (sim.degree[ii])))
}
n=100
power.degree = array(0, c(l.2s,length(alpha)))

sim.h0 = sim.degree
hist(sim.h0,30)

for (i in 1:l.2s){
  coef.h1 = c(-2, coef.2s[i], 0.01)
  g.sim.data <- simulate(model0, nsim=n,
                         coef=coef.h1,control=control.simulate(MCMC.burnin=1000, MCMC.interval=50))
  pvalue = matrix(0,n)
  N = length(sim.h0)
  for (j in 1:n){
    X.sim = g.sim.data[[j]][,]
    stats = generate.degree.stats(t_fun, X.sim, idx, mean=mean.degree)[[1]]
    pvalue[j] = mean(rep(stats, N)<sim.h0)
  }
  for (c in 1:length(alpha)){
    power.degree[i,c] = mean(pvalue<alpha[c])
    print(power.degree[i,c])
  }
}

power.degree


##graphical tests based on degree

generate.method = generate.graphical.stats.degree

#simluate the null distribution
N = 1000

sim.deg <- matrix(1, N)

idx = sample.int(d^2, size = B, replace = TRUE)
for (ii in 1:N){
  method.out=generate.method(coef.h0, g.sim[[ii]],idx)
  sim.deg[ii] = method.out[['stats.value']]
  if(ii%%100==0)print(paste("Iter",ii, (sim[ii])))
}
test.stats.total.deg = array(0, c(l.2s,n))
power.deg = array(0, c(l.2s,length(alpha)))

sim.h0 = matrix(sim.deg)
for (i in 1:l.2s){
  coef.h1 = c(-2, coef.2s[i], 0.01)
  res = perform.test.graphical(generate.method, model0, coef.h1, n, B, sim.h0)
  pvalue = res[['pvalue']]
  test.stats.total[j,i,] = res[['stats']]
  for (c in 1:length(alpha)){
    power.deg[i,c] = mean(pvalue<alpha[c])
    print(power.deg[i,c])
  }
}


##graphical tests based on espart

generate.method = generate.graphical.stats.espart
N = 100
sim.esp <- matrix(1, N)

idx = sample.int(d^2, size = B, replace = TRUE)
for (ii in 1:N){
  method.out=generate.method(coef.h0, g.sim[[ii]],idx)
  sim.esp[ii] = method.out[['stats.value']]
  if(ii%%100==0)print(paste("Iter",ii, (sim.esp[ii])))
}


test.stats.total.esp = array(0, c(l.2s,n))
power.esp = array(0, c(l.2s,length(alpha)))

sim.h0 = matrix(sim.esp)
n=50
for (i in 1:l.2s){
  coef.h1 = c(-2, coef.2s[i], 0.01)
  res = perform.test.graphical(generate.method, model0, coef.h1, n, B, sim.h0)
  pvalue = res[['pvalue']]
  test.stats.total.esp[i,] = res[['stats']]
  for (c in 1:length(alpha)){
    power.esp[i,c] = mean(pvalue<alpha[c])
    print(power.esp[i,c])
  }
}




## mahalanobis distance based tests

N = 1000
sim.md <- matrix(1, N)

n=200
test.stats.total.md = array(0, c(l.2s,n))
power.md = array(0, c(l.2s,length(alpha)))

#compute the null mean/cov
X = g.sim[[1]]
gest = ergm(X~ edges + kstar(2) + triangles, estimate = "MPLE")
gnull <- gof(gest, GOF=~degree, coef=coef.h0)
sim.null = gnull$sim.deg
mu = colMeans(sim.null)
Sig = cov(sim.null)

for (i in 1:l.2s){
  coef.h1 = c(-2, coef.2s[i], 0.01)
  g.sim.data <- simulate(model0, nsim=n,
                         coef=coef.h1,control=control.simulate(MCMC.burnin=1000, MCMC.interval=500))
  pvalue = matrix(0, n)
  for (j in 1:n){
    sim.deg = degreedist(g.sim.data[[j]],print=FALSE)
    l = length(sim.deg)
    Ax = matrix(0, d)
    Ax[1:l] = sim.deg
    stat.md = t(Ax-mu) %*% ginv(Sig) %*%(Ax-mu) /10
    test.stats.total.md[i,j] = stat.md
    pvalue[j] =1- pchisq(stat.md, df=1)
    }
  
  for (c in 1:length(alpha)){
    power.md[i,c] = mean(pvalue<alpha[c])
    print(power.md[i,c])
  }
}


