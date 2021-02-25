library(ergm)
library(network)
library("xtable")
library("ggplot2")
library(igraph)
library(graphkernels)

# re-direct to the desired directory
setwd("~/Workspace/gksd/mc_gksd/code") #for ubuntu working directory
#setwd("C:/Users/user/Documents/GitHub/gksd/gksd_r/") #for windows working directory
source("Base_fun.R")

# Figure path and result path
path.fig = "~/Workspace/gksd/mc_gksd/fig"
path.res = "~/Workspace/gksd/mc_gksd/res"

#########################################
### Lazega laywer connection network ####
#########################################

library(Bergm)
g.laywer <- lazega
plot(g.laywer)

summary(g.laywer)
X = g.laywer[,]

#estimate E2ST
gest = ergm(g.laywer~edges + kstar(2) + triangles, directed=FALSE, estimate = "MPLE")


p = 2*sum(X)/(35*36)
a.star = logit(p)


# coef.h0 = c(a.star,0,0) #a* for E2ST

coef.h0 = c(-2.8547114,  -0.0002635,   0.6882067  )#E2ST

# coef.h0 = c(-2.75,  0.00, 0.0) #ER

t_fun = function(x){
  Del = coef.h0[1] + coef.h0[2]*del_t2s(x) + coef.h0[3] * del_tri(x)
  S = sigmoid(Del)
}

#simulate the null
d = 36
un=network(d, directed = FALSE)
model0<- un ~ edges+ kstar(2) + triangles

generate.method = generate.one.GKSS.sampled

#sample index 
idx = sample.int(d^2, size = 200, replace = TRUE)

#simulated statistics from the null
N = 100
sim <- matrix(1,N)
g.sim <- simulate(model0, nsim=N,
                  coef=coef.h0,control=control.simulate(MCMC.burnin=1000+10*d, MCMC.interval=d))

for (ii in 1:N){
  X.sim = g.sim[[ii]][,]
  method.out=generate.method(t_fun, X.sim,idx)
  sim[ii] = method.out[['stats.value']]
  if(ii%%100==0)print(paste("Iter",ii, (sim[ii])))
}

res = generate.method(t_fun, X, idx)
stats.laywer = res[['stats.value']]
K.laywer = res[['J.kernel']]

mean(rep(stats.laywer,N) < sim[1:N])

pval = mean(rep(stats.laywer, N)<sim.teen)
pval 




#####################################
### Teenager friendship network #####
#####################################
library("RSiena")
g.teen <- s501
X = g.teen
p.teen=sum(X)/(25*49)
a.star.teen = logit(p.teen)
p.teen
a.star.teen

d=50
N=200
un.teen=network(d, directed = FALSE)
model0.teen<- un.teen ~ edges
sim.teen <- matrix(1,N)
g.sim.teen <- simulate(model0.teen, nsim=N,
                  coef=a.star.teen,control=control.simulate(MCMC.burnin=1000+10*d, MCMC.interval=d))

coef.h0 = c(a.star.teen,0,0)
idx =  sample.int(d^2, size = 300, replace = TRUE)
for (ii in 1:N){
  X.sim = g.sim.teen[[ii]][,]
  method.out=generate.one.GKSS.sampled(t_fun, X.sim, idx)
  sim.teen[ii] = method.out[['stats.value']]
  if(ii%%20==0)print(paste("Iter",ii, (sim[ii])))
}

res.teen = generate.one.GKSS.sampled(t_fun, X,idx)
stats.teen = res.teen[['stats.value']]#/sqrt(50)
K.teen = res.teen[['J.kernel']]
stats.teen

mean(rep(stats.teen,N) < sim.teen[1:N])


####################################
##### Senate Cosponsor Network #####
####################################

dat.senate = read.delim("./108_senmatrix.txt", sep = ",", header=FALSE)
dim(dat.senate)
dat.senate
Co.ind.mat = data.matrix(dat.senate[,1:7804] ==2)
Co.ind = (colSums(Co.ind.mat)>1)
sum(Co.ind)
senate.mat = data.matrix(dat.senate[,Co.ind])
d = dim(senate.mat)[2]
random.draw = 0.5
senate.cor =  cor(senate.mat)
senate.adjacency = 1*(senate.cor>random.draw)
diag(senate.adjacency) = 0
X = senate.adjacency

a.star = logit(sum(X)/(d*(d-1)))
a.star
coef.h0 = c(a.star,0,0)

dat.csv = read.csv("./senate.csv")
party.108 = dat.csv[dat.csv[,1] == 108,6]
length(party.108)

del_k.alternate = function(x, lamb=0.4975){
  n = dim(x)[1]
  coef = (1-1/lamb)^seq(1,n-1)
  X2 = x%*%t(x)
  diag(X2) <- 0
  X2
}

t_fun = function(x){
  Del = coef.h0[1] + coef.h0[2]*del_t2s(x) + coef.h0[3] * del_tri(x)
  S = sigmoid(Del)
}


idx = sample.int(d^2, size = 200, replace = TRUE)
res = generate.one.GKSS.sampled(t_fun, X, idx)
stats.cs = res[['stats.value']]



d = dim(X)[1]
un=network(d, directed = FALSE)
model0<- un ~ edges + kstar(2) + triangles

#simluate the null distribution
N = 500
sim.cs <- matrix(1, N)
a.star
g.sim <- simulate(model0, nsim=N,
                  coef=coef.h0,control=control.simulate(MCMC.burnin=1100, MCMC.interval=d))

for (ii in 1:N){
  X.sim = g.sim[[ii]][,]
  method.out=generate.one.GKSS.linear.sampled(t_fun, X.sim, idx)
  sim.cs[ii] = method.out[['stats.value']]
  if(ii%%10==0)print(paste("Iter",ii, (sim.cs[ii])))
}

pval = mean(rep(stats.cs,N) < sim.cs[1:N])


#########################################
###### Florentine Marriage Network ######
#########################################
data(package='ergm') # tells us the datasets in our packages
data(florentine)
X = flomarriage[,]


g.flo = network(X, directed = FALSE)
gest = ergm(g.flo~edges + kstar(2) + triangles, directed=FALSE, estimate = "MPLE")
gest = ergm(g.flo~edges, directed=FALSE)
gest

p = sum(X)/(15*16)
a.star = logit(p)
a.star 

# coef.h0 = c(a.star,0,0)#ER
coef.h0 = c(-1.62319,  -0.01884,   0.24593) #E2ST

d = 16
un=network(d, directed = FALSE)
model0<- un ~ edges+ kstar(2) + triangles

N = 100
sim <- matrix(1,N)
g.sim <- simulate(model0, nsim=N,
                  coef=coef.h0,control=control.simulate(MCMC.burnin=1000+10*d, MCMC.interval=d))

ed=rep(0,N)
for (i in 1:N){X = g.sim[[i]][,]
ed[i] = sum(X)/(d*(d-1))
}


idx = sample.int(d^2, size = 200, replace = TRUE)

t_fun = function(x){
  # Del = (-2 + 0.01*del_tri(x))
  Del = coef.h0[1] + coef.h0[2]*del_t2s(x) + coef.h0[3] * del_tri(x)
  S = sigmoid(Del)
}


generate.method = generate.one.GKSS.sampled

idx = sample.int(d^2, size = 50, replace = TRUE)
# number of simulated null samples
N = 100
sim <- matrix(1,N)
g.sim <- simulate(model0, nsim=N,
                  coef=coef.h0,control=control.simulate(MCMC.burnin=1000+10*d, MCMC.interval=d))

for (ii in 1:N){
  X.sim = g.sim[[ii]][,]
  method.out=generate.method(t_fun, X.sim,idx)
  sim[ii] = method.out[['stats.value']]
  if(ii%%100==0)print(paste("Iter",ii, (sim[ii])))
}

res = generate.method(t_fun, X, idx)
stats.flo = res[['stats.value']]
K.flo = res[['J.kernel']]

mean(rep(stats.flo,N) < sim[1:N])
















