library(grid)
library(png)
library(MASS)
library(mvtnorm)
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

# read in data file 
d = read.csv("./BrunoData/training_ss_151_blocks.csv")
means = read.csv("/Users/caseynold/Desktop/STAT610/Unsupervised/kmeans/kmeans.csv")
r = read.csv("/Users/caseynold/Desktop/STAT610/Unsupervised/kmeans/ri.csv")
r = r[-1]
means = means[-1]
#remove the index from the image
d = as.matrix(d[,-1])
#number of samples
n = dim(d)[1]
#sample dimension
p = dim(d)[2]
# number of clusers
k = 3
# create responsibility matrix
ri = matrix(0,nrow=n,ncol=k)
#r = matrix(unlist(r),ncol=k,byrow=TRUE)

# parameter matrices 
mu = matrix(unlist(means),ncol=p,byrow=TRUE)
sigma1 = diag(rnorm(1, mean=0.5,sd=0.5),nrow=p,ncol=p ) 
sigma2 = diag(rnorm(1, mean=0.5,sd=0.5),nrow=p,ncol=p ) 
sigma3 = diag(rnorm(1, mean=0.5,sd=0.5),nrow=p,ncol=p ) 

SIGMA = list(sigma1,sigma2,sigma3)
pi = matrix(c(.33,.34,.33), nrow = k,ncol=1 ) # was p

m = matrix(0,nrow=k,ncol=p)
############# Start Initialization ###################

rk = colSums(r)
# update pi
pi = rk/n  #notes  02/FEB/2017
mu_old = mu
# Initialize mu

for(i in (1:k)){ 
  dv_mu = 0
  for(j in (1:n)){
    dv_mu %+=% ((r[j,i]) %*% d[j,])
  }
  mu[i,] = (dv_mu/rk[i])
}

dv_sig = matrix(0,p,p)
#Initialize sigma
for(i in 1:k){
  for(j in 1:n){
    dv_sig = dv_sig + (r[j,i]* ((d[j,] - mu_old[i,])%*%t(d[j,] - mu_old[i,])))
  }
  SIGMA[[i]] = (dv_sig/rk[i])
}
####### End Initialization ########

#Repeat
step = 0
repeat{
  step = step+1
  
  #expectation: 1.) compute responsibilities
 
  ri_old = ri
  
  for(i in 1:k){
    for(j in 1:n){
      sumProb =(pi[1]*dmvnorm(d[j,],mu[1,],SIGMA[[1]]) +
                pi[2]*dmvnorm(d[j,],mu[2,],SIGMA[[2]]) +
                pi[3]*dmvnorm(d[j,],mu[3,],SIGMA[[3]]) )
      if(sumProb == 0){
        ri[j,i] = ri_old[j,i]
      }else{
        ri[j,i] = (pi[i] * dmvnorm(x =d[i,],mean=mu[i,],sigma=SIGMA[[i]]))/sumProb
      }
    }
  }

  #ri = ri/rowSums(ri)
  #print(ri)

  #maximazation
  # weighted average of cluster membership
  rk = colSums(ri)
  # update pi
  pi = rk/n  #notes  02/FEB/2017
  
  mu_old = mu
  # update mu
  for(i in (1:k)){ 
    dv_mu = 0
    for(j in (1:n)){
      dv_mu %+=% ((ri[j,i]) * d[j,])
    }
   if(rk[i]==0){
      mu[i,] = dv_mu
    }else{
      mu[i,] = (dv_mu/rk[i])
    }
  }
  
  #update sigma
  for(i in 1:k){
    dv_sig = matrix(0,p,p)
    for(j in 1:n){
      dv_sig %+=% ri[j,i]* ((d[j,] - mu[i,])%*%t(d[j,] - mu[i,]))
    }
    if(rk[i] == 0){
      SIGMA[[i]] = 0
    }else{
      SIGMA[[i]] = (dv_sig/rk[i])
    }
  }
  
  if(step == 5){
    break
  }
}

r.new = rep(0,times=n)
for(i in 1:n){
  r.new[i] = which.max(ri[i,])
}

for (j in (1:k)){
  if (sum(r.new==j)>0)
    m[j,]=as.vector(colMeans(d[which(r.new==j),]))
}
 
write.csv(m,file="/Users/caseynold/Desktop/STAT610/Unsupervised/kmeans/em_means.csv")
