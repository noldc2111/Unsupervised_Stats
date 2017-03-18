library(grid)
library(png)

`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

# read in data file 
d = read.csv("/Users/caseynold/Desktop/STAT610/Unsupervised/kmeans/training_ss_151_blocks.csv")
#remove the index from the image
d = as.matrix(d[,-1])
#number of samples
n = dim(d)[1]
#sample dimension
p = dim(d)[2]
# number of clusers
k = 3
# create responsibility matrix
# may need to start with some initial values
ri = matrix(rnorm(k*p, mean=0.5,sd=0.5),nrow=n,ncol=k)
# parameter matrices 

mu = matrix(rnorm(k*p, mean=0.5,sd=0.5),nrow=k,ncol=p ) 
sigma = matrix(rnorm(k*p, mean=0.5,sd=0.5),nrow=k,ncol=p ) 
pi = matrix(rnorm(k*p,mean=0.5,sd=0.5), nrow = k,ncol=p ) 

m = matrix(0,nrow=k,ncol=p)
#stochastic row matrix
pi = pi/rowSums(pi)


#Repeat
step = 0
repeat{
  step = step+1
  
  #expectation: 1.) compute responsibilities
  for(i in 1:k){
    #ri[,i] = pi[i] %*% dnorm(d[i,],mu[i],sd=sqrt(abs(sigma[i])),log=F)
    ri[,i] = pi[i] %*% rnorm(d[i,],mu[i],sd=sigma[i],log=F)
  }

  ri = ri/rowSums(ri)
  #print(ri)

  #maximazation
  # weighted average of cluster membership
  rk = colSums(ri)
  # update pi
  pi = rk/n  #notes  02/FEB/2017
  
  m_t = matrix(0,nrow=n,ncol=p)
   # update mu
  dv_mu = 0
  for(i in (1:k)){ 
    for(j in (1:n)){
      dv_mu %+=% (ri[j,i]) * d[j,]
    }
    mu[i] = (dv_mu/rk[i])
   }
  
  dv_sig = 0
  #update sigma
  for(i in 1:k){
    for(j in 1:n){
      dv_sig %+=% ri[j,i]* (d[j,] - mu[i,])^2
    }
    sigma[i] = (dv_sig/rk[i])
  }
  
  if(step == 50){
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
