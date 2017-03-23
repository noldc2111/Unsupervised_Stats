library(grid)
library(png)

boltzman = function(x,t){
  exp(x/t)/sum(exp(x/t))
}

`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

# read in data file 
d = read.csv("/Users/caseynold/Desktop/STAT610/Unsupervised/kmeans/data/training_ss_151_blocks.csv")
#remove the index from the image
d = as.matrix(d[,-1])
#number of samples
n = dim(d)[1]
#sample dimension
p = dim(d)[2]
# number of clusers
k = 3
T=1000
# create responsibility matrix
ri = matrix(rnorm(k*p, mean=0.5,sd=0.5),nrow=n,ncol=k)
ri_hat = matrix(rnorm(k*p, mean=0.5,sd=0.5),nrow=n,ncol=p)
# parameter matrices 
mu = matrix(rnorm(k*p, mean=0.5,sd=0.5),nrow=k,ncol=p ) 
sigma = matrix(rnorm(k*p, mean=0.5,sd=0.5),nrow=k,ncol=p ) 
pi = matrix(rnorm(k*p,mean=0.5,sd=0.5), nrow = k,ncol=p ) 

#output matrix
m = matrix(0,nrow=k,ncol=p)

#stochastic row matrix
pi = pi/rowSums(pi)


#Repeat
step = 0
repeat{
  step = step+1
  
  temp = ri
  
  #expectation: 1.) compute responsibilities using boltzman function
  for(i in 1:k){
    temp_val = pi[i] * boltzman(rnorm(d*k,mu[i],sd=sqrt(abs(sigma[i]))),T)
    temp[,i] = temp_val/rowSums(temp_val)
  }
  
  # Get stochastic row matrix n x k
  temp = temp/rowSums(temp)
  
  for(i in 1:n){
    ## Check if the new probability is lower than the old
    ## if so use the new probability
    if(all( ri[i,] > temp[i,])){
       ri[i,] = temp[i,]
    }else if(T <= 0.001){
       g = sample(0:1,1)
       if(g > 0.82){break}
      }
  }
  
  # Update Temperature
  T = T*0.8
  print(T)
  
  # Maximazation Step as Usual
  
  # weighted average of cluster membership
  rk = colSums(ri)
 
   # update pi
  pi = rk/n  #notes  02/FEB/2017
  
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
}

r.new = rep(0,times=n)
for(i in 1:n){
  r.new[i] = which.max(ri[i,])
}

for (j in (1:k)){
  if (sum(r.new==j)>0)
    m[j,]=as.vector(colMeans(d[which(r.new==j),]))
}

write.csv(m,file="/Users/caseynold/Desktop/STAT610/Unsupervised/kmeans/annealed_em.csv")

