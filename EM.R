library(grid)
library(png)

# read in data file 
d = read.csv("/Users/caseynold/Desktop/STAT610/Unsupervised/kmeans/training_ss_151_blocks.csv")
#remove the index from the image
d = as.matrix(d[,-1])
#number of samples
n = dim(d)[1]
#sample dimension
p = dim(d)[2]
# number of clusers
k = 2
# create responsibility matrix
ri = matrix(0,nrow=n,ncol=k)
# parameter matrices 
mu = matrix(rnorm(n*1, mean=0.5,sd=0.5),nrow=n,ncol=k ) 
sigma = matrix(rnorm(n*1, mean=0.5,sd=0.5),nrow=n,ncol=k ) 
pi = matrix(rnorm(1*k,mean=0.5,sd=0.5), nrow = n,ncol=k ) 

#stochastic row matrix
#pi = pi/rowSums(pi)
#Repeat
step = 0
repeat{
  step = step+1
  
  #expectation: 1.) compute responsibilities
  for(i in 1:k){
      ri[,i] = pi[,i] %*% dnorm(d,mu[i],sd=sqrt(sigma[i]),log=F)
  }

  ri = ri/rowSums(ri)
  print(r)

  #maximazation
  # update pi
  rk = colSums(ri)
  pi = rk/n  #notes  02/FEB/2017
 
   # update mu
   for(i in (1:n)){
     for(j in 1:p){
       for(m in (1:k)){
         mu[i,m] = (ri[i,m]/rk[m]) * d[i,k]
       }
     }
   }
  
  #update sigma
  for(i in 1:n){
    for(j in 1:p){
      for(m in 1:k){
      sigma[i,m] = (ri[i,m]/rk[m])* (d[i,j] - mu[i,m])^2
      }
    }
  }
  
  if(step == 3){
    break
  }
  write.csv(mu,file="/Users/caseynold/Desktop/STAT610/Unsupervised/kmeans/em_means.csv")
}