library(png)
library(grid)
library(mvtnorm)
library(MASS)
rm(list=ls())


##########

im <- readPNG("./LisaData/001_4_orig.png")
#grid.raster(im)
### getting the picture
im.seg=readPNG("./LisaData/001_4_seg.png")

white = c(1,1,1)
yellow = c(1,1,0) # white matter lesions
red = c(1,0,0) # csf
green = c(0,1,0) # white matter
blue = c(0,0,1) # grey matter
black = c(0,0,0)
im.back=im
im.csf=im
im.grey=im
im.white=im
for (i in (1:dim(im)[1])){
  for (j in (1:dim(im)[2])){
    if(all(im.seg[i,j,]==red)){
      for(c in seq(3)){im.csf[i,j,c]=white[c]}
    }
    else{
      for(c in seq(3)){im.csf[i,j,c]=black[c]}
    }
    if(all(im.seg[i,j,]== yellow) || all(im.seg[i,j,]== green)){
      for(c in seq(3)){im.white[i,j,c]=white[c]}
    }
    else{
      for(c in seq(3)){im.white[i,j,c]=black[c]}
    }
    if(all(im.seg[i,j,]==blue)){
      for(c in seq(3)){im.grey[i,j,c]=white[c]}
    }
    else{
      for(c in seq(3)){im.grey[i,j,c]=black[c]}
    }
  }
}
# grid.raster(im.back)
# grid.raster(im.csf)
# grid.raster(im.grey)
# grid.raster(im.white)


im.total <- im.csf+im.grey+im.white
for (i in (1:dim(im)[1])){
  for (j in (1:dim(im)[2])){
    im.back[i,j,]=ifelse(im.total[i,j,1]==0,c(1,1,1),c(0,0,0)) # just not the others
  }
}
# grid.raster(im.back)
# # visualize the white matter only on the original image
# grid.raster(im*im.white)
##########################
## extract the data for the clustering
##########################

im <- im

n = sum(!im.back[,,1])
coords=matrix(nrow=n, ncol=2)
index = 1
d=matrix(nrow=n,ncol=25)
for (i in (1:dim(im)[1])){
  for (j in (1:dim(im)[2])){
    if (!im.back[i,j,1]){
      d[index,]=as.vector(im[(i-2):(i+2),(j-2):(j+2),1])
      coords[index,]=as.vector(c(i,j))
      index=index+1
    }
  }
}


data <- data.frame(d)
# integer values
data <- data*255
data$membership <- 0
dnum <- as.matrix(data[,1:25])
# initialize means
# picking 3 based on low, med, high
M1 <- rep(1,25)
M2 <- rep(128,25)
M3 <- rep(255,25)


# data is a row of data for use in apply
mem <- function(drow,M1,M2,M3){
  dist<- c(0,0,0)
  dist[1] <- (drow-M1)%*%(drow-M1)
  dist[2] <- (drow-M2)%*%(drow-M2)
  dist[3] <- (drow-M3)%*%(drow-M3)
  return(which.min(dist))
}

# couldn't figure apply

epsilon <- 1200
while(epsilon > 1000){
  # assign membership
  for(i in 1:dim(data)[1]){
    data$membership[i] <- mem(dnum[i,],M1,M2,M3)
  }
  # update centroid
  M1_up <- apply(subset(data,data$membership==1),2,mean)[-26]
  M2_up <- apply(subset(data,data$membership==2),2,mean)[-26]
  M3_up <- apply(subset(data,data$membership==3),2,mean)[-26]
  
  epsilon <- sum(c((M1_up-M1)^2,(M2_up-M2)^2,(M3_up-M3)^2))
  M1 <- M1_up
  M2 <- M2_up
  M3 <- M3_up
}
# two ways 
# b = apply(a,2,mean)
#  c = colMeans(a)

clst <- cbind(coords,data$membership)
clst <- as.data.frame(clst)
names(clst) <- c('x','y','membership')
clst$value <- lapply(clst$membership,function(x){
  if(x==1){
    mean(M1)/255
  } else if(x==2){
    mean(M2)/255
  } else {
    mean(M3)/255
  }
})
im3 <- im
a = 0
for(i in 1:dim(clst)[1]){
  im3[clst$x[i],clst$y[i],] <- clst$value[[i]]
}
#grid.raster(im3)
writePNG(im3,'./Data/001_4_kmeans.png')
#write.csv(clst,'./Data/001_4_kmeans_clst.csv')
Data_001_4_kmeans_clst <- clst
Data_001_4_kmeans_means <- list(M1,M2,M3)

saveRDS(Data_001_4_kmeans_clst,'./Data/001_4_kmeans_clst')
saveRDS(Data_001_4_kmeans_means,'./Data/001_4_kmeans_means')


n <- dim(data)[1]
p <- 25
K <- 3


pi <- matrix(1, ncol = K)/K # max H
mu <- matrix(0, nrow = K, ncol = p)
sig_dim = matrix(0, nrow = p, ncol = p)
sig <- list(sig_dim, sig_dim, sig_dim)

# from k-means and data is d*255
for(k in seq(K)){
  mu[k,] = apply(subset(data,data$membership==k)[,-26],2,mean)
}

for(k in seq(K)){
  sig[[k]] = cov(subset(data,data$membership==k)[,-26])
}

# initialize r (nxk) matrix of responsibilities
r <- matrix(0, nrow = n, ncol = K)



# initial E for initial likelihood calculation
# E-step
for(i in seq(n)){
  denom <- 0
  for(k in seq(K)){
    denom <- denom + pi[k]*dmvnorm(dnum[i,],mu[k,],sig[[k]])
  }
  for(k in seq(K)){
    r[i,k] <- (pi[k]*dmvnorm(dnum[i,],mu[k,],sig[[k]]))/denom
  }
}

# r check
R <- colSums(r)
n == sum(R)
head(rowSums(r))

# log likelihood
# log likelihood of pi mixture
lpi <- function(r, pi){
  l = 0
  for(i in seq(dim(r)[1])){
    for(k in seq(dim(r)[2])){
      l = l + r[i,k]*log(pi[k]) 
    }
  }
  return(l)
}
l_pi <- lpi(r, pi)
# log likelihood of theta
ltheta <- function(r,x,mu,sig){
  l=0
  for(i in seq(dim(r)[1])){
    for(k in seq(dim(r)[2])){
      l = l + r[i,k]*dmvnorm(x[i,], mu[k,], sig[[k]],log = T) 
    }
  }
  return(l)
}

l_th <- ltheta(r,dnum,mu,sig)

likelihood <- function(lpi,ltheta){
  return(sum(ltheta,lpi))
}

l_tot <- likelihood(l_th,l_pi)
epsilon <- 1000
# done with initializations and initial checks

# EM
repeat{
  l_tot_pre <- l_tot
  # E-step
  for(i in seq(n)){
    denom <- 0
    for(k in seq(K)){
      denom <- denom + pi[k]*dmvnorm(dnum[i,],mu[k,],sig[[k]])
    }
    for(k in seq(K)){
      r[i,k] <- (pi[k]*dmvnorm(dnum[i,],mu[k,],sig[[k]]))/denom
    }
  }
  
  # M-step
  # pi update
  R <- colSums(r)
  pi <- R/n
  
  # mu update
  for(k in seq(K)){
    if(R[k] != 0){
      mu[k,] <- colSums(r[,k]*dnum)/R[k]
    }
  }
  
  # sigma update
  for(i in seq(n)){
    for(k in seq(K)){
      sig[[k]] <- sig[[k]] + r[i,k]*dnum[i,]%*%t(dnum[i,])
    }
  } 
  for(k in seq(K)){
    if(R[k] != 0){
      sig[[k]] <- sig[[k]]/R[k] - mu[k,]%*%t(mu[k,])
    }
  }
  l_tot <- likelihood(ltheta(r = r,x = dnum,mu = mu,sig = sig),lpi(r = r,pi = pi))
  print(l_tot)
  if((l_tot-l_tot_pre)^2 < epsilon){break}
}


