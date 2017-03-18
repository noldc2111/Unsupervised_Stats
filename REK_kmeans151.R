# K means on 151.csv from GetClusteringData.R

# cleanup
#rm(list=ls())
library(png)
library(grid)
setwd("/Users/kramerPro/Google Drive/Graphical models stat II/code")

# data
data <- read.csv("training_ss_151_blocks.csv")

# first dim is a repeat of index
data <- data.frame(data[,2:26])
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

epsilon <- 120
while(epsilon > 10){
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

### getting the picture
im=readPNG("training_ss_151.png")
grid.raster(im)

## im has 3 color channels. Since the original image is a grey level image
## all 3 channels are identical
## the first dimension is vertical, bottom to top
## the second dimension is horizontal, left to right

## read the associated segmented image
im.seg=readPNG("training_seg_151.png")
grid.raster(im.seg)
im.seg.int=im.seg*255 # get integer values
## visualize the segmentation on top of the original image
# a trick to get the numeric values associated with each label
v=as.numeric(names(table(im.seg.int)))
# create a "background" image
# a "csf" image
# a "white matter" image and a "grey matter image"
im.back=im
im.csf=im
im.grey=im
im.white=im
for (i in (1:dim(im)[1])){
  for (j in (1:dim(im)[2])){
    im.back[i,j,]=ifelse(im.seg.int[i,j,1]==v[1],c(1,1,1),c(0,0,0))
    im.csf[i,j,]=ifelse(im.seg.int[i,j,1]==v[2],c(1,1,1),c(0,0,0))
    im.grey[i,j,]=ifelse(im.seg.int[i,j,1]==v[3],c(1,1,1),c(0,0,0))
    im.white[i,j,]=ifelse(im.seg.int[i,j,1]==v[4],c(1,1,1),c(0,0,0))
  }
}
#grid.raster(im.back)
#grid.raster(im.csf)
#grid.raster(im.grey)
#grid.raster(im.white)
## visualize the white matter only on the original image
#grid.raster(im*im.white)
##########################
## extract the data for the clustering
##########################

n = sum(!im.back[,,1])
coords=matrix(nrow=n, ncol=2)
index = 1
d=matrix(nrow=n,ncol=25)
for (i in (1:dim(im)[1])){
  for (j in (1:dim(im)[2])){
    if (!im.back[i,j,1]){
      if (index==1){
        print(i)
        print(j)
      }
      d[index,]=as.vector(im[(i-2):(i+2),(j-2):(j+2),1])
      coords[index,]=as.vector(c(i,j)) 
      index=index+1
    }
  }
}
# ## visualize the blocks
# 
im2=im
im2[(61-2):(61+2),(109-2):(109+2),1]*255
im2[(61-2):(61+2),(109-2):(109+2),]=c(0,1,1)
grid.raster(im2)

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
grid.raster(im3)
writePNG(im3,'REK_kmeans151.png')
