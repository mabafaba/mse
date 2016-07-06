
library("sp")
library("rgdal")
library("maptools")
library("ggplot2")
library("UScensus2010")
library("geospacom")
library("pdist")
library("mgcv")
library("abind")
library('raster')
library(doParallel)
setwd('~/multiscale-entropy')
#read spatial data source and transform to UTM coordinate system

L1875<-readOGR("all_LU.shp", layer="lu_ 1875_all")
L1895<-readOGR("all_LU.shp", layer="lu_1895_all")
L1915<-readOGR("all_LU.shp", layer="lu_1915_all")
L1935<-readOGR("all_LU.shp", layer="lu_1935_all")
L1960<-readOGR("all_LU.shp", layer="lu_1960_all")
L1985<-readOGR("all_LU.shp", layer="lu_1985_all")
L2005<-readOGR("all_LU.shp", layer="lu_2005_all")
L1875<-spTransform(L1875, CRS("+proj=utm +zone=30 ellps=WGS84")) # the utm  zone must be specified
L1895<-spTransform(L1895, CRS("+proj=utm +zone=30 ellps=WGS84")) # the utm  zone must be specified
L1915<-spTransform(L1915, CRS("+proj=utm +zone=30 ellps=WGS84")) # the utm  zone must be specified
L1935<-spTransform(L1935, CRS("+proj=utm +zone=30 ellps=WGS84")) # the utm  zone must be specified
L1960<-spTransform(L1960, CRS("+proj=utm +zone=30 ellps=WGS84")) # the utm  zone must be specified
L1985<-spTransform(L1985, CRS("+proj=utm +zone=30 ellps=WGS84")) # the utm  zone must be specified
L2005<-spTransform(L2005, CRS("+proj=utm +zone=30 ellps=WGS84")) # the utm  zone must be specified


years<-c("1875","1895","1915","1935","1960","1985","2005")
spdfNames<-paste("L",years,sep="")


#translate categories to integers

catints<-c(1:32)
names(catints)<-levels(L2005@data$LU)

L1875@data$LU<-catints[L1875@data$LU]
L1895@data$LU<-catints[L1895@data$LU]
L1915@data$LU<-catints[L1915@data$LU]
L1935@data$LU<-catints[L1935@data$LU]
L1960@data$LU<-catints[L1960@data$LU]
L1985@data$LU<-catints[L1985@data$LU]
L2005@data$LU<-catints[L2005@data$LU]

# rasterize spatial Polygons

#rasterise with 4 times the resolution used in analysis
#this is necessary because the rasterisation only looks at which polygon covers the pixels centroid.
#rasterising at high resolution and then aggregating later for expensive analysis makes much more precise.

L1875<-rasterize(L1875, y = raster(x = L1875,ncol=68*12,nrow=32*12),field="LU")
L1895<-rasterize(L1895, y = raster(x = L1895,ncol=68*12,nrow=32*12),field="LU")
L1915<-rasterize(L1915, y = raster(x = L1915,ncol=68*12,nrow=32*12),field="LU")
L1935<-rasterize(L1935, y = raster(x = L1935,ncol=68*12,nrow=32*12),field="LU")
L1960<-rasterize(L1960, y = raster(x = L1960,ncol=68*12,nrow=32*12),field="LU")
L1985<-rasterize(L1985, y = raster(x = L1985,ncol=68*12,nrow=32*12),field="LU")
L2005<-rasterize(L2005, y = raster(x = L2005,ncol=68*12,nrow=32*12),field="LU")

#transform raster objects back into matrices. for visualisation, matrix is transposed and rotated to original direction.
#the same code from category to int (catints) was used for all matrices. this is given as a parameter to use the lowest unused number for empty fields.
#make matrix and rotate 90:
raster2matrix<-function(rasterobj,catints){
print('raster2matrix: ...')
m<-t(matrix(data = rasterobj@data@values,nrow = rasterobj@nrows,ncol = rasterobj@ncols,byrow=T))
m<-m[,ncol(m):1]
m[which(is.na(m))]<-max(catints)+1
return(m)
}


hist(L2005m,143)
L1875m<-raster2matrix(L1875,catints)
L1895m<-raster2matrix(L1895,catints)
L1915m<-raster2matrix(L1915,catints)
L1935m<-raster2matrix(L1935,catints)
L1960m<-raster2matrix(L1960,catints)
L1985m<-raster2matrix(L1985,catints)
L2005m<-raster2matrix(L2005,catints)

image(L2005m)
#######################################################
#######################################################
# 16:38
#######################################################
#######################################################
randompatchmatrixGeneric(L1895m,patchsize = 1)

#because rasterised in over resolution, need to overwrite the singlecatrep function to include 4x aggregation 
singleCatRep<-function(r,id,fw){
  # if(i%%5 ==0) {print(paste('singleCatRep: moving window of category ',id, 'for one fw'))}
  
  rcat<-r
  rcat@data@values[which(rcat@data@values!=id)]<-0
  rcat@data@values[which(rcat@data@values==id)]<-1
  # plot(rcat)
  ################ this line is extra:
   # rcat<-aggregate(x = rcat,fact=6)
  ####################################
  #   plot(rcat)
  focaled <- focal(x = rcat,w=fw,fun=sum,na.rm=T,pad=T,padValue=0)
  focaled <- fixEdgeEffects(focaled,fw)
  return(as.vector(focaled))
}


#now calculate the entropy for each year:
L2005msmall<-raster(L2005m)
length(as.matrix(aggregate(L2005msmall,6)))

matrix(L2005m)

image(randompatchmatrixGeneric(originalMatrix = as.matrix(L2005msmall),patchsize = 128))

#parameters: (missing bin size)
d<-c(1,2,4,8,16,32,64,128)
require(doParallel)
# sum(table(L2005@data@values))/min(table(L2005@data@values))

#parallel computation:
L_all_m<-list(L1875m,L1895m,L1915m,L1935m,L1960m,L1985m,L2005m)

stopCluster(cl)
cl <- makeCluster(2)
registerDoParallel(cl)

Sys.time()
ereal<- foreach(i = 1:7, .packages=c("raster","entropy"),.export=c(
  'fixEdgeEffects', 
  'singleCatRep', 
  'dummyreps', 
  'H', 
  'binmethod', 
  'randomuniformmatrix', 
  'orderedMatrix', 
  'randomcascade', 
  'addLevel', 
  'replist',
  'batch_pca',
  'batch_ents',
  'batch_mlist_rl_pca_ents')) %dopar% batch_mlist_rl_pca_ents(list(L_all_m[[i]]))
Sys.time()

stopCluster(cl)
cl <- makeCluster(4)
registerDoParallel(cl)
epatchesrandom<-list()
#do for each patch size..
Sys.time()
epatchesrandomthis<-list()
# for(patchsize in patchsizes){
#and for all years parallel:
  epatchesrandom2005<- foreach(patchsize=patchsizes, .packages=c("raster","entropy"),.export=c(
  'fixEdgeEffects', 
  'singleCatRep', 
  'dummyreps', 
  'H', 
  'binmethod', 
  'randomuniformmatrix', 
  'orderedMatrix', 
  'randomcascade', 
  'addLevel', 
  'replist',
  'batch_pca',
  'batch_ents',
  'batch_mlist_rl_pca_ents')) %dopar% batch_mlist_rl_pca_ents(list(
                  randompatchmatrixGeneric(originalMatrix = L2005m,patchsize = patchsize)),returnpcas = TRUE)
# }
#epatchesrandom structure: epatchesrandom -> patchsizes -> years -> radii
Sys.time()
stopCluster(cl)



epatchesrealthis
lapply(epatchesrealthis,lines,function(x){
  lines
})


#trying for a good bin size.
#take 2005, and randomise once for each patch size
#then 
randomms<-list()
image(randomms[[8]])
for(i in c(1:length(d))){
randomms<-c(randomms,list(randompatchmatrixGeneric(originalMatrix = L2005m,patchsize = d[i])))
}

!!!
rl<-replist(mlist = randomms,ds = d)
rlpca<-batch_pca(rl)

rlents10
rlents5 <-batch_ents(pcas = rlpca,bins = 5)

rlents20 <-batch_ents(pcas = rlpca,bins = 20)
rlents30 <-batch_ents(pcas = rlpca,bins = 30)
rlents100 <-batch_ents(pcas = rlpca,bins = 100)
rlents200 <-batch_ents(pcas = rlpca,bins = 200)
rlents400 <-batch_ents(pcas = rlpca,bins = 400)
par(mfrow=c(3,3))
for(i in c(1:8)){
  this<-i
plot(d,rlents20[,this]/log(20),type='n',log='x',ylim=c(0,20),main=d[this])
abline(v=d[this])
lines(d,rlents5[,this]/log(5),type='l',col='blue')
lines(d,rlents10[,this]/log(10),type='l',col='blue')
lines(d,rlents20[,this]/log(20),type='l',col='blue')
lines(d,rlents100[,this]/log(100),type='l',col='blue')
lines(d,rlents200[,this]/log(200),type='l',col='blue')
lines(d,rlents400[,this]/log(400),type='l',col='blue')
}
d
d
par(mfrow=c(1,1))
binmethod
!!!



plot(d,epatchesrealthis[[4]],type="l",log='x')
lines(d,epatchesrealthis[[1]],type="l")
lines(d,epatchesrealthis[[2]],type="l")
lines(d,epatchesrealthis[[3]],type="l")
lines(d,epatchesrealthis[[4]],type="l")
lines(d,epatchesrealthis[[5]],type="l")
lines(d,epatchesrealthis[[6]],type="l")
lines(d,epatchesrealthis[[7]],type="l")
lines(d,epatchesrealthis[[8]],type="l")

lines(d,ereal[[7]],type="l")
plot(d,ereal[[7]],type="l",log='x')
ereal[[1]]
ereal
epatchesrealthis


plot(d,ereal[[7]],type='l',log='x')
l2005
lines(d,ereal[[3]],type='l',col='red')

lines(d,epatchesrealthis[[8]],type='l')
str(rlpca$rotateddata)
l_random_m<-lapply(L_all_m,function(m){m<-matrix(sample(m,size = length(m),replace = F),nrow(m),ncol(m))})
image(l_random_m[[3]])
str(rlpca$rotateddata[[1]][[1]])
par(mfrow=c(8,1))
hist(rlpca$rotateddata[[8]][[7]],main='',xlab='')
hist(rlpca$rotateddata[[2]][[7]],xlim=c(-0.3,0.3),main='',xlab='')

hist(rlpca$rotateddata[[1]][[2]])

lines(density(rlpca$rotateddata[[1]][[8]]))



cl <- makeCluster(4)
registerDoParallel(cl)
Sys.time()
erandom2005many<- foreach(i = 1:100, .packages=c("raster","entropy"),.export=c(
  'fixEdgeEffects', 
  'singleCatRep', 
  'dummyreps', 
  'H', 
  'binmethod', 
  'randomuniformmatrix', 
  'orderedMatrix', 
  'randomcascade', 
  'addLevel', 
  'replist',
  'batch_pca',
  'batch_ents',
  'batch_mlist_rl_pca_ents')) %dopar% batch_mlist_rl_pca_ents(list(matrix(sample(L2005m,size = length(L2005m),replace = F),nrow(L2005m),ncol(L2005m))))
Sys.time()
saveRDS(erandom2005many,file = "erandom2005many.ROBJ")



saveRDS(ereal,file = "ereal.ROBJ")
saveRDS(erandom,file = "erandom.ROBJ")
xval<-log(d)
plot(xval,ereal[[3]],ylim=c(0,40))
lines(xval,ereal[[1]])
lines(xval,erandom[[1]])

# sequential:
real_data_entropy_results_all_years <- list()
for(i in 1:7){
  real_data_entropy_results_all_years <- c(real_data_entropy_results_all_years,list(batch_mlist_rl_pca_ents(list(L_all_m[[i]]))))
}
erandom2005manyM <- do.call(cbind, erandom2005many)

plot(d,erandom2005manyM[,1],type='l',log="x")
for(i in c(1:ncol(erandom2005manyM))){
  lines(d,erandom2005manyM[,i])
  
  }
install.packages("psych")
library("psych")
error.bars(erandom2005manyM[1,])
hist(erandom2005manyM[2,],breaks=30)
rm(ps)
mtest<-cbind(randomuniformmatrix(128),randomuniformmatrix(128))
sort(unique(as.vector(L2005m)))
image(L2005m)
batch_mlist_rl_pca_ents(list(mtest))
boxplot(L2005m[1,])
batch_mlist_rl_pca_ents(list(L2005m))
ids
a<-apply(erandom2005manyM,1,mean)
s<-apply(erandom2005manyM,1,sd)
n<-ncol(erandom2005manyM)
error <- qnorm(0.975)*s/sqrt(n)
left<-a-error
right<-a+error
plot(d,a,type="l",log="x")
lines(d,left)
lines(d,right)

df<-data.frame(a=a,d=log(d),u=left,l=right,real=)
require(ggplot2)
ggplot(df, aes(x = d, y = a)) 
   +geom_line(size = 0.1) 
   +geom_point(size = 4,) 
   +geom_errorbar(aes(ymax = u, ymin = l))
   +geom_line(aes) +
df<-data.frame(x=log(d),y=ereal[[7]])

plot(data = df, aes(x=x,y=y),color = "red")
plot(log(d),ereal[[7]],type='l')
lines(log(d),ereal[[2]],type='l')

str(ereal)

?geom_point
str(mtest)
table(L2005m)
table(mtest)
#JUNK
L2005@data@attributes[[1]]$LU
L2005
str(L2005@data@attributes[[1]])
matrix(as.numeric(L2005@data@attributes[[1]]$LU
)L2005@nrows,L2005@ncols)

hist(as.numeric(L2005@data@attributes[[1]]$LU
)
which(as.numeric(L2005@data@attributes[[1]]$LU)==)

length(L2005@data@attributes[[1]]$LU)
L2005@data@values


?rasterize
#JUNK
?rasterize
lo<-L2005
test<-rasterize(lo,y = raster(x = lo,ncol=68,nrow=32),field="LU")
plot(test)
test@data@values[which(is.na(test@data@values))]<-33    
lo@data$LU
types<-(unique(as.integer(levels(lo@data$LU)))
types<-as.integer(as.factor(levels(lo@data$LU)))        
names(types)<-levels(lo@data$LU)
types
types[lo@data$LU]
types[lo@data$LU[test@data@values]]

test@data@values<-lo@data$LU[test@data@values]
str(test)
?rasterize
plot(test)
test2<-aggregate(x = test,fact=4)

a = as.factor(c("A", "B", "A", "C"))
b = as.integer(factor(a))


