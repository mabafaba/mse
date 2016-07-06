require('raster')
require('entropy')
require(factoextra)
#small module functions
makefws <- function(d){
  print("makefws: creating a circle weight matrix for focal moving window for each radius")
  fws<-list()
  for(i in c(1:length(d))){
    fw<-focalWeight(x = raster(matrix(1,d[i],d[i])), d=1, type=c('circle'))
    fw<-fw*length(which(fw!=0))
    
    fws<-c(fws,list(fw))
  }
  
  return(fws)
}

fixEdgeEffects <- function(focaled,fw){
  # print('fixEdgeEffects: normalising window sum by valid number of cells in window')
  count<-focaled
  # print(0)
  count@data@values<-1
  # print(1)
  fcount<-focal(x = count,w=fw,fun=sum,na.rm=T,pad=T,padValue=0)
  # print(2)
  focaled@data@values<-focaled@data@values/fcount@data@values
  # print(3)
  focaled@data@values[is.na(focaled@data@values)]<-0
  # print(4)
  return(focaled)
}
singleCatRep<-function(r,id,fw){
  # print(paste('singleCatRep: moving window of category ',id, ' at radius', r))
  rcat<-r
  rcat@data@values[which(rcat@data@values!=id)]<-0
  rcat@data@values[which(rcat@data@values==id)]<-1
  #plot(rcat)
  # rcat<-aggregate(x = rcat,fact=2)
  #  plot(rcat)
  focaled <- focal(x = rcat,w=fw,fun=sum,na.rm=T,pad=T,padValue=0)
  focaled <- fixEdgeEffects(focaled,fw)
  is.na(focaled[1])
  return(as.vector(focaled))
}
mtoents<-function(m,d){
  print("mtoents: binning reps; for each column of all matrices in a list using function discretize; then function entropy")
  reps<-dummyreps(m,d)
  counts<-lapply(reps,FUN = function(x){apply(x,MARGIN = 2,FUN = function(x){discretize(x,numBins = 20,range(0,1))})})
  ents<-unlist(lapply(counts,entropy))
  #plot(d,ents,type="l")
  return(ents)
}

replistrange<-function(rl){
print('replistrange: get minimum and maximum values from entire replist')
    d<-1
  for(d in length(rl)){
    apply(do.call("rbind",lapply(rl, function(x) x[[d]])),2,range)
  }
}


dummyreps<-function(m,ds){
  # print('dummyreps: for all d, for all categories, singleCatRep')
  print('go')
  r<-raster(m)
  fws<-makefws(ds)
  reps<-list()
  ids<-unique(r)
  for (d in c(1:length(ds))){
    pixelsXids<-matrix(NA,nrow = length(r@data@values),ncol=length(ids))
    # aggregation thingy
    # pixelsXids<-matrix(NA,nrow = length(aggregate(r,6)@data@values),ncol=length(ids))
        print(paste('radius:',d,'now for all cats'))

    for (i in c(1:length(ids))){
      thisDandIDcatRep<-singleCatRep(r,ids[i],fws[[d]])
      pixelsXids[,i]<-thisDandIDcatRep
      }
    reps<-c(reps,list(pixelsXids))
  }
  return(reps)
}




H<-function(counts){
  # print('H: sum -p log(p)')
  probs = counts/sum(counts)
  probs<-probs[which(probs>0)]
  #probdens norm by amin
  return(
    -sum(probs*log2(probs))
  )
  
  #prob, discrete enropy
  #-sum(probs*log2(probs))
  
  #probdens, relative entropy 1/6
  #sum(probs*log2((probs/Ai)/fullspace))
  
}
binmethod <-function(counts,bins){
  # print('binmethod: discretize')
  return(discretize(counts-min(counts),numBins=bins,r=c(-1,1)))
}

# make dummy data functions
randomuniformmatrix <-function(matrixsize,ps){
  print('randomuniformmatrix: fills ncol=nrow matrix uniformly. depends on global ids.')
  #   m<-matrix(sample(ids,size = ((matrixsize/4)^2),prob = ps,replace = T),matrixsize/4,matrixsize/4)
  #   m<-raster(m)
  #   m<-disaggregate(m,4)
  #   m<-as.matrix(m)
    m<-matrix(sample(ids,size = ((matrixsize)^2),prob = ps,replace = T),matrixsize,matrixsize)
  
    return(m)
}

orderedMatrix <-function(matrixsize){
  cumsumPS<-c(0,cumsum(ps))
  print('orderedMatrix: ...')
  
  m2<-matrix(0,matrixsize,matrixsize)
  m2[,c(0:(ncol(m2)*ps[1]))]<-1
  for(i in 1:length(ps)){
    m2[,c((ncol(m2)*cumsumPS[i]):(ncol(m2)*cumsumPS[i+1]))]<-i
  }
  m2[which(m2==0)]<-1
  return(m2)
}
randomcascade<-function(matrixsize){
print('randomcascade: ...')
    cascade <- matrix(1,ncol=1,nrow=1,byrow = T)
  probs <- c(1.1,1.2,1.3,1.4)
  p <- matrix(probs,ncol=2,nrow=2,byrow = T)
  #add cascade level:
  for(i in c(1:log2(matrixsize))){
    cascade <- addLevel(cascade, p)
  }
  
   cas<-cascade
#   reds<-which(cas>median(cas))
#   yellows<-which(cas<=median(cas))
#   cas[yellows]<-0
#   cas[reds]<-1
#   cascade<-cas
cas<- (cas - min(cas))/diff(range(cas))
cumsumPS<-c(0,cumsum(ps))
quantiles<-quantile(cas,cumsumPS)
quantcas<-cas
for(i in c(1:(length(quantiles)-1))){
  quantcas[which(cas>=quantiles[i] & cas<=quantiles[i+1] )]<-i
}
return(quantcas)
}
randompatches<-function(matrixsize,patchsize){
  print('randompatches: ...')
  m<-matrix(sample(ids,size = ((matrixsize/patchsize)^2),prob = ps,replace = T),matrixsize/patchsize,matrixsize/patchsize)
  m<-raster(m)
  m<-disaggregate(m,patchsize)
  m<-as.matrix(m)
  return(m)
}


addLevel <- function(m,p){
  nc<-ncol(m)
  nr <-nrow(m)
  ncNew<-nc*2
  nrNew<-nr*2
  #make a matrix twice as big
  m2 <- matrix(1, ncol=ncNew,nrow=nrNew)
  for (r in c(1:nr)){
    for (c in c(1:nc)){
      #probs <- rnorm(4000, mean = 1, sd=0.001)
      #p <- matrix(probs,ncol=2,nrow=2,byrow = T)
      #multiply each cell of original matrix with the probability matrix (shuffeld), then put values into new matrix
      m2[((((r-1)*2))+1):((((r-1)*2))+2),((((c-1)*2))+1):((((c-1)*2))+2)]<- m[r,c]*p[sample.int(2),sample.int(2)]
    }
  }
  return(m2)
}


#batch functions
replist<-function(mlist,ds){
  print('replist: apply making reps to list of matrices')
  print(ds)
  lapply(mlist,FUN = dummyreps,ds=ds)
}
batch_pca <- function(rl){
  print('batch_pca: list of reps, for each rep, for each radius, make the pcas, rotate data')
  pcas<-list()
  rotated_values<-list()
  for(case in c(1:length(rl))){
    pcas[[case]]<-list()
    rotated_values[[case]]<-list()
    for(radius in c(1:length(rl[[case]]))){
      print(paste('case',case,'radius',radius))
      pcas[[case]][[radius]]  <- prcomp(rl[[case]][[radius]], scale. = F, center = T)
      rotated_values[[case]][[radius]] <- predict(pcas[[case]][[radius]],rl[[case]][[radius]])
    }
  }
  result <- list()
  result$pcas <- pcas
  result$rotateddata <- rotated_values
  return(result)
}
batch_ents<-function(pcas,bins=100){
  counts <- list()
  entropies <- list()
  for(case in c(1:length(pcas$rotateddata))){
    counts[[case]]<-list()
    entropies[[case]]<-list()
    for(radius in c(1:length(pcas$rotateddata[[case]]))){
      
      counts[[case]][[radius]]<-apply(X = pcas$rotateddata[[case]][[radius]],MARGIN = 2,FUN = binmethod,bins=bins)
      
      # counts[[case]][[radius]]<-apply(rl[[case]][[radius]],2,binmethod)
      entropies[[case]][[radius]]<-apply(counts[[case]][[radius]],2,H)
    }
  }
  
  sumentropies<-lapply(entropies,function(x){lapply(x,sum)})
  sumentropies<-matrix(unlist(sumentropies),ncol = length(entropies),byrow = F)
  return(sumentropies)
}
batch_mlist_rl_pca_ents<-function(mlist,returnpcas=FALSE,ds){
  print('batch_mlist_rl_pca_ents: making reps, pcas, and ents for all')
  rl<-replist(mlist,ds)
  pcas<-batch_pca(rl)
  ents <-batch_ents(pcas)
  if(returnpcas){
  ep<-list()
  ep$ents<-ents
  ep$pcas<-pcas
  return(ep)}
  else{
    return(ents)
    }
}

ranents<-function(matrixsize,matrixmaker){
  paste('ranents')
  e<-batch_mlist_rl_pca_ents(list(matrixmaker(matrixsize)))
  return(e)
}

ranentspatch<-function(matrixsize,patchsize,matrixmaker){
  print('ranentspatch')
  e<-batch_mlist_rl_pca_ents(list(matrixmaker(matrixsize,patchsize)),returnpcas = TRUE)
  return(e)
}

#parameters
# ps<-c(1/2,1/4,1/4)
ps<-rep(1/2,2)
ids<-c(1:length(ps))
d<-c(1,2,3,4,5,6,10,20,40,80,128)
d<-c(1,2,4,8,16,32,64,128)

matrixsize<-128
3
library(doParallel)
stopCluster(cl)
cl <- makeCluster(3)
registerDoParallel(cl)
matrixmakerscollection<-c(randomuniformmatrix,orderedMatrix,randomcascade)
getwd()
Sys.time()
entsFE<- foreach(i=1:3, .packages=c("raster","entropy"),.export=c(
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
  'batch_mlist_rl_pca_ents')) %dopar% ranents(matrixsize,randomuniformmatrix)
Sys.time()
entsFEordered<- foreach(i=1:3, .packages=c("raster","entropy"),.export=c(
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
  'batch_mlist_rl_pca_ents')) %dopar% ranents(matrixsize,orderedMatrix)
saveRDS(entsFEpatch,'entsFEpatch.rds')
Sys.time()
entsFEcasc<- foreach(i=1:3, .packages=c("raster","entropy"),.export=c(
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
  'batch_mlist_rl_pca_ents')) %dopar% ranents(matrixsize,randomcascade)



