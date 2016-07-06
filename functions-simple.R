# randompatchmatrixGeneric
# randomuniformmatrix
# fixEdgeEffects
# makefws
# singleCatRep
# batch_ents
# batch_pca

makefwsrect <- function(d){
  print("makefws: creating a circle weight matrix for focal moving window for each radius")
  fws<-list()
  for(i in c(1:length(d))){
    fw<-focalWeight(x = raster(matrix(1,d[i],d[i])), d=1, type=c('rectangle'))
    fw<-fw*length(which(fw!=0))
    
    fws<-c(fws,list(fw))
  }
  
  return(fws)
}








batch_ents_comparebins<-function(pcas,binlist=seq(10,500,10)){
allsumentropies<-vector()
for (bins in binlist) {
  counts <- list()
  entropies <- list()
  for(case in cbind(1:length(pcas$rotateddata))){
    counts[[case]]<-list()
    entropies[[case]]<-list()
    for(radius in c(1:length(pcas$rotateddata[[case]]))){

      counts[[case]][[radius]]<-apply(X = pcas$rotateddata[[case]][[radius]],MARGIN = 2,FUN = binmethod,bins=bins)
      hist(counts[[case]][[radius]],breaks=1000)
      
      # counts[[case]][[radius]]<-apply(rl[[case]][[radius]],2,binmethod)
      entropies[[case]][[radius]]<-apply(counts[[case]][[radius]],2,H)
    print('this case this radius done')

    }
  }
  
  sumentropies<-lapply(entropies,function(x){lapply(x,sum)})
  sumentropies<-matrix(unlist(sumentropies),ncol = length(entropies),byrow = F)

allsumentropies<-cbind(allsumentropies,sumentropies)
}
  return(allsumentropies)
}

