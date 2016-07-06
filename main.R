setwd('~/multiscale-entropy')

# categories
ids = c(1:2)
radii<-c(3,5,7,9,11,13,15,17,19,23,29,31,33,49,63,127)
radii<-(radii-1)/2
patchsizes<-c(1,2,4,8,16,32,64)

binnum<-20

# random matrix with ents:
randompatchmatrixGeneric

randomCase<-function(matrixsize,patchsize,radii,bins){
    m<-randompatchmatrixGeneric(randomuniformmatrix(matrixsize,rep(1/max(ids),max(ids))),patchsize)
    # radii 

    #representation
    reps<-dummyreps(m,radii)

    pca<-batch_pca(list(reps))
    me<-list()
    me$entropy<-batch_ents_comparebins(pca)
    me$patchsize<-patchsize
    me$radii<-radii
    me$matrixsize<-matrixsize
    return(me)
}

randomCasesFromReal<-function(originalMatrix,originalMatrixName,patchsizes,radii){
    all<-list()
    for(patchsize in patchsizes){
    print(paste('patchsize:',patchsize))
    m<-randompatchmatrixGeneric(originalMatrix,patchsize)
    # radii 

    #representation
    reps<-dummyreps(m,radii)

    pca<-batch_pca(list(reps))
    me<-list()
    me$entropy<-batch_ents_comparebins(pca)
    me$patchsize<-patchsize
    me$radii<-radii
    me$original<-originalMatrixName
    all<-c(all,list(me))
    }
    return(all)
}


cascadeCase<-function(matrixsize,radii,bins,rectangle=false){
    casc<-randomcascade(matrixsize)
    if(rectangle){


makefwsoriginal<-makefws
makefws<-makefwsrect
    }
    reps<-dummyreps(casc,radii)
    pca<-batch_pca(list(reps))
    me<-list()
    me$entropy<-batch_ents(pca,bins=bins)
    me$patchsize<-patchsize
    me$radii<-radii
    me$matrixsize<-matrixsize

if(rectangle){
makefws<-makefwsoriginal

}
    return(me)
}
  
m<-matrix(floor(runif(10000)*5)+1,100,100)
m<-L1875m
singleCatRep 
realcase<-function(m){

    # radii 

    #representation
    reps<-dummyreps(m,radii)

    pca<-batch_pca(list(reps))
    me<-list()
    me$entropy<-batch_ents_comparebins(pca)
    # me$patchsize<-patchsize
    me$radii<-radii
    me$matrixsize<-c(nrow(m),ncol(m))
    return(me)   
}


require(doParallel)
stopCluster(cl)
cl <- makeCluster(3)
registerDoParallel(cl)
singleCatRep
test<-realcase(m)

singleCatRep
yearmatrixnames<-paste('L',years,'m',sep='')

ent1875real<-realcase(L1875m)
saveRDS(ent1875real,'ent1875real.rds')
ent1895real<-realcase(L1895m)
saveRDS(ent1895real,'ent1875real.rds')
ent1915real<-realcase(L1915m)
saveRDS(ent1915real,'ent1915real.rds')
ent1935real<-realcase(L1935m)
saveRDS(ent1935real,'ent1935real.rds')
ent1960real<-realcase(L1960m)
saveRDS(ent1960real,'ent1960real.rds')
ent1985real<-realcase(L1985m)
saveRDS(ent1985real,'ent1985real.rds')
ent2005real<-realcase(L2005m)
saveRDS(ent2005real,'ent2005real.rds')



few_patchsizes<-c(1,4,8,16,32,64)

Sys.time()
entPatches1875 <-randomCasesFromReal(L1875m,'L1875m',few_patchsizes,radii[-16],bins)
saveRDS(entPatches1875,'entPatches1875.rds')
Sys.time()
entPatches1895 <-randomCasesFromReal(L1895m,few_patchsizes,radii[-16],bins)
saveRDS(entPatches1895,'entPatches1895.rds')
Sys.time()
entPatches1915 <-randomCasesFromReal(L1915m,few_patchsizes,radii[-16],bins)
saveRDS(entPatches1915,'entPatches1915.rds')
Sys.time()
entPatches1935 <-randomCasesFromReal(L1935m,few_patchsizes,radii[-16],bins)
entPatches1960 <-randomCasesFromReal(L1960m,few_patchsizes,radii[-16],bins)
entPatches1985 <-randomCasesFromReal(L1985m,few_patchsizes,radii[-16],bins)
entPatches2005 <-randomCasesFromReal(L2005m,few_patchsizes,radii[-16],bins)



randompatchmatrixGeneric(L1875m,10)










ctest2d<- foreach(i=yearmatrixnames, .packages=c("raster","entropy"),.export=c(
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
  'batch_mlist_rl_pca_ents',
  'randomCase',
  'batch_ents_comparebins'
  ,'randompatchmatrixGeneric'
  )) %dopar% realcase()
Sys.time()



ctest2d<- foreach(i=patchsizes, .packages=c("raster","entropy"),.export=c(
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
  'batch_mlist_rl_pca_ents',
  'randomCase',
  'batch_ents_comparebins'
  ,'randompatchmatrixGeneric'
  )) %dopar% random_batch_bins(i)
Sys.time()

asdfasfd



random_batch_bins<-function(x){
print(paste('patchsize:::',x))
thisone<-list()
thisone$patchsize<-x
thisone$me<-randomCase(256,x,radii,binnum)
return(thisone)
}





ids<-c(1,2)
ctest2d<-lapply(patchsizes,function(x){
print(paste('patchsize:::',x))
thisone<-list()
thisone$patchsize<-x
thisone$me<-randomCase(256,x,radii,binnum)
return(thisone)
    })


ids<-c(1:8)
ctest8d<-lapply(patchsizes,function(x){
print(paste('patchsize:::',x))
thisone<-list()
thisone$patchsize<-x
thisone$me<-randomCase(256,x,radii,binnum)
return(thisone)
    })

ids<-c(1:16)
ctest16d<-lapply(patchsizes,function(x){
print(paste('patchsize:::',x))
thisone<-list()
thisone$patchsize<-x
thisone$me<-randomCase(256,x,radii,binnum)
return(thisone)
    })


ids<-c(1:32)
ctest32d<-lapply(patchsizes,function(x){
print(paste('patchsize:::',x))
thisone<-list()
thisone$patchsize<-x
thisone$me<-randomCase(256,x,radii,binnum)
return(thisone)
    })



recttest<-lapply(patchsizes,function(x){
thisone<-list()
thisone$patchsize<-x
thisone$me<-randomCase(256,x,radii,binnum)
return(thisone)
    })

ids


makefwsoriginal<-makefws
makefws<-makefwsrect
recttest<-lapply(patchsizes,function(x){
thisone<-list()
thisone$patchsize<-x
thisone$me<-randomCase(200,x,radii,binnum)
return(thisone)
    })
makefws<-makefwsoriginal

cascRect<-cascadeCase(256,radii,binnum,TRUE)






randomCasesFromReal




