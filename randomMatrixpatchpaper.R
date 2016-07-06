randompatchmatrixGeneric<-function(originalMatrix,patchsize){
    rows = nrow(originalMatrix)
    cols = ncol(originalMatrix)
    ps = table(originalMatrix)/sum(table(originalMatrix))

    #empty matrix
    m<-matrix(0,rows,cols)
    #set proportion of categories (including "empty" cells)
    cats<-c(1:length(ps))
    #set size of randomly distributed patches

    #calculate how many individual cells of each category
    totals<-length(m)*ps/patchsize
    totals<-ceiling(totals)

    arr<-array()
    #make an array of patches
    for(i in c(1:length(cats))){
    fullpatches<-rep(patchsize,floor(totals[i]/patchsize)) # add columns for full patches
    restpatch<-ceiling((totals[i]/patchsize - floor(totals[i]/patchsize))*patchsize)
    patches<-rbind(cats[i],c(fullpatches,restpatch))
    arr<-cbind(arr, patches)#repeat for each category

    }


    #shuffle array
    arr<-arr[,sample(x = c(1:ncol(arr)),size = ncol(arr),replace = F)]
    arr<-arr[,which(!is.na(arr[1,]))]
    #expand array of patches to array of cells
    arr<-unlist(apply(arr,2,function(x){rep(x[1],x[2])}))
    #account for rounding errors
    # arr<-arr[1:length(m)]
    #expand to matrix

    m2<-matrix(100,nrow =floor(nrow(m)/patchsize) #number of rows: rows*cols
     ,ncol = ncol(m),byrow=FALSE)


    m2<-t(m2)
    m2[,]<-arr[1:length(m2)]
    m2<-t(m2)


    # m2<-matrix(m2,,nrow =floor(nrow(m)/patchsize) #number of rows: rows*cols
    #            ,ncol = ncol(m),byrow=TRUE)

    m2<-(m2[sort(rep(c(1:nrow(m2)),patchsize)), ])
    return(m2)
}

p<-rep(1/4,4)
testm<-matrix(ps,100,100)
uni<-randomuniformmatrix(matrixsize = 256,ps = c(0.5,0.5))

test<-randompatchmatrixGeneric(originalMatrix = uni,patchsize = 25)
testms<-lapply(c(1,2,4,8,16,64,128,254),function(x){
    randompatchmatrixGeneric(uni,x)
    })
testms[[1]]
testes<-lapply(c(1,2,4,8,16,64,128,254),function(x){
  randompatchmatrixGeneric(uni,x)
  })

length(testms)

ds<-c(1,2,4,8,16,64,128)
teste$ents
teste<-batch_mlist_rl_pca_ents(testms,returnpcas = TRUE,ds=ds)
str(teste)

dev.off()
par(mfrow=c(8,2))
par(mfrow=c(1,1))

plot(ds*6,teste$ents[,1],type='n',main=1,ylim=c(0,6),xlim=c(1,128*6),log='x')
dev.off()
for (i in c(1:8)) {
    image(testms[[i]],col=c('black','white'))
    plot(ds*6,teste$ents[,i],type='l',main=c(1,2,4,8,16,64,128,254)[i],ylim=c(0,6),xlim=c(0,800))

  # lines(ds*6,teste$ents[,i],main=i,ylim=c(0,6),xlim=c(0,130))
}
text(ds[5]*6,teste$ents[5,],c(1,2,4,8,16,64,128,254))

?text
plot(teste$ents[,1],ds,type='l')

library(entropy)

(table(test)/sum(table(test)))-ps

L2005ms<-L2005m
L2005ms<-raster(L2005m)
L2005ms<-aggregate(L2005ms,4)
image(L2005ms)
L2005ms<-as.matrix(L2005ms)
L2005m32<-randompatchmatrixGeneric(originalMatrix = L2005m,patchsize = 32)
L2005m128<-randompatchmatrixGeneric(originalMatrix = L2005m,patchsize = 128)
L2005m1<-randompatchmatrixGeneric(originalMatrix = L2005m,patchsize = 1)
image(L2005m128)
image(L2005m1)

table(L2005ms)
singleCatRep
L2005m

binmethod

L2005msEnt<-batch_mlist_rl_pca_ents(list(L2005m),returnpcas = TRUE)
L2005mEnt1<-batch_mlist_rl_pca_ents(list(L2005m1),returnpcas = TRUE)
L2005mEnt128<-batch_mlist_rl_pca_ents(list(L2005m128),returnpcas = TRUE)
L2005mEnt32<-batch_mlist_rl_pca_ents(list(L2005m32),returnpcas = TRUE)


plot(d,L2005msEnt$ents,type='l',ylim=c(0,100),log='x')
lines(d,L2005mEnt1$ents)
lines(d,L2005mEnt128$ents)

saveRDS(L2005msEnt,file = 'L2005msEnt.rds')
saveRDS(L2005mEnt1,file = 'L2005mEnt1.rds')
saveRDS(L2005mEnt128,file = 'L2005mEnt128.rds')
dev.off()
par(mfrow=c(4,4))
principleComponent<-2
scale<-1
for(scale in c(1:8)){
    plot(density(L2005mEnt128$pcas$rotateddata[[1]][[scale]][,principleComponent]),type='l',ylim=c(0,100),main=paste('scale:', d[scale]),col='red')
# lines(density(L2005mEnt1$pcas$rotateddata[[1]][[scale]][,principleComponent]),xlim=c(-1,1),col='blue')
lines(density(L2005msEnt$pcas$rotateddata[[1]][[scale]][,principleComponent]),xlim=c(-1,1),col='black')
}

# plot(d,L2005msEnt$ents,type='l',ylim=c(0,100),log='x',col='black')
# lines(d,L2005mEnt1$ents,col='blue')
# lines(d,L2005mEnt128$ents,col='red')

for(scale in c(1:8)){
  hist(L2005mEnt128$pcas$rotateddata[[1]][[scale]][,principleComponent],main=paste('scale:', d[scale]),col=rgb(1,0,0,0.2),breaks=seq(from=-1,to=1,by=0.05))
  # hist(L2005mEnt1$pcas$rotateddata[[1]][[scale]][,principleComponent],main=paste('scale:', d[scale]),col=rgb(0,0,1,0.2),breaks=seq(from=-1,to=1,by=0.05),add=T)
  hist(L2005msEnt$pcas$rotateddata[[1]][[scale]][,principleComponent],main=paste('scale:', d[scale]),col=rgb(0,0,0,0.2),breaks=seq(from=-1,to=1,by=0.05),add=T)
  
}

image(L2005ms)

image(L2005m1)
image(L2005m128)
image(fws[[4]])
L2005m128circ<-L2005m128

L2005m128circ[which(fws[[4]]==1)] <-100
image(L2005m128circ)

fws<-makefws(d)

enttest<- function(x,numBins=100,r=c(-1,1)){ 
    entropy(discretize(x[,1],numBins = numBins,r = r))
}

enttest(L2005mEnt128)

entropy(L2005mEnt128$pcas$rotateddata[[1]][[scale]][,1])
entropy(L2005mEnt1$pcas$rotateddata[[1]][[scale]][,1])


counts<-lapply(L2005msEnt$rotateddata[[1]][[1]],FUN = function(x){apply(x,MARGIN = 2,FUN = function(x){discretize(x,numBins = 20,range(-1,1))})})
ents<-unlist(lapply(counts,entropy))


L2005msEnt
L2005msEnt$pcas[[1]][[1]][[1]]
plot(L2005msEnt$ents)
attributes(L2005msEnt$pcas[[1]][[1]][[1]]$rotation)
L2005ms
