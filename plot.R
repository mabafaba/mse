
for(i in 1:length(test)){ 
    
    pdf(sprintf("myplot%d.pdf",i)) 
    
    plot(   
        ((2*test[[i]]$me$radii)+1)
        ,test[[i]]$me$entropy
        ,type='l'
        ,ylim=c(0,5)
        ,xlim=c(1,max(test[[i]]$me$radii))
        ,log='x'
        ,xlab='radius'
        ,ylab='entropy'
        ,main=paste('patch size:',test[[i]]$patchsize)
    ) 

        # abline(v=unlist(lapply(test,function(x){(2*x$me$radii)+1})),col='grey')
    abline(v=(test[[i]]$patchsize))

    image(
        randompatchmatrixGeneric(
            randomuniformmatrix(test[[i]]$me$matrixsize,c(0.5,0.5))
            ,test[[i]]$patchsize)
        ,col=c('black','white')
        ,xaxt='n'
        ,yaxt='n'
    )

    dev.off() 
}



pdf("myplotcascRect.pdf") 
    plot(   
        ((2*cascRect$radii)+1)
        ,cascRect$entropy
        ,type='l'
        ,ylim=c(0,5)
        ,xlim=c(1,max(cascRect$radii*2)+1)
        ,log='x'
        ,xlab='radius'
        ,ylab='entropy'
        ,main='multiplicative cascade rect')
    ) 

        # abline(v=unlist(lapply(test,function(x){(2*x$me$radii)+1})),col='grey')

    image(randomcascade(256)
        ,col=c('black','white')
        ,xaxt='n'
        ,yaxt='n'
    )

    dev.off() 
}
makefws
casc
    for(i in 1:length(test)){ 
        
        pdf(sprintf("2d-b20-log%d.pdf",i)) 
        
        plot(   
            ((2*recttest[[i]]$me$radii)+1)
            ,recttest[[i]]$me$entropy
            ,type='l'
            ,ylim=c(0,5)
            ,xlim=c(1,256)#max(recttest[[i]]$me$radii))
            ,log='x'
            ,xlab='radius'
            ,ylab='entropy'
            ,main=paste('patch size:',recttest[[i]]$patchsize)
        ) 
    plotcasc()

            # abline(v=unlist(lapply(test,function(x){(2*x$me$radii)+1})),col='grey')
        abline(v=(recttest[[i]]$patchsize))
        abline(h=(log2(20)),col='red')
    abline(h=seq(1,4,0.5),lwd=0.1,col='grey')
        image(
            randompatchmatrixGeneric(
                randomuniformmatrix(recttest[[i]]$me$matrixsize,c(0.5,0.5))
                ,recttest[[i]]$patchsize)
            ,col=c('black','white')
            ,xaxt='n'
            ,yaxt='n'
        )

        dev.off() 
    }



ids<-c(1:32)
testplot<-test32d
for(i in 1:length(testplot)){ 
    
    pdf(sprintf("32d-b20-patch%d.pdf",i)) 
    
    plot(   
        ((2*testplot[[i]]$me$radii)+1)
        ,testplot[[i]]$me$entropy
        ,type='l'
        ,ylim=c(0,log2(binnum)*(max(ids)-1)+1)
        ,xlim=c(1,1+2*max(testplot[[i]]$me$radii))
        ,log='x'
        ,xlab='radius'
        ,ylab='entropy'
        ,main=paste('patch size:',testplot[[i]]$patchsize)
        ,bty="n"
    ) 

        # abline(v=unlist(lapply(test,function(x){(2*x$me$radii)+1})),col='grey')
    abline(v=(test[[i]]$patchsize))
        abline(h=(log2(binnum)*(max(ids)-1)),col='red')
            abline(h=seq(1,floor((log2(binnum)*max(ids))),0.5),lwd=0.1,col='grey')
axis(1, col = "white", tcl = 0)


    image(
        randompatchmatrixGeneric(
            randomuniformmatrix(testplot[[i]]$me$matrixsize,rep(1/max(ids),max(ids)))
            ,test[[i]]$patchsize)
        # ,col=c('black','white')
        ,xaxt='n'
        ,yaxt='n'
    )

    dev.off() 
}



image(randompatchmatrixGeneric(randomuniformmatrix(test3d[[i]]$me$matrixsize,c(1/3,1/3,1/3)),10))


i
i<-1




plotcasc<-function(){
lines(   (2*casc$radii)+1
        ,casc$entropy
        ,col='blue'
        )

    }











#plot different binnings

testplot<-ent1875real
for(i in 1:length(testplot)){ 
    
    pdf(sprintf("testest.pdf",i)) 
    
    plot(   
        ((2*testplot[[i]]$me$radii)+1)
        ,testplot[[i]]$me$entropy
        ,type='l'
        ,ylim=c(0,log2(binnum)*(max(ids)-1)+1)
        ,xlim=c(1,1+2*max(testplot[[i]]$me$radii))
        ,log='x'
        ,xlab='radius'
        ,ylab='entropy'
        ,main=paste('patch size:',testplot[[i]]$patchsize)
        ,bty="n"
    ) 

        # abline(v=unlist(lapply(test,function(x){(2*x$me$radii)+1})),col='grey')
    abline(v=(test[[i]]$patchsize))
        abline(h=(log2(binnum)*(max(ids)-1)),col='red')
            abline(h=seq(1,floor((log2(binnum)*max(ids))),0.5),lwd=0.1,col='grey')
axis(1, col = "white", tcl = 0)


    image(
        randompatchmatrixGeneric(
            randomuniformmatrix(testplot[[i]]$me$matrixsize,rep(1/max(ids),max(ids)))
            ,test[[i]]$patchsize)
        # ,col=c('black','white')
        ,xaxt='n'
        ,yaxt='n'
    )

    dev.off() 
}




str(ctest2d)
testplot<-ent1875real
for(i in 1:length(testplot)){ 
    
    pdf(sprintf("c2d-patch%d.pdf",i)) 
    
    plot(   
        ((2*testplot[[i]]$me$radii)+1)
        ,matrix(testplot[[i]]$me$entropy,ncol=50,nrow=15)[,1]/seq(10,500,10)[1]
        ,type='l'
        ,ylim=c(0,1)
        ,xlim=c(1,1+2*max(testplot[[i]]$me$radii))
        ,log='x'
        ,xlab='radius'
        ,ylab='entropy'
        ,main=paste('patch size:',testplot[[i]]$patchsize)
        ,bty="n"
    ) 

    for (bin in c(1:50)) {
    lines((2*testplot[[i]]$me$radii)+1
    ,matrix(testplot[[i]]$me$entropy,ncol=50,nrow=15)[,bin]/log2(seq(10,500,10)[bin])
    ,col=rainbow(50)[bin])
    }

        # abline(v=unlist(lapply(test,function(x){(2*x$me$radii)+1})),col='grey')
    abline(v=(test[[i]]$patchsize))
        abline(h=(log2(binnum)*(max(ids)-1)),col='red')
            abline(h=seq(1,floor((log2(binnum)*max(ids))),0.5),lwd=0.1,col='grey')
axis(1, col = "white", tcl = 0)


    image(
        randompatchmatrixGeneric(
            randomuniformmatrix(testplot[[i]]$me$matrixsize,rep(1/max(ids),max(ids)))
            ,test[[i]]$patchsize)
        # ,col=c('black','white')
        ,xaxt='n'
        ,yaxt='n'
    )

    dev.off() 
}


plot(   
        ((2*testplot$radii)+1)
        ,testplot$entropy[,1]
        ,type='l'
        ,ylim=c(0,log2(binnum)*(max(ids)-1)+1)
        ,xlim=c(1+2*min(testplot$radii),1+2*max(testplot$radii))
        ,log='x'
        ,xlab='radius'
        ,ylab='entropy'
        ,main=paste('year')
        ,bty="n"
    ) 

years


matrix(ctest2d[[1]]$me$entropy,ncol=50,nrow=15)
par(mfrow=c(1,1))
par(new=FALSE)
image(L1915m)
image(L)
pdf('years.pdf')
plot( ent2005real$radii,ent2005real$entropy[,10],type='n', xlim=c(1,60),ylim=c(0,100),log='x')
lines(ent1875real$radii,ent1875real$entropy[,10],lty=1, pch=1 ,type='b',cex=0.5)
lines(ent1895real$radii,ent1895real$entropy[,10],lty=1, pch=2 ,type='b',cex=0.5)
lines(ent1915real$radii,ent1915real$entropy[,10],lty=1, pch=3 ,type='b',cex=0.5)
lines(ent1935real$radii,ent1935real$entropy[,10],lty=1, pch=4 ,type='b',cex=0.5)
lines(ent1960real$radii,ent1960real$entropy[,10],lty=1, pch=5 ,type='b',cex=0.5)
lines(ent1985real$radii,ent1985real$entropy[,10],lty=1, pch=6 ,type='b',cex=0.5)
lines(ent2005real$radii,ent2005real$entropy[,10],lty=1, pch=7 ,type='b',cex=0.5)

dev.off()
str()
