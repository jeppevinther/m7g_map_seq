## getFreq2000 (Anders Albrechtsen, 2017)
## Function to parse mpileup file, count mutations for each position and perform statistical tests


getFreq <- function(inFile, #input multisample mpileup file (samtools mpileup -b bam.list -f seq.fasta)
                    
                    CC = c(0,0,0,1,1,1), #vector defining the samples (here 3x control, 3x treated)
                    
                    CSVout = "out.txt", #output file

                    minFrac = 0.001, #minimum fraction of mutatation for analysis

                    minCounts = 2, #position only analysed if totalCount>minCounts
                    
                    pvalThres = -1, #position only analysed if SNPPval>pvalThres ( for test of all sites even if non-variable use -1).

                    nCores = 50, #number of cores used for analysis

                    nSites = 500, #number of positions parsed before analysis

                    maxSites = 1e6 #maxnumber of positions analysed
                     ) {

    
resAll<-list()
siteList<-list()
siteNum<-1

con <-file(inFile,"r")
write.table( rbind(c("chr","pos","ref","totalCounts","relFreq","CCPval","indPval","CC0Pval","CC1Pval","SNPPval","countsCC0","countsCC1","typeCountCC0","typeCountCC1","totalCountsCC0","totalCountsCC1","AltFreq","AltFreqCC0","AltFreqCC1","SNPfreq","CC0freq","CC1freq","counts","indCounts","indFreq")),file=CSVout,row=F,qu=F,col=F,sep="\t")

for(iSites in 1:maxSites){
                                        #r<-readLines(con,n=1)
    cat("\r", iSites)
    r<-scan(con,nl=1,what="df",quiet=TRUE,sep="\t",quote=NULL)
    {if(length(r)>0)   { ## parse data
        siteList[[siteNum]] <- r
        siteNum <- siteNum +1
    }
    else{

        cat("last chunk\n")
    }}
    
    if(siteNum%%nSites==0 | length(r)==0 | iSites==maxSites){#variant analysis
      
       resAll<-parallel::mclapply(siteList,oneSiteNoBaq,minFrac=minFrac,minCounts=minCounts,pvalThres=pvalThres,mc.cores=nCores,CC=CC)
    # resAllmat <-list()
    # for(h in 1:length(siteList) )
    #     resAllmat[[h]] <- oneSiteNoBaq(siteList[[h]],minFrac=minFrac,minCounts=minCounts,pvalThres=pvalThres,CC=CC)

        resAllmat <- do.call(rbind,resAll)

        write.table(resAllmat,file=CSVout,append=TRUE,col=FALSE,row=FALSE,sep="\t",qu=FALSE)
        siteNum<-1
         if(length(r)==0 | iSites==maxSites)
           break
       resAll <- list()
       siteList<-list()
    }
    
}

close(con)
}



parseMpileupNoBaq<-function(string,Q,minFrac,minCounts){
# string=r[2,1+i];Q=r[3,1+i] ;minFrac=minFrac;minCounts=minCounts;pvalThres=pvalThres


    
    string <- gsub("\\^.","",string)
    string <- gsub("\\$","",string) 
    string <- gsub("\\,",".",string)
    x <- unlist(strsplit(string,""))
   x <- toupper(x) ## bases covert to capital
 

    
    w<- which(x%in%c("+","-"))
    ## if integer larger than 9
    while(any(x[w+2]%in%c(0:9))){
        wh <- w[which(x[w+2]%in%c(0:9))][1]
        x[wh+1]<-paste0(x[wh+1:2],collapse="")
        x <- x[-(wh+2)]
        w<- which(x%in%c("+","-"))
    }
    ## entries for indels c(-/+, LENGHT,TYPE)
    rmC <- unlist(lapply(w,function(y) y+1+c(-1:as.integer(x[y+1]))))
   
    ##type of indel
    type <- sapply(w,function(y) paste(x[y+1+c(1:as.integer(x[y+1]))],collapse=""))
    typeDel<-c()
    typeIns<-c()
    if(length(type)>0){
        for(i in 1:length(type)){   
            if(x[w[i]]=="-")
                typeDel<- c(typeDel, tolower(type[i]))
            else
                typeIns<-c(typeIns,  toupper(type[i]))
        }
        x <- x[-rmC]
    }
    
    star<-c()
    starProb<-c()
   
    prob <- QtoProb(Q)


    if(length(wAster<-which(x=="*"))>0){
        star<-x[wAster]
        starProb<-prob[wAster]
        prob<-prob[-wAster]
        x<-x[-wAster]
    }
    highI<-list()

    typeAll<-c(typeIns,star)
    depthX<-length(x) + length(typeAll)
    highB<-table(x)>minFrac*depthX  & table(x)>minCounts
    highB <- names(highB)[highB]
     if(length(typeAll)>0){
        highI<-table(typeAll)>minFrac*depthX & table(typeAll)>minCounts
        highI <- names(highI)[highI]
    }

    
    
    x[x==","] <- "."
    if(any(x=="*"))
        cat("* in file\n")

    if(length(x)!=length(prob))
        stop("wrong length Q!=b")

 
    
    list(base=x,type=typeAll,Q=prob,depthX=depthX,keepBase=highB,keepIndel=highI,nextType=typeDel)
}

oneSiteNoBaq<-function(r,minFrac,minCounts,pvalThres,CC){
    ##r <- siteList[[h]] ;minFrac=minFrac;minCounts=minCounts;pvalThres=pvalThres
    if(length(CC)+1 != length(r)/3)
        stop("wrong size of CC vector\n")
    
    dim(r)<-c(3,length(CC)+1)
    pos <- r[,1]
    names(pos)=c("chr","pos","ref")
    site<-list()
    for(i in 2:ncol(r)-1)
        site[[i]] <- parseMpileupNoBaq(r[2,1+i],r[3,1+i],minFrac=minFrac,minCounts=minCounts)
    site[[1]]$pos<-pos
     
    
    nIndel <- length(unlist( sapply(site,function(x) x$type) ))
    pos <- site[[1]]$pos
    nSNP <- sapply(site,function(x) length(x$base))
    
   
    if(all(sapply(site,function(x) length(x$Q))!=sapply(site,function(x) length(x$base))))
        stop("Q not equal depth")
    
    rAll<-analyse(site,multi=TRUE,pos)
    if(rAll$pval0_all>pvalThres){
        r0<-analyse(site[CC==0],multi=TRUE,ref=pos,types<-rAll$types)
        r1<-analyse(site[CC==1],multi=TRUE,ref=pos,types=rAll$types)
        pvalCC <- pval(2*(r0$like+r1$like - rAll$like),length(rAll$freq)-1)
        rMulti<-  multiAnalyse(site,pos,types=rAll$types)
        pval0 <- pval(2*(sum(sapply(rMulti,function(x) x$like)[CC==0])-r0$like),length(r0$freq)-1)
        pval1 <- pval(2*(sum(sapply(rMulti,function(x) x$like)[CC==1])-r1$like),length(r1$freq)-1)
        
        results<-c(pos,totalCounts=sum(rAll$counts),
                   relFreq=round( (1-r1$freq[[1]]) -(1-r0$freq[[1]]) ,4),
                   CCpval=pvalCC,
                   indPval=paste(sapply(rMulti,function(x) x$pval),collapse=";"),
                   CC0Pval=pval0,
                   CC0Pval=pval1,
                   SNPpval=rAll$pval0_all,
                   countsCC0=compressFreq(r0$counts),
                   countsCC1=compressFreq(r1$counts),
                   typeCountsCC0=length(r0$counts),
                   typeCountsCC1=length(r1$counts),
                   totalCountsCC0=sum(r0$counts),
                   totalCountsCC1=sum(r1$counts),         
                   altFreq=round(1-rAll$freq[[1]],4),
                   altFreqCC0=round(1-r0$freq[[1]],4),
                   altFreqCC1=round(1-r1$freq[[1]],4),
                   SNPfreq=compressFreq(rAll$freq),          
                   CC0freq=compressFreq(r0$freq),
                   CC1freq=compressFreq(r1$freq),
                   counts=compressFreq(rAll$counts),
                   indCounts=paste(sapply(rMulti,function(x) compressFreq(x$counts)),collapse=";"),
                   indFreq=paste(sapply(rMulti,function(x) compressFreq(x$freq)),collapse=";")                 
                   )
    }
    if(rAll$pval0_all<=pvalThres)
        results<-c(pos,totalCounts=sum(rAll$counts),counts=compressFreq(rAll$counts),SNPpval=rAll$pval0_all,rep(NA,16) )
    results
    
}

oneSite<-function(site,pvalThres,CC){

    nIndel <- length(unlist( sapply(site,function(x) x$type) ))
    pos <- site[[1]]$pos
    nSNP <- sapply(site,function(x) length(x$base))
    
   
    if(all(sapply(site,function(x) length(x$Q))!=sapply(site,function(x) length(x$base))))
        stop("Q not equal depth")
    
    rAll<-analyse(site,multi=TRUE,pos)
    if(rAll$pval0_all>pvalThres){
        r0<-analyse(site[CC==0],multi=TRUE,ref=pos,types<-rAll$types)
        r1<-analyse(site[CC==1],multi=TRUE,ref=pos,types=rAll$types)
        pvalCC <- pval(2*(r0$like+r1$like - rAll$like),length(rAll$freq)-1)
        rMulti<-  multiAnalyse(site,pos,types=rAll$types)
        pval0 <- pval(2*(sum(sapply(rMulti,function(x) x$like)[CC==0])-r0$like),length(r0$freq)-1)
        pval1 <- pval(2*(sum(sapply(rMulti,function(x) x$like)[CC==1])-r1$like),length(r1$freq)-1)
        
        results<-c(pos,totalCounts=sum(rAll$counts),counts=compressFreq(rAll$counts),
                   SNPpval=rAll$pval0_all,
                   altFreq=round(1-rAll$freq[[1]],4),
                   SNPfreq=compressFreq(rAll$freq),
                   CCpval=pvalCC,CC0freq=compressFreq(r0$freq),CC1freq=compressFreq(r1$freq),
                   controlpval=pval0,casesPval=pval1,
                   altFreqCC0=round(1-r0$freq[[1]],4),
                   altFreqCC1=round(1-r1$freq[[1]],4),
                   totalCountsCC0=sum(r0$counts),
                   totalCountsCC1=sum(r1$counts),
                   countsCC0=compressFreq(r0$counts),
                   countsCC1=compressFreq(r1$counts),
                   indPval=paste(sapply(rMulti,function(x) x$pval),collapse=";"),
                   indFreq=paste(sapply(rMulti,function(x) compressFreq(x$freq)),collapse=";"),
                   indCounts=paste(sapply(rMulti,function(x) compressFreq(x$counts)),collapse=";")
                   )
    }
    if(rAll$pval0_all<=pvalThres)
             results<-c(pos,counts=compressFreq(rAll$counts),SNPpval=rAll$pval0_all,
                        altFreq=round(1-rAll$freq[[1]],4),
                        SNPfreq=compressFreq(rAll$freq),
                        rep(NA,15) )

    results
    
}


len<-function(x){
    print(table(strsplit(x,"")[[1]]))
    length(strsplit(x,"")[[1]])
}

pval<-function(x,df=1){
    if(df==0)
        return(0)
    if(any(x< -10)) #happends for ultra high depth
        stop("negative test statistics")
    min(999,round(-10*log10(1 - pchisq(x,df))))
}
compressFreq <- function(x){
    paste(paste0(names(x),"=",round(x,4)),collapse=",")
}

#function that seperates indels and SNPs
parseMpileup<-function(string,Q,preS,minFrac,minCounts){


    
    string <- gsub("\\^.","",string)
    string <- gsub("\\$","",string) 
    string <- gsub("\\,",".",string)
    x <- unlist(strsplit(string,""))
    x <- toupper(x) ## bases covert to capital


    
    w<- which(x%in%c("+","-"))
    ## if integer larger than 9
    while(any(x[w+2]%in%c(0:9))){
        wh <- w[which(x[w+2]%in%c(0:9))][1]
        x[wh+1]<-paste0(x[wh+1:2],collapse="")
        x <- x[-(wh+2)]
        w<- which(x%in%c("+","-"))
    }
    ## entries for indels c(-/+, LENGHT,TYPE)
    rmC <- unlist(lapply(w,function(y) y+1+c(-1:as.integer(x[y+1]))))
   
    ##type of indel
    type <- sapply(w,function(y) paste(x[y+1+c(1:as.integer(x[y+1]))],collapse=""))
    typeDel<-c()
    typeIns<-c()
    if(length(type)>0){
        for(i in 1:length(type)){   
            if(x[w[i]]=="-")
                typeDel<- c(typeDel, tolower(type[i]))
            else
                typeIns<-c(typeIns,  toupper(type[i]))
        }
        x <- x[-rmC]
    }
    
    
   
    prob <- QtoProb(Q)


    if(length(wAster<-which(x=="*"))>0){
        prob<-prob[-wAster]
        x<-x[-wAster]
    }
    highI<-list()

    typeAll<-c(typeIns,preS$nextType)
    depthX<-length(x) + length(typeAll)
    highB<-table(x)>minFrac*depthX  & table(x)>minCounts
    highB <- names(highB)[highB]
     if(length(typeAll)>0){
        highI<-table(typeAll)>minFrac*depthX & table(typeAll)>minCounts
        highI <- names(highI)[highI]
    }

    
    
    x[x==","] <- "."
    
    if(length(x)!=length(prob))
        stop("wrong length Q!=b")

    if(any(x=="*"))
        stop("* in file")

    
    list(base=x,type=typeAll,Q=prob,depthX=depthX,keepBase=highB,keepIndel=highI,nextType=typeDel)
}
 QtoProb<-function(x){
   y <- utf8ToInt(x) ## ascii to integer
   Q <- y - 33 #offset
   10^(-Q/10)
}
normF<-function(x) x/sum(x)
if(FALSE){
    q<-Q[[1]]
    ref<-pos
    b<-base[[1]][[1]]
    i<-base[[1]][[2]]

    

    microbenchmark(t(t(GL) * fn), 
                   GL %*% diag(fn),
                   GL * rep(fn, rep.int(nrow(GL),length(fn))), 
                   GL * rep(fn, rep(nrow(GL),length(fn))), 
                   GL * rep(fn, each = nrow(GL)),
                   GL*fn
                   )

}


getLike<-function(GL,fn){

    Qz <-  GL * rep(fn, rep(nrow(GL), length(fn)))
    sum(log(rowSums(Qz)))

}

EM <- function(GL){

    #start guess
    fnPlus1 <- rep(1,ncol(GL))/ncol(GL)
    fn <- 100

    #EM
    iter<-1
    while(any(abs(fn-fnPlus1)>0.00001)){ #stop criteria
        fn <- fnPlus1
        ##Q function:
        Qz <-  GL * rep(fn, rep(nrow(GL), length(fn)))
        norm <- rowSums(Qz)
	#M step
        fnPlus1 <- normF(colSums(Qz/norm))
        iter <- iter +1 
    }
    names(fnPlus1) <- colnames(GL)
    like<-sum(log(rowSums(Qz)))
    w<-which.max(fnPlus1)

    f0 <- rep(0,length(fnPlus1))
    f0[w]<-1


    like0<-getLike(GL,f0)
    if(fnPlus1[w]>0.999)
        pval<-1
    else
        pval<-pval(-2*(like0-like),length(fn)-1)
    return(list(freq=fnPlus1,like=like,major=names(w),like0=like0,pval0_all=pval))
}

addProtect<-function(x)
    log(sum(exp(x-max(x))))+max(x)

EMunknownMinor <- function(GL,major){

    #start guess
    fnPlus1 <- c(0.2,0.8)
 
    fn <- 100

    nonMajor<-(1:ncol(GL))[-major]
    #EM
    iter<-1
    QzL<-list()
    while(any(abs(fn-fnPlus1)>0.000001)){ #stop criteria
        fn <- fnPlus1
        ##Q function:
        for(i in 1:length(nonMajor))
            QzL[[i]] <-  GL[,c(major,nonMajor[i])] * rep(fn, rep(nrow(GL), length(fn)))
        s<-sapply(QzL,function(x) sum(log(rowSums(x))))
        post<-normF(exp( s-max(s)))
        Qz<-QzL[[1]]/rowSums(QzL[[1]])*post[1]
        for(i in 2:length(nonMajor))
            Qz <- Qz + QzL[[i]]/rowSums(QzL[[i]])*post[i]
       print( addProtect(sapply(QzL,function(x)(sum(log(rowSums(x))))))-log(length(nonMajor)))
    	#M step
        fnPlus1 <- normF(colSums(Qz))
        iter <- iter +1 
    }
    names(fnPlus1) <- c("major","unknownMinor")
    like <- addProtect(sapply(QzL,function(x)(sum(log(rowSums(x))))))-log(length(nonMajor))
      like0<-sum(log(GL[,major]))

    pval<-pval(2*(like0-like),length(fn)-1)
    return(list(freq=fnPlus1,like=like,major=names(w),like0=like0,pval0_all=pval))
}

analyse<-function(theSite,multi=FALSE,ref,types=NULL,indelErr=0.005){

# b<-site[[i]];multi=FALSE;ref<-pos;types=types;indelErr<-0.00
    if(multi){
        if(is.null(types)){
            types<-list(
                B=unique(unlist(lapply(theSite,function(x)x$keepBase))),
                I=unique(unlist(lapply(theSite,function(x)x$keepIndel)))
                )
        }
        i<-unlist(lapply(theSite,function(x)x[[2]]))
        q<-unlist(lapply(theSite,function(x)x$Q))
        b<-unlist(lapply(theSite,function(x)x[[1]]))
    }
    else{
        i<-theSite$type
        q<-theSite$Q
        b<-theSite$base
    }
    if(!is.null(types)){
        rmB <- !b%in%types$B
        b <- b[!rmB]
        q <- q[!rmB]
        if(length(types$I)>0)
            i<-i[i%in%types$I]
        else
            i<-c()
    }
    if(length(i)==0)
         { #no indels
          ty <- types$B
          if(length(ty)==0 | length(b)==0)
		res<- list(pval0_all=0,like=1,like0=1,freq=1)
          else if(length(ty)==1){
              lik <- sum(log(1-q))
              res<-list(pval0_all=0,like=lik,like0=lik,freq=1)
          }
          else{
              GL<-matrix(NA,nrow=length(b),ncol=length(ty))
              colnames(GL) <- ty
              for( ti  in ty)
                  GL[,ti]<-ifelse(b==ti,1-q,q/3)

              res<-EM(GL)
          }
          res[["counts"]]<-table(b)
          res[["totalCounts"]]<-length(b)
      }
    else{
        ##with indels
        namesI<-paste0("i",types$I)
        GL<-matrix(NA,nrow=length(b)+length(i),ncol=length(types$B)+length(types$I))
        colnames(GL) <- c(types$B,namesI)

        for( ti  in types$B)
            GL[,ti]<-c( ifelse(b==ti,1-q,q/3) , rep(indelErr/4,length(i)))

        for(ii in 1:length(types$I))
            GL[,namesI[ii]]<-c(  rep(0 , length(b)) , ifelse(i==types$I[ii],1-indelErr,indelErr/length(types$I)) )
                
        res<-EM(GL)
        res[["counts"]]<-table(c(b,paste0("i",i)))
        res[["totalCounts"]]<-length(b)+length(i)

    }
          

    res[["types"]] <- types

    
    
    res
}


multiAnalyse<-function(site,pos,types){

    res<-list()
    nInd<-length(site)
    for(i in 1:nInd)
        res[[i]]<-analyse(site[[i]],multi=FALSE,pos,types=types)

    res
}
