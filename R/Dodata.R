
##########################################################
#' To perform method
#' @param dir the path of the output
#' @param calculatedata the inputdata
#' @param chr chromosome
#' @param color1 the color of the chromosome
#' @param CLO the numbers of CPU
#' @return list
#' @export
#' @examples
#' data(BSA)
#' dir<- tempdir()
#' data.calculatedata<-Readdata1(BSA)
#' Dodata(dir,calculatedata=data.calculatedata,chr="all",color1="blue",CLO=1)

Dodata<-function(dir,calculatedata,chr,color1,CLO){
  
  # dir<-dir1
  # calculatedata<-data.calculatedata
  # chr="all";color1="blue"
  
  Fileformat<-calculatedata$FileFormat; population<- calculatedata$Population
  windowsize<-as.numeric(calculatedata$windowsize);
  np<-as.numeric(calculatedata$np)/100;
  Smooth_method<-calculatedata$Smooth_method;
  repetition<-as.numeric(calculatedata$Repetition)
  DrawPlot<-calculatedata$DrawPlot
  Resolution<-calculatedata$Resolution
  Plotformat<-calculatedata$Plotformat
  nmarker<-500;qq = matrix(nmarker/2,byrow = TRUE);h2<-0.20;v1<-10;vg<-v1*h2/(1-h2);b0<-100
  nsam<-as.numeric(calculatedata$nindividual)
  #nsam<-26
  chrom<-chr;color<-color1
  if(!((Fileformat=="BSA")|(Fileformat=="Extreme individuals")|(Fileformat=="CIM")|(Fileformat=="ICIM")|(Fileformat== "GCIM")
  )) stop ("Please input the correct Data-file-format")
  if(!((population=="F2")|(population=="BC")|(population=="DH")|(population=="RIL")
  )) stop ("Please input the correct Population-type")
  if(!is.numeric(np)) stop("Please input the correct Sample-size")
  if(!((Smooth_method=="None")|(Smooth_method=="Default")|(Smooth_method=="AIC")|(Smooth_method=="Window size")|(Smooth_method=="Block")|(length(strsplit(Smooth_method,split = " ")[[1]])==3)
       | (length(strsplit(Smooth_method,split = " ")[[1]])==2))) stop ("Please input the correct Smooth_method")
  if(!is.numeric(repetition)) stop("Please input the correct No. of permutations")
  if(!((DrawPlot=="TRUE")|(DrawPlot=="FALSE")
  )) stop ("Please input the correct Figure")
  if(!((Resolution=="Low")|(Resolution=="High")
  )) stop ("Please input the correct Figure-resolution")
  if(!((Plotformat=="jpeg")|(Plotformat=="png")|(Plotformat=="tiff")|(Plotformat=="pdf")
  )) stop ("Please input the correct Figure-file-format")
  
  
  
  
  #####no-smooth########
  None<-function(result.sta){
    result.smooth<-result.sta
    return(result.smooth)
  }
  ######AIC#########
  AIC.f<-function(G){
    #G<-result
    chr = as.numeric(G$Chromosome);cc<-unique(chr)
    ff<-numeric();
    for(i in 1:length(cc)){
      #i<-1
      r<-cc[i]
      marker<-G[which(chr==r),1]
      pos<-as.numeric(G[which(chr==r),3])
      yy=as.matrix(G[which(chr==r),4:dim(G)[2]])
      numeric.change<-function(x){
        x<-as.numeric(x)
        return(x)
      }
      yy<-apply(yy,2,numeric.change)
      
      loess.f<-function(y){
        x=as.matrix(as.numeric(pos)) #chromosomal coordinate
        y=as.matrix(y)
        opt.span<-function(model,criterion = c("aicc","gcv"),span.range = c(0.05,0.95) ){
          as.crit <-function(x){
            span <- x$pars$span;
            traceL <- x$trace.hat;
            sigma2 <-sum(x$residuals^2)/(x$n-1);
            aicc <-log(sigma2)+1+2*(2*(traceL+1))/(x$n-traceL-2);
            gcv<-x$n*sigma2/(x$n-traceL)^2;
            result<-list(span=span,aicc=aicc,gcv=gcv);
            return(result)
          }
          critertion <- match.arg(criterion);
          fn <- function(span){
            mod <- update(model,span =span)
            as.crit(mod)[[criterion]]
          }
          
          result<-optimize(fn,span.range);
          return(list(span = result$minimum,criterion =result$objective));
        }
        locfit_by_loess <- function(x,y){
          fit0 <- loess(y~x,degree = 2)
          span1 <-opt.span(fit0,criterion = "aicc")$span
          fit1 <- loess(y~x,degree = 2,span = span1)
          #predict
          plo <- predict(fit1,x,se=TRUE)
          return(plo)
        }
        smoothdata<-locfit_by_loess(x,y)
        smoothg<-smoothdata$fit
        return(smoothg)
      }
      result.loess<-apply(yy, 2, loess.f)
      ff<-rbind(ff,result.loess)
      
    }
    ff1<-cbind(G[,1:3],ff)
    return(ff1)
  }
  
  #####Windowsize#########
  Windowsize.f<-function(G,Windowsize){
    #G<-result;Windowsize<-1000
    chr = as.numeric(G$Chromosome);cc<-unique(chr)
    ff<-numeric();
    for(i in 1:length(cc)){
      #i<-1
      r<-cc[i]
      marker<-G[which(chr==r),1]
      pos<-as.numeric(G[which(chr==r),3])
      yy=as.matrix(G[which(chr==r),4:dim(G)[2]])
      numeric.change<-function(x){
        x<-as.numeric(x)
        return(x)
      }
      yy<-apply(yy,2,numeric.change)
      
      
      winsize.f<-function(y){
        #y<-yy[,1]
        x=as.matrix(as.numeric(pos)) #chromosomal coordinate
        y=as.matrix(y)
        w=Windowsize #3375000       # width
        olen = length(x)
        
        #beginning intervals
        xfirst = x[1,]      #start site
        startband = x <= xfirst + w
        xstart = xfirst - (x[startband] - xfirst)
        ystart = y[startband]
        nstart = length(xstart)-1
        
        #end of intervals
        xlast=x[length(x),]#last site
        endband = x >= xlast - w
        xend = xlast + (xlast - x[endband])
        yend = y[endband]
        nend = length(xend) - 1
        
        #at the head and tail of the sequence reflect the right and left (respectively)
        a=rev(xstart)   #reversal
        a=a[-length(a)] #remove the last site of b
        b=rev(xend)
        b=b[-1]         #remove the first site of b
        x<-c(a,x,b)
        c=rev(ystart)
        c=c[-length(c)] #remove the last site of b
        d=rev(yend)
        d=d[-1]
        y<-c(c,y,d)
        lx = length(x);ly= length(y)
        
        if(is.null(CLO)==TRUE){
          cl.cores <- detectCores()
          if((cl.cores<=2)){
            cl.cores<-1
          }else if(cl.cores>2){
            if(cl.cores>10){
              cl.cores<-10
            }else {
              cl.cores <- detectCores()-1
            }
          }else{
            if(CLO>10){
              cl.cores<-10
            }else{
              cl.cores<-CLO
            }
          }
          
        }
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        #Calculates moving average of y at coordinate x_i by weighted averaging
        #over all points in range (x_i-w, x_i+w)
        yw=foreach(i=1:lx)%dopar%
          {
            c = x[i]
            inband = ( x >= (c-w)&x <= (c+w))
            xfrac = (abs(x[inband] - c))/w
            xwt = (1.0 - abs(xfrac)^3)^3    #Weights
            xwt[abs(xfrac) >= 1] = 0
            ywin = sum(y[inband]*xwt)/sum(xwt)
          }
        stopCluster(cl)
        #trim off the reflected right/left bandwidths
        g2=cbind(y,yw)
        h1=length(a)+1
        h2=length(a)+olen
        G2=g2[h1:h2,]
        smooth_G<-unlist(G2[,2])
        return(smooth_G)
      }
      result.winsize<-apply(yy, 2, winsize.f)
      ff<-rbind(ff,result.winsize)
      
    }
    ff1<-cbind(G[,1:3],ff)
    return(ff1)
  }
  
  ######Block###############
  Block.f<-function(data,blocknumber){
    #data<-result;blocknumber<-10
    chr = as.numeric(data$Chromosome);cc<-unique(chr)
    ff<-numeric();
    for(i in 1:length(cc)){
      #i<-6
      r<-cc[i];
      marker<-data[which(chr==r),1]
      pos<-as.numeric(data[which(chr==r),3])
      yy=as.matrix(data[which(chr==r),4:dim(data)[2]])
      numeric.change<-function(x){
        x<-as.numeric(x)
        return(x)
      }
      yy<-apply(yy,2,numeric.change)
      
      block.f<-function(y){
        #y<-yy[,1]
        block<-blocknumber
        x=as.matrix(as.numeric(pos)) #chromosomal coordinate
        y=as.matrix(y)
        endband = length(x) - block+1
        yend = y[endband:length(x)]
        d.yend<-rev(yend)
        newy<-c(y,d.yend)
        
        if(is.null(CLO)==TRUE){
          cl.cores <- detectCores()
          if((cl.cores<=2)){
            cl.cores<-1
          }else if(cl.cores>2){
            if(cl.cores>10){
              cl.cores<-10
            }else {
              cl.cores <- detectCores()-1
            }
          }else{
            if(CLO>10){
              cl.cores<-10
            }else{
              cl.cores<-CLO
            }
          }
          
        }
        
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        
        y.mean<-foreach(i =1:length(x))%dopar%
          {
            start1<-i
            end<-start1+block-1
            block.y<-newy[start1:end]
            block.mean<-mean(block.y)
          }
        stopCluster(cl)
        
        by<-unlist(y.mean)
        return(by)
      }
      result.block<-apply(yy, 2, block.f)
      ff<-rbind(ff,result.block)
    }
    ff1<-cbind(data[,1:3],ff)
    return(ff1)
  }
  
  #####significance##########
  Significance<-function(data,threshold){
    
    #data<-result.all
    nc<-(dim(data)[2]-3)/2
    result.sig<-list(NULL)
    for(nn in 1:nc){
      #nn<-1
      result<-data[,c(1,2,3,(2*nn-1)+3,2*nn+3)]
      chr = as.numeric(result$Chromosome);cc<-unique(chr);aa<-numeric()
      for(i in 1:length(cc)){
        l<-cc[i]
        bb<-result[which(result[,2]==l),]
        bb<-na.omit(bb)
        bb1<-as.numeric(bb[,5])
        
        if(bb1[1] > bb1[2]){
          aa<-rbind(aa,bb[1,])
        }
        for (i in 1:(dim(bb)[1]-2)){
          if ((bb1[i+1]>=bb1[i]) & (bb1[i+1]>=bb1[i+2])){
            aa<-rbind(aa,bb[i+1,])
          }
        }
        if((bb1[dim(bb)[1]])>(bb1[dim(bb)[1]-1])){
          aa<-rbind(aa,bb[dim(bb)[1],])
        }
      }
      res<-aa[which(abs(aa[,5])>=threshold[nn]),]
      if(dim(res)[1]==0){
        result.sig[[nn]]<-NA
      }else{
        Threshold<-rep(NA,dim(res)[1])
        Threshold[1]<-threshold[nn]
        result.sig[[nn]]<-cbind(res,Threshold)
      }
    }
    return(result.sig)
  }
  
  
  ##########Permutation###########
  gen.f <- function(map.data,p.type,nsam){
    sim.cross(map = map.data, model = NULL, n.ind = nsam, type=p.type,
              error.prob=1e-10, missing.prob=0, partial.missing.prob=0,
              keep.qtlgeno=FALSE, keep.errorind=FALSE, m=0, p=0)
  }
  
  F2.data<-function(map.data,p.type,nsam,vg,qq,b0,v1 ){
    F2.result<-list(NULL)
    xq <- matrix(c(1, 0, -1, 0, 1, 0), 2, 3, byrow = TRUE);
    rd<-1;a<-sqrt(4*vg/(rd^2+2));d<-a*rd;bq <- matrix(c(a,d), 1, 2)#Additive and dominant effect values  20%
    sim.f2.gen.data<-gen.f(map.data,p.type,nsam)
    genotype <- sim.f2.gen.data$geno[[1]]$data
    geno<-t(genotype)-2
    f2<-function(x){
      gt<-x;eff<-matrix(xq[,2 - gt[qq[1]]] ,nrow =2)
      yq= bq[1,] %*% eff
      w1 = b0 + yq + rnorm(1) * sqrt(v1)
      s1 <- c(w1, gt)
      return(s1)
    }
    genotype.f2<-apply(geno,2,f2)
    F2.result[[1]]<- genotype.f2
    F2.result[[2]]<-sim.f2.gen.data$geno[[1]]$map
    return(F2.result)
  }
  
  BC.data<-function(map.data,p.type,nsam,vg,qq,b0,v1 ){
    BC.result<-list(NULL)
    rd<-0;a<-sqrt(4*vg/(rd^2+2));
    #sim.bc.gen.data<-gen.f("bc")
    sim.bc.gen.data<-gen.f(map.data,p.type,nsam)
    genotype <- sim.bc.gen.data$geno[[1]]$data
    geno<-t(genotype)
    geno[which(geno=="2")]<--1
    bc<-function(x){
      gt<-x
      yq= a %*% gt[qq[1]]
      w1 = b0 + yq + rnorm(1) * sqrt(v1)
      s1 <- c(w1, gt)
      return(s1)
    }
    genotype.bc<-apply(geno,2,bc)
    BC.result[[1]]<- genotype.bc
    BC.result[[2]]<- sim.bc.gen.data$geno[[1]]$map
    return(BC.result)
  }
  
  DH.data<-function(map.data,p.type,nsam,vg,qq,b0,v1 ){
    DH.result<-list(NULL)
    rd<-0;a<-sqrt(4*vg/(rd^2+2));
    # sim.dh.gen.data<-gen.f("bc")
    sim.dh.gen.data<-gen.f(map.data,p.type,nsam)
    genotype <- sim.dh.gen.data$geno[[1]]$data
    geno<-t(genotype)
    geno[which(geno=="2")]<--1
    dh<-function(x){
      gt<-x
      yq= a %*% gt[qq[1]]
      w1 = b0 + yq + rnorm(1) * sqrt(v1)
      s1 <- c(w1, gt)
      return(s1)
    }
    genotype.dh<-apply(geno,2,dh)
    DH.result[[1]]<-genotype.dh
    DH.result[[2]]<-sim.dh.gen.data$geno[[1]]$map
    return(DH.result)
  }
  
  
  RIL.data<-function(map.data,p.type,nsam,vg,qq,b0,v1 ){
    RIL.result<-list(NULL)
    rd<-0;a<-sqrt(4*vg/(rd^2+2));
    #sim.ril.gen.data<-gen.f("risib")
    sim.ril.gen.data<-gen.f(map.data,p.type,nsam)
    genotype <- sim.ril.gen.data$geno[[1]]$data
    geno<-t(genotype)
    geno[which(geno=="2")]<--1
    ril<-function(x){
      gt<-x
      yq= a %*% gt[qq[1]]
      w1 = b0 + yq + rnorm(1) * sqrt(v1)
      s1 <- c(w1, gt)
      return(s1)
    }
    genotype.ril<-apply(geno,2,ril)
    RIL.result[[1]]<- genotype.ril
    RIL.result[[2]]<- sim.ril.gen.data$geno[[1]]$map
    return( RIL.result)
  }
  
  dodata1<-function(data.result,nsam,np){
    ss<-data.result[[1]]
    ott1<-ss[,order(ss[1,])];ss<-ott1
    phe<-as.numeric(ss[1,]);sx=1:nsam
    resam<-sample(x=sx);newphe<-phe[resam]
    ss1<-rbind(newphe,ss[-1,])
    ott<-ss1[,order(ss1[1,])]
    l12<-ott[,1:ceiling(nsam*np)];
    h12<-ott[,(dim(ott)[2]-(ceiling(nsam*np)-1)):dim(ott)[2]]
    frq<-function(m){
      nrow=dim(ott)[1]-1
      f<-matrix(0,ncol=3,nrow=nrow)
      for(i in 1:nrow){
        f[i,1]<-sum(m[i,]==1)
        f[i,2]<-sum(m[i,]==0)
        f[i,3]<-sum(m[i,]==-1)
      }
      fA<-f[,1]*2+f[,2]
      fa<-f[,2]+f[,3]*2
      f1<-cbind(fA,fa)
      return(f1)
    }
    lfrq12<-frq(l12[-1,]);hfrq12<-frq(h12[-1,])
    res12<-cbind(lfrq12,hfrq12)
    result<-cbind(data.result[[2]],rep(1,500),data.result[[2]],res12)
    colnames(result)<-c("Maker","Chromosome","Position","AL","aL","AH","aH")
    return(result)
  }
  
  dodata2<-function(data.result,nsam,np){
    ss<-data.result[[1]]
    ott1<-ss[,order(ss[1,])];ss<-ott1
    phe<-as.numeric(ss[1,]);sx=1:nsam
    resam<-sample(x=sx);newphe<-phe[resam]
    ss1<-rbind(newphe,ss[-1,])
    ott<-ss1[,order(ss1[1,])]
    l12<-ott[,1:ceiling(nsam*np)];
    h12<-ott[,(dim(ott)[2]-(ceiling(nsam*np)-1)):dim(ott)[2]]
    frq<-function(m){
      nrow=dim(ott)[1]-1
      f<-matrix(0,ncol=2,nrow=nrow)
      for(i in 1:nrow){
        f[i,1]<-sum(m[i,]==1)
        f[i,2]<-sum(m[i,]==-1)
      }
      fA<-f[,1]*2+f[,2]
      fa<-f[,2]
      f1<-cbind(fA,fa)
      return(f1)
    }
    lfrq12<-frq(l12[-1,]);hfrq12<-frq(h12[-1,])
    res12<-cbind(lfrq12,hfrq12)
    result<-cbind(data.result[[2]],rep(1,500),data.result[[2]],res12)
    colnames(result)<-c("Maker","Chromosome","Position","AL","aL","AH","aH")
    return(result)
  }
  
  
  dodata3<-function(data.result,nsam,np){
    ss<-data.result[[1]]
    ott1<-ss[,order(ss[1,])];ss<-ott1
    phe<-as.numeric(ss[1,]);sx=1:nsam
    resam<-sample(x=sx);newphe<-phe[resam]
    ss1<-rbind(newphe,ss[-1,])
    ott<-ss1[,order(ss1[1,])]
    l12<-ott[,1:ceiling(nsam*np)];
    h12<-ott[,(dim(ott)[2]-(ceiling(nsam*np)-1)):dim(ott)[2]]
    frq<-function(m){
      nrow=dim(ott)[1]-1
      f<-matrix(0,ncol=2,nrow=nrow)
      for(i in 1:nrow){
        f[i,1]<-sum(m[i,]==1)
        f[i,2]<-sum(m[i,]==-1)
      }
      fA<-f[,1]*2
      fa<-f[,2]*2
      f1<-cbind(fA,fa)
      return(f1)
    }
    lfrq12<-frq(l12[-1,]);hfrq12<-frq(h12[-1,])
    res12<-cbind(lfrq12,hfrq12)
    result<-cbind(data.result[[2]],rep(1,500),data.result[[2]],res12)
    colnames(result)<-c("Maker","Chromosome","Position","AL","aL","AH","aH")
    return(result)
  }
  
  dodata4<-function(data.result,nsam,np){
    #data.result<-F2.data("f2")
    ss<-data.result[[1]]
    ott1<-ss[,order(ss[1,])];ss<-ott1
    phe<-as.numeric(ss[1,]);sx=1:nsam
    resam<-sample(x=sx);newphe<-phe[resam]
    ss1<-rbind(newphe,ss[-1,])
    ott<-ss1[,order(ss1[1,])]
    l12<-ott[,1:ceiling(nsam*np)];
    h12<-ott[,(dim(ott)[2]-(ceiling(nsam*np)-1)):dim(ott)[2]]
    frq<-function(m){
      nrow=dim(ott)[1]-1
      f<-matrix(0,ncol=3,nrow=nrow)
      for(i in 1:nrow){
        f[i,1]<-sum(m[i,]==1)
        f[i,2]<-sum(m[i,]==0)
        f[i,3]<-sum(m[i,]==-1)
      }
      fA<-f[,1]*2+f[,2]
      fa<-f[,2]+f[,3]*2
      f1<-cbind(f,fA,fa)
      return(f1)
    }
    lfrq12<-frq(l12[-1,]);hfrq12<-frq(h12[-1,])
    res12<-cbind(lfrq12,hfrq12)
    result<-cbind(data.result[[2]],rep(1,500),data.result[[2]],res12[,4:5],res12[,9:10],res12[,1:3],res12[,6:8])
    colnames(result)<-c("Maker","Chromosome","Position","AL","aL","AH","aH","AAL","aaL","AaL","AAH","aaH","AaH")
    return(result)
  }
  
  dQTGseq2.p<-function(x){
    #x<-list(dodata4(F2.data(map.data,"f2",nsam)))
    result<-dQTGseq2(x)
    f.result<-result[order(as.numeric(result))][0.95*500]
    return(f.result)
  }
  
  G.p<- function(x) {
    result<-G.function(x)
    f.result<-result[order(as.numeric(result))][500]
    return(f.result)
  }
  
  SNP.p<- function(x) {
    result<-SNP.function(x)
    f.result<-result[order(as.numeric(result))][500]
    return(f.result)
  }
  
  
  
  ED.p<- function(x) {
    result<-ED.function(x)
    f.result<-result[order(as.numeric(result))][500]
    return(f.result)
  }
  
  SmoothLOD.p<- function(x) {
    result<-SmoothLOD.function(x)
    f.result<-result[order(as.numeric(result))][500]
    return(f.result)
  }
  
  draw1<-function(result.sta,Fileformat,population,color,Smooth_method,threshold){
    
    names(result.sta)<-c("Maker","Chromosome","Position","GW","Smooth_Gw","LOD","Smooth_LOD","G","Smooth_G","deltaSNP","Smooth_deltaSNP","ED","Smooth_ED")
    margin_space<-0.55;mgp4<-c(3,1,1);mgp2<-c(3,1.15,1);line1<-2.3;line2<-2.15
    lwd=0.8;pch=15;cexm3<-0.6;
    #color<-c("blue")
    col1<-color
    layout(matrix(c(rep(1,7),rep(2,7),rep(3,7),rep(4,7),rep(5,7),rep(6,1)),ncol = 1))
    #col1<-c("blue","red")
    d1 <- result.sta[order(as.numeric(result.sta$Chromosome),as.numeric(result.sta$Position)), ]
    d1<-d1[,-1]
    ff<-function(x){
      y<-as.numeric(x)
      return(y)
    }
    d<-apply(d1, 2, ff)
    d<-as.data.frame(d)
    index<-c("index")
    d$pos=NA
    d$index=NA
    ind = 0
    for (i in unique(d$Chromosome)){
      ind = ind + 1
      d[d$Chromosome==i,]$index = ind
    }
    nchr = length(unique(d$Chromosome))
    if (nchr==1) { ## For a single chromosome
      d$pos=as.numeric(d$Position)/1000000
      xlabel = paste('Chromosome',unique(d$Chromosome),'(Mb)')
    } else { ## For multiple chromosomes
      lastbase=0
      ticks=NULL; ticks1=NULL;
      for (i in unique(d$index)) {
        if (i==1) {
          d[d$index==i, ]$pos = d[d$index==i, ]$Position
        } else {
          lastbase=lastbase+tail(as.numeric(subset(d,index==i-1)$Position), 1)
          d[d$index==i, ]$pos=d[d$index==i, ]$Position+lastbase
        }
        ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
        ticks1 = c(ticks1,  max(d[d$index == i,]$pos ))
      }
      labs <- unique(d$Chromosome)
    }
    xmax = ceiling(max(d$pos))
    xmin = floor(min(d$pos))
    ymax<-ceiling(max(d$GW))
    margin_space<-margin_space
    
    
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=pch,
                                      xlim=c(xmin,xmax),ylim=c(0,ymax),
                                      xlab="", cex.lab=0.6,ylab=""))
    dotargs <- as.list(match.call())[-1L]
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      with(d,  points(d$pos, GW, col="LightGray", pch=pch,cex=0.2))
    } else {
      icol=1
      for (i in unique(d$index)) {
        #i<-2
        with(d[d$index==unique(d$index)[i], ], points(pos, GW, col="LightGray", pch=pch,cex=0.25))
        icol=icol+1
      }
    }
    suppressWarnings(axis(4,ylim=c(0,ceiling(max(d$GW))),cex.lab=0.6,mgp=mgp4,tcl=-0.15,cex.axis=0.7,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray",lwd=lwd))
    mtext((expression(G[w])),side=4,line=line1,cex=0.5,font = 1,col="DarkGray")
    par(new=TRUE,mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=pch, yaxt="n" ,
                                      xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$Smooth_Gw))),
                                      xlab="", ylab="",cex.lab=0.6,font = 1))#line=0.78,
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (nchr==1) {
      suppressWarnings(axis(1,cex.lab=0.6,tcl=-0.2,cex.axis=0.4))
    } else {
      suppressWarnings(axis(1, at=c(0,ticks1), labels=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
      suppressWarnings(axis(1, at=ticks, labels= labs,tick=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
    }
    #box(which ="plot",bty="l",col ='black',lwd = 0.8)
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      if(Smooth_method=="None"){
        with(d, lines(d$pos,d$Smooth_Gw, type = "l", cex=0.2,col=col[1]))
      }else{
        with(d, points(d$pos,d$Smooth_Gw, pch=pch, cex=0.2,col=col[1]))
      }
    } else {
      icol=1
      for (i in unique(d$index)) {
        if(Smooth_method=="None"){
          with(d[d$index==unique(d$index)[i], ], lines(pos,Smooth_Gw, col=col[icol],type = "l",cex=0.25))
        }else{
          with(d[d$index==unique(d$index)[i], ], points(pos,Smooth_Gw, col=col[icol],pch=pch,cex=0.25))
        }
        icol=icol+1
      }
    }
    abline(h=threshold[1], lty=3,lwd = 2.5,col="black")
    if(nchr!=1){
      abline(v=ticks1,lty=3,lwd = 2.5,col="black")
    }
    
    axis(2,ylim=c(0,ceiling(max(d$Smooth_Gw))),cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    mtext((expression(Smooth_G[w])),side=2,line=line2,cex=0.5,font = 1,col = col[1],at=max(d$Smooth_Gw)/2)
    #mtext("Smooth",side=2,line=line2,cex=0.5,font = 1,col = col[1],at=ymax/6)
    if(Fileformat=="BSA" & population=="F2"){
      mtext("A  dQTG-seq1",side=3,col="black",cex=cexm3,font=2,at=xmin+(xmax-xmin)/20)
    } else {
      mtext("A  dQTG-seq2",side=3,col="black",cex=cexm3,font=2,at=xmin+(xmax-xmin)/20)
    }

    
    
    
    xmax = ceiling(max(d$pos))
    xmin = floor(min(d$pos))
    ymax<-ceiling(max(d$LOD))
    margin_space<-margin_space
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=pch,
                                      xlim=c(xmin,xmax),ylim=c(0,ymax),
                                      xlab="", cex.lab=0.6,ylab=""))
    dotargs <- as.list(match.call())[-1L]
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      with(d,  points(pos, d$LOD, col="LightGray", pch=pch,cex=0.2))
    } else {
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, LOD, col="LightGray", pch=pch,cex=0.25))
        icol=icol+1
      }
    }
    suppressWarnings(axis(4,ylim=c(0,ceiling(max(d$LOD))),cex.lab=0.6,mgp=mgp4,tcl=-0.15,cex.axis=0.7,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray",lwd=lwd))
    mtext("LOD",side=4,line=line1,cex=0.5,font = 1,col="DarkGray")
    par(new=TRUE,mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=pch, yaxt="n" ,
                                      xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$Smooth_LOD))),
                                      xlab="", ylab="",cex.lab=0.6,font = 1))#line=0.78,
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (nchr==1) {
      suppressWarnings(axis(1,cex.lab=0.6,tcl=-0.2,cex.axis=0.4))
    } else {
      suppressWarnings(axis(1, at=c(0,ticks1), labels=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
      suppressWarnings(axis(1, at=ticks, labels= labs,tick=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
    }
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      if(Smooth_method=="None"){
        with(d, lines(pos,d$Smooth_LOD, type="l", cex=0.2,col=col[1]))
      }else{
        with(d, points(pos,d$Smooth_LOD, pch=pch, cex=0.2,col=col[1]))
      }
    } else {
      icol=1
      for (i in unique(d$index)) {
        if(Smooth_method=="None"){
          with(d[d$index==unique(d$index)[i], ], lines(pos, Smooth_LOD, col=col[icol],type="l",cex=0.25))
        }else{
          with(d[d$index==unique(d$index)[i], ], points(pos, Smooth_LOD, col=col[icol],pch=pch,cex=0.25))
        }
        icol=icol+1
      }
    }
    abline(h=threshold[2], lty=3,lwd = 2.5,col="black")
    if(nchr!=1){
      abline(v=ticks1,lty=3,lwd = 2.5,col="black")
    }
    axis(2,ylim=c(0,ceiling(max(d$Smooth_LOD))),cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    mtext("Smooth_LOD",side=2,line=line2,cex=0.5,font = 1,col = col[1])
    mtext("B  Smooth_LOD",col="black",cex=cexm3,font=2,side = 3,at=xmin+(xmax-xmin)/18)

    
    
    xmax = ceiling(max(d$pos))
    xmin = floor(min(d$pos))
    margin_space<-margin_space
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=pch,
                                      xlim=c(xmin,xmax),ylim=c(0,ceiling(max(d$G))),
                                      xlab="", cex.lab=0.6,ylab=""))
    dotargs <- as.list(match.call())[-1L]
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    col=rep(col1, max(d$Chromosome))
    # Add points to the plot
    if (nchr==1) {
      with(d,  points(pos, d$G, col="LightGray", pch=pch,cex=0.2))
    } else {
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, G,col="LightGray", pch=pch,cex=0.25))
        icol=icol+1
      }
    }
    suppressWarnings(axis(4,ylim=c(0,ceiling(max(d$G))),cex.lab=0.6,mgp=mgp4,tcl=-0.15,cex.axis=0.7,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray",lwd=lwd))
    mtext("G",side=4,line=line1,cex=0.5,font = 1,col="DarkGray")
    par(new=TRUE,mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=pch, yaxt="n" ,
                                      xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$Smooth_G))),
                                      xlab="", ylab="",cex.lab=0.6,font = 1))#line=0.78,
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (nchr==1) {
      suppressWarnings(axis(1,cex.lab=0.6,tcl=-0.2,cex.axis=0.4))
    } else {
      suppressWarnings(axis(1, at=c(0,ticks1), labels=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
      suppressWarnings(axis(1, at=ticks, labels= labs,tick=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
    }
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      if(Smooth_method=="None"){
        with(d, lines(pos, Smooth_G, type="l", cex=0.1,col=col[1]))
      }else{
        with(d, points(pos, Smooth_G, pch=pch, cex=0.1,col=col[1]))
      }
    } else {
      icol=1
      for (i in unique(d$index)) {
        if(Smooth_method=="None"){
          with(d[d$index==unique(d$index)[i], ], lines(pos,Smooth_G, col=col[icol],type="l",cex=0.25))
        }else{
          with(d[d$index==unique(d$index)[i], ], points(pos,Smooth_G, col=col[icol],pch=pch,cex=0.25))
        }
        icol=icol+1
      }
    }
    abline(h=threshold[3], lty=3,lwd = 2.5,col="black")
    if(nchr!=1){
      abline(v=ticks1,lty=3,lwd = 2.5,col="black")
    }
    axis(2,ylim=c(0,ceiling(max(d$Smooth_G))),cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    mtext("G'",side=2,line=line2,cex=0.5,font = 1,col = col[1])
    mtext("C  G'",col="black",cex=cexm3,font=2,side = 3,at=xmin+(xmax-xmin)/36)
    #on.exit(par(old.par))
    
    
    
    
    xmax = ceiling(max(d$pos))
    xmin = floor(min(d$pos))
    margin_space<-margin_space
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=pch,
                                      xlim=c(xmin,xmax),ylim=c(min(d$deltaSNP),max(d$deltaSNP)),
                                      xlab="", cex.lab=0.6,ylab=""))
    dotargs <- as.list(match.call())[-1L]
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    col=rep(col1, max(d$Chromosome))
    # Add points to the plot
    if (nchr==1) {
      with(d,  points(pos,d$deltaSNP, col="LightGray", pch=pch,cex=0.2))
    } else {
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, deltaSNP, col="LightGray", pch=pch,cex=0.25))
        icol=icol+1
      }
    }
    suppressWarnings(axis(4,at=c(round(min(d$deltaSNP),1),round(min(d$deltaSNP),1)/2,0,round(max(d$deltaSNP),1)/2,round(max(d$deltaSNP),1)),
                          cex.lab=0.6,mgp=mgp4,tcl=-0.15,cex.axis=0.7,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray",lwd=lwd))
    mtext("deltaSNP",side=4,line=line1,cex=0.5,font = 1,col="DarkGray")
    par(new=TRUE,mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=pch, yaxt="n" ,
                                      xlim=c(xmin,xmax), ylim=c(min(d$Smooth_deltaSNP),max(d$Smooth_deltaSNP)),
                                      xlab="", ylab="",cex.lab=0.6,font = 1))#line=0.78,
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (nchr==1) {
      suppressWarnings(axis(1,cex.lab=0.6,tcl=-0.2,cex.axis=0.4))
    } else {
      suppressWarnings(axis(1, at=c(0,ticks1), labels=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
      suppressWarnings(axis(1, at=ticks, labels= labs,tick=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
    }
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      if(Smooth_method=="None"){
        with(d, lines(pos,d$Smooth_deltaSNP, type="l", cex=0.1,col=col[1]))
      }else{
        with(d, points(pos,d$Smooth_deltaSNP, pch=pch, cex=0.1,col=col[1]))
      }
    } else {
      icol=1
      for (i in unique(d$index)) {
        if(Smooth_method=="None"){
          with(d[d$index==unique(d$index)[i], ], lines(pos, Smooth_deltaSNP, col=col[icol],type="l",cex=0.25))
        }else{
          with(d[d$index==unique(d$index)[i], ], points(pos, Smooth_deltaSNP, col=col[icol],pch=pch,cex=0.25))
        }
        icol=icol+1
      }
    }
    abline(h=threshold[4], lty=3,lwd = 2.5,col="black")
    
    if(nchr!=1){
      abline(v=ticks1,lty=3,lwd = 2.5,col="black")
    }
    # axis(2,ylim=c(round(min(d$deltaSNP),1),round(max(d$deltaSNP),1)),cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    axis(2,at=c(round(min(d$Smooth_deltaSNP),1),round(min(d$Smooth_deltaSNP),1)/2,0,round(max(d$Smooth_deltaSNP),1)/2,round(max(d$Smooth_deltaSNP),1)),
         cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    
    mtext("Smooth_deltaSNP",side=2,line=line2,cex=0.5,font = 1,col = col[1])
    mtext("D  deltaSNP",col="black",cex=cexm3,font=2,side = 3,at=xmin+(xmax-xmin)/21)
    #on.exit(par(old.par))
    
    
    
    xmax = ceiling(max(d$pos))
    xmin = floor(min(d$pos))
    margin_space<-margin_space
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=pch,
                                      xlim=c(xmin,xmax),ylim=c(0,max(d$ED)),
                                      xlab="", cex.lab=0.6,ylab=""))
    dotargs <- as.list(match.call())[-1L]
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    col=rep(col1, max(d$Chromosome))
    # Add points to the plot
    if (nchr==1) {
      with(d,  points(pos, d$ED, col="LightGray", pch=pch,cex=0.2))
    } else {
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, ED, col="LightGray", pch=pch,cex=0.25))
        icol=icol+1
      }
    }
    suppressWarnings(axis(4,ylim=c(0,max(d$ED)),cex.lab=0.6,mgp=mgp4,tcl=-0.15,cex.axis=0.7,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray",lwd=lwd))
    mtext("ED",side=4,line=line1,cex=0.5,font = 1,col="DarkGray")
    par(new=TRUE,mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=pch, yaxt="n" ,
                                      xlim=c(xmin,xmax), ylim=c(0,max(d$Smooth_ED)),
                                      xlab="", ylab="",cex.lab=0.6,font = 1))#line=0.78,
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (nchr==1) {
      suppressWarnings(axis(1,cex.lab=0.6,tcl=-0.2,cex.axis=0.4))
    } else {
      suppressWarnings(axis(1, at=c(0,ticks1), labels=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
      suppressWarnings(axis(1, at=ticks, labels= labs,tick=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
    }
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      if(Smooth_method=="None"){
        with(d, lines(d$pos,d$Smooth_ED, type="l", cex=0.1,col=col[1]))
      }else{
        with(d, points(d$pos,d$Smooth_ED, pch=pch, cex=0.1,col=col[1]))
      }
    } else {
      icol=1
      for (i in unique(d$index)) {
        if(Smooth_method=="None"){
          with(d[d$index==unique(d$index)[i], ], lines(pos, Smooth_ED, col=col[icol],type="l",cex=0.25))
        }else{
          with(d[d$index==unique(d$index)[i], ], points(pos, Smooth_ED, col=col[icol],pch=pch,cex=0.25))
        }
        icol=icol+1
      }
    }
    abline(h=threshold[5], lty=3,lwd = 2.5,col="black")
    if(nchr!=1){
      abline(v=ticks1,lty=3,lwd = 2.5,col="black")
    }
    axis(2,ylim=c(0,max(d$Smooth_ED)),cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    mtext("Smooth_ED",side=2,line=line2,cex=0.5,font = 1,col = col[1])
    mtext("E   ED",col="black",cex=cexm3,font=2,side = 3,at=xmin+(xmax-xmin)/30)
    #on.exit(par(old.par))
    
    
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3*margin_space,margin_space,3*margin_space)+0.8*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    if (nchr==1) {
      mtext(xlabel,side=1,line=1.4,cex=0.6,font=2)#linepos
    } else {
      mtext("Chromosome",side=1,line=1.4,cex=0.6,font=2)
    }
    #on.exit(par(old.par))
    
    
    
  }
  
  draw2<-function(result.sta,population,color,Smooth_method,threshold){
    names(result.sta)<-c("Maker","Chromosome","Position","LOD","Smooth_LOD","G","Smooth_G","deltaSNP","Smooth_deltaSNP","ED","Smooth_ED")
    margin_space<-0.55;mgp4<-c(3,1,1);mgp2<-c(3,1.15,1);line1<-2.3;line2<-2.15
    lwd=0.8;pch=15;cexm3<-0.6;
    #color<-c("blue")
    col1<-color
    layout(matrix(c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,1)),ncol = 1))
    
    
    d1 <- result.sta[order(as.numeric(result.sta$Chromosome),as.numeric(result.sta$Position)), ]
    d1<-d1[,-1]
    fff<-function(x){
      y<-as.numeric(x)
      return(y)
    }
    d<-apply(d1, 2, fff)
    d<-as.data.frame(d)
    index=c("index")
    d$pos=NA
    d$index=NA
    ind = 0
    for (i in unique(d$Chromosome)){
      ind = ind + 1
      d[d$Chromosome==i,]$index = ind
    }
    nchr = length(unique(d$Chromosome))
    if (nchr==1) { ## For a single chromosome
      d$pos=as.numeric(d$Position)/1000000
      xlabel = paste('Chromosome',unique(d$Chromosome),'(Mb)')
    } else { ## For multiple chromosomes
      lastbase=0
      ticks=NULL; ticks1=NULL;
      for (i in unique(d$index)) {
        if (i==1) {
          d[d$index==i, ]$pos = d[d$index==i, ]$Position
        } else {
          lastbase=lastbase+tail(as.numeric(subset(d,index==i-1)$Position), 1)
          d[d$index==i, ]$pos=d[d$index==i, ]$Position+lastbase
        }
        ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
        ticks1 = c(ticks1,  max(d[d$index == i,]$pos ))
        
      }
      labs <- unique(d$Chromosome)
    }
    
    
    
    xmax = ceiling(max(d$pos))
    xmin = floor(min(d$pos))
    ymax<-ceiling(max(d$LOD))
    margin_space<-margin_space
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=pch,
                                      xlim=c(xmin,xmax),ylim=c(0,ymax),
                                      xlab="", cex.lab=0.6,ylab=""))
    dotargs <- as.list(match.call())[-1L]
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      
      with(d,  points(pos, d$LOD, col="LightGray", pch=pch,cex=0.2))
      
    } else {
      icol=1
      for (i in unique(d$index)) {
        
        with(d[d$index==unique(d$index)[i], ], points(pos, LOD, col="LightGray", pch=pch,cex=0.25))
        
        
        icol=icol+1
      }
    }
    suppressWarnings(axis(4,ylim=c(0,ceiling(max(d$LOD))),cex.lab=0.6,mgp=mgp4,tcl=-0.15,cex.axis=0.7,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray",lwd=lwd))
    mtext("LOD",side=4,line=line1,cex=0.5,font = 1,col="DarkGray")
    
    par(new=TRUE,mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=pch, yaxt="n" ,
                                      xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$Smooth_LOD))),
                                      xlab="", ylab="",cex.lab=0.6,font = 1))#line=0.78,
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (nchr==1) {
      suppressWarnings(axis(1,cex.lab=0.6,tcl=-0.2,cex.axis=0.4))
    } else {
      suppressWarnings(axis(1, at=c(0,ticks1), labels=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
      suppressWarnings(axis(1, at=ticks, labels= labs,tick=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
    }
    box(which ="plot",bty="l",col ='black',lwd = 0.8)
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      if(Smooth_method=="None"){
        with(d, lines(pos,d$Smooth_LOD, type="l", cex=0.2,col=col[1]))
      }else{
        with(d, points(pos,d$Smooth_LOD, pch=pch, cex=0.2,col=col[1]))
      }
    } else {
      icol=1
      for (i in unique(d$index)) {
        if(Smooth_method=="None"){
          
          with(d[d$index==unique(d$index)[i], ], lines(pos, Smooth_LOD, col=col[icol],type="l",cex=0.25))
        }else{
          with(d[d$index==unique(d$index)[i], ], points(pos, Smooth_LOD, col=col[icol],pch=pch,cex=0.25))
          
        }
        icol=icol+1
      }
    }
    abline(h=threshold[1], lty=3,lwd = 2.5,col="black")
    if(nchr!=1){
      abline(v=ticks1,lty=3,lwd = 2.5,col="black")
    }
    axis(2,ylim=c(0,ceiling(max(d$Smooth_LOD))),cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    mtext("Smooth_LOD",side=2,line=line2,cex=0.5,font = 1,col = col[1])
    mtext("A  Smooth_LOD",col="black",cex=cexm3,font=2,side = 3,at=xmin+(xmax-xmin)/18)
    #on.exit(par(old.par))
    
    xmax = ceiling(max(d$pos))
    xmin = floor(min(d$pos))
    margin_space<-margin_space
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=pch,
                                      xlim=c(xmin,xmax),ylim=c(0,ceiling(max(d$G))),
                                      xlab="", cex.lab=0.6,ylab=""))
    dotargs <- as.list(match.call())[-1L]
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    col=rep(col1, max(d$Chromosome))
    
    # Add points to the plot
    if (nchr==1) {
      
      with(d,  points(pos, d$G, col="LightGray", pch=pch,cex=0.2))
      
    } else {
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, G,col="LightGray", pch=pch,cex=0.25))
        icol=icol+1
      }
    }
    suppressWarnings(axis(4,ylim=c(0,ceiling(max(d$G))),cex.lab=0.6,mgp=mgp4,tcl=-0.15,cex.axis=0.7,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray",lwd=lwd))
    mtext("G",side=4,line=line1,cex=0.5,font = 1,col="DarkGray")
    par(new=TRUE,mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=pch, yaxt="n" ,
                                      xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$Smooth_G))),
                                      xlab="", ylab="",cex.lab=0.6,font = 1))#line=0.78,
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (nchr==1) {
      suppressWarnings(axis(1,cex.lab=0.6,tcl=-0.2,cex.axis=0.4))
    } else {
      suppressWarnings(axis(1, at=c(0,ticks1), labels=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
      suppressWarnings(axis(1, at=ticks, labels= labs,tick=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
    }
    box(which ="plot",bty="l",col ='black',lwd = 0.8)
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      if(Smooth_method=="None"){
        with(d, lines(pos, Smooth_G, type="l", cex=0.1,col=col[1]))
      }else{
        with(d, points(pos, Smooth_G, pch=pch, cex=0.1,col=col[1]))
      }
    } else {
      icol=1
      for (i in unique(d$index)) {
        if(Smooth_method=="None"){
          with(d[d$index==unique(d$index)[i], ], lines(pos,Smooth_G, col=col[icol],type="l",cex=0.25))
        }else{
          with(d[d$index==unique(d$index)[i], ], points(pos,Smooth_G, col=col[icol],pch=pch,cex=0.25))
        }
        icol=icol+1
      }
    }
    abline(h=threshold[2], lty=3,lwd = 2.5,col="black")
    if(nchr!=1){
      abline(v=ticks1,lty=3,lwd = 2.5,col="black")
    }
    axis(2,ylim=c(0,ceiling(max(d$Smooth_G))),cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    mtext("G'",side=2,line=line2,cex=0.5,font = 1,col = col[1])
    mtext("B  G'",col="black",cex=cexm3,font=2,side = 3,at=xmin+(xmax-xmin)/36)
    #on.exit(par(old.par))
    
    
    xmax = ceiling(max(d$pos))
    xmin = floor(min(d$pos))
    margin_space<-margin_space
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=pch,
                                      xlim=c(xmin,xmax),ylim=c(min(d$deltaSNP),max(d$deltaSNP)),
                                      xlab="", cex.lab=0.6,ylab=""))
    dotargs <- as.list(match.call())[-1L]
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    col=rep(col1, max(d$Chromosome))
    # Add points to the plot
    if (nchr==1) {
      with(d,  points(pos,d$deltaSNP, col="LightGray", pch=pch,cex=0.2))
    } else {
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, deltaSNP, col="LightGray", pch=pch,cex=0.25))
        icol=icol+1
      }
    }
    suppressWarnings(axis(4,at=c(round(min(d$deltaSNP),1),round(min(d$deltaSNP),1)/2,0,round(max(d$deltaSNP),1)/2,round(max(d$deltaSNP),1)),
                          cex.lab=0.6,mgp=mgp4,tcl=-0.15,cex.axis=0.7,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray",lwd=lwd))
    mtext("deltaSNP",side=4,line=line1,cex=0.5,font = 1,col="DarkGray")
    
    par(new=TRUE,mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=pch, yaxt="n" ,
                                      xlim=c(xmin,xmax), ylim=c(min(d$Smooth_deltaSNP),max(d$Smooth_deltaSNP)),
                                      xlab="", ylab="",cex.lab=0.6,font = 1))#line=0.78,
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (nchr==1) {
      suppressWarnings(axis(1,cex.lab=0.6,tcl=-0.2,cex.axis=0.4))
    } else {
      suppressWarnings(axis(1, at=c(0,ticks1), labels=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
      suppressWarnings(axis(1, at=ticks, labels= labs,tick=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
    }
    box(which ="plot",bty="l",col ='black',lwd = 0.8)
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      if(Smooth_method=="None"){
        with(d, lines(pos,d$Smooth_deltaSNP, type="l", cex=0.1,col=col[1]))
      }else{
        with(d, points(pos,d$Smooth_deltaSNP, pch=pch, cex=0.1,col=col[1]))
      }
    } else {
      icol=1
      for (i in unique(d$index)) {
        if(Smooth_method=="None"){
          with(d[d$index==unique(d$index)[i], ], lines(pos, Smooth_deltaSNP, col=col[icol],type="l",cex=0.25))
        }else{
          with(d[d$index==unique(d$index)[i], ], points(pos, Smooth_deltaSNP, col=col[icol],pch=pch,cex=0.25))
        }
        icol=icol+1
      }
    }
    abline(h=threshold[3], lty=3,lwd = 2.5,col="black")
    if(nchr!=1){
      abline(v=ticks1,lty=3,lwd = 2.5,col="black")
    }
    axis(2,at=c(round(min(d$Smooth_deltaSNP),1),round(min(d$Smooth_deltaSNP),1)/2,0,round(max(d$Smooth_deltaSNP),1)/2,round(max(d$Smooth_deltaSNP),1)),
         cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    mtext("Smooth_deltaSNP",side=2,line=line2,cex=0.5,font = 1,col = col[1])
    mtext("C  deltaSNP",col="black",cex=cexm3,font=2,side = 3,at=xmin+(xmax-xmin)/21)
    #on.exit(par(old.par))
    
    
    
    xmax = ceiling(max(d$pos))
    xmin = floor(min(d$pos))
    margin_space<-margin_space
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=pch,
                                      xlim=c(xmin,xmax),ylim=c(0,max(d$ED)),
                                      xlab="", cex.lab=0.6,ylab=""))
    dotargs <- as.list(match.call())[-1L]
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    col=rep(col1, max(d$Chromosome))
    # Add points to the plot
    if (nchr==1) {
      with(d,  points(pos, d$ED, col="LightGray", pch=pch,cex=0.2))
    } else {
      icol=1
      for (i in unique(d$index)) {
        with(d[d$index==unique(d$index)[i], ], points(pos, ED, col="LightGray", pch=pch,cex=0.25))
        icol=icol+1
      }
    }
    suppressWarnings(axis(4,ylim=c(0,max(d$ED)),cex.lab=0.6,mgp=mgp4,tcl=-0.15,cex.axis=0.7,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray",lwd=lwd))
    mtext("ED",side=4,line=line1,cex=0.5,font = 1,col="DarkGray")
    
    par(new=TRUE,mar=c(3*margin_space,3.5*margin_space,margin_space,5*margin_space)+2*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=pch, yaxt="n" ,
                                      xlim=c(xmin,xmax), ylim=c(0,max(d$Smooth_ED)),
                                      xlab="", ylab="",cex.lab=0.6,font = 1))#line=0.78,
    dotargs <- list()
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (nchr==1) {
      suppressWarnings(axis(1,cex.lab=0.6,tcl=-0.2,cex.axis=0.4))
    } else {
      suppressWarnings(axis(1, at=c(0,ticks1), labels=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
      suppressWarnings(axis(1, at=ticks, labels= labs,tick=FALSE,cex.lab=0.6,tcl=-0.15,cex.axis=0.7,lwd=lwd,mgp=c(3,0.2,0.4)))#第三个调坐标轴
    }
    box(which ="plot",bty="l",col ='black',lwd = 0.8)
    col=rep(col1, max(d$Chromosome))
    if (nchr==1) {
      if(Smooth_method=="None"){
        with(d, lines(d$pos,d$Smooth_ED, type="l", cex=0.1,col=col[1]))
      }else{
        with(d, points(d$pos,d$Smooth_ED, pch=pch, cex=0.1,col=col[1]))
      }
    } else {
      
      icol=1
      for (i in unique(d$index)) {
        if(Smooth_method=="None"){
          with(d[d$index==unique(d$index)[i], ], lines(pos, Smooth_ED, col=col[icol],type="l",cex=0.25))
        }else{
          with(d[d$index==unique(d$index)[i], ], points(pos, Smooth_ED, col=col[icol],pch=pch,cex=0.25))
        }
        icol=icol+1
      }
    }
    abline(h=threshold[4], lty=3,lwd = 2.5,col="black")
    if(nchr!=1){
      abline(v=ticks1,lty=3,lwd = 2.5,col="black")
    }
    axis(2,ylim=c(0,max(d$Smooth_ED)),cex.lab=0.6,mgp=mgp2,tcl=-0.15,cex.axis=0.7,col.axis = col[1],col.ticks=col[1],col =col[1],lwd=lwd )
    mtext("Smooth_ED",side=2,line=line2,cex=0.5,font = 1,col = col[1])
    mtext("D   ED",col="black",cex=cexm3,font=2,side = 3,at=xmin+(xmax-xmin)/30)
    #on.exit(par(old.par))
    
    
    
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar=c(3*margin_space,3*margin_space,margin_space,3*margin_space)+0.8*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
    if (nchr==1) {
      mtext(xlabel,side=1,line=1.4,cex=0.6,font=2)#linepos
    } else {
      mtext("Chromosome",side=1,line=1.4,cex=0.6,font=2)
    }
  }
  #on.exit(par(old.par))
  
  G.function <- function(data.G){
    #data.G<-extreme.result
    data.G<-as.data.frame(data.G)
    Al= as.numeric(data.G$AL); al= as.numeric(data.G$aL);
    Ah= as.numeric(data.G$AH); ah= as.numeric(data.G$aH);
    exp <- c(
      (Al + Ah) * (Al + al) / (Al + Ah + al + ah),
      (Al + Ah) * (Ah + ah) / (Al + Ah + al + ah),
      (Al + al) * (al + ah) / (Al + Ah + al + ah),
      (al + ah) * (Ah + ah) / (Al + Ah + al + ah)
    )
    obs <- c(Al, Ah, al, ah)
    G <-2 * (rowSums((obs+0.00000001) * log(
      matrix((obs+0.00000001), ncol = 4) / matrix((exp+0.00000001), ncol = 4)
    )))
    result.G<-G
    return(result.G)
  }
  
  SNP.function <- function(data.SNP){
    #data.SNP<-extreme.result
    data.SNP<-as.data.frame(data.SNP)
    Al= as.numeric(data.SNP$AL); al= as.numeric(data.SNP$aL);
    Ah= as.numeric(data.SNP$AH); ah= as.numeric(data.SNP$aH);
    deltaSNP<-Al/(Al+al)-Ah/(Ah+ah)
    result.SNP<-deltaSNP
    return(result.SNP)
  }
  
  ED.function <- function(data.ED){
    data.ED<-as.data.frame(data.ED)
    Al= as.numeric(data.ED$AL); al= as.numeric(data.ED$aL);
    Ah= as.numeric(data.ED$AH); ah= as.numeric(data.ED$aH);
    
    ED<-function(a,b,c,d){
      ED<-sqrt((a/(a+b)-c/(c+d))^2+(b/(a+b)-d/(c+d))^2)^2
      return (ED)
    }
    ED<-ED(Al,al,Ah,ah)
    result.ED<-ED
    return(result.ED)
  }
  
  SmoothLOD.function <- function(data.SmoothLOD){
    #data.SmoothLOD<-allel.f
    data.SmoothLOD<-as.data.frame(data.SmoothLOD)
    Al= as.numeric(data.SmoothLOD$AL); al= as.numeric(data.SmoothLOD$aL);
    Ah= as.numeric(data.SmoothLOD$AH); ah= as.numeric(data.SmoothLOD$aH);
    
    al[which(al==0)]<-0.00000001
    Al[which(Al==0)]<-0.00000001
    ah[which(ah==0)]<-0.00000001
    Ah[which(Ah==0)]<-0.00000001
    P_QTGL<-Al/(Al+al);P_QTGH<-Ah/(Ah+ah)
    P_QTGL<-matrix(P_QTGL,ncol = 1);P_QTGH<-matrix(P_QTGH,ncol = 1)
    nn<-dim(P_QTGH)[1]
    nQTG<-cbind(Al,Ah)
    LOD<-matrix(NA,nrow = nn,ncol = 1)
    for(i in 1:nn){
      #i<-10
      c1<-Al+al;c2<-Ah+ah
      LOD[i]<-nQTG[i,1]*log10(P_QTGL[i]/0.5)+(c1[i]-nQTG[i,1])*log10((1-P_QTGL[i])/0.5)+
        nQTG[i,2]*log10(P_QTGH[i]/0.5)+(c2[i]-nQTG[i,2])*log10((1-P_QTGH[i])/0.5)
    }
    result.SmoothLOD<-LOD
    return(result.SmoothLOD)
  }
  
  getG1 <- function(allel){
    Al= as.numeric(allel$AL); al=as.numeric( allel$aL);
    Ah= as.numeric(allel$AH); ah= as.numeric(allel$aH);
    exp <- c(
      (Al + Ah) * (Al + al) / (Al + Ah + al + ah),
      (Al + Ah) * (Ah + ah) / (Al + Ah + al + ah),
      (Al + al) * (al + ah) / (Al + Ah + al + ah),
      (al + ah) * (Ah + ah) / (Al + Ah + al + ah)
    )
    obs <- c(Al, Ah, al, ah)
    
    G <-2 * (rowSums((obs+0.00000001) * log(
      matrix(obs+0.00000001, ncol = 4) / matrix((exp+0.00000001), ncol = 4)
    )))
    return(G)
  }
  
  getG2 <- function(genotype){
    nAl<-as.numeric(genotype[1]);nBl<-as.numeric(genotype[2]);nHl<-as.numeric(genotype[3])
    nAh<-as.numeric(genotype[4]);nBh<-as.numeric(genotype[5]);nHh<-as.numeric(genotype[6])
    if((nAl+nBl+nHl+nAh+nBh+nHh)==0){
      G<-0
    }else{
      exp <- c(
        (nAl + nAh) * (nAl+nBl+nHl) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nAl + nAh) * (nAh+nBh+nHh) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nBl + nBh) * (nAl+nBl+nHl) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nBl + nBh) * (nAh+nBh+nHh) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nHl + nHh) * (nAl+nBl+nHl) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nHl + nHh) * (nAh+nBh+nHh) / (nAl+nBl+nHl+nAh+nBh+nHh)
      )
      obs <- c(nAl,nAh,nBl,nBh,nHl,nHh)
      
      G <-2 * (rowSums((obs+0.00000001) * log(
        matrix(obs+0.00000001, ncol = 6) / matrix((exp+0.00000001), ncol = 6)
      )))
    }
    return(G)
  }
  
  dQTGseq1<-function(resallel,np,windowsize){
    chr = as.numeric(resallel[,2]);mm<-unique(chr);allmap<-as.numeric(resallel[,3])/windowsize;
    resallel<-cbind(resallel[,1:3],allmap,resallel[4:7])
    estimate_genotype<-list(NULL)
    
    for(j in 1:length(mm)){
      #j<-1
      r<-mm[j]
      marker<-resallel[which(chr==r),1]
      pos<-as.numeric(resallel[which(chr==r),3])
      map<-as.numeric(resallel[which(chr==r),4])
      AL= as.numeric(resallel[which(chr==r),5]);aL=as.numeric(resallel[which(chr==r),6]);AH=as.numeric(resallel[which(chr==r),7]);aH=as.numeric(resallel[which(chr==r),8])
      ####genotype#####
      L<-AL/(AL+aL)
      H<-AH/(AH+aH)
      absA<-abs(L-H)
      qtlpos<-which.max(absA)
      nn<-length(L)
      n<-matrix(NA,nrow = nn,ncol =11)
      rec = matrix(NA, nrow = nn, ncol = 1)
      startx<-matrix(NA,nrow = nn,ncol = 5)
      ff1<-function(a0,d0,pL,pH,b,x){
        lis<-list(NULL);length(lis)<-1
        STEP  <- 1;
        THETA <- 1E-10;
        cal_selection_proportion <- function(b, a, d, x0){
          p <- b*pnorm(x0-a) + (1-2*b)*pnorm(x0-d) + b*pnorm(x0+a);
          return(p);
        }
        bisearch_cut_point <- function(b, a, d, p0){
          L <- -2;
          R <- 2;
          pl <- cal_selection_proportion(b, a, d, L);
          while( (pl-p0) > 0 ){
            L <- L-STEP;
            pl<- cal_selection_proportion(b, a, d, L);
          }
          pr <- cal_selection_proportion(b, a, d, R);
          while( (pr-p0) < 0 ){
            R <- R+STEP;
            pr<- cal_selection_proportion(b, a, d, R);
          }
          x0<- (L+R)/2;
          p <- cal_selection_proportion(b, a, d, x0);
          while( abs(p-p0)>THETA*p0 ){
            if( (p-p0)<0 ) L<-x0;
            if( (p-p0)>0 ) R<-x0;
            x0 <- (L+R)/2;
            p <- cal_selection_proportion(b, a, d, x0);
          }
          return(x0);
        }
        cal_y1 <- function(x, a0, d0){
          b/(sqrt(2*pi)*exp((x+a0)^2/2) )
        }
        cal_y2 <- function(x, a0, d0){
          (1-2*b)/(sqrt(2*pi)*exp((x-d0)^2/2) )
        }
        cal_y3 <- function(x, a0, d0){
          b/(sqrt(2*pi)*exp((x-a0)^2/2) )
        }
        cal_fx <- function(x, a0, d0){
          cal_y1(x, a0, d0) + cal_y2(x, a0, d0) + cal_y3(x, a0, d0)
        }
        xL <- bisearch_cut_point(b, a0, d0, pL)
        xH <- bisearch_cut_point(b, a0, d0, pH)
        point<-c(xL,xH);
        lis[[1]]<-point;
        return(lis)
      }
      
      a0<-0;d0<-0; b<-0.25;x <- seq(-2, 2, by = 0.001);pL <- np
      XL0<-ff1(a0,d0,pL,1-pL,b,x)[[1]][1]
      XH0<-ff1(a0,d0,pL,1-pL,b,x)[[1]][2]
      for (i in 1:nn) {
        #i<-6848
        rec[i] = 0.5 * (1 - exp(-2 * abs((map[i] - map[qtlpos])) / 100))#Haldane function
        n[i,]<-tryCatch ({
          estimate_by_full_mode <- function( fAH, fAL, pH, pL){
            #fAH=H[i];fAL= L[i]; pH=np;pL=np
            full_mode <- function(x){
              F <- numeric(length(x));
              F[1] <- (0.25)*pnorm(x[3]-x[1])*(1-(x[5])^2)^2+(0.5)*pnorm(x[3]-x[2])*(x[5])^2*(1-(x[5])^2)+(0.25)*pnorm(x[3]+x[1])*(x[5])^2^2+(0.25)*pnorm(x[3]-x[1])*(1-(x[5])^2)*2*(x[5])^2+(0.5)*pnorm(x[3]-x[2])*(1-2*(x[5])^2+2*(x[5])^2^2)+(0.25)*pnorm(x[3]+x[1])*2*(x[5])^2*(1-(x[5])^2)+(0.25)*pnorm(x[3]-x[1])*(x[5])^2^2+(0.5)*pnorm(x[3]-x[2])*(x[5])^2*(1-(x[5])^2)+(0.25)*pnorm(x[3]+x[1])*(1-(x[5])^2)^2- pL;
              F[2] <- (0.25)*(1-pnorm(x[4]-x[1]))*(1-(x[5])^2)^2+(0.5)*(1-pnorm(x[4]-x[2]))*(x[5])^2*(1-(x[5])^2)+(0.25)*(1-pnorm(x[4]+x[1]))*(x[5])^2^2+(0.25)*(1-pnorm(x[4]-x[1]))*(1-(x[5])^2)*2*(x[5])^2+(0.5)*(1-pnorm(x[4]-x[2]))*(1-2*(x[5])^2+2*(x[5])^2^2)+(0.25)*(1-pnorm(x[4]+x[1]))*2*(x[5])^2*(1-(x[5])^2)+(0.25)*(1-pnorm(x[4]-x[1]))*(x[5])^2^2+(0.5)*(1-pnorm(x[4]-x[2]))*(x[5])^2*(1-(x[5])^2)+(0.25)*(1-pnorm(x[4]+x[1]))*(1-(x[5])^2)^2- pH;
              F[3] <- (((0.25)*pnorm(x[3]-x[1])*(1-(x[5])^2)^2+(0.5)*pnorm(x[3]-x[2])*(x[5])^2*(1-(x[5])^2)+(0.25)*pnorm(x[3]+x[1])*(x[5])^2^2))/pL+0.5*(((0.25)*pnorm(x[3]-x[1])*(1-(x[5])^2)*2*(x[5])^2+(0.5)*pnorm(x[3]-x[2])*(1-2*(x[5])^2+2*(x[5])^2^2)+(0.25)*pnorm(x[3]+x[1])*2*(x[5])^2*(1-(x[5])^2)))/pL-fAL;
              F[4] <- (((0.25)*(1-pnorm(x[4]-x[1]))*(1-(x[5])^2)^2+(0.5)*(1-pnorm(x[4]-x[2]))*(x[5])^2*(1-(x[5])^2)+(0.25)*(1-pnorm(x[4]+x[1]))*(x[5])^2^2))/pH+0.5*(((0.25)*(1-pnorm(x[4]-x[1]))*(1-(x[5])^2)*2*(x[5])^2+(0.5)*(1-pnorm(x[4]-x[2]))*(1-2*(x[5])^2+2*(x[5])^2^2)+(0.25)*(1-pnorm(x[4]+x[1]))*2*(x[5])^2*(1-(x[5])^2)))/pH-fAH;
              F[5] <- abs(fAL-fAH)-(1-2*(x[5])^2)*abs(L[qtlpos]-H[qtlpos]);
              F;
            }
            startx <- matrix(c(0, 0, XL0, XH0, rec[i]),ncol = 5) #else startx <- matrix(c(n[i-1,7], n[i-1,8], n[i-1,9], n[i-1,10],sqrt(n[i-1,11])),ncol = 5)
            result <- dfsane(startx, full_mode,control=list(maxit=5000,trace = FALSE))
            theta  <- result$par;
            
            a <- theta[1];
            d <- theta[2];
            XL<- theta[3];
            XH<- theta[4];
            r <- (theta[5])^2;
            PAAL <- (0.25)*pnorm(XL-a)/pL*(1-r)^2+(0.5)*pnorm(XL-d)/pL*r*(1-r)+(0.25)*pnorm(XL+a)/pL*r^2;
            PAaL<-(0.25)*pnorm(XL-a)/pL*2*(1-r)*r+(0.5)*pnorm(XL-d)/pL*(1-2*r+2*r^2)+(0.25)*pnorm(XL+a)/pL*2*r*(1-r);
            PaaL<-(0.25)*pnorm(XL-a)/pL*r^2+(0.5)*pnorm(XL-d)/pL*r*(1-r)+(0.25)*pnorm(XL+a)/pL*(1-r)^2;
            nAAL<- PAAL*(AL+aL)[i]/2
            nAaL<- PAaL*(AL+aL)[i]/2
            naaL<- PaaL*(AL+aL)[i]/2
            
            PAAH<-(0.25)*(1-pnorm(XH-a))/pH*(1-r)^2+(0.5)*(1-pnorm(XH-d))/pH*r*(1-r)+(0.25)*(1-pnorm(XH+a))/pH*r^2;
            PAaH<-(0.25)*(1-pnorm(XH-a))/pH*2*(1-r)*r+(0.5)*(1-pnorm(XH-d))/pH*(1-2*r+2*r^2)+(0.25)*(1-pnorm(XH+a))/pH*2*r*(1-r);
            PaaH<-(0.25)*(1-pnorm(XH-a))/pH*r^2+(0.5)*(1-pnorm(XH-d))/pH*r*(1-r)+(0.25)*(1-pnorm(XH+a))/pH*(1-r)^2;
            nAAH<- PAAH*(AH+aH)[i]/2
            nAaH<- PAaH*(AH+aH)[i]/2
            naaH<- PaaH*(AH+aH)[i]/2
            return(c( nAAL=nAAL, naaL=naaL,nAaL=nAaL,nAAH=nAAH, naaH=naaH,nAaH=nAaH,a=a,d=d,XL=XL,XH=XH,r=r) );
          }
          h<- estimate_by_full_mode( fAH=H[i],fAL= L[i], pH=np,pL=np)
          n[i,]<-c(h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8],h[9],h[10],h[11])
        }, warning = function(w) {
          n[i,]<-c(0,0,0,0,0,0,0,0,0,0,0)
        })
        
      }
      estimate_genotype[[r]]<-n
    }
    final_es_genotype=Reduce("rbind",estimate_genotype)
    
    allel<-resallel[,5:8]
    G1<-getG1(allel);
    genotype<-as.data.frame(final_es_genotype[,1:6])
    colnames(genotype)<-c("AAL","aaL","AaL","AAH","aaH","AaH")
    es_G2<-apply(genotype,1,getG2)
    es_resG<-rep(NA,times=length(G1) )
    for(i in 1:length(G1)){
      if((G1[i]+es_G2[i])==0){
        es_resG[i]<-0
      }else{
        es_resG[i]<-(G1[i]+0.00000001)*(G1[i]+0.00000001)/((G1[i]+0.00000001)+(es_G2[i]+0.00000001))+(es_G2[i]+0.00000001)*(es_G2[i]+0.00000001)/((G1[i]+0.00000001)+(es_G2[i]+0.00000001))
      }
    }
    result.dQTGseq1<-es_resG
    return(result.dQTGseq1)
  }
  
  dQTGseq2<-function(extreme.result){
    extreme.result<-as.data.frame(extreme.result)
    allel<-extreme.result[,4:7]
    genotype<-extreme.result[,8:13]
    
    numeric.change<-function(x){
      x<-as.numeric(x)
      return(x)
    }
    allel.f<-apply(allel,2,numeric.change)
    genotype.f<-apply(genotype,2,numeric.change)
    getG1 <- function(allel)
    {
      Al=allel[1];al=allel[2];
      Ah=allel[3];ah=allel[4];
      exp <- c(
        (Al + Ah) * (Al + al) / (Al + Ah + al + ah),
        (Al + Ah) * (Ah + ah) / (Al + Ah + al + ah),
        (Al + al) * (al + ah) / (Al + Ah + al + ah),
        (al + ah) * (Ah + ah) / (Al + Ah + al + ah)
      )
      obs <- c(Al, Ah, al, ah)
      
      G <-2 * (rowSums((obs+0.00000001) * log(
        matrix((obs+0.00000001), ncol = 4) / matrix((exp+0.00000001), ncol = 4)
        #matrix((obs+0.00000001), ncol = 4) / matrix((exp), ncol = 4)
      )))
      return(G)
    }
    
    getG2 <- function(genotype)
    {
      # genotype<-genotype.f[1,]
      
      nAl<-genotype[1];nBl<-genotype[2];nHl<-genotype[3]
      nAh<-genotype[4];nBh<-genotype[5];nHh<-genotype[6]
      exp <- c(
        (nAl + nAh) * (nAl+nBl+nHl) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nAl + nAh) * (nAh+nBh+nHh) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nBl + nBh) * (nAl+nBl+nHl) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nBl + nBh) * (nAh+nBh+nHh) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nHl + nHh) * (nAl+nBl+nHl) / (nAl+nBl+nHl+nAh+nBh+nHh),
        (nHl + nHh) * (nAh+nBh+nHh) / (nAl+nBl+nHl+nAh+nBh+nHh)
      )
      obs <- c(nAl,nAh,nBl,nBh,nHl,nHh)
      
      G <-2*(rowSums((obs+0.00000001)*log(
        matrix(obs+0.00000001, ncol = 6) / matrix((exp+0.00000001), ncol = 6)
        #matrix(obs+0.00000001, ncol = 6) / matrix((exp), ncol = 6)
      )))
      return(G)
    }
    
    G1<-apply(allel.f, 1,getG1 );
    G2<-apply(genotype.f, 1,getG2 )
    resG<-(G1+0.00000001)*(G1+0.00000001)/((G1+0.00000001)+(G2+0.00000001))+((G2+0.00000001)*(G2+0.00000001))/((G1+0.00000001)+(G2+0.00000001))
    result.dQTGseq2<-resG
    #}
    return(result.dQTGseq2)
  }
  
  
  
  
  
  BSA<-function(calculatedata){
    result<-calculatedata$gendata
    result<-result[-1,]
    colnames(result)<-c("Maker","Chromosome","Position","AL","aL","AH","aH")
    result1<-result
    extreme.result<-result1
    return(extreme.result)
  }
  
  Extreme.individuals<-function(calculatedata){
    #calculatedata
    #phenotype<-calculatedata$Phenotype
    # if(phenotype==TRUE){
    extreme.result<-list(NULL)
    gen<-calculatedata$gendata
    gen1<-gen[-c(1,dim(gen)[1]),]
    phe<-gen[dim(gen)[1],]
    alldata<-rbind(phe[,-c(1:3)],gen1[,-c(1:3)])
    alldata<-alldata[,order(as.numeric(alldata[1,]))]
    np<-as.numeric(calculatedata$np)/100
    pp<-c(np/4,(np/4)*2,(np/4)*3,(np/4)*4)
    individuals<-as.numeric(calculatedata$nindividual)
    
    for(i in 1:4){
      #i<-1
      lgen<-as.matrix(alldata[-1,1:floor(pp[i]*individuals)])
      hgen<-as.matrix(alldata[-1,(dim(alldata)[2]-floor(pp[i]*individuals)+1):dim(alldata)[2]])
      
      num<-function(x){
        nA<-length(which(x=="2"))
        nB<-length(which(x=="0"))
        nH<-length(which(x=="1"))
        nA<-matrix(nA,1)
        nB<-matrix(nB,1)
        nH<-matrix(nH,1)
        NN<-rbind(nA,nB,nH)
        return(NN)
      }
      
      
      ###A,B,H##########
      nl<-apply(lgen,1,num);nl<-t(nl)
      nh<-apply(hgen,1,num);nh<-t(nh)
      genotype<-cbind(nl,nh)
      
      
      f1<-function(x){
        A<-2*x[1]+x[3]
        a<-2*x[2]+x[3]
        return(c(A,a))
      }
      
      
      ###allele#####
      allel_h<-apply(nh, 1, f1);allel_h<-t(allel_h)
      allel_l<-apply(nl, 1, f1);allel_l<-t(allel_l)
      allel<-cbind(allel_l,allel_h)
      result<-cbind(gen1[,1:3],allel,genotype)
      colnames(result)<-c("Maker","Chromosome","Position","AL","aL","AH","aH","AAL","aaL","AaL","AAH","aaH","AaH")
      #result1<-Filterdata(result)
      extreme.result[[i]]<-result1
    }
    #}
    return(extreme.result)
  }
  
  CIM<-function(calculatedata){
    #calculatedata<-data.calculatedata
    extreme.result<-list(NULL)
    gen<-calculatedata$gendata
    caldata1<-gen[,1:(dim(gen)[2]-1)]
    caldata2<-gen[,c(1,dim(gen)[2])]
    caldata1.1<-caldata1[-4,]
    caldata2.1<-caldata2[-c(1,2,3,4),]
    phe<-caldata2.1[,c(1,2)]
    phe1<-phe[-c(which(phe[,2]==".")),]
    phe2<-phe1[order(as.numeric(phe1[,2])),]
    nindividual<-dim(phe2)[1]
    np<-as.numeric(calculatedata$np)/100
    pp<-c(np/4,(np/4)*2,(np/4)*3,(np/4)*4)
    num<-function(x){
      nA<-length(which(x=="2"))
      nB<-length(which(x=="0"))
      nH<-length(which(x=="1"))
      nA<-matrix(nA,1)
      nB<-matrix(nB,1)
      nH<-matrix(nH,1)
      NN<-rbind(nA,nB,nH)
      return(NN)
    }
    
    
    f1<-function(x){
      A<-2*x[1]+x[3]
      a<-2*x[2]+x[3]
      return(c(A,a))
    }
    
    
    for(i in 1:4){
      sam<-ceiling(dim(phe2)[1]*pp[i])
      l<-phe2[1:sam,1]
      h<-phe2[(dim(phe2)[1]-sam+1):dim(phe2)[1],1]
      headnames<-caldata1[,1]
      lloc<-match(l,headnames);
      hloc<-match(h,headnames);
      lgen<-caldata1.1[lloc,-1];
      hgen<-caldata1.1[hloc,-1];
      nl<-apply(lgen,2,num);
      nh<-apply(hgen,2,num);
      genotype<-rbind(nl,nh)
      allel_h<-apply(nh,2, f1);
      allel_l<-apply(nl,2, f1);
      allel<-rbind(allel_l,allel_h)
      caldata3<-t(rbind(caldata1.1[1:3,-1],allel,genotype))
      result<-cbind(caldata3[,1],caldata3[,3],caldata3[,2],caldata3[,4:13])
      colnames(result)<-c("Maker","Chromosome","Position","AL","aL","AH","aH","AAL","aaL","AaL","AAH","aaH","AaH")
      extreme.result[[i]]<- result
    }
    extreme.result[[5]]<-dim(phe2)[1]
    # Filterdata<-function(x){
    #   #x<-result
    #   x<-na.omit(x)
    #   AL<-x$AL;aL<-x$aL;AH<-x$AH;aH<-x$aH
    #   x1<-x[-c(which(as.numeric(x$AL)==0) & which(as.numeric(x$aL)==0)),]
    #   x2<-x1[-c(which( as.numeric(x1$AH)==0) & which(as.numeric(x1$aH)==0)),]
    #
    #   return(x2)
    #
    # }
    # extreme.result<-Filterdata(extreme.result)
    return(extreme.result)
  }
  
  GCIM<-function(calculatedata){
    extreme.result<-list(NULL)
    gen<-calculatedata$gendata
    phe<-gen[c(1,dim(gen)[1]),-c(1,2,3)]
    gen<-gen[-dim(gen)[1],]
    phe1<-phe[,order(as.numeric(phe[2,]))]
    phe2<-phe1[,colSums(is.na(phe1))==0]
    np<-as.numeric(calculatedata$np)/100
    pp<-c(np/4,(np/4)*2,(np/4)*3,(np/4)*4)
    num<-function(x){
      nA<-length(which(x=="2"))
      nB<-length(which(x=="0"))
      nH<-length(which(x=="1"))
      nA<-matrix(nA,1)
      nB<-matrix(nB,1)
      nH<-matrix(nH,1)
      NN<-rbind(nA,nB,nH)
      return(NN)
    }
    
    f1<-function(x){
      A<-2*x[1]+x[3]
      a<-2*x[2]+x[3]
      return(c(A,a))
    }
    
    for(i in 1:4){
      #i<-4
      sam<-ceiling(dim(phe2)[2]*pp[i])
      l<-phe2[1,1:sam]
      h<-phe2[1,(dim(phe2)[2]-sam+1):dim(phe2)[2]]
      headnames<-gen[1,]
      lloc<-match(l,headnames);
      hloc<-match(h,headnames);
      lgen<-gen[-1,lloc];
      hgen<-gen[-1,hloc];
      nl<-apply(lgen,1,num);nl<-t(nl)
      nh<-apply(hgen,1,num);nh<-t(nh)
      genotype<-cbind(nl,nh)
      allel_h<-apply(nh, 1, f1); allel_h<-t(allel_h)
      allel_l<-apply(nl, 1, f1); allel_l<-t(allel_l)
      allel<-cbind(allel_l,allel_h)
      result<-cbind(gen[-1,1:3],allel,genotype)
      colnames(result)<-c("Maker","Chromosome","Position","AL","aL","AH","aH","AAL","aaL","AaL","AAH","aaH","AaH")
      extreme.result[[i]]<- result
    }
    extreme.result[[5]]<-dim(phe2)[2]
    # Filterdata<-function(x){
    #   #x<-result
    #   x<-na.omit(x)
    #   AL<-x$AL;aL<-x$aL;AH<-x$AH;aH<-x$aH
    #   x1<-x[-c(which(as.numeric(x$AL)==0) & which(as.numeric(x$aL)==0)),]
    #   x2<-x1[-c(which( as.numeric(x1$AH)==0) & which(as.numeric(x1$aH)==0)),]
    #
    #   return(x2)
    #
    # }
    #extreme.result<-Filterdata(extreme.result)
    return(extreme.result)
  }
  
  ICIM<-function(calculatedata){
    
    extreme.result<-list(NULL)
    phe<-calculatedata$phe
    gen<-calculatedata$gen
    pos<-calculatedata$pos
    
    dataall1<-rbind(phe[,-1],gen[,-1])
    dataall2<-dataall1[,-c(which(dataall1[1,]=="NA"))]
    dataall3<-dataall2[,order(as.numeric(dataall2[1,]))]
    np<-as.numeric(calculatedata$np)/100
    pp<-c(np/4,(np/4)*2,(np/4)*3,(np/4)*4)
    
    num<-function(x){
      nA<-length(which(x=="2"))
      nB<-length(which(x=="0"))
      nH<-length(which(x=="1"))
      nA<-matrix(nA,1)
      nB<-matrix(nB,1)
      nH<-matrix(nH,1)
      NN<-rbind(nA,nB,nH)
      return(NN)
    }
    
    f1<-function(x){
      A<-2*x[1]+x[3]
      a<-2*x[2]+x[3]
      return(c(A,a))
    }
    for(i in 1:4){
      #i<-1
      sam<-ceiling(dim(dataall3)[2]*pp[i])
      lgen<- dataall3[-1,1:sam];
      hgen<- dataall3[-1,(dim(dataall3)[2]-sam+1):dim(dataall3)[2]];
      nl<-apply(lgen,1,num);
      nh<-apply(hgen,1,num);
      genotype<-cbind(t(nl),t(nh))
      allel_h<-t(apply(nh, 2, f1))
      allel_l<-t(apply(nl, 2, f1))
      allel<-cbind(allel_l,allel_h)
      result<- cbind(pos[,1:3],allel,genotype)
      colnames(result)<-c("Maker","Chromosome","Position","AL","aL","AH","aH","AAL","aaL","AaL","AAH","aaH","AaH")
      extreme.result[[i]]<- result
    }
    extreme.result[[5]]<-dim(dataall3)[2]
    # extreme.result<-Filterdata(extreme.result)
    return(extreme.result)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if(Fileformat=="BSA" & population=="F2"){
    extreme.result<-BSA(calculatedata)
    result.sta<-cbind(extreme.result[,1:3],dQTGseq1(extreme.result,np,windowsize),SmoothLOD.function(extreme.result),G.function(extreme.result),
                      SNP.function(extreme.result),ED.function(extreme.result)
    )
    colnames(result.sta)<-c("Maker","Chromosome","Position","dQTGseq1","SmoothLOD","G'","deltaSNP","ED")
    if(Smooth_method=="None"){
      result.smooth<-None(result.sta)
    }
    if(Smooth_method=="AIC"){
      result.smooth<-AIC.f(result.sta)
    }
    if(Smooth_method == "Default"){
      smoothdQTGseq1<-Windowsize.f(result.sta[,1:6],windowsize)
      smoothSNP<-Block.f(cbind(result.sta[,1:3],result.sta[,7]),20)
      smoothED<-AIC.f(cbind(result.sta[,1:3],result.sta[,8]))
      result.smooth<-cbind(smoothdQTGseq1,smoothSNP[,-c(1:3)],smoothED[,-c(1:3)])
    }
    if(str_detect(c(Smooth_method), "Window size")==TRUE){
      result.smooth<-Windowsize.f(result.sta,windowsize)
    }
    
    if(str_detect(c(Smooth_method), "Block")==TRUE){
      if(length(strsplit(Smooth_method,split = " ")[[1]])==2){
        blocknumber<-as.numeric(strsplit(Smooth_method,split = " ")[[1]][2])
      }else{
        blocknumber<-20
      }
      result.smooth<-Block.f(result.sta,blocknumber)
    }
    
    
    map.data <- sim.map(len = c(999), n.mar = nmarker,
                        anchor.tel = TRUE, include.x = FALSE, sex.sp = FALSE, eq.spacing = TRUE)
    f.result1<-NULL;f.result2<-NULL;f.result3<-NULL;f.result4<-NULL;f.result5<-NULL;
    for(i in 1:repetition){
      result1<-dQTGseq2.p(list(dodata4(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np)))
      f.result1<-rbind(f.result1,result1)
      result2<-G.p(as.data.frame(dodata1(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np)))
      f.result2<-rbind(f.result2,result2)
      result3<-ED.p(as.data.frame(dodata1(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np)))
      f.result3<-rbind(f.result3,result3)
      result4<-SNP.p(as.data.frame(dodata1(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np)))
      f.result4<-rbind(f.result4,result4)
      result5<-SmoothLOD.p(as.data.frame(dodata1(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np)))
      f.result5<-rbind(f.result5,result5)
    }
    p.1<-f.result1[order(f.result1)][0.95*repetition]
    p.2<-f.result2[order(f.result2)][0.95*repetition]
    p.3<-f.result3[order(f.result3)][0.95*repetition]
    p.4<-f.result4[order(f.result4)][0.95*repetition]
    p.5<-f.result5[order(f.result5)][0.95*repetition]
    threshold <- c(p.1,p.5,p.2,p.4,p.3)
    
    #threshold <- c(100,100,100,100,100)
    
    allresult<-cbind(result.sta[,1:3],result.sta[,4],result.smooth[,4],result.sta[,5],result.smooth[,5],result.sta[,6],result.smooth[,6],result.sta[,7],result.smooth[,7],result.sta[,8],result.smooth[,8])
    names(allresult)<-c("Marker","Chromosome","Position","Gw","Smooth_Gw","LOD","Smooth_LOD","G","G'","deltaSNP","Smooth_deltaSNP","ED","Smooth_ED")
    #write_xlsx(allresult,path = paste(dir,"Allresult.xlsx",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    
    write.csv(allresult,paste(dir,"/","Allresult.csv",sep=""))
    
    #write.csv(allresult,path = paste(dir,"Allresult.csv",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    
    
    
    result.significance<-Significance(allresult,threshold)
    names(result.significance)<-c("dQTGseq1","Smooth_LOD","G","deltaSNP","ED")
    tryCatch({
      write_xlsx(result.significance,path = paste(dir,"/","Significanceresult.xlsx",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    },error=function(w){
      cc<-which(is.na(result.significance))
      vec<-rep(1:length(result.significance))
      for(jj in 1:length(cc)){
        #jj<-1
        cat(paste("No significant markers are detected in the",names(cc[jj]),"method!  "))
      }
      result.significance1<-result.significance[vec[-cc]]
      write_xlsx(result.significance1,path = paste(dir,"/","Significanceresult.xlsx",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    })
    
    if(DrawPlot==TRUE){
      if(Resolution=="Low"){
        manwidth<-960;manhei<-960;manwordre<-20;manfigurere<-72
      }else if(Resolution=="High"){
        #width=6500;height=4500;res=1000;compression = c("zip")
        manwidth=20000; manhei=20000;units= "px";manwordre =30;manfigurere=1000;compression = c("zip")
      }
      if(Plotformat=="png"){
        png(paste(dir,"/","LOD.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat=="tiff"){
        tiff(paste(dir,"/","LOD.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat=="jpeg"){
        jpeg(paste(dir,"/","LOD.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat=="pdf"){
        pdf(paste(dir,"/","LOD.pdf",sep=""),width=15,fonts = "sans")
      }
      if("all"%in%chrom){
        draw1(allresult,Fileformat,population,color,Smooth_method,threshold)
      }
      
      
      if("all"%in%chrom == FALSE){
        data1<-allresult
        data3<-NULL
        for (i in 1:length(unique(chrom))){
          data3[[i]]<-data1[which(as.numeric(data1[,2])==chrom[i]),]
        }
        data3<-do.call(rbind,data3)
        draw1(data3,Fileformat,population,color,Smooth_method,threshold)
      }
      dev.off()
    }
    
  }
  
  
  #######################################################################
  if(Fileformat=="BSA" & (population=="BC" | population=="DH" | population=="RIL")){
    
    extreme.result<-BSA(calculatedata)
    result.sta<-cbind(extreme.result[,1:3],SmoothLOD.function(extreme.result),G.function(extreme.result),
                      SNP.function(extreme.result),ED.function(extreme.result)
    )
    colnames(result.sta)<-c("Maker","Chromosome","Position","SmoothLOD","G'","deltaSNP","ED")
    
    
    if(Smooth_method=="None"){
      result.smooth<-None(result.sta)
    }
    
    if(Smooth_method=="AIC"){
      result.smooth<-AIC.f(result.sta)
    }
    
    if(Smooth_method == "Default"){
      smoothdQTGseq1<-Windowsize.f(result.sta[,1:5],windowsize)
      smoothSNP<-Block.f(cbind(result.sta[,1:3],result.sta[,6]),20)
      smoothED<-AIC.f(cbind(result.sta[,1:3],result.sta[,7]))
      result.smooth<-cbind(smoothdQTGseq1,smoothSNP[,-c(1:3)],smoothED[,-c(1:3)])
    }
    
    if(str_detect(c(Smooth_method), "Window size")==TRUE){
      result.smooth<-Windowsize.f(result.sta,windowsize)
    }
    
    if(str_detect(c(Smooth_method), "Block")==TRUE){
      if(length(strsplit(Smooth_method,split = " ")[[1]])==2){
        blocknumber<-as.numeric(strsplit(Smooth_method,split = " ")[[1]][2])
      }else{
        blocknumber<-20
      }
      result.smooth<-Block.f(result.sta,blocknumber)
    }
    map.data <- sim.map(len = c(999), n.mar = nmarker,
                        anchor.tel = TRUE, include.x = FALSE, sex.sp = FALSE, eq.spacing = TRUE)
    if(population=="BC"){
      f.result1<-NULL;f.result2<-NULL;f.result3<-NULL;f.result4<-NULL;f.result5<-NULL;
      for(i in 1:repetition){
        result2<-G.p(as.data.frame(dodata2(BC.data(map.data,"bc",nsam,vg,qq,b0,v1),nsam,np)))
        f.result2<-rbind(f.result2,result2)
        result3<-ED.p(as.data.frame(dodata2(BC.data(map.data,"bc",nsam,vg,qq,b0,v1),nsam,np)))
        f.result3<-rbind(f.result3,result3)
        result4<-SNP.p(as.data.frame(dodata2(BC.data(map.data,"bc",nsam,vg,qq,b0,v1),nsam,np)))
        f.result4<-rbind(f.result4,result4)
        result5<-SmoothLOD.p(as.data.frame(dodata2(BC.data(map.data,"bc",nsam,vg,qq,b0,v1),nsam,np)))
        f.result5<-rbind(f.result5,result5)
      }
      p.2<-f.result2[order(f.result2),][0.95*repetition]
      p.3<-f.result3[order(f.result3),][0.95*repetition]
      p.4<-f.result4[order(f.result4),][0.95*repetition]
      p.5<-f.result5[order(f.result5),][0.95*repetition]
      threshold<- c(p.5,p.2,p.4,p.3)
    }
    
    if(population=="DH"){
      f.result1<-NULL;f.result2<-NULL;f.result3<-NULL;f.result4<-NULL;f.result5<-NULL;
      for(i in 1:repetition){
        result2<-G.p(as.data.frame(dodata3(DH.data(map.data,"bc",nsam,vg,qq,b0,v1),nsam,np)))
        f.result2<-rbind(f.result2,result2)
        result3<-ED.p(as.data.frame(dodata3(DH.data(map.data,"bc",nsam,vg,qq,b0,v1),nsam,np)))
        f.result3<-rbind(f.result3,result3)
        result4<-SNP.p(as.data.frame(dodata3(DH.data(map.data,"bc",nsam,vg,qq,b0,v1),nsam,np)))
        f.result4<-rbind(f.result4,result4)
        result5<-SmoothLOD.p(as.data.frame(dodata3(DH.data(map.data,"bc",nsam,vg,qq,b0,v1),nsam,np)))
        f.result5<-rbind(f.result5,result5)
      }
      p.2<-f.result2[order(f.result2),][0.95*repetition]
      p.3<-f.result3[order(f.result3),][0.95*repetition]
      p.4<-f.result4[order(f.result4),][0.95*repetition]
      p.5<-f.result5[order(f.result5),][0.95*repetition]
      threshold<- c(p.5,p.2,p.4,p.3)
    }
    
    
    if(population=="RIL"){
      f.result1<-NULL;f.result2<-NULL;f.result3<-NULL;f.result4<-NULL;f.result5<-NULL;
      for(i in 1:repetition){
        result2<-G.p(as.data.frame(dodata3(RIL.data(map.data,"risib",nsam,vg,qq,b0,v1),nsam,np)))
        f.result2<-rbind(f.result2,result2)
        result3<-ED.p(as.data.frame(dodata3(RIL.data(map.data,"risib",nsam,vg,qq,b0,v1),nsam,np)))
        f.result3<-rbind(f.result3,result3)
        result4<-SNP.p(as.data.frame(dodata3(RIL.data(map.data,"risib",nsam,vg,qq,b0,v1),nsam,np)))
        f.result4<-rbind(f.result4,result4)
        result5<-SmoothLOD.p(as.data.frame(dodata3(RIL.data(map.data,"risib",nsam,vg,qq,b0,v1),nsam,np)))
        f.result5<-rbind(f.result5,result5)
      }
      p.2<-f.result2[order(f.result2),][0.95*repetition]
      p.3<-f.result3[order(f.result3),][0.95*repetition]
      p.4<-f.result4[order(f.result4),][0.95*repetition]
      p.5<-f.result5[order(f.result5),][0.95*repetition]
      #p.5<-1
      threshold<- c(p.5,p.2,p.4,p.3)
    }
    
    allresult<-cbind(result.sta[,1:3],result.sta[,4],result.smooth[,4],result.sta[,5],result.smooth[,5],result.sta[,6],result.smooth[,6],result.sta[,7],result.smooth[,7])
    names(allresult)<-c("Maker","Chromosome","Position","LOD","Smooth_LOD","G","G'","deltaSNP","Smooth_deltaSNP","ED","Smooth_ED")
    #write_xlsx(allresult,path = paste(dir,"allresult.xlsx",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    #write.csv(allresult,path = paste(dir,"Allresult.csv",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    write.csv(allresult,paste(dir,"/","Allresult.csv",sep=""))
    
    
    result.significance<-Significance(allresult,threshold)
    names(result.significance)<-c("Smooth_LOD","G","deltaSNP","ED")
    tryCatch({
      write_xlsx(result.significance,path = paste(dir,"/","Significanceresult.xlsx",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    },error=function(w){
      cc<-which(is.na(result.significance))
      vec<-rep(1:length(result.significance))
      for(jj in 1:length(cc)){
        cat(paste("No significant markers are detected in the",names(cc[jj]),"method!  "))
      }
      result.significance1<-result.significance[vec[-cc]]
      write_xlsx(result.significance1,path = paste(dir,"/","Significanceresult.xlsx",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    })
    
    if(DrawPlot==TRUE){
      if(Resolution=="Low"){
        manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
      }else if(Resolution=="High"){
        manwidth=20000; manhei=20000;units= "px";manwordre =30;manfigurere=1000;compression = c("zip")
      }
      if(Plotformat=="png"){
        png(paste(dir,"/","LOD.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat=="tiff"){
        tiff(paste(dir,"/","LOD.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat=="jpeg"){
        jpeg(paste(dir,"/","LOD.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat=="pdf"){
        
        pdf(paste(dir,"/","LOD.pdf",sep=""),width=15,fonts = "sans")
      }
      if("all"%in%chrom){
        draw2(allresult,population,color,Smooth_method,threshold)
      }
      if("all"%in%chrom == FALSE){
        data1<-allresult
        data3<-NULL
        for (i in 1:length(unique(chrom))){
          data3[[i]]<-data1[which(as.numeric(data1[,2])==chrom[i]),]
        }
        data3<-do.call(rbind,data3)
        draw2(data3,population,color,Smooth_method,threshold)
      }
      dev.off()
    }
    
  }
  
  
  ##################################################################################
  if((Fileformat=="Extreme individuals")|(Fileformat=="CIM")|(Fileformat=="ICIM")|(Fileformat=="GCIM")){
    
    if(Fileformat =="Extreme individuals"){
      extreme.result<-Extreme.individuals(calculatedata)
    }
    
    if(Fileformat =="CIM"){
      extreme.result<-CIM(calculatedata)
      nsam<-as.numeric(extreme.result[[5]])
      extreme.result[[5]]<-NULL
    }
    if(Fileformat =="ICIM"){
      extreme.result<-ICIM(calculatedata)
      nsam<-as.numeric(extreme.result[[5]])
      extreme.result[[5]]<-NULL
    }
    if(Fileformat =="GCIM"){
      extreme.result<-GCIM(calculatedata)
      nsam<-as.numeric(extreme.result[[5]])
      extreme.result[[5]]<-NULL
    }
    
    result.sta <- list(NULL);result.significance<-list(NULL)
    map.data <- sim.map(len = c(999), n.mar =nmarker,anchor.tel = TRUE,include.x = FALSE, sex.sp = FALSE,eq.spacing = TRUE)
    f.result1<-NULL;f.result2<-NULL;f.result3<-NULL;f.result4<-NULL;f.result5<-NULL;
    for(i in 1:repetition){
      #for(i in 1:1000){
      result1<-dQTGseq2.p(dodata4(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np))
      f.result1<-rbind(f.result1,result1)
      result2<-G.p(as.data.frame(dodata1(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np)))
      f.result2<-rbind(f.result2,result2)
      result3<-ED.p(as.data.frame(dodata1(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np)))
      f.result3<-rbind(f.result3,result3)
      result4<-SNP.p(as.data.frame(dodata1(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np)))
      f.result4<-rbind(f.result4,result4)
      result5<-SmoothLOD.p(as.data.frame(dodata1(F2.data(map.data,"f2",nsam,vg,qq,b0,v1),nsam,np)))
      f.result5<-rbind(f.result5,result5)
    }
    p.1<-f.result1[order(f.result1)][0.95*repetition]
    #p.1<-f.result1[order(f.result1)][0.95*4000]
    p.2<-f.result2[order(f.result2)][0.95*repetition]
    p.3<-f.result3[order(f.result3)][0.95*repetition]
    p.4<-f.result4[order(f.result4)][0.95*repetition]
    p.5<-f.result5[order(f.result5)][0.95*repetition]
    threshold <- c(p.1,p.5,p.2,p.4,p.3)
    
    allresult<-NULL;allresult.smooth<-NULL;result.smooth<-NULL;allresult.ori<-NULL
    for(j in 1:length(extreme.result)){
      allel<-extreme.result[[j]][,4:7]
      numeric.change<-function(x){
        x<-as.numeric(x)
        return(x)
      }
      allel.f<-apply(allel,2,numeric.change)
      colnames(allel.f)<-c("AL","aL","AH","aH")
      allel.f<-as.data.frame(allel.f)
      round.f<-function(x){
        x<-round(x,5)
        return(x)
      }
      result.sta<-cbind(extreme.result[[1]][,1:3],round.f(dQTGseq2(extreme.result[[j]])),
                        round.f(SmoothLOD.function(allel.f)),round.f(G.function(allel.f)),
                        round.f(SNP.function(allel.f)),round.f(ED.function(allel.f)))
      
      result.sta<-as.data.frame(result.sta)
      colnames(result.sta)<-c("Maker","Chromosome","Position","dQTGseq2","SmoothLOD","G'","deltaSNP","ED")
      
      if(Smooth_method=="None"){
        result.smooth<-None(result.sta)
      }
      
      
      if(Smooth_method == "Default"){
        smoothdQTGseq1<-Windowsize.f(result.sta[,1:6],windowsize)
        smoothSNP<-Block.f(cbind(result.sta[,1:3],result.sta[,7]),20)
        smoothED<-AIC.f(cbind(result.sta[,1:3],result.sta[,8]))
        result.smooth<-cbind(smoothdQTGseq1,smoothSNP[,-c(1:3)],smoothED[,-c(1:3)])
      }
      
      if(Smooth_method=="AIC"){
        result.smooth<-AIC.f(result.sta)
      }
      if(str_detect(c(Smooth_method), "Window size")==TRUE){
        result.smooth<-Windowsize.f(result.sta,windowsize)
      }
      
      if(str_detect(c(Smooth_method), "Block")==TRUE){
        if(length(strsplit(Smooth_method,split = " ")[[1]])==2){
          blocknumber<-as.numeric(strsplit(Smooth_method,split = " ")[[1]][2])
        }else{
          blocknumber<-20
        }
        result.smooth<-Block.f(result.sta,blocknumber)
      }
      
      
      allresult.smooth[[j]]<-result.smooth
      allresult.ori[[j]]<-result.sta
      result.all<-cbind(result.sta[,1:3],result.sta[,4],result.smooth[,4],result.sta[,5],result.smooth[,5],result.sta[,6],result.smooth[,6],result.sta[,7],result.smooth[,7],result.sta[,8],result.smooth[,8])
      names(result.all)<-c("Maker","Chromosome","Position","Gw","Smooth_Gw","LOD","Smooth_LOD","G","G'","deltaSNP","Smooth_deltaSNP","ED","Smooth_ED")
      allresult[[j]]<-result.all;
      newthreshold<-c(0,0,0,0,0)
      result.significance[[j]]<-Significance(result.all,newthreshold)
    }
    
    names(allresult)<-c("0.05","0.10","0.15","0.20")
    
    #write_xlsx(allresult,path = paste(dir,"Allresult.xlsx",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    
    # write.csv(allresult[[1]],path = paste(dir,"0.05result.csv",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    # write.csv(allresult[[2]],path = paste(dir,"0.10result.csv",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    # write.csv(allresult[[3]],path = paste(dir,"0.15result.csv",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    # write.csv(allresult[[4]],path = paste(dir,"0.20result.csv",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    
    
    write.csv(allresult[[1]],paste(dir,"/","0.05result.csv",sep=""))
    write.csv(allresult[[2]],paste(dir,"/","0.10result.csv",sep=""))
    write.csv(allresult[[3]],paste(dir,"/","0.15result.csv",sep=""))
    write.csv(allresult[[4]],paste(dir,"/","0.20result.csv",sep=""))
    
    
    
    
    
    ds<-1000000;
    ff.1<-function(ff,ds){
      ff.sig<-list(NULL)
      for(n in 1:length(ff)){
        caldata<-ff[[n]]
        if(dim(caldata)[2]==1){
          ff.sig[[n]]<-NA
        }else{
          cc<-numeric()
          chr<-na.omit(unique(as.numeric(caldata[,2])))
          for(rr in 1:length(chr)){
            caldata1<-caldata[as.numeric(caldata[,2])==chr[rr],1:5]
            caldata1<-caldata1[order(as.numeric(caldata1$Position)),]
            caldata1<-na.omit(caldata1)
            qq<-numeric()
            aa<-caldata1
            
            
            #threshold[1]<-7.9
            for(j in 1:length(unique(aa[,5]))){
              repaa<-aa[which(aa[,5]==unique(aa[,5])[j]),]
              if (dim(repaa)[1]==1) {
                aa1<-repaa
              }else{
                ord<-which(aa[,5]==unique(aa[,5])[j])
                index<-diff(ord)
                index2<-numeric()
                if(length(unique(index)) == 1){
                  index2<-repaa[ceiling(dim(repaa)[1]/2),]
                }else{
                  for(k in 1:((length(which(index!=1)))+1)){
                    index5<-length(which(index!=1))+1
                    if(k!=index5){
                      index3<-aa[ord[1:which(index!=1)[k]],]
                      index4<-index3[ceiling(dim(index3)[1]/2),]
                      index2<-rbind(index2,index4)
                    }else{
                      index3<-aa[ord[(which(index!=1)[k-1]+1):(length(index)+1)],]
                      index4<-index3[ceiling(dim(index3)[1]/2),]
                      index2<-rbind(index2,index4)
                    }
                  }
                }
                aa1<-index2
              }
              qq<-rbind(qq,aa1)
            }
            
            
            qq5<-qq[order(as.numeric(qq$Position)),]
            qq5<-qq5[which(abs(qq5[,5])>=threshold[n]),]
            qq5[,3]<-as.numeric(qq5[,3])
            qq6<-qq5
            if((dim(qq5)[1]==1) |(dim(qq5)[1]==0) ){
              qq6<-qq5
            }else{
              for (i in 1:(dim(qq5)[1]-1)){
                distance<-qq5[i+1,3]-qq5[i,3]
                if(distance<=ds){
                  qq6[which(qq5[,5]==min(qq5[i+1,5],qq5[i,5])),4:5]<-NA
                }else{
                  qq6[i,]<-qq6[i,]
                }
              }
              if((qq5[dim(qq5)[1],3]-qq5[dim(qq5)[1]-1,3])<=ds){
                qq6[which(qq5[,5]==min(qq5[dim(qq5)[1],5],qq5[dim(qq5)[1]-1,5])),4:5]<-NA
              }else{
                qq6[dim(qq5)[1],]<-qq5[dim(qq5)[1],]
              }
            }
            qq6<-na.omit(qq6)
            cc<-rbind(cc,qq6)
          }
          
          if(dim(cc)[1]==0){
            cc<-NA
          }else{
            Threshold<-rep(NA,dim(cc)[1])
            Threshold[1]<-threshold[n]
            cc<-cbind(cc,Threshold)
          }
          ff.sig[[n]]<-cc
        }
      }
      return(ff.sig)
    }
    
    pp1<-ff.1(result.significance[[1]],ds)
    pp2<-ff.1(result.significance[[2]],ds)
    pp3<-ff.1(result.significance[[3]],ds)
    pp4<-ff.1(result.significance[[4]],ds)
    
    fff<-list(NULL)
    for(a in 1:5){
      fff[[a]]<-rbind(pp1[[a]],pp2[[a]],pp3[[a]],pp4[[a]])
    }
    
    ff.sig1<-list(NULL); ff.sig2<-list(NULL)
    for(ii in 1:length(fff)){
      if(dim(fff[[ii]])[2]==1){
        ff.sig2[[ii]]<-NA
      }else{
        data1<-fff[[ii]][order(as.numeric(fff[[ii]]$Chromosome),as.numeric(fff[[ii]]$Position)),]
        data2<-data1[,1:5]
        sigaa<-numeric()
        ff.ds<-function(xx,ds){
          #xx<-data2
          data<-xx
          chrome<-as.numeric(data[,2]);
          cc<-unique(chrome);
          for(ll in 1:length(cc)){
            qq5<-data[which(chrome==cc[ll]),]
            qq6<-qq5
            if((dim(qq5)[1]==1) |(dim(qq5)[1]==0) ){
              qq6<-qq5
            }else{
              for (i in 1:(dim(qq5)[1]-1)){
                distance<-qq5[i+1,3]-qq5[i,3]
                if((distance<=ds) | (distance =0)){
                  qq6[which(qq5[,5]==min(qq5[i+1,5],qq5[i,5])),4:5]<-NA
                }else{
                  qq6[i,]<-qq6[i,]
                }
                
              }
              if((qq5[dim(qq5)[1],3]-qq5[dim(qq5)[1]-1,3])<=ds){
                qq6[which(qq5[,5]==min(qq5[dim(qq5)[1],5],qq5[dim(qq5)[1]-1,5])),4:5]<-NA
              }else{
                qq6[dim(qq5)[1],]<-qq5[dim(qq5)[1],]
              }
            }
            sigaa<-rbind(sigaa,qq6)
          }
          sigaa1<-na.omit(sigaa)
          return(sigaa1)
        }
        
        ff.sig1[[ii]]<-ff.ds(data2,ds)
        ff.sig2[[ii]]<-ff.ds(ff.sig1[[ii]],ds)
        
        if(dim(ff.sig2[[ii]])[1]==0){
          ff.sig2[[ii]]<-NA
        }else{
          Threshold<-rep(NA,dim(ff.sig2[[ii]])[1])
          Threshold[1]<-threshold[ii]
          ff.sig2[[ii]]<-cbind(ff.sig2[[ii]],Threshold)
        }
      }
      
    }
    
    
    
    
    
    names(ff.sig2)<-c("dQTGseq2","Smooth_LOD","G","deltaSNP","ED")
    tryCatch({
      write_xlsx(ff.sig2,path = paste(dir,"/","Significanceresult.xlsx",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    },error=function(w){
      vec.na<-which(is.na(ff.sig2))
      vec<-rep(1:length(ff.sig2))
      for(jj in 1:length(vec.na)){
        cat(paste("No significant markers are detected in the",names(vec.na[jj]),"method!  "))
      }
      result.significance1<-ff.sig2[vec[-vec.na]]
      write_xlsx(result.significance1,path = paste(dir,"/","Significanceresult.xlsx",sep=""),col_names = TRUE,format_headers = TRUE,use_zip64 = FALSE)
    })
    
    
    
    drawresult.smooth1<-matrix(NA,ncol = 5,nrow = dim(allresult.smooth[[1]])[1])
    drawresult.ori1<-matrix(NA,ncol = 5,nrow = dim(allresult.ori[[1]])[1])
    
    for(jj in 1:5){
      #jj<-1
      result1<-cbind(allresult.smooth[[1]][,jj+3],allresult.smooth[[2]][,jj+3],allresult.smooth[[3]][,jj+3],allresult.smooth[[4]][,jj+3])
      result2<-cbind(as.numeric(allresult.ori[[1]][,jj+3]),as.numeric(allresult.ori[[2]][,jj+3]),as.numeric(allresult.ori[[3]][,jj+3]),as.numeric(allresult.ori[[4]][,jj+3]))
      
      ff.num<-function(x){
        y<-max(as.numeric(x))
        return(y)
      }
      drawresult.smooth1[,jj]<-matrix(apply(result1,1,ff.num),ncol = 1)
      drawresult.ori1[,jj]<-matrix(apply(result2,1,ff.num),ncol = 1)
    }
    
    drawresult.smooth<-cbind(allresult.smooth[[1]][,1:3],drawresult.smooth1)
    drawresult.ori<-cbind(allresult.ori[[1]][,1:3],drawresult.ori1)
    all<-cbind(drawresult.ori[,1:3],drawresult.ori[,4],drawresult.smooth[,4],drawresult.ori[,5],drawresult.smooth[,5],drawresult.ori[,6],drawresult.smooth[,6],
               drawresult.ori[,7],drawresult.smooth[,7],drawresult.ori[,8],drawresult.smooth[,8])
    names(all)<-c("Marker","Chromosome","Position","Gw","Smooth_Gw","LOD","Smooth_LOD","G","Smooth_G","deltaSNP","Smooth_deltaSNP","ED","Smooth_ED")
    all1<-na.omit(all)
    #write.csv(all1,"all1.csv")
    
    
    if(DrawPlot==TRUE){
      if(Resolution=="Low"){
        manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
      }else if(Resolution=="High"){
        manwidth=20000; manhei=20000;units= "px";manwordre =30;manfigurere=1000;compression = c("zip")
      }
      if(Plotformat=="png"){
        png(paste(dir,"/","LOD.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat=="tiff"){
        tiff(paste(dir,"/","LOD.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat=="jpeg"){
        jpeg(paste(dir,"/","LOD.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat=="pdf"){
        
        pdf(paste(dir,"/","LOD.pdf",sep=""),width=15,fonts = "sans")
      }
      if("all"%in%chrom){
        draw1(all1,Fileformat,population,color,Smooth_method,threshold)
        
      }
      if("all"%in%chrom == FALSE){
        data1<-all1
        data3<-NULL
        for (i in 1:length(unique(chrom))){
          data3[[i]]<-data1[which(as.numeric(data1[,2])==chrom[i]),]
        }
        data3<-do.call(rbind,data3)
        draw1(data3,Fileformat,population,color,Smooth_method,threshold)
      }
      dev.off()
    }
    
    
    
  }
  
}


