####readdata#####
#' Title readdata function
#' @param File  the input file
#'
#' @return list
#' @export
#' @examples
#' data(BSA)
#' Readdata1(BSA)
Readdata1<-function(File){
  if(is.character(File)==TRUE){
    data1<-fread(File,header = FALSE,stringsAsFactors=FALSE)
    #titlenameGen<-colnames(genRaw)[1:3]
  }else{
    data1<-File
  }



  result<-as.data.frame(data1)
  Species.windowsize<-function(Species,Smooth_method){
    #Species<-"Maize";Smooth_method<-"None"
    spe<- c("Arabidopsis","Cucumber","Maize","Brassica juncea","Brassica napus","Rice","Tobacco","Tomato","Wheat",
            "Yeast","Glycine soja","Glycine max","Gossypium hirsutum L","Gossypium barbadense","Brassica pekinensis"
    )
    index<-c(0.2083,0.2529,1.1208,0.3797,0.3380,0.1373,0.8257,0.8180,6.1151,0.0024,0.3580,0.4358,0.5672,0.7774,0.3134)
    cc<-cbind(spe,index)
    #Species="Arabidopsis"
    if(Smooth_method=="None"|Smooth_method=="AIC"|Smooth_method == "Block"|Smooth_method=="Default"){
      windowsize<-1000000
    }else{
      #Smooth_method<-"windown size 1"
      if(length(strsplit(Smooth_method,split = " ")[[1]])==3){
        windowsize<-as.numeric(strsplit(Smooth_method,split = " ")[[1]][3])*1000000
      }else{
        if (Species %in% spe){
          windowsize<-as.numeric(cc[which(cc[,1]==Species),2])*1000000
        }else{
          windowsize<-1000000
        }
      }
    }
    return(windowsize)

  }
  if(result[3,1]=="Sample-size"){
    gendata<-result[-c(1:10),]
    calculatedata<-list(Species.windowsize(result[1,2],result[6,2]),result[2,2],result[3,2],result[4,2],result[5,2],result[6,2],result[7,2],
                        result[8,2],result[9,2],result[10,2],gendata)
    names(calculatedata)<-c("windowsize","FileFormat","nindividual","Population","np","Smooth_method","Repetition",
                            "DrawPlot","Resolution","Plotformat","gendata")
  }else{
    gendata<-result[-c(1:9),]
    calculatedata<-list(Species.windowsize(result[1,2],result[5,2]),result[2,2],result[3,2],result[4,2],result[5,2],result[6,2],result[7,2],
                        result[8,2],result[9,2],gendata)
    names(calculatedata)<-c("windowsize","FileFormat","Population","np","Smooth_method","Repetition",
                            "DrawPlot","Resolution","Plotformat","gendata")
  }
  return(calculatedata)
}




