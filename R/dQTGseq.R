# library("stringr")
# library("doParallel")
# library("data.table")
# library("qtl")
# library("writexl")
# library("BB")
# library("openxlsx")


#' Title The function of dQTG.seq
#' @param dir the path of the output
#' @param filegen the input data
#' @param chr the chromosome
#' @param color the color
#' @param CLO the numbers of CPU
#' @return list
#' @export
#' @examples
#' data(BSA)
#' dQTGseq(dir=tempdir(),filegen=BSA,chr="all",color="blue",CLO=1)
dQTGseq<-function(dir,filegen,chr,color,CLO){



  Readdata2<-function(File){
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
    GeneralInfo<-read.xlsx(File,sheet = "GeneralInfo",colNames = FALSE,na.strings=TRUE)
    gen<-read.xlsx(File,sheet = "Genotype",colNames = FALSE,na.strings=TRUE)
    phe<-read.xlsx(File,sheet = "Phenotype",colNames = FALSE,na.strings=TRUE)
    pos<-read.xlsx(File,sheet = "LinkageMap",colNames = FALSE,na.strings=TRUE)

    result<-as.data.frame(GeneralInfo)
    calculatedata<-list(Species.windowsize(result[1,2],result[5,2]),result[2,2],result[3,2],result[4,2],result[5,2],result[6,2],result[7,2],
                        result[8,2],result[9,2],gen,phe,pos)
    names(calculatedata)<-c("windowsize","FileFormat","Population","np","Smooth_method","Repetition",
                            "DrawPlot","Resolution","Plotformat","gen","phe","pos")
    return(calculatedata)
  }

  if(is.character(filegen)==TRUE){

    if(str_detect(c(filegen), ".xlsx")==TRUE){
      data.calculatedata<-Readdata2(filegen)
      Dodata(dir,data.calculatedata,chr,color,CLO)
    }else{
      data.calculatedata<-Readdata1(filegen)
      Dodata(dir,data.calculatedata,chr,color,CLO)
    }
  }else{
      data.calculatedata<-Readdata1(filegen)
      Dodata(dir,data.calculatedata,chr,color,CLO)

    }

  }









