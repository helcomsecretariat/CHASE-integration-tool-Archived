library("dplyr")
library("tidyr")

#===============================================================================
# function Assessment
Assessment<- function(assessmentdata,summarylevel=1){

  # Confidence Penalty Criteria (minimum numbers of indicators)
  CountHM<-2
  CountOrg<-3
  PenaltyHM<-0.5
  PenaltyOrg<-0.5

  # Get column names from the input data
  cnames<-names(assessmentdata)

  # If CR is not included, then Threshold and Status are required columns
  if('CR' %in% toupper(cnames)){
    UserCR=TRUE
    requiredcols <- c("Matrix","Substance")
  }else{
    UserCR=FALSE
    requiredcols <- c("Matrix","Substance","Threshold","Status")
  }

  extracols <- c("Waterbody","Response","Type","Confidence")

  #Check column names in the imported data
  nimp = ncol(assessmentdata)
  nreq = length(requiredcols)
  nextra = length(extracols)

  ok <- rep(0, nreq)
  okextra <- rep(0, nextra)
  foundresponse=FALSE


  for (i in 1:nimp){
    for (j in 1:nreq){
      if(toupper(requiredcols[j])==toupper(cnames[i])){
        names(assessmentdata)[i] <- requiredcols[j]
        ok[j]=1
      }
    }
    for (j in 1:nextra){
      if(toupper(extracols[j])==toupper(cnames[i])){
        names(assessmentdata)[i] <- extracols[j]
        okextra[j]=1
      }
    }
  }

  for(j in 1:nextra){
    if(okextra[j]==0){
      if(regexpr("Confidence",extracols[j],ignore.case=TRUE)>0){
        assessmentdata[[extracols[j]]]<-'Low'
      }else{
        assessmentdata[[extracols[j]]]<-1
      }
    }
  }

  n<-sum(ok, na.rm = TRUE)

  if(n<nreq){
    # The required columns were not found in the input data
    message("Error in CHASE Assessment. Required column(s) were not found in the input data:")
    for (j in 1:nreq){
      if(ok[j]!=1){
        cat(paste0("    ",requiredcols[j]),"\n")
      }
    }
    return(NA)
  }else{
    # The required columns are present - do the assessment


    # Change order of matrices factors
    mat1<-data.frame(unique(assessmentdata$Matrix))
    names(mat1)[1] <- 'Matrix'
    mat1$char<-as.character(mat1$Matrix)
    mat1$len<-nchar(mat1$char)
    mat1<-arrange(mat1,len)

    assessmentdata$Matrix <- factor(assessmentdata$Matrix, levels = mat1$char)

    # All combinations of matrices and waterbodies
    # This is used to ensure that a NA is returned where the combinations are missing
    waterbodies<-unique(assessmentdata$Waterbody)
    matrices<-unique(assessmentdata$Matrix)
    matrices<-expand.grid(waterbodies, matrices)
    names(matrices)[1] <- 'Waterbody'
    names(matrices)[2] <- 'Matrix'

    if(UserCR==FALSE){
      assessmentdata$CR<-ContaminationRatio(assessmentdata$Threshold,assessmentdata$Status,assessmentdata$Response)
    }
    # assessmentdata$ConfScore<-mapply(ConfValue,assessmentdata$Confidence)
    assessmentdata$ConfScore<-Confidence(assessmentdata$ConfStatus,assessmentdata$ConfThresh)


    # print(head(assessmentdata$ConfScore))

    # Add columns specifying if the inidcator is an organic or heavy metal. Used in confidence penalty calculations
    assessmentdata$HM<-mapply(IsHeavyMetal,assessmentdata$Type)
    assessmentdata$Org<-mapply(IsOrganic,assessmentdata$Type)


    MatchListBioEffects<-c("bioeffect","bioeffects","bio effects","bio effect",
                           "biological effects","biological effect")
    MatchListBiota<-c("biota","biot")
    MatchListSed<-c("sediment","sed","sedi")

    # Count the number of distinct indicators in each category (type) HM and Org
    DistinctIndicators<-assessmentdata %>%
      filter(!is.na(Type),!is.na(CR)) %>%
      group_by(Waterbody,Matrix,Type,Substance) %>%
      summarise(Count=n()) %>%
      ungroup
    DistinctIndicators$HM<-mapply(IsHeavyMetal,DistinctIndicators$Type)
    DistinctIndicators$Org<-mapply(IsOrganic,DistinctIndicators$Type)

    QEtypeCount<-DistinctIndicators %>%
      group_by(Waterbody) %>%
      summarise(HM=sum(HM,na.rm=TRUE),Org=sum(Org,na.rm=TRUE))


    QEtypeCount$HM_OK<-ifelse(QEtypeCount$HM<CountHM,FALSE,TRUE)
    QEtypeCount$Org_OK<-ifelse(QEtypeCount$Org<CountOrg,FALSE,TRUE)
    QEtypeCount$Penalty<-(1-((1-PenaltyHM*(!QEtypeCount$HM_OK))*(1-PenaltyOrg*(!QEtypeCount$Org_OK))))
    QEtypeCount$HM_OK<-NULL
    QEtypeCount$Org_OK<-NULL


    QEdata<-assessmentdata %>%
      filter(!is.na(CR)) %>%
      group_by(Waterbody,Matrix) %>%
      summarise(sumCR=sum(CR,na.rm = TRUE),
                sumConf=sum(ConfScore,na.rm = TRUE),
                countHM=sum(HM,na.rm = TRUE),
                countOrg=sum(Org,na.rm = TRUE),
                Count=n()) %>%
      ungroup()


    #QEdata$ConSum<-QEdata$sumCR/sqrt(QEdata$Count)
    QEdata$IsBio<-ifelse(tolower(QEdata$Matrix) %in% MatchListBioEffects,TRUE,FALSE)
    QEdata$IsBiota<-ifelse(tolower(QEdata$Matrix) %in% MatchListBiota,TRUE,FALSE)
    QEdata$IsSed<-ifelse(tolower(QEdata$Matrix) %in% MatchListSed,TRUE,FALSE)
    QEdata$ConSum<-QEdata$sumCR/ifelse(QEdata$IsBio,QEdata$Count,sqrt(QEdata$Count))

    QEdata$ConfScore<-QEdata$sumConf/QEdata$Count


    QEdata$Confidence<-mapply(ConfidenceStatus,QEdata$ConfScore,TRUE)
    QEdata$MultiplierHM<-NULL
    QEdata$MultiplierOrg<-NULL
    QEdata$sumConf<-NULL
    QEdata$IsBio<-NULL
    QEdata$sumCR <- NULL
    QEdata$Count <- NULL



    QEspr<-QEdata %>%
      select(Waterbody,Matrix,ConSum) %>%
      spread(Matrix,ConSum)

    QEspr1<-QEspr

    QEdata$QEStatus<-CHASEStatus(QEdata$ConSum,2)
    QEdata<-left_join(matrices,QEdata,c('Waterbody','Matrix'))
    QEdata<-arrange(QEdata,Waterbody,Matrix)
    QEdataOut<-QEdata %>%
      select(Waterbody,Matrix,ConSum,QEStatus,ConfScore,Confidence)

    CHASE<-QEdata %>%
      group_by(Waterbody) %>%
      summarise(ConSum=max(ConSum, na.rm = TRUE),Sed=sum(IsSed, na.rm = TRUE),
                Biota=max(IsBiota, na.rm = TRUE),ConfScore=mean(ConfScore, na.rm = TRUE))

    #CHASE$Waterbody<-NULL
    CHASEQE<-QEdata %>%
      select(Waterbody,Matrix,ConSum,QEStatus) %>%
      inner_join(CHASE, by=c("Waterbody"="Waterbody","ConSum"="ConSum"))

    CHASEQE<-rename(CHASEQE,Status=QEStatus,Worst=Matrix)

    CHASEQE$Confidence <- NA

    CHASEQE<-CHASEQE %>%
      select(Waterbody,Worst,ConSum,Status,ConfScore,Confidence) %>%
      left_join(QEtypeCount,by=c("Waterbody"="Waterbody"))

  #  cat(nrow(CHASEQE))
  #  print(CHASEQE$Penalty)
  #  print(CHASEQE$ConfScore)

    CHASEQE$ConfScore<-CHASEQE$ConfScore*(1-CHASEQE$Penalty)
    CHASEQE$Confidence<-mapply(ConfidenceStatus,CHASEQE$ConfScore,TRUE)
    CHASEQE$Penalty<-scales::percent(CHASEQE$Penalty)

   # print(CHASEQE$ConfScore)

    assessmentdata$HM<-NULL
    assessmentdata$Org<-NULL
    if(!'RESPONSE' %in% toupper(cnames)){
      assessmentdata$Response<-NULL
    }

    assessmentdata<-assessmentdata %>%
      left_join(rename(QEdataOut,QEConfidence=Confidence,QEConfScore=ConfScore),c('Waterbody','Matrix'))

    QEspr<-inner_join(QEspr, CHASEQE, 'Waterbody')

    if(summarylevel==2){
      return(QEspr)
    }else if(summarylevel==3){
      return(QEdataOut)
    }else if(summarylevel==4){
      return(CHASEQE)
    }else{
      return(assessmentdata)
      #return()
    }
    #
  }
}

#===============================================================================
# function ContaminationRatio
ContaminationRatio<- function(threshold, status, response=1){
  # If response is not specified, it will be assumed to be positive
  # i.e. ContaminationRatio increases (worsens) with increasing status value
  if (missing(response)){
    response=1
  }
  response<-ifelse(is.na(response), 1, response)

  # ContaminationRatio calculated depending on Response direction
  cr<-ifelse(response>0, status/threshold, threshold/status)
  return(cr)
}

#===============================================================================
#Function CHASEStatus
CHASEStatus<-function(CRsum,nCat=5){
  if(nCat==5){
    status<-ifelse(CRsum>0.5, "Good", "High")
    status<-ifelse(CRsum>1, "Moderate", status)
    status<-ifelse(CRsum>5, "Poor", status)
    status<-ifelse(CRsum>10, "Bad",status )
  }else{
    status<-ifelse(CRsum>1, "Not good","Good")
  }
  return(status)
}

CHASEStatus1<-function(CRsum){
  status<-ifelse(CRsum>0.5, "Good", "High")
  status<-ifelse(CRsum>1, "Moderate", status)
  status<-ifelse(CRsum>5, "Poor", status)
  status<-ifelse(CRsum>10, "Bad",status )
  return(status)
}
# Colours associated with Status classes - used by Shiny App
AddColours<-function(CRsum){
  co<-ifelse(CRsum>0.5, '#66FF66', '#3399FF')
  co<-ifelse(CRsum>1, '#FFFF66', co)
  co<-ifelse(CRsum>5, '#FF9933', co)
  co<-ifelse(CRsum>10, '#FF6600',co)
  return(co)
}


#===============================================================================
# function Threshold
Threshold <- function(substance){
  substances<- read.csv('./data/substances.csv', header = TRUE,  sep=";")
  return(value)
}

# Function to calculate numeric confidence from string --------------------------------------------------------------------
# The function will always return a numeric value between 0 and 1, depending on the argument sConf
# Given a numeric argument between 0 and 1, the function returns the same value
# e.g. ConfValue(0.37) returns a value of 0.37
# Passing a numeric argument less than 0, the function returns a value of 0
# Passing a numeric argument greater than 1, the function returns a value of 1
# The function recognizes the following words and returns the respective values:
#    High = 1.0
#    Intermediate = 0.5
#    Medium = 0.5
#    Moderate = 0.5
#    Low = 0.0
# The function is case-insensitive.
# Starting from the leftmost characte, the function recognizes any part of the key words
# e.g. "H", "hi", "hig" will all result in a value of 1.0
#      "med", "m", "int", "I", "in" will all return values of 0.5
#      "lo", "l" will all give a value of 0.0
#
# Any other argument passed to the function will give a result equal to the argument NAvalue (default=0)

ConfValue<-function(sConf,NAvalue=NA){
  if(is.numeric(sConf)){
    return(sConf)
  }else{
    sConf<-tolower(sConf)
    l<-nchar(sConf)
    if(l<1){
      return(NAvalue)
    }else{
      desc<-c("low          ","intermediate","medium          ","moderate          ","high          ")
      value<-c(0,0.5,0.5,0.5,1)
      df<-data.frame(desc,value)
      df$desc<-substr(df$desc, 1, l)
      if(sConf %in% df$desc){
        n<-df$value[df$desc==sConf][1]
        return(n)
      }else{
        n<-suppressWarnings((as.numeric(sConf)))
        if(is.na(n)){
          n<-NAvalue
        }else{
          n<-min(1,max(0,n))
        }
        return(n)
      }
    }
  }
}

# Confidence<-function(Confidence1,Confidence2,Confidence3){
#   n1<-lapply(Confidence1,ConfValue)
#   if(missing(Confidence2)){
#
#     n<-n1
#   }else{
#
#     n2<-lapply(Confidence2,ConfValue)
#     if(missing(Confidence3)){
#       n<-mapply(ConfAvg, n1, n2)
#     }else{
#       n3<-lapply(Confidence3,ConfValue)
#       n<-mapply(ConfAvg, n1, n2, n3)
#     }
#   }
#   return(n)
# }

Confidence<-function(Confidence1,Confidence2,Confidence3){
  n1=lapply(Confidence1,ConfValue)
  n2=lapply(Confidence2,ConfValue)
  # print(head(n1))
  # print(head(n2))
  if(missing(Confidence3)){
    n=mapply(ConfAvg, n1, n2)
  }else{
    n3=lapply(Confidence3,ConfValue)
    n=mapply(ConfAvg, n1, n2, n3)
  }
  return(n)
}


ConfAvg<-function(C1,C2,C3){
  if(missing(C3)){
    n=(C1+C2)/2
  }else{
    n=(C1+C2+C3)/3
  }
  return(n)
}


ConfidenceS<-function(arrayConf){
  arrayConf<-sapply(arrayConf, ConfValue)
  return(arrayConf)
}

ConfidenceStatus<-function(Score,Roman=FALSE){
  Status<-ifelse(Roman==TRUE,
                 ifelse(Score<0.5,"Class III",
                        ifelse(Score<0.75,"Class II","Class I")),
                 ifelse(Score<0.5,"Low",
                        ifelse(Score<0.75,"Medium","High"))
  )
  return(Status)
}



IsHeavyMetal<-function(sType,NAvalue=NA){
  n<-NAvalue
  if(regexpr("hm",sType,ignore.case=TRUE)>0){
    n<-1
  }
  if(regexpr("heavy",sType,ignore.case=TRUE)>0){
    n<-1
  }
  return(n)
}

IsOrganic<-function(sType,NAvalue=NA){
  n<-NAvalue
  if(regexpr("organ",sType,ignore.case=TRUE)>0){
    n<-1
  }
  if(regexpr("org",sType,ignore.case=TRUE)>0){
    n<-1
  }
  return(n)
}
