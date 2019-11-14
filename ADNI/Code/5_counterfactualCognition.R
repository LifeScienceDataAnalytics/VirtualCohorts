#########################contrafactual for cognition####################################
load("~/bnCreation.RData")
#load adni data#
library(h2o)
library(randomForest)
library(data.table)
library(bnlearn)
library(arules)
library(missForest)
library(parallel)
library(caret)
library(ggplot2)
library(gtools)
library(plyr)
library(stringi)
library(stringr)
library(doMC)
library(gRain)
library(lattice)
library(binst)

##find all the patients who remain normal throughout##
dataAlz$dementiaOrNot <- apply(dataAlz[,diag], 1, function(x)all(x == "Dementia"))

cogAtBL <- grep("(CDRSB|ADAS11|ADAS13|MMSE|FAQ|RAVLT.*)\\.bl", colnames(impData),value = TRUE)
diag <- grep("^DX", colnames(impData), value = TRUE)
aggDiagnosis <- ddply(adnimerge, .(PTID), summarize, DX = toString(DX))
aggDiagnosisMonth <- ddply(adnimerge, .(PTID), summarize, Month = toString(Month))
nlPatientID <- c()
length <- c()
for(i in 1:nrow(aggDiagnosis)){
  testDf <- aggDiagnosis[i,2]
  test <- unlist(strsplit(testDf, split= " "))
  length <- c(length, length(test))
  if(testDf != "NA"){
    if(sum(str_count(sapply(test, grepl, "NL,|NL|NA,|NA"), "TRUE")) == length(test)){
    nlPatientID <- c(nlPatientID, aggDiagnosis[i,1])
    }
  }
}

aggDiagnosisMerge <- merge(aggDiagnosis, aggDiagnosisMonth)
aggDiagnosisMerge$TotalNumberOfMonths <- length
##aggregate for the normal group##
nl_pat <- subset(aggDiagnosisMerge, aggDiagnosisMerge$PTID %in% nlPatientID)
##patients who remain NL throughout##
controlGroup <- subset(adnimerge, adnimerge$PTID %in% nlPatientID)
controlGroupCog <- controlGroup[,c("PTID",cogAtBL)]
controlGroupCog <- controlGroupCog[complete.cases(controlGroupCog), ]
controlGroupCog <- controlGroupCog[!duplicated(controlGroupCog), ]
##Find the median of all the normal subjects##
medianCog <- as.data.frame(sapply(controlGroup[,cogAtBL],function(x)(median(x, na.rm = TRUE))))
##make rows as first column ##
library(data.table)
medianCog <- setDT(medianCog, keep.rownames = TRUE)[]
##change the column names##
names(medianCog)[1] <- "CogTest"
names(medianCog)[2] <- "Median"

##find all the dementia subjects##
cog.dementia <- subset(impData[, c(cogAtBL, diag)], impData$DX.bl == "Dementia")
allDemids <- c()
for(i in 1:nrow(cog.dementia)){
  if((cog.dementia[i, "DX.bl"] == "Dementia") & (cog.dementia[i, "DX.m06"] == "Dementia") & (cog.dementia[i, "DX.m12"] == "Dementia") & 
     (cog.dementia[i, "DX.m24"] == "Dementia")){
    allDemids <- c(allDemids, i)
  }
}

##find all subjects  remaining in MCI and MCI convertors##
allMCIbl <- which(impData$DX.bl == "MCI")
allDementiaMCI <- c(allDemids, allMCIbl)
GetCog <- impData[,cogAtBL]
GetCog[allDementiaMCI,"MMSE.bl"] = medianCog[which(medianCog$CogTest == "MMSE.bl"), "Median"]
GetCog[allDementiaMCI,"CDRSB.bl"] = medianCog[which(medianCog$CogTest == "CDRSB.bl"), "Median"]
GetCog[allDementiaMCI,"ADAS11.bl"] = medianCog[which(medianCog$CogTest == "ADAS11.bl"), "Median"]
GetCog[allDementiaMCI,"ADAS13.bl"] = medianCog[which(medianCog$CogTest == "ADAS13.bl"), "Median"]
GetCog[allDementiaMCI,"FAQ.bl"] = medianCog[which(medianCog$CogTest == "FAQ.bl"), "Median"]
GetCog[allDementiaMCI,"RAVLT.forgetting.bl"] = medianCog[which(medianCog$CogTest == "RAVLT.forgetting.bl"), "Median"]
GetCog[allDementiaMCI,"RAVLT.immediate.bl"] = medianCog[which(medianCog$CogTest == "RAVLT.immediate.bl"), "Median"]
GetCog[allDementiaMCI,"RAVLT.learning.bl"] = medianCog[which(medianCog$CogTest == "RAVLT.learning.bl"), "Median"]
GetCog[allDementiaMCI,"RAVLT.perc.forgetting.bl"] = medianCog[which(medianCog$CogTest == "RAVLT.perc.forgetting.bl"), "Median"]



ModMMSECDRSB <- GetCog
colnames(ModMMSECDRSB) <- cogAtBL

h2o.init()
ModMMSECDRSB = as.h2o(ModMMSECDRSB)
cogModel.bl <- h2o.loadModel(paste(getwd(), "/Auto_Model_New/", meta_feature_cog_list[[1]][[4]], sep =""))


ModMMSECDRSB <- as.data.frame(h2o.deepfeatures(cogModel.bl,ModMMSECDRSB, meta_feature_cog_list[[1]][[3]]))
colnames(ModMMSECDRSB) = paste("Cog", "bl",sep= ".")


getCogscore <- as.data.frame(impData[,cogAtBL])

cogOriginalAutoen <- as.vector(meta_feature_cog_list[[1]][[2]]) 
cogModifiedAutoen <- as.vector(ModMMSECDRSB)
DiscModCog <- DiscAutoEnADNI
DiscModCog$Cog.bl <- NULL
DiscModCog$Cog.bl<- cogModifiedAutoen$Cog.bl

dt_breaks <- create_breaks(cogOriginalAutoen$Cog.bl, DiscAutoEnADNI$DX.bl, method="dt")
DiscModCog$Cog.bl <- create_bins(cogModifiedAutoen$Cog.bl, dt_breaks)

##Map the binned data to the factors in r##
bins <- c("<-0.06151640", 
          "-0.06151639 to  0.00978962", 
          "0.00978963 to  0.03233920", 
          "0.03233921 to 0.11123562",
          "0.11123563 to 0.38520119",
          ">0.38520119")
facLevels <- c(1,2,3,4,5,6)
mapbinsToLevels <- cbind.data.frame(bins, facLevels)

DiscModCog <- as.data.frame(lapply(DiscModCog, factor))
finaltestBN <- finalBN
finaltestBN <- drop.arc(finaltestBN, "SNP.bl", "Cog.bl")
fittedNew = bn.fit(finaltestBN, DiscAutoEnADNI, method = "bayes")
fittedGrain <- as.grain(fittedNew)

set_node_in_model <- setFinding(fittedGrain, nodes = c("Cog.bl","allbrainvol.bl.aux","csf.bl.aux",
                                                       "snp.aux", "brain68.aux", "allbrainvol.m06.aux",
                                                       "CogScore.m06.aux", "allbrainvol.m12.aux", "CogScore.m12.aux",
                                                       "allbrainvol.m24.aux", "CogScore.m24.aux"), 
                                                        c("6", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"))

set.seed("12345")
sim_based_on_fixed_edges = simulate(set_node_in_model, n= nrow(DiscAutoEnADNI))

cog_pred <- data.frame(1:689)
cog_pred = cbind(sim_based_on_fixed_edges[,diag] ,DiscModCog[,diag])
colnames(cog_pred) = c(paste("pred_", colnames(cog_pred[,1:4]), sep="") , colnames(cog_pred[,5:8]))

GetDxfromCog <- grep("^(pred_DX|DX).*(?<!aux)$",colnames(cog_pred),value=TRUE,perl=TRUE)

PredObs <- function(giveFeatures, DXs){
  dfAll <- data.frame()
  for (i in 1:(length(colnames(giveFeatures))/2)){
    print(i)
    grepDX <- grep(DXs[i],colnames(giveFeatures),value=TRUE)  
    print(grepDX)
    dfPred <- as.data.frame(giveFeatures[,grepDX[1]])
    names(dfPred) <- "Values"
    dfPred$Category = 'Simulated'
    dfPred$time = DXs[i]
    #print(dfPred$time)
    dfObs <- as.data.frame(giveFeatures[,grepDX[2]])
    names(dfObs) <- "Values"
    dfObs$time = DXs[i]
    #print(dfObs$time)
    dfObs$Category = 'Real'
    AddDFs <- rbind.data.frame(dfPred,dfObs)
    #print(AddDFs)
    dfAll <- rbind.data.frame(dfAll,AddDFs)
  }
  return(dfAll)
}


dfDiag <- PredObs(cog_pred[,GetDxfromCog], c("bl", "m06", "m12", "m24"))  

png(file = "~/output/contraCogDX.png",width = 8, height = 6, units = 'in', res = 300)
ggplot(dfDiag,aes(x=Values, group=Category,fill = Category)) + ylab("Fraction of subjects") + xlab("Stages") +
geom_histogram(aes(y=stat(count/689)),position="dodge", binwidth=0.25, stat="count") +  facet_wrap(~dfDiag$time, scales= "free")  + theme(axis.text = element_text(size=10),
                                                                                                                      axis.text.x = element_text(size = 10, angle=90, hjust=1),
                                                                                                                      axis.title = element_text(size = 10, face = "bold"),
                                                                                                                      strip.text = element_text(size = 10))

dev.off()
save.image("~/counterfactualCognition.RData")

