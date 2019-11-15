library(magrittr)
library(ggpubr)
library(gridExtra)
load("~/conservative.RData")
realPat <- subset(allPep, allPep$typeOfPatient == "1")
virtPat <- subset(allPep, allPep$typeOfPatient == "0")
colnames(virtPat) <- paste("Virt", colnames(virtPat), sep = "_")
totalPat <- cbind.data.frame(realPat, virtPat)

GetDX <- grep("DX",colnames(totalPat),value=TRUE)
GetBrain <- grep("^(Virt_brain|brain).*(?<!aux)$",colnames(totalPat),value=TRUE, perl= T)
GetBrain <- setdiff(GetBrain, c("Virt_brain68.bl", "brain68.bl"))
GetCSF <- grep("CSF",colnames(totalPat),value=TRUE)
GetFDG <- grep("FDG",colnames(totalPat),value=TRUE)
GetCog <- grep("^(Virt_Cog|Cog).*(?<!aux)$",colnames(totalPat),value=TRUE,perl=TRUE)
GetDemogs <- grep("^(Virt_PT|PT|Virt_AGE|AGE)", colnames(totalPat), value = TRUE, perl = TRUE)
GetSNP <- grep("^(Virt_SNP|SNP)", colnames(totalPat), value = TRUE, perl = TRUE)
GetCortical <- grep("^(Virt_brain68|brain68).*(?<!aux)$", colnames(totalPat), value = TRUE, perl = TRUE)
DXs <- c("bl","m06","m12","m24")

PredObs <- function(giveFeatures){
  dfAll <- data.frame()
  for(i in 1:(length(colnames(giveFeatures))/2)){
    print(i)
    grepDX <- grep(DXs[i],colnames(giveFeatures),value=TRUE)  
    dfPred <- as.data.frame(giveFeatures[,grepDX[1]])
    print(dfPred)
    names(dfPred) <- "Values"
    dfPred$Category = 'Simulated'
    dfPred$time = DXs[i]
    #print(dfPred$time)
    dfObs <- as.data.frame(giveFeatures[,grepDX[2]])
    names(dfObs) <- "Values"
    dfObs$time = DXs[i]
    print(dfObs$time)
    dfObs$Category = 'Real'
    AddDFs <- rbind.data.frame(dfPred,dfObs)
    #print(AddDFs)
    dfAll <- rbind.data.frame(dfAll,AddDFs)
  }
  return(dfAll)
}

######dx######
dfDXpred <- totalPat[,GetDX]
dfDx <- PredObs(dfDXpred)  
##BarPlot for DX
dxBar <- dfDx
##For density plot
dfDx$Values <- as.numeric(dfDx$Values)

#####cognition#####
dfCogpred <- totalPat[,GetCog]
dfCog <- PredObs(dfCogpred)
##BarPlot for DX
cogBar <- dfCog
##For density plot                                                                                                                         
dfCog$Values <- as.numeric(dfCog$Values)                                                                                                                     
                                         
###brain####
dfBrainpred <- totalPat[,GetBrain]
dfBrain <- PredObs(dfBrainpred)
brainBar <- dfBrain
##For density plot   
dfBrain$Values <- as.numeric(dfBrain$Values)

###brain 68, cortical####
dfCorticalPred <- totalPat[,GetCortical]
dfCort <- PredObs(dfCorticalPred)
cortBar <- dfCort
##For density plot   
dfCort$Values <- as.numeric(dfCort$Values)

##csf##
dfCSFpred <- totalPat[,GetCSF]
dfcsf <- PredObs(dfCSFpred)
csfBar <- dfcsf
##For density plot   
dfcsf$Values <- as.numeric(dfcsf$Values)

##fdg##
dfFDGpred <- totalPat[,GetFDG]
dfFDG <- PredObs(dfFDGpred)
fdgBar <- dfFDG
##For density plot
dfFDG$Values <- as.numeric(dfFDG$Values)

##education##
dfDemogEducatpred <- totalPat[,c("PTEDUCAT.bl","Virt_PTEDUCAT.bl")]
dfdemogEducat <- PredObs(dfDemogEducatpred)
eduBar <- dfdemogEducat
##For Density plot
dfdemogEducat$Values <- as.numeric(dfdemogEducat$Values)

##gender##
dfDemogGenderpred <- totalPat[,c("PTGENDER.bl","Virt_PTGENDER.bl")]
dfdemogGender <- PredObs(dfDemogGenderpred)
genBar <- dfdemogGender
##For Density plot
dfdemogGender$Values <- as.numeric(dfdemogGender$Values)

##ethnicity##
dfDemogEthcatpred <- totalPat[,c("PTETHCAT.bl","Virt_PTETHCAT.bl")]
dfdemogEthcat <- PredObs(dfDemogEthcatpred)
ethBar <- dfdemogEthcat
##For density plot
dfdemogEthcat$Values <- as.numeric(dfdemogEthcat$Values)

##marriage##
dfDemogMarrypred <- totalPat[,c("PTMARRY.bl","Virt_PTMARRY.bl")]
dfdemogMarry <- PredObs(dfDemogMarrypred)
marBar <- dfdemogMarry
##For density plot
dfdemogMarry$Values <- as.numeric(dfdemogMarry$Values)

##race##
dfDemogRacepred <- totalPat[,c("PTRACCAT.bl","Virt_PTRACCAT.bl")]
dfdemogRace <- PredObs(dfDemogRacepred)
raceBar <- dfdemogRace
##For density plot
dfdemogRace$Values <- as.numeric(dfdemogRace$Values)

##age##
dfDemgAgepred <- totalPat[,c("AGE.bl","Virt_AGE.bl")]
dfdemogAge <- PredObs(dfDemgAgepred)
ageBar <- dfdemogAge
##For density plot
dfdemogAge$Values <- as.numeric(dfdemogAge$Values)

##snp##
dfSNPpred <- totalPat[,GetSNP]
dfSNP <- PredObs(dfSNPpred)
snpBar <- dfSNP
##For density plot
dfSNP$Values <- as.numeric(dfSNP$Values)

##changing the values in time column according to the group##
dfDx$time = paste("DX", dfDx$time, sep = ".")
dfCog$time = paste("Cog", dfCog$time, sep = ".")
dfBrain$time = paste("brain", dfBrain$time, sep = ".")
dfCort$time = paste("brain68", dfCort$time, sep = ".")
dfcsf$time = paste("csf", dfcsf$time, sep = ".")
dfFDG$time = paste("FDG", dfFDG$time, sep = ".")
dfdemogEducat$time = "Education"
dfdemogGender$time = "Gender"
dfdemogEthcat$time = "Ethnicity"
dfdemogMarry$time = "Marital Status"
dfdemogRace$time = "Race"
dfdemogAge$time = "Age"
dfSNP$time = "SNP"

###for denisty##
dfAll <- rbind.data.frame(dfDx, dfCog, dfBrain, dfCort, dfcsf, dfFDG, dfdemogEducat, dfdemogGender, dfdemogEthcat,
                          dfdemogMarry, dfdemogRace, dfdemogAge, dfSNP)
dfAll$Values <- as.numeric(dfAll$Values) 

names(dfAll)[3] <- "variable"
names(dfAll)[1] <- "value"
png(file = "/output/densityADNI.png",
    width = 12, height = 8, units = 'in', res = 300)
ggplot(data = dfAll, aes(value,colour=Category)) +  facet_wrap(~ variable, scales = "free") +theme(strip.text.x = element_text(size=4)) + geom_density()+theme_classic() + theme(axis.text = element_text(size=10),
                                                                                                                                                                                 axis.text.x = element_text(size = 10),
                                                                                                                                                                                 axis.title = element_text(size = 10, face = "bold"),
                                                                                                                                                                                 strip.text = element_text(size = 10))  
dev.off()


##mutual informartion between real and virtual subjects##
library(infotheo)
auxVar <- grep("aux", colnames(real), value = TRUE)
features <- setdiff(colnames(real), auxVar)
miDf <- data.frame()
for(name in features){
  mi <- mutinformation(real[,name], simulated[,name], method="emp")
  miDf<-  rbind(miDf, data.frame(name, mi))
}



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("maigesPack")
library(maigesPack)


pValueDf = data.frame()
for(name in features){
  pVal <- bootstrapMI(as.numeric(real[,name]), as.numeric(simulated[,name]), bRep=1000, ret="p-value")
  pValueDf = rbind(pValueDf, data.frame(name, pVal))
}
pValADNI <- pValueDf

##for histogram
##changing the values in time column according to the group##
dxBar$time = paste("DX", dxBar$time, sep = ".")
cogBar$time = paste("Cog", cogBar$time, sep = ".")
brainBar$time = paste("brain", brainBar$time, sep = ".")
cortBar$time = paste("brain68", cortBar$time, sep = ".")
csfBar$time = paste("csf", csfBar$time, sep = ".")
fdgBar$time = paste("FDG", fdgBar$time, sep = ".")
eduBar$time = "Education"
genBar$time = "Gender"
ethBar$time = "Ethnicity"
marBar$time = "Marital Status"
raceBar$time = "Race"
ageBar$time = "Age"
snpBar$time = "SNP"
dfAllBar <- rbind.data.frame(dxBar, cogBar, brainBar, cortBar, csfBar, fdgBar, eduBar, genBar, ethBar,
                          marBar, raceBar, ageBar, snpBar)
names(dfAllBar)[3] <- "variable"
dfAllBar$variable <- as.factor(dfAll$variable)
names(dfAllBar)[1] <- "value"

names(pValADNI)[2] = "variable"
pValADNI$X1 <- NULL
pValADNI$variable <- as.factor(pValADNI$variable )
pValADNI$variable = levels(dfAllBar$variable)
pValADNICp <- pValADNI
pValADNICp$variable = as.character(pValADNICp$variable)
pValADNICp$pVal <- paste(pValADNICp$variable, ", ", "pval=", pValADNICp$pVal, sep = "")
colnames(pValADNICp) <- pValADNICp[1,]
pValADNICp <- t(pValADNICp)
pValADNICp <- as.data.frame(pValADNICp)
pValADNICp <- pValADNICp[-1,]
pValADNICp <- as.list(pValADNICp)
variable_labeller <- function(variable,value){
  return(pValADNICp[value])
}
png(file = "/output/barPlotADNI.png",
    width = 15, height = 13, units = 'in', res = 300)
ggplot(data = dfAllBar, aes(value,fill=Category, colour=Category)) +  facet_wrap(~ variable, scales = "free", labeller = variable_labeller) +theme(strip.text.x = element_text(size=4)) + geom_histogram(position="dodge",stat="count")+theme_classic() + theme(axis.text = element_text(size=10),
                                                                                                                                                                                 axis.text.x = element_text(size = 10,  angle=90, hjust = 1),
                                                                                                                                                                                 axis.title = element_text(size = 10, face = "bold"),
                                                                                                                                                                                 strip.text = element_text(size = 10))



dev.off()

