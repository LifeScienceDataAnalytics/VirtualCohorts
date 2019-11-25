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
png(file = "output/densityADNI.png",
    width = 12, height = 8, units = 'in', res = 300)
ggplot(data = dfAll, aes(value,colour=Category)) +  facet_wrap(~ variable, scales = "free") +theme(strip.text.x = element_text(size=4)) + geom_density()+theme_classic() + theme(axis.text = element_text(size=10),
                                                                                                                                                                                 axis.text.x = element_text(size = 10),
                                                                                                                                                                                 axis.title = element_text(size = 10, face = "bold"),
                                                                                                                                                                                 strip.text = element_text(size = 10))  
dev.off()

##for histogram
## first perform chi square test, and calculate adjusted p values using bonferroni correction
auxVar <- grep("aux", colnames(real), value = TRUE)
features <- setdiff(colnames(real), auxVar)
chiDf <- data.frame()
library(entropy)
for(name in features){
  print(name)
  df <- rbind(table(real[,name]), table(simulated[,name]))
  chiSQ =  chisq.test(df)
  chiDf <-  rbind(chiDf, data.frame(name, p.value = chiSQ$p.value))
  df <- NULL
}

chiDf$adj.p = p.adjust(chiDf$p.value, method = "bonferroni", n = length(chiDf$p.value))


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
dfAllBar$variable <- as.factor(dfAllBar$variable)
names(dfAllBar)[1] <- "value"

names(pValADNI)[2] = "variable"

##trace back the ranges of factor values for each variable
decimalnumcount<-function(x){
  x<-gsub("(.*)(\\.)|([0]*$)","",x)
  nchar(x)
}
disc <- meta_rf
names(disc)[12] <- "AGE.bl"
names(disc)[13] <- "FDG.bl"
dfMerge <- data.frame()
for(names in colnames(disc)){
  new_dat <- data.frame()
  #print(names)
  sr = sort(create_breaks(disc[,names], input_df[,"DX.bl"], method="dt"), decreasing = FALSE)
  #print(sr)
  mn <- data.frame()
  if(length(sr) == 1){
    mn[1,"value"] <- 1
    mn[1, "valueRange"] <-  paste("<=", round(sr, digits = 4), sep = "")
    mn[1, "variable"] <- names
    mn[2,"value"] <- 2
    mn[2, "valueRange"] <-  paste(">", round(sr, digits = 4), sep = "")
    mn[2, "variable"] <- names
    new_dat <- mn
  }

  if(length(sr) > 1) {
    print(names)
    disc[,names] <- create_bins(disc[,names], sr)
    print(sr)
    mn <- data.frame()
    for(j in 1:length(sr)){
      print(j)
      mn[j,"value"] <- j+1
      mn[j, "valueRange"] <- paste((round(sr[j], digits = 4)) + (as.numeric(paste(0,".", strrep(0, decimalnumcount(round(sr[j], digits = 4)-1)), "1", sep =""))), "-", round(sr[j+1], digits = 4), sep="")
      mn[j, "variable"] <- names
  }
  new_dat <- rbind(c(1, paste("<=", round(sr[1],digits = 4), sep ="")), mn)
  new_dat[1, "variable"] <- names
  new_dat[length(table(disc[,names])),"valueRange"] <-  paste(">=", round(sr[length(sr)],digits = 4), sep ="")
  print(new_dat[, "valueRange"])
  }
  new_dat$valueRange <- factor(new_dat$valueRange, levels = new_dat$valueRange[order(new_dat$value)])
  dfMerge <- rbind.data.frame(dfMerge, new_dat)
} 

sameValues <- setdiff(unique(dfAllBar$variable), unique(dfMerge$variable))
dfAllNoCateog <- dfAllBar[!dfAllBar$variable %in% sameValues,]
dfAllCateog <- subset(dfAllBar, dfAllBar$variable %in% sameValues)
dfAllCateog$valueRange <- dfAllCateog$value
dfMergeAll1 <- join(dfMerge, dfAllNoCateog)
dfMergeAll1 <- dfMergeAll1[names(dfAllCateog)]
dfMergeAll2 <- rbind.data.frame(dfMergeAll1, dfAllCateog)
dfMergeAll2$value <- as.factor(dfMergeAll2$value)

setDT(dfAllBar)[setDT(dfMergeAll2), valueRange := i.valueRange, on=c("value", "Category", "variable")]
dfAllBar$variable = factor(dfAllBar$variable, levels= unique(dfAllBar$variable))

chiDf$adj.p <- round(chiDf$adj.p, digits = 4)
chiDf$adj.p[chiDf$adj.p<0.001] <- "<0.001"


chiDf[chiDf=="PTEDUCAT.bl"]<-"Education"
chiDf[chiDf=="PTGENDER.bl"]<-"Gender"
chiDf[chiDf=="PTETHCAT.bl"]<-"Ethnicity"
chiDf[chiDf=="PTMARRY.bl"]<-"Marital Status"
chiDf[chiDf=="PTRACCAT.bl"]<-"Race"
chiDf$name <- factor(chiDf$name, levels =unique(as.character(dfAllBar$variable)))
chiDf <- chiDf[order(chiDf$name),]
adjpValADNI <- chiDf
adjpValADNI$p.value <- NULL
names(adjpValADNI)[2] = "variable"
adjpValADNI$X1 <- NULL
adjpValADNI$variable <- factor(adjpValADNI$variable)
adjpValADNICp <- adjpValADNI
adjpValADNICp$variable = as.character(adjpValADNICp$variable)
adjpValADNICp$adj.p <- paste(adjpValADNICp$variable, ", ", "adj.p=", adjpValADNICp$adj.p, sep = "")
adjpValADNICp <- t(adjpValADNICp)
colnames(adjpValADNICp) <- adjpValADNICp[1,]
adjpValADNICp <- as.data.frame(adjpValADNICp)
adjpValADNICp <- adjpValADNICp[-1,]
adjpValADNICp <- as.list(adjpValADNICp)
variable_labeller <- function(variable,value){
  return(adjpValADNICp[value])
}

require(data.table)

png(file = "/Users/Meems/Documents/PhDWork/BayesianNetworkAD/CodeAndResultsforPaper/PaperRevisions/barPlotADNI_Vs5.png",
    width = 20, height = 20, units = 'in', res = 300)
p <- ggplot(data = dfAllBar, aes(valueRange, fill = Category, colour=Category)) + geom_histogram(aes(y = (..count..)/sum(..count..)),position="dodge",stat="count") +theme_classic() + theme(legend.title = element_text(face = "bold", size = 14),
                                                                                                                                                                                 legend.text=element_text(size=12),
                                                                                                                                                                                 axis.text = element_text(size=14),
                                                                                                                                                                                 axis.text.x = element_text(size = 14,  angle=90, hjust = 1),
                                                                                                                                                                                 axis.text.y = element_text(size = 14),
                                                                                                                                                                                 axis.title = element_text(size = 14, face = "bold"),
                                                                                                                                                                                 strip.text = element_text(size = 14))
p1 <- p + facet_wrap(~ variable, scales = "free", labeller=variable_labeller) +theme(strip.text.x = element_text(size=14, face = "bold"),
                                                                               strip.background = element_rect(
                                                                                 color="black", size=1.5, linetype="solid"
                                                                               )) + theme(panel.spacing = unit(1, "lines"))

p1 + labs(x = "Values", y = "Relative Frequencies") + theme(axis.text=element_text(size=14)) 


dev.off()



