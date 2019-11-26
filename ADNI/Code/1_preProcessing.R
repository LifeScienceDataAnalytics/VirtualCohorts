
##Script name: 1_preProcessing.R
##Purpose of Script: preprocesses ADNI data, deals with missing values, creation of auxilliary variables and imputation. 
##Author: Meemansa Sood
##Date Created:October 2018

##load all the required libraries
library(data.table)
library(bnlearn)
library(arules)
library(missForest)
library(parallel)
library(caret)
library(ggplot2)
library(gtools)
library(plyr)
library(stringr)
library(doMC)
library(readxl)
library(dplyr)

##load adni dataset adnimerge (requires ADNI data access)##

##extract all the convertors and the de-novo subjects##
aggDiagnosis <- ddply(adnimerge, .(PTID), summarize, DX = toString(DX))
patientid <- c()
for(i in 1:nrow(aggDiagnosis)){
  if(grepl("Dementia", aggDiagnosis[i,2])){
    patientid <- c(patientid, aggDiagnosis[i,1])
  }
}

##dataFrame with these convertors and denovo ids##
convDenovoDf <- subset(adnimerge, adnimerge$PTID %in% patientid)

##generalized function to extract demographic columns from a data frame##
func.demog.id <- function(input_data){
  ids <- grep(".*ID$", names(input_data), value=TRUE)
  if(length(ids)>1){
    id_idx1 <- which(colnames(input_data) %in% ids[1])
    id_idx2 <- which(colnames(input_data) %in% ids[2])
  }
  time_period <- c("M", "Month", "month", "Time", "time", "Years")
  time <-  colnames(input_data)[which(names(input_data) %in% time_period)]
  if(length(time)>1){
    time_idx <- which(colnames(input_data) %in% time[1])
  }
  visit_id <- c("VISCODE", "visit", "VISIT")
  visit <-  colnames(input_data)[which(names(input_data) %in% visit_id)]
  visit_idx <- which(colnames(input_data) %in% visit)
  prm_key <- c("primary_key", "primary_key1", "primary_key2")
  pk <- colnames(input_data)[which(names(input_data) %in% prm_key)]
  pk_id <- which(colnames(input_data) %in% pk)
  age <- c("Age", "age", "AGE")
  age_col <- colnames(input_data)[which(names(input_data) %in% age)]
  age_idx <- which(colnames(input_data) %in% age_col)
  diagnosis_bl <- c("DX.bl", "Diagnosis.bl", "Stage.bl")
  diag_bl <-  colnames(input_data)[which(names(input_data) %in% diagnosis_bl)]
  diag_bl_idx <- which(colnames(input_data) %in% diag_bl)
  diagnosis <- c("DX", "Diagnosis", "Stage")
  diag <-  colnames(input_data)[which(names(input_data) %in% diagnosis)]
  diag_idx <- which(colnames(input_data) %in% diag)
  education <- c("Education", "EDUCATION", "PTEDUCAT")
  edu <-  colnames(input_data)[which(names(input_data) %in% education)]
  edu_idx <- which(colnames(input_data) %in% edu)
  gender <- c("GENDER", "PTGENDER", "Gender") 
  gen <-  colnames(input_data)[which(names(input_data) %in% gender)]
  gen_idx <- which(colnames(input_data) %in% gen)
  ethnicity <- c("PTETHCAT", "ethinicity", "ethics")
  eth <-  colnames(input_data)[which(names(input_data) %in% ethnicity)]
  eth_idx <- which(colnames(input_data) %in% eth)
  race <- c("PTRACCAT", "RACE", "race", "Race")
  race_col <-  colnames(input_data)[which(names(input_data) %in% race)]
  race_idx <- which(colnames(input_data) %in% race_col)
  marital_status <- c("PTMARRY", "MARRIED", "Married", "married", "MARRY")
  marital <-  colnames(input_data)[which(names(input_data) %in% marital_status)]
  marital_idx <- which(colnames(input_data) %in% marital)
  demog_cols <- c(ids, time,visit, pk, age_col, diag_bl, diag, edu, gen, eth, race_col, marital)
  if(length(pk_id == 2)){
    pk_upd <- c(pk_id[1], pk_id[2])
  }
  else{
    pk_upd <- pk_id
  }
  demog_col_id <- c(id_idx1, id_idx2, time_idx,visit_idx,pk_upd, age_idx, 
                    diag_bl_idx, diag_idx, edu_idx, gen_idx, eth_idx, race_idx, marital_idx)
  return(list(demog_col_id, demog_cols))
}

demog_list <- func.demog.id(convDenovoDf)  
demog_cols <- unlist(demog_list[2])

##Vector for each type of biomarkers##
mri <- c("FDG", "AV45","Ventricles", "Hippocampus", "WholeBrain", "Entorhinal", "Fusiform", "MidTemp", "ICV")
ravlt <- grep("^RAVLT.*[^bl]$", colnames(adnimerge), value = TRUE)
ecog <- grep("^Ecog.*[^bl]$",  colnames(adnimerge), value = TRUE)
cognitionFunctional <- c("CDRSB", "MMSE", "MOCA", "ADAS11", "ADAS13", "FAQ")
csf <- c("ABETA", "TAU", "PTAU")
apoe <- c("APOE4")
mostRelevantCols <- c(mri, ravlt, ecog, cognitionFunctional, csf, apoe)

##subset the data frame according to the most relevant columns you need##
selectedFeatures <- convDenovoDf[,c(demog_cols, mostRelevantCols)]

# CSF biomarkers values are represented as characters due to presence of "<" or ">" symbols
#convert them to numeric
print("converting to Numeric")
selectedFeatures$ABETA <- gsub('(<|>)','',selectedFeatures$ABETA)
selectedFeatures$PTAU <- gsub('(<|>)','',selectedFeatures$PTAU)
selectedFeatures$TAU <- gsub('(<|>)','',selectedFeatures$TAU)

selectedFeatures$ABETA <- as.numeric(selectedFeatures$ABETA)
selectedFeatures$PTAU <- as.numeric(selectedFeatures$PTAU)
selectedFeatures$TAU <- as.numeric(selectedFeatures$TAU)

##Removing rows with no diagnostic information
rem.rid <- c()
for(i in 1:nrow(selectedFeatures)){
  if(is.na(selectedFeatures[i,"DX"]) == TRUE){
    if(sum(is.na(selectedFeatures[i,])) != (ncol(selectedFeatures)-10)){
      rem.rid <- c(rem.rid, i)
    }
  }
}
selectedFeatures <- selectedFeatures[-rem.rid,] 

##removing baseline columns as they are repititive##
blCols <- grep("bl", colnames(selectedFeatures), value = TRUE)
selectedFeatures[,blCols] <- NULL

##extracting demographic features##
demogs <- grep("^PT[^ID|AU].*", colnames(selectedFeatures), value = TRUE)
demogsAndOthers <- c(demogs, "AGE", "M", "PTID", "Month", "APOE4")

##converting the data frame to wider format##
##extracting the columns that need to be converted to wider format##
selectedCols <- setdiff(colnames(selectedFeatures), demogsAndOthers)
convDenovoSelectedDf <- selectedFeatures[,selectedCols]
reshapeDf = reshape(convDenovoSelectedDf,idvar='RID',timevar="VISCODE",dir='w')

# baseline features with less than 50% missing value
baselineFeatures <- reshapeDf[,grep("bl", colnames(reshapeDf), value = TRUE)]
baselineFeaturesUpd = baselineFeatures[ , -which(colMeans(is.na(baselineFeatures)) > 0.5)]

##removing the columns having more than 50% missing data from the main data frame##
remove50PercMissingCols <- setdiff(colnames(baselineFeatures),colnames(baselineFeaturesUpd))
toMatch <- (unlist(strsplit(remove50PercMissingCols, "\\.bl")))
matches <- unique(grep(paste(toMatch,collapse="|"), colnames(reshapeDf), value=TRUE))
reshapeDf[,matches] <- NULL

##csf data frame##
csfDf <- reshapeDf[,grep("ABETA|TAU|PTAU", colnames(reshapeDf), value = TRUE)]

##Volumetric data frame##
volumeDf <- reshapeDf[,grep("Entorhinal|Fusiform|MidTemp|Hippocampus|ICV|Ventricles|WholeBrain", colnames(reshapeDf),
                            value = TRUE)]

##Cognitive test data frame##
cogTestDf <- reshapeDf[,grep("CDRSB|ADAS11|ADAS13|MMSE|FAQ|RAVLT", colnames(reshapeDf), value = TRUE)]

##fdg##
fdgDf <- reshapeDf[,grep("FDG", colnames(reshapeDf), value = TRUE)]

##Diagnostic data frame##
dxDf <- reshapeDf[,grep("PTID|DX", colnames(reshapeDf), value = TRUE)]

# Remove columns with less than 50% missing value
csfDf = csfDf[ , -which(colMeans(is.na(csfDf)) > 0.5)]
volumeDf = volumeDf[ , -which(colMeans(is.na(volumeDf)) > 0.5)]
cogTestDf = cogTestDf[ , -which(colMeans(is.na(cogTestDf)) > 0.5)]
fdgDf = fdgDf[ , -which(colMeans(is.na(fdgDf)) > 0.5)]
dxDf = dxDf[ , -which(colMeans(is.na(dxDf)) > 0.5)]

##extract data related to above features##
demogsAndOthersDf <- convDenovoDf[rownames(dxDf),demogsAndOthers]


##Add bl extension to the demogsAndOthers features##
colnames(demogsAndOthersDf) <- paste(colnames(demogsAndOthersDf) , ".bl", sep = "")

##Combine all the data frames##
allFeatures <- cbind.data.frame(demogsAndOthersDf, csfDf, volumeDf, cogTestDf, fdgDf, dxDf)

##change the column name of fdgdf##
names(allFeatures)[names(allFeatures) =="fdgDf"] <- "FDG.bl"

## Add group name to columns
colnames(csfDf) = paste0("csf_",colnames(csfDf))
colnames(volumeDf) = paste0("brain_",colnames(volumeDf))
colnames(cogTestDf) = paste0("Cog_",colnames(cogTestDf))


# Auxiliary variables keep track of visit-wise and group-wise patient dropout. 
# Measurements of features are marked by value missing not at random (MNAR).
# MNAR  results from a systematic absence of subject data for a measurement type (feature group). 
get_aux_all_groups = function(cohortdata){
  mysample = cohortdata
  timepoint = str_extract(colnames(mysample), "m[0-9][0-9]")[1]
  #print(timepoint)
  csf = select(mysample,grep( "csf",colnames(mysample),value=TRUE)) 
  #print(csf)
  volumes = select(mysample,grep( "brain",colnames(mysample),value=TRUE)) 
  #print(volumes)
  cogTest = select(mysample,grep( "Cog",colnames(mysample),value=TRUE)) 
  #print(cogTest)
  snp = select(mysample,grep( "snp",colnames(mysample),value=TRUE)) 
  #select fdg
  fdg = select(mysample,grep( "fdg",colnames(mysample),value=TRUE)) 
  #path = select(mysample,grep( "path",colnames(mysample),value=TRUE)) 
  cort = select(mysample,grep( "cortical",colnames(mysample),value=TRUE)) 
  output_aux = function(mysubsample){
    
    #return_df = data.frame()
    #groupname = deparse(substitute(a))
    print("nsbsn")
    if(dim(mysubsample)[2] != 0 ){
      if(dim(mysubsample)[2]== 1){
        print("memem")
        mysubsample = mysubsample
      }else{
        print("jdjdcd")
        # Add a new column for AUX
        new = "new"
        in.loop = mysubsample
        mysubsample[new] <- 0
        
        # Get rownames where all value is NA 
        mysubsample_NA = which(apply(in.loop, 1, function(x) all(is.na(x))))
        mysubsample_pat = names(mysubsample_NA)
        print(mysubsample_pat)
        if(length(mysubsample_pat) !=0 ){
          mysubsample[which(rownames(mysubsample) %in% mysubsample_pat ),]$new <- 1
          print(mysubsample)
          # Annoate aux column with group name and visit number 
          groupName = sub("_.*$", "", colnames(mysubsample)[1])
          new2 = paste(groupName,"aux",timepoint, sep = "_")
          colnames(mysubsample)[which(names(mysubsample) == "new")] <- new2
          print(paste0("Aux available for ", groupName, "at" ,timepoint ))
        }else{ 
          mysubsample$new  = NULL
          groupName = sub("_.*$", "", colnames(mysubsample)[1])
          print(paste0("Aux unavailable for ", groupName, "at" ,timepoint ))}
      }
      
    }else{  print(paste0("Missing group at",timepoint ))}
    
    return(mysubsample)
  }
  
  csf_aux = output_aux(mysubsample = csf)
  volume_aux = output_aux(mysubsample = volumes)
  cogTest_aux = output_aux(mysubsample = cogTest)
  snp_aux = output_aux(mysubsample = snp)
  fdg_aux = output_aux(mysubsample = fdg)
  #path_aux = output_aux(mysubsample = path)
  cort_aux = output_aux(mysubsample = cort)
  #diagnostics_aux = output_aux(mysubsample = diagnostics)
  
  outputdf <- data.frame(matrix("removelater", ncol = 1, nrow = nrow(mysample)))
  names(outputdf)[1]<- "toremove"
  if(dim(csf)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , csf_aux))
  } else{
    print(paste0("CSF data unavailable for visit", timepoint))
  }
  if(dim(volumes)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , volume_aux))
  } else{
    print(paste0("Volumetric data unavailable for visit", timepoint))
  }
  if(dim(cogTest)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , cogTest_aux))
  } else{
    print(paste0("Cognition test unavailable for visit", timepoint))
  }
  if(dim(snp)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , snp_aux))
  } else{
    print(paste0("Snp unavailable for visit", timepoint))
  }
  if(dim(fdg)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , fdg_aux))
  } else{
    print(paste0("fdg unavailable for visit", timepoint))
  }
  # if(dim(path)[2] != 0 ){
  #   outputdf = as.data.frame(cbind(outputdf , path_aux))
  # } else{
  #   print(paste0("Pathway unavailable for visit", timepoint))
  # }
  if(dim(cort)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , cort_aux))
  } else{
    print(paste0("Cortical brain region unavailable for visit", timepoint))
  }
  outputdf$toremove = NULL
  return(outputdf)
}

##add PTID to reshapeDf##
allFeatures$PTID <- demogsAndOthersDf$PTID.bl
##add snp and pathway data ##
#load("~/dat_full.rda")
                                     
snps.dat.full = c(colnames(dat.full)[grep("(PTID|rs[0-9]+)", colnames(dat.full))])
pathways = c(colnames(dat.full)[grep("PTID|^.*Homo.sapiens|\\.[a-z]+$|\\.[A-Z][a-z]+$|^[A-Z][a-z]+$",colnames(dat.full))])
snp.pathway.data <- dat.full[,c(snps.dat.full, pathways)]
snp.pathway.data <- setDT(snp.pathway.data, keep.rownames = TRUE)[]
names(snp.pathway.data)[1] <- "PTID"
                                     
##load cortical brain volume baseline data##                                     
#EMC_ADNI_FS60_Phenotypes_Desikan_20180219 <- read_excel("~/EMC_ADNI_FS60_Phenotypes_Desikan_20180219.xlsx")
volumeRegion <- as.data.frame(EMC_ADNI_FS60_Phenotypes_Desikan_20180219)
names(volumeRegion)[1] <- "PTID"
brainRegion <- c(colnames(volumeRegion)[grep("Left|Right", colnames(volumeRegion))])
volumeRegion <- volumeRegion[,c("PTID",brainRegion)]
volumeRegion <- setDT(volumeRegion, keep.rownames = TRUE)[]

##snp data for 689 patients##
allFeatures <- merge(allFeatures, snp.pathway.data, all.x = TRUE)

##cortical desikan volume data for 689 patients##
allFeatures <- merge(allFeatures, volumeRegion , all.x = TRUE)

##snp data frame##
snpDf <- allFeatures[,snps.dat.full]


##Pathway data frame##
#pathDf <- reshapeDf[,pathways]

##cortical##
corticalDesikanDf <- allFeatures[,brainRegion]

colnames(snpDf) = paste0("snp_",colnames(snpDf), ".bl")
#colnames(pathDf) = paste0("path_",colnames(pathDf), ".bl")
colnames(corticalDesikanDf) = paste0("cortical_",colnames(corticalDesikanDf), ".bl")

##all visit Data##
allData = cbind.data.frame(csfDf, volumeDf, cogTestDf, snpDf, corticalDesikanDf, allFeatures$FDG.bl)
names(allData)[names(allData) == "allFeatures$FDG.bl"] <- "fdg_FDG.bl"

allData[-snp_idx] = lapply(allData[-snp_idx], as.numeric)
allData <- as.data.frame(allData)
snps <- grep("snp", colnames(allData), value = TRUE)
snp_idx <- grep("snp", colnames(allData))


allData <- cbind.data.frame(allData, dxDf)
                
visitData = list("visitbl" = allData[, grep(".bl", colnames(allData), value = TRUE)],
                 "visit1" = allData[, grep(".m06", colnames(allData), value = TRUE)],
                 "visit2" = allData[, grep(".m12", colnames(allData), value = TRUE)],
                 "visit3" = allData[, grep(".m24", colnames(allData), value = TRUE)])

colnames(visitData$visitbl) = sub("bl", "m00", colnames(visitData$visitbl)) 
                                     
# Create auxillary columns - group-wise and visit-wise
visitData_aux = sapply(visitData, get_aux_all_groups)

#Impute value (visit wise)
set.seed(123)
visitData_imputed = imputedData = sapply(visitData_aux, function(x)missForest::missForest(x, ntree = 500)[1])

set.seed(123)
##Make it as one data frame##
imputedData <- cbind.data.frame(demogsAndOthersDf, visitData_imputed$visitbl.ximp, visitData_imputed$visit1.ximp,
                            visitData_imputed$visit2.ximp, visitData_imputed$visit3.ximp, dxDf)
colClean <- function(x){ colnames(x) <- gsub("m00$", "bl", colnames(x)); x } 
imputedData <- colClean(imputedData)
colClean <- function(x){ colnames(x) <- gsub("m00.bl", "bl", colnames(x)); x } 
imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("csf_aux", "csf.aux", colnames(x)); x } 
imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("brain_aux", "brain.aux", colnames(x)); x } 
imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("Cog_aux", "Cog.aux", colnames(x)); x } 
imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("csf_", "", colnames(x)); x } 
imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("brain_", "", colnames(x)); x } 
imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("Cog_", "", colnames(x)); x } 
imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("snp_snp_", "", colnames(x)); x } 
imputedData <- colClean(imputedData)
colClean <- function(x){ colnames(x) <- gsub("snp_", "snp.", colnames(x)); x } 
imputedData <- colClean(imputedData) 
#colClean <- function(x){ colnames(x) <- gsub("path_", "", colnames(x)); x } 
#imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("cortical_aux", "cortical.aux", colnames(x)); x } 
imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("cortical_", "", colnames(x)); x } 
imputedData <- colClean(imputedData) 
colClean <- function(x){ colnames(x) <- gsub("_", ".", colnames(x)); x } 
imputedData <- colClean(imputedData) 
setnames(imputedData, old = grep("rs[0-9]+", colnames(imputedData), value= TRUE), new = gsub("\\.bl","", grep("rs[0-9]+", colnames(imputedData), value= TRUE)))
setnames(imputedData, old = grep("Left|Right", colnames(imputedData), value= TRUE), new = gsub("\\.bl","", grep("Left|Right", colnames(imputedData), value= TRUE)))
names(imputedData)[names(imputedData) == "fdg.FDG.bl"] <- "FDG.bl"
imputedData$M.bl <- NULL
visitData_OBB = sapply(visitData_aux, function(x) missForest::missForest(x, ntree = 100)[2])


visitData_complete <- cbind.data.frame(demogsAndOthersDf, visitData_aux$visitbl, visitData_aux$visit1,
                                       visitData_aux$visit2, visitData_aux$visit3, dxDf)

colClean <- function(x){ colnames(x) <- gsub("m00", "bl", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("m00.bl", "bl", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("csf_aux", "csf.aux", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("brain_aux", "brain.aux", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("Cog_aux", "Cog.aux", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("csf_", "", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("brain_", "", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("Cog_", "", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("snp_snp_", "", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete)
colClean <- function(x){ colnames(x) <- gsub("snp_", "snp.", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
#colClean <- function(x){ colnames(x) <- gsub("path_", "", colnames(x)); x } 
#visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("cortical_aux", "cortical.aux", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("cortical_", "", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 
colClean <- function(x){ colnames(x) <- gsub("_", ".", colnames(x)); x } 
visitData_complete <- colClean(visitData_complete) 

library(data.table)
setnames(visitData_complete, old = grep("rs[0-9]+", colnames(visitData_complete), value= TRUE), new = gsub("\\.bl.bl","", grep("rs[0-9]+", colnames(visitData_complete), value= TRUE)))
setnames(visitData_complete, old = grep("Left|Right", colnames(visitData_complete), value= TRUE), new = gsub("\\.bl","", grep("Left|Right", colnames(visitData_complete), value= TRUE)))
names(visitData_complete)[names(visitData_complete) == "fdg.FDG.bl"] <- "FDG.bl"
names(visitData_complete)[names(visitData_complete) == "cortical.aux.bl"] <- "brain68.aux"

visitData_complete$M.bl <- NULL
visitData_complete <- as.data.frame(visitData_complete)

imputedData[,grep("DX", colnames(imputedData), value = TRUE)] <- lapply(imputedData[,grep("DX", colnames(imputedData), value = TRUE)], factor)

levels(imputedData$DX.bl) <- c(levels(imputedData$DX.bl), "Unknown")
levels(imputedData$DX.m06) <- c(levels(imputedData$DX.m06), "Unknown")
levels(imputedData$DX.m12) <- c(levels(imputedData$DX.m12), "Unknown")
levels(imputedData$DX.m24) <- c(levels(imputedData$DX.m24), "Unknown")
imputedData[is.na(imputedData)] <- "Unknown"
names(imputedData)[names(imputedData) == "cortical.aux.bl"] <- "brain68.aux"


write.csv(visitData_complete, file = "visitData_complete.csv")
saveRDS(visitData_complete, file = "visitData_complete.rds")
dxInVisitData <- grep("DX", colnames(visitData_complete), value = TRUE)
write.csv(imputedData, file = "imputedData.csv")
saveRDS(imputedData, file = "imputedData.rds")

save.image("~/preprocessing.RData")
