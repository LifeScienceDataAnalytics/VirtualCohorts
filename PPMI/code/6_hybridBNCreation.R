##Script name: 6_mixedBNppmi.R
##Purpose of Script: generation of hybrid bayesian network and comparison with the one generated using discrete data
##Author: Meemansa Sood
##Date Created:October 2018

##load allData_meta from 1_getData.R
load("~/allData_meta.RData")
meta_ppmi <- allData_meta

names(meta_ppmi)[names(meta_ppmi) == 'Patient_Demographic_V00'] <- 'Patient_Demographic'
names(meta_ppmi)[names(meta_ppmi) == 'Patient_PDhistory_V00'] <- 'Patient_PDhistory'
load("~/Documents/PhDWork/BayesianNetworkAD/CodeAndResultsforPaper/PPMI/Workspace/stable_network.RData")
rm( list = setdiff(ls() , c("meta_ppmi", "dt_bl", "dt_wl" ) ))
library(randomForest)
library(bnlearn)

library(gdata)
real <- meta_ppmi
realCopy = real
simulate_VPs = function(real, n_VP=NROW(real), iterative=TRUE){
  # n_VP : number of virtual patients to simulate
  # estimate the structure and parameters of the Bayesian network.
  res = tabu(real, maxp=5, blacklist=dt_bl,whitelist=dt_wl,  score="bic-cg") # assuming tabu was the best structure learning approach (currently it seems like for PPMI hc is better!)
  fitted = bn.fit(res,real)
  VP = c()
  iter = 1
  while(NROW(VP) <= n_VP){
    print(NROW(VP))
    cat("iteration = ", iter, "\n")
    generatedDF = rbn(fitted, n = n_VP) # draw virtual patients
    if(iterative){
      #y = factor(c(rep("original", NROW(real)), rep("generated", n_VP)))
      y = factor(c(rep("original", NROW(realCopy)), rep("generated", NROW(generatedDF) + NROW(real) - NROW(realCopy))))
      df = data.frame(y=y, x=rbind(real, generatedDF))
      fit = rfsrc(y ~ ., data=df, case.wt=c(rep(1,sum(y=="original")), rep(0.2*sum(y=="original")/sum(y=="generated"), sum(y=="generated"))))
      print(fit)
      
      DGz = predict(fit)$predicted[(NROW(df) - NROW(generatedDF) + 1):NROW(df), "original"]
      DGz = (DGz >0.5)*1
      acceptedVPs = generatedDF[DGz == 1,]
    }
    else
      acceptedVPs = generatedDF
    VP = rbind.data.frame(VP, acceptedVPs)
    iter = iter + 1
  }
  return(VP)
}

set.seed(51)
simulatedVP <- simulate_VPs(real)
simulated <- simulatedVP[1:nrow(real),]
#simulated <- simulatedVP[(nrow(simulatedVP)-nrow(real)+1):nrow(simulatedVP),]
intersect(simulated, real)
simulated[] <- lapply(simulated, factor)

#remove DX and  aux 
to.remove <- grep("aux",colnames(real),value=TRUE)
realCp <- real

original <- realCp[,-(which(colnames(realCp) %in% to.remove))]
virtual <- simulated[,-(which(colnames(simulated) %in% to.remove))]

intersect(virtual, original)

RealPep = real
VirtPep = simulated 

rownames(RealPep) = paste0("R", 1:nrow(RealPep))
rownames(VirtPep) = paste0("V", 1:nrow(VirtPep))

#Add category variable to each dataframe 
VirtPep$typeOfPatient = "VP"
RealPep$typeOfPatient = "RP"
allPep = rbind(VirtPep, RealPep)

# Chnage outcome to factor
allPep$typeOfPatient = as.factor(allPep$typeOfPatient)
levels(allPep$typeOfPatient)[levels(allPep$typeOfPatient)=="RP"] <- "1"
levels(allPep$typeOfPatient)[levels(allPep$typeOfPatient)=="VP"] <- "0"

library(caret)
# Data split 
rf_pred_df <- setNames(data.frame(matrix(ncol=4, nrow=0)), c("Rownames", "1", "0", "fold"))
selected_index3 = createMultiFolds(allPep[,"typeOfPatient"], k=10, times = 10) 
rf_auc_df = data.frame("auc" = 100)
rf_feature_list = list()
for(i in 1:length(selected_index3)){
  
  f_index = unlist(selected_index3[[i]])
  rf_train = allPep[f_index,]
  rf_test = allPep[-f_index,]
  
  print(i)
  
  
  #Model fit
  rf_model = rfsrc(typeOfPatient ~. , data = rf_train , ntree = 500, na.action = "na.impute", forest = TRUE , importance = TRUE)
  
  
  # Prediction 
  rf_pred <- predict(rf_model, newdata = rf_test , na.action = "na.impute", type ="prob")
  #print("memem")
  rf_pred_prob = as.data.frame(rf_pred$predicted)
  #print(rf_pred_prob)
  rf_pred_prob$Rownames <- rownames(rf_test)
  rf_pred_prob$fold <- rep(i, nrow(rf_test))
  rf_pred_df <- rbind.data.frame(rf_pred_df,rf_pred_prob)
  # Calcualte AUC
  rf_auc= auc(as.numeric(as.character(rf_test$typeOfPatient)), rf_pred_prob$`1`, partial.auc.focus="sens",partial.auc=c(0.9,1), partial.auc.correct=TRUE)
  #print(rf_auc)
  #rf_auc= auc(as.numeric(as.character(rf_test$typeOfPatient)), rf_pred_prob$`1`, partial.auc.focus="sens",partial.auc=c(0.9,1), partial.auc.correct=TRUE)
  
  if(is.na(rf_auc) == FALSE){
    rf_auc_df = cbind(rf_auc_df ,rf_auc)
    #Rename columns
    last_index = ncol(rf_auc_df)
    fold_name = paste("Fold" , as.character(i))
    colnames(rf_auc_df)[last_index] <- fold_name
  }
  # #Save features
  # fold_feature =   rf_model$importance
  # rf_feature_list = append(rf_feature_list , list(fold_feature),0)
  
}

rf_auc_df[1,1] = 0
rf_auc_df$auc = NULL

rf_auc_df_mean =  mean(unlist(rf_auc_df)) 
bx = boxplot(unlist(rf_auc_df) , main = "AUC" , col = c("ivory4"))
rf_auc_df_rec <- rf_auc_df
boxplot(unlist(rf_auc_df_rec), col = "#00BFC4" , main = "AUC of cross-validation of classifier (virtual vs real patients)")

meta_ppmi <- allData_meta

names(meta_ppmi)[names(meta_ppmi) == 'Patient_Demographic_V00'] <- 'Patient_Demographic'
names(meta_ppmi)[names(meta_ppmi) == 'Patient_PDhistory_V00'] <- 'Patient_PDhistory'
rm( list = setdiff( ls() , c("meta_ppmi", "dt_wl", "dt_bl")))
finalBN = tabu(meta_ppmi, maxp=5, blacklist=dt_bl,  whitelist=dt_wl, restart=50, score="bic-cg")
library(randomForestSRC)
library(pROC)
#load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/bitcluster/ClusteringandBN/stable_network_cluster.RData")

## Simulate virtual patientsm
simulated = rbn(finalBN, n = 362, meta_ppmi, debug = FALSE)
finalSIM = tabu(simulated, maxp=5, blacklist=dt_bl,whitelist=dt_wl, score="bic-cg") 


#remove DX and  aux 
to.remove <- grep("aux|Status|DX",colnames(meta_ppmi),value=TRUE)

# real : real data,  virtual:  simulated data
real <- meta_ppmi[,-(which(colnames(meta_ppmi) %in% to.remove))]
virtual <- simulated[,-(which(colnames(simulated) %in% to.remove))]


# Save result 
rf_auc_df = data.frame("auc" = 100)
test_d = setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("D","D.str","M1","M2"))

# Copy data
RealPep = real 
VirtPep = virtual
rownames(RealPep) = paste0("R", 1:nrow(RealPep))
rownames(VirtPep) = paste0("V", 1:nrow(VirtPep))

#Add category variable to each dataframe 
RealPep$typeOfPatient = "RP"
VirtPep$typeOfPatient = "VP"
allPep = rbind(RealPep, VirtPep)

# Chnage outcome to factor
allPep$typeOfPatient = as.factor(allPep$typeOfPatient)
levels(allPep$typeOfPatient)[levels(allPep$typeOfPatient)=="RP"] <- "1"
levels(allPep$typeOfPatient)[levels(allPep$typeOfPatient)=="VP"] <- "0"

selected_index3 = createMultiFolds(allPep[,"typeOfPatient"], k=10, times = 10) 

for(i in 1:length(selected_index3)){
  
  f_index = unlist(selected_index3[[i]])
  rf_train = allPep[f_index,]
  rf_test = allPep[-f_index,]
  
  print(i)
  
  
  #Model fit
  rf_model = rfsrc(typeOfPatient ~. , data = rf_train , ntree = 500, na.action = "na.impute", forest = TRUE , importance = TRUE)
  
  
  # Prediction 
  rf_pred <- predict(rf_model, newdata = rf_test , na.action = "na.impute", type ="prob")
  rf_pred_prob = as.data.frame(rf_pred$predicted)
  
  # Calcualte AUC
  rf_auc= auc(as.numeric(as.character(rf_test$typeOfPatient)), rf_pred_prob$`1`, partial.auc.focus="sens",partial.auc=c(0.9,1), partial.auc.correct=TRUE)
  rf_auc_df = cbind(rf_auc_df ,rf_auc)
  
  #Rename columns
  last_index = ncol(rf_auc_df)
  fold_name = paste("Fold" , as.character(i))
  colnames(rf_auc_df)[last_index] <- fold_name
  
  # #Save date for ROC curve
  D.ex <- as.numeric(as.character(rf_test$typeOfPatient))
  M1 <- rf_pred_prob$`1`
  M2 <- rf_pred_prob$`0`
  test_d_loop <- data.frame(D = D.ex, D.str = c("Virtual", "Real")[D.ex + 1],M1 = M1, M2 = M2, stringsAsFactors = FALSE ) 
  test_d = rbind(test_d,test_d_loop)
}

rf_auc_df[1,1] = 0
rf_auc_df$auc = NULL

rf_auc_df_mean =  mean(unlist(rf_auc_df) )   # Mean auc  0.6340365 !!!
bx = boxplot(unlist(rf_auc_df) , main = "AUC" , col = c("ivory4"))
bx = boxplot(unlist(rf_auc_df), col = "#00BFC4" , main = "AUC of cross-validation of classifier (virtual vs real patients) " )

intersect(simulated, real)
save.image("~/mixedBNPPMI.RData")
