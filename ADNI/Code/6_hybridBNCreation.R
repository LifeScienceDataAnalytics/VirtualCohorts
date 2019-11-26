##Script name: "6_hybridBNCreation.R"
##Purpose of Script: generation of hybrid bayesian network and comparison with the one generated using dicrete data
##Author: Meemansa Sood
##Date Created: October 2019

##load the data (mixedMeta_rf) before discretization from AutoencodingAndBNCreation.Rmd##
library(bnlearn)
library(h2o)
library(randomForest)
library(parallel)
library(doMC)
mixedMeta_rf <- meta_rf
names(mixedMeta_rf)[12] <- "AGE.bl"
names(mixedMeta_rf)[13] <- "FDG.bl"

mixed_meta_rf_merge <- cbind(mixedMeta_rf, input_df[,diagnosis], input_df[,unique(c(demogs,auxBl,auxm06, auxm12, auxm24))])

set.seed(1234)
cvres1Gn = bn.cv(mixed_meta_rf_merge, "rsmax2", runs=10,  loss="logl-cg",algorithm.args = list(blacklist=blNoPath, whitelist=wlNoPath), cluster= clGn)
print(cvres1Gn)

cvres2Gn = bn.cv(mixed_meta_rf_merge, "mmhc", runs=10,  loss="logl-cg", algorithm.args = list(blacklist=blNoPath,  whitelist=wlNoPath), cluster=clGn)
print(cvres2Gn)

cvres3Gn = bn.cv(mixed_meta_rf_merge, "hc", runs=10, loss="logl-cg", algorithm.args = list(maxp=5, blacklist=blNoPath,  whitelist=wlNoPath, restart=10, score="bic-cg", debug = TRUE), cluster=clGn)

cvres4Gn = bn.cv(mixed_meta_rf_merge, "tabu", runs=10, loss="logl-cg", algorithm.args = list(maxp=5, blacklist=blNoPath,  whitelist=wlNoPath, restart=10, score="bic-cg"), cluster=clGn)
print(cvres4Gn)

cvres5Gn = bn.cv(mixed_meta_rf_merge, "si.hiton.pc", runs=10, loss="logl-cg", algorithm.args = list(blacklist=blNoPath,  whitelist=wlNoPath, undirected=FALSE),cluster=clGn)
print(cvres5Gn)

cvres6Gn = bn.cv(mixed_meta_rf_merge, "mmpc", runs=10, loss="logl-cg", algorithm.args = list(blacklist=blNoPath, whitelist=wlNoPath, undirected=FALSE), cluster=clGn)
print(cvres6Gn)


plot(cvres1Gn, cvres2Gn, cvres3Gn, cvres4Gn, xlab=c("rsmax2", "mmhc", "hc", "tabu"))

##Constrained based algorithms are inferior to score based algorithms, mmhc has very low loss as compared to
##other score based algorithms, considering all these issues we select tabu score based algorithm
rm( list = setdiff(ls() , c("mixed_meta_rf_merge", "blNoPath", "wlNoPath" ) ))
clGn =makeCluster(10)
mix.boot.stren = boot.strength(mixed_meta_rf_merge, algorithm="tabu", R=1000, algorithm.args = list(maxp=5, blacklist=blNoPath, whitelist=wlNoPath, restart=50, score="bic-cg"), cluster=clGn)
mixedBN = tabu(mixed_meta_rf_merge, maxp=5, blacklist=blNoPath,  whitelist=wlNoPath, score="bic-cg")
mix.finalBNarcs <- as.data.frame(mixedBN$arcs)
stopCluster(clGn)
thresh = 0.1
bn.av = averaged.network(mix.boot.stren, threshold = thresh)

boot.strenwithThreshold = mix.boot.stren[mix.boot.stren$strength >= 0.5 & mix.boot.stren$direction >= 0.5, ]
############################Non-Conservative#######################################################
## Simulate virtual patients
simulated = rbn(mixedBN, n = nrow(mixed_meta_rf_merge), mixed_meta_rf_merge, debug = FALSE)
mixedSIM = tabu(simulated, blacklist=blNoPath, whitelist=wlNoPath, score="bic-cg")
adni_sim_bayes = bnlearn::BF(mixedBN,mixedSIM, mixed_meta_rf_merge)

logLikADNI <- logLik(mixedBN,mixed_meta_rf_merge)
logLikSIM <- logLik(mixedSIM,mixed_meta_rf_merge)

# #remove DX and  aux
to.remove <- grep("aux",colnames(mixed_meta_rf_merge),value=TRUE)


realCp <- mixed_meta_rf_merge
# real : real data,  virtual:  simulated data
real <- realCp[,-(which(colnames(realCp) %in% to.remove))]
virtual <- simulated[,-(which(colnames(simulated) %in% to.remove))]


library(randomForestSRC)
library(pROC)

RealPep = real
VirtPep = virtual

#Add category variable to each dataframe
VirtPep$typeOfPatient = "VP"
RealPep$typeOfPatient = "RP"
allPep = rbind(VirtPep, RealPep)


# Chnage outcome to factor
allPep$typeOfPatient = as.factor(allPep$typeOfPatient)
levels(allPep$typeOfPatient)[levels(allPep$typeOfPatient)=="RP"] <- "1"
levels(allPep$typeOfPatient)[levels(allPep$typeOfPatient)=="VP"] <- "0"
rf_auc_df = data.frame("auc" = 100)

library(caret)
# Data split
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
  print(rf_auc)
  if(is.na(rf_auc) == FALSE){
    rf_auc_df = cbind(rf_auc_df ,rf_auc)
    #Rename columns
    last_index = ncol(rf_auc_df)
    fold_name = paste("Fold" , as.character(i))
    colnames(rf_auc_df)[last_index] <- fold_name
  }
}

rf_auc_df_not_iter <- rf_auc_df
rf_auc_df_not_iter[1,1] = 0
rf_auc_df_not_iter$auc = NULL

rf_auc_df_mean =  mean(unlist(rf_auc_df_not_iter) )

bx2 = boxplot(unlist(rf_auc_df_not_iter) , main = "AUC" , col = c("ivory4"))
boxplot(unlist(rf_auc_df_not_iter), col = "#00BFC4" , main = "Partial AUC of cross-validation of classifier (virtual vs real patients)")
save.image("/mixedBNADNI.RData")



############################Conservative#######################################################
library(gdata)
real <- mixed_meta_rf_merge 

simulate_VPs = function(real, n_VP=NROW(real), iterative=TRUE){
  # n_VP : number of virtual patients to simulate
  # estimate the structure and parameters of the Bayesian network.
  res = tabu(real, maxp=5, blacklist=blNoPath,whitelist=wlNoPath,  score="bic-cg") # assuming tabu was the best structure learning approach (currently it seems like for PPMI hc is better!)
  fitted = bn.fit(res,real)
  VP = c()
  iter = 1
  while(NROW(VP) <= n_VP){
    print(NROW(VP))
    cat("iteration = ", iter, "\n")
    generatedDF = rbn(fitted, n = n_VP) # draw virtual patients
    if(iterative){
      y = factor(c(rep("original", NROW(real)), rep("generated", n_VP)))
      df = data.frame(y=y, x=rbind(real, generatedDF))
      vec <- c(rep(1,sum(y=="original")), rep(0.2*sum(y=="original")/sum(y=="generated"), sum(y=="generated")))
      print(vec)
      print(typeof(vec))
      fit = rfsrc(y ~ ., data=df, case.wt=vec)
    
      print(c(rep(1,sum(y=="original")), rep(0.2*sum(y=="original")/sum(y=="generated"), sum(y=="generated"))))
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


set.seed(123)
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

