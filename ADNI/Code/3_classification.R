##Script name: "3_classification.R"
##Purpose of Script: classifier to classify virtual and real patients using coonservative and non-conservative methods
##Author: Meemansa Sood
##Date Created: October 2018

##Partial AUC based on real patients and all virtual patients#
#simulate patients
library(bnlearn)
library(cluster)
library(dplyr)
library(gdata)
library(factoextra)
library(FactoMineR)

##load the saved workspace from 2_AutoencodingAndBNCreation.Rmd
load("/bnCreation.RData")

###################################Non-Conservative approach to simulate virtual patients###################################
## Simulate virtual patients 
simulated = rbn(finalBN, n = nrow(DiscAutoEnADNI), DiscAutoEnADNI, fit = "bayes", debug = FALSE)
simulated[] <- lapply(simulated, factor)


finalSIM = tabu(simulated, maxp=5, blacklist=blNoPath, whitelist=wlNoPath, score="bic") 
adni_sim_bayes = bnlearn::BF(finalBN,finalSIM, DiscAutoEnADNI)

logLikADNI <- logLik(finalBN,DiscAutoEnADNI)
logLikSIM <- logLik(finalSIM,simulated) 

#remove DX and  aux 
to.remove <- grep("aux",colnames(DiscAutoEnADNI),value=TRUE)


realCp <- DiscAutoEnADNI
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


# change outcome to factor
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

rf_auc_df[1,1] = 0
rf_auc_df$auc = NULL

rf_auc_df_non_cons <- rf_auc_df
rf_auc_df_mean =  mean(unlist(rf_auc_df) )   # Mean auc  0.568!!!

bx2 = boxplot(unlist(rf_auc_df) , main = "AUC" , col = c("ivory4"))
boxplot(unlist(rf_auc_df), col = "#00BFC4" , main = "Partial AUC of cross-validation of classifier (virtual vs real patients)")
save.image("~/non-conservative.RData")

###################################Conservative approach to simulate virtual patients###################################
real <- DiscAutoEnADNI

simulate_VPs = function(real, n_VP=NROW(real), iterative=TRUE){
  # n_VP : number of virtual patients to simulate
  # estimate the structure and parameters of the Bayesian network.
  res = tabu(real, maxp=5, blacklist=blNoPath,whitelist=wlNoPath,  score="bic") # assuming tabu was the best structure learning approach (currently it seems like for PPMI hc is better!)
  fitted = bn.fit(res,real, method="bayes")
  VP = c()
  iter = 1
  while(NROW(VP) <= n_VP){
    print(NROW(VP))
    cat("iteration = ", iter, "\n")
    generatedDF = rbn(fitted, n = n_VP) # draw virtual patients
    if(iterative){
      y = factor(c(rep("original", NROW(real)), rep("generated", n_VP)))
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
rf_auc_df_cons <- rf_auc_df
boxplot(unlist(rf_auc_df_rec), col = "#00BFC4" , main = "AUC of cross-validation of classifier (virtual vs real patients)")

save.image("~/conservative.RData")

#png(file = "/output/aucCompared.png",
#    width = 12, height = 8, units = 'in', res = 300)
par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
bx2 <- boxplot(unlist(rf_auc_df_non_cons), col = "#00BFC4", main = "Virtual patients generated by non-conservative method") + theme(axis.text = element_text(size=10),
                                                                                                                        axis.text.x = element_text(size = 10),
                                                                                                                        axis.title = element_text(size = 10, face = "bold"),
                                                                                                                        strip.text = element_text(size = 10))
bx <- boxplot(unlist(rf_auc_df_cons), col = "#00BFC4", main = "Virtual patients generated by conservative method") + theme(axis.text = element_text(size=10),
                                                                                                                       axis.text.x = element_text(size = 10),
                                                                                                                       axis.title = element_text(size = 10, face = "bold"),
                                                                                                                       strip.text = element_text(size = 10))
mtext("Partial AUC of cross-validation of classifier (virtual vs real patients)", outer = TRUE, cex = 1.5)
#dev.off()


######MCA plots for real and virtual patients (use data generated from conservative approach)##

##probabilities of virtual patient##
rf_pred_df_virtual <- subset(rf_pred_df, rf_pred_df$Rownames %in% (grep("V", rf_pred_df$Rownames,value = TRUE )))
rf_pred_df_virtual <- rf_pred_df_virtual[,c(1,3,4)]
rf_pred_df_virtual_agg <- aggregate(. ~Rownames, data=rf_pred_df_virtual, mean, na.rm=TRUE)
names(rf_pred_df_virtual_agg)[names(rf_pred_df_virtual_agg) == "1"] <- "Mean_Probability"
rf_pred_df_virtual_agg$fold <- NULL


rf_pred_df_real <- subset(rf_pred_df, rf_pred_df$Rownames %in% (grep("R", rf_pred_df$Rownames,value = TRUE )))
rf_pred_df_real <- rf_pred_df_real[,c(1,3,4)]
rf_pred_df_real_agg <- aggregate(. ~Rownames, data=rf_pred_df_real, mean, na.rm=TRUE)
names(rf_pred_df_real_agg)[names(rf_pred_df_real_agg) == "1"] <- "Mean_Probability"
rf_pred_df_real_agg$fold <- NULL
rf_pred_df_real_agg$Mean_Probability <- 100

rf_pred_all <- rbind.data.frame(rf_pred_df_real_agg, rf_pred_df_virtual_agg)
rownames(rf_pred_all) <- rf_pred_all[,1]
rf_pred_all$Rownames <- NULL

real_copy <- RealPep
real_copy$typeOfPatient <- NULL
vir_copy <- VirtPep
vir_copy$typeOfPatient <- NULL
real_copy$TypeOfSubject = "Real"
vir_copy$TypeOfSubject = "Virtual"

all_subject = rbind.data.frame(real_copy, vir_copy)
all_subject_prob <- merge(all_subject, rf_pred_all, by=0, all = TRUE)

rownames(all_subject_prob) <- all_subject_prob[,1]
all_subject_prob$Row.names <- NULL
#all_subject$TypeOfSubject = as.factor(all_subject$TypeOfSubject)
library("FactoMineR")
library("factoextra")
all_subject_prob_cp <- all_subject_prob
all_subject_prob$Mean_Probability <- NULL


res.mca =  MCA(all_subject_prob, graph = FALSE)
##to calculate the variance explained##
fviz_mca_biplot(res.mca,
                repel = FALSE, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal(),habillage = allPep$typeOfPatient,palette = c("#00AFBB", "#E7B800","#ff0000"), invisible=c("var"))
mca1_obs_df = data.frame(res.mca$ind$coord)
mca1_obs_df_two <- mca1_obs_df[,1:2]

names(all_subject_prob_cp)[names(all_subject_prob_cp)== "Mean_Probability"] <- "Probability"
mca_obs_merge <- merge(x=mca1_obs_df_two, y=all_subject_prob_cp[,"Probability", drop=FALSE], by.x = "row.names", by.y = "row.names")

##divide real and virtual patient##
mca_obs_merge_rp <- subset(mca_obs_merge, mca_obs_merge$Row.names %in% grep("R", mca_obs_merge$Row.names, value = TRUE))
mca_obs_merge_rp$Category = "Real"
mca_obs_merge_rp$Category = as.factor(mca_obs_merge_rp$Category)
mca_obs_merge_vp <- subset(mca_obs_merge, mca_obs_merge$Row.names %in% grep("V", mca_obs_merge$Row.names, value = TRUE))
mca_obs_merge_vp$Category = "Virtual"
mca_obs_merge_vp$Category = as.factor(mca_obs_merge_vp$Category)

#png(file = "/output/allMCAupdNov.png",
#    width = 8, height = 6, units = 'in', res = 300)

ggplot(mca_obs_merge_rp, aes(x = Dim.1, y = Dim.2)) +   xlab("Dim1 (9.8%)") + ylab("Dim2 (5.5%)") +
  geom_hline(yintercept = 0, colour = "black")+
  geom_vline(xintercept = 0, colour = "black")+
  geom_point(colour= "#619CFF", aes(shape=Category)) +  
  geom_point(data = mca_obs_merge_vp, aes(x=Dim.1, y=Dim.2, colour=Probability, shape = Category)) + 
  scale_colour_gradient(low = "Yellow", high = "Red") + theme_classic()

#dev.off()

save.image("/virtPartients.RData")
