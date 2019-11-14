load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/stable_network.RData")
rm( list = setdiff(ls() , c("disc_meta", "dt_bl", "dt_wl","dics_data" ) ))
library(randomForest)
library(bnlearn)

################# Use random forest to classify real and virtual patients  ######################################

real = dics_data   # x
realCopy = real  # xorig
rownames(realCopy) = paste0("R",rownames(realCopy))
converged = FALSE
thresh = 0.8
iter = 1

while(!converged & iter <= 100){
  cat("iteration = ", iter, "\n")
  # learn the network structure - toy example.
  res = tabu(real, maxp=5, blacklist=dt_bl,whitelist=dt_wl,  score="bic")
  # estimate the parameters of the Bayesian network.
  #fitted = bn.fit(res,real)
  
  generatedDF = rbn(res, n = nrow(realCopy), real, fit = "bayes", debug = FALSE) # my data
  
  y = factor(c(rep("original", NROW(realCopy)), rep("generated", NROW(generatedDF) + NROW(real) - NROW(realCopy))))
  df = data.frame(y=y, x=rbind(real, generatedDF))
  #fit = randomForest(y ~ ., data=df, classwt=c("original" = 1000, "generated" = 1), ntree = 500)
  #fit = tuneRF(df,y, stepFactor=1.5, improve=1e-5, ntree=500,doBest = TRUE)
  #fit <- rfsrc(y ~ ., data = df, case.wt = c(rep(1,sum(y =="original")), rep(0.2, sum(y== "generated")))   )
  fit = rfsrc(y ~ ., data=df, case.wt=c(rep(1,sum(y=="original")), rep(0.2*sum(y=="original")/sum(y=="generated"), sum(y=="generated"))))
  print(fit)
  
  DGz = predict(fit)$predicted[(NROW(df) - NROW(generatedDF) + 1):NROW(df), "original"]
  DGz = (DGz >0.5)*1
  converged = mean(DGz == 1) >= thresh
  real = rbind.data.frame(real, generatedDF[DGz == 1,]) 
  print(converged)
  iter = iter + 1
}


################# Use random forest to classify real and virtual patients  ######################################

library(randomForestSRC)
library(pROC)
rf_auc_df = data.frame("auc" = 100)
rf_feature_list = list()
test_d = setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("D","D.str","M1","M2"))
# 
RealPep = real[1:362,] 
VirtPep = real[362:nrow(real),]  

# RealPep = as.data.frame(sapply(real, as.numeric))
# VirtPep  =  as.data.frame(sapply(virtual , as.numeric))


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
#rf_auc_df$auc = NULL

rf_auc_df_mean =  mean(unlist(rf_auc_df),na.rm = TRUE)   # Mean auc  0.6340365 !!!
bx = boxplot(unlist(rf_auc_df) , main = "AUC" , col = c("ivory4"))
bx = boxplot(unlist(rf_auc_df), col = "#00BFC4" , main = "AUC of cross-validation of classifier (virtual vs real patients) " )

library(plotROC)
set.seed(2529)
aucplot <- ggplot(test_d, aes(d = D, m = M1)) + geom_roc() + style_roc()

# workspace : 3_VirtualPatient_120.RData


# MCA

levels(allPep$typeOfPatient)[levels(allPep$typeOfPatient)=="1"] <- "Real"
levels(allPep$typeOfPatient)[levels(allPep$typeOfPatient)=="0"] <- "Simulated"

res.mca =  MCA(allPep, graph = FALSE)
mca_plot = fviz_mca_biplot(res.mca, 
                repel = FALSE, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal(),habillage = allPep$typeOfPatient,palette = c("#00AFBB", "#E7B800","#ff0000"), invisible=c("var"))

plot.MCA(res.mca,invisible=c("var"))

#https://cran.r-project.org/web/packages/plotROC/vignettes/examples.html      
################################################################################################################
################# Simulate using rbn() with adding more patients ######################################
load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/stable_network.RData")
rm( list = setdiff( ls() , c("disc_meta", "dics_data" , "disc_bnlearn", "get_bl_wl","finalBN" , "discPCA2", "dt_wl", "dt_wl" ,"bnoutput_mean")))
library(randomForestSRC)
library(pROC)
#load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/bitcluster/ClusteringandBN/stable_network_cluster.RData")

## Simulate virtual patients 
simulated = rbn(finalBN, n = 362, discPCA2, fit = "bayes", debug = FALSE)
finalSIM = tabu(simulated, maxp=5, blacklist=dt_wl,whitelist=dt_wl, score="bic") 


#remove DX and  aux 
to.remove <- grep("aux|Status|DX",colnames(discPCA2),value=TRUE)

# real : real data,  virtual:  simulated data
real <- discPCA2[,-(which(colnames(discPCA2) %in% to.remove))]
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

library(plotROC)
set.seed(2529)
aucplot <- ggplot(test_d, aes(d = D, m = M1)) + geom_roc() + style_roc()

# workspace 
#save.image("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/3_VirtualPatient_211.RData")
#========= MCA plot =======================

# Apply MCA

real_copy = RealPep
vir_copy = VirtPep
real_copy$typeOfPatient = "Real"
vir_copy$typeOfPatient = "Simulated"

all_subject = as.data.frame(rbind(real_copy,vir_copy))
all_subject$typeOfPatient = as.factor(all_subject$typeOfPatient)

res.mca =  MCA(all_subject, graph = FALSE)
# visualize plots to see outliers 
#The function fviz_mca_biplot() [factoextra package] is used to draw the biplot of individuals and variable categories:

fviz_mca_biplot(res.mca, 
                repel = TRUE, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal())


fviz_mca_biplot(res.mca, 
                repel = FALSE, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal(),habillage = all_subject$typeOfPatient,palette = c("#00AFBB", "#E7B800"))

mca_plot = fviz_mca_biplot(res.mca, 
                repel = FALSE, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal(),habillage = all_subject$typeOfPatient,palette = c("#00AFBB", "#E7B800","#ff0000"), invisible=c("var"))

plot.MCA(res.mca,invisible=c("var"))

################################################################################################################
#                           Plot density of real data and simulated data
#################################################################################################################
library(reshape2)
# Generative method
#load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/3_VirtualPatient_120.RData")
load("~/Documents/PhDWork/BayesianNetworkAD/CodeAndResultsforPaper/PPMI/Workspace/3_VirtualPatient_120.RData")
rm( list = setdiff(ls() , c("dics_data", "dt_bl", "dt_wl","real" ) ))
xorig = real[1:362,]
xsimulated = real[362:nrow(real),]

# remove aux
remove_row = grep("aux", colnames(xorig) , value = TRUE)
xorig <- xorig[, ! colnames(xorig) %in% remove_row ]
##(Meemansa)for bar plot
xorigBar <- xorig
remove_row2 = grep("aux", colnames(xsimulated) , value = TRUE)
xsimulated <- xsimulated[, ! colnames(xsimulated) %in% remove_row2 ]
##(Meemansa)for bar plot
xsimulatedBar <- xsimulated
print(setdiff(colnames(xorig), colnames(xsimulated)))

# Reshape
xorig = as.data.frame(sapply(xorig, as.numeric))
real_melt = melt(xorig)
real_melt$Category = "Real"
##(Meemansa)for bar plot
real_melt_Bar = reshape(xorigBar, direction = "long", varying = list(colnames(xorigBar)),  v.names = "value")
##change colids with colnames##
real_melt_Bar[["time"]] <-rep(colnames(xorigBar), each = nrow(xorigBar))
real_melt_Bar$Category = "Real"
real_melt_Bar$id = NULL

xsimulated = as.data.frame(sapply(xsimulated, as.numeric))
simulated_melt = melt(xsimulated)
simulated_melt$Category = "Simulated"
##(Meemansa)for bar plot
simulated_melt_Bar =  reshape(xsimulatedBar, direction = "long", varying = list(colnames(xsimulatedBar)),  v.names = "value")
##change colids with colnames##
simulated_melt_Bar[["time"]] <-rep(colnames(xsimulatedBar), each = nrow(xsimulatedBar)) 
simulated_melt_Bar$Category = "Simulated"
simulated_melt_Bar$id = NULL

x_xorig = rbind(real_melt,simulated_melt)
##(Meemansa)for bar plot
x_xorig_Bar = rbind(real_melt_Bar,simulated_melt_Bar)
names(x_xorig_Bar)[1] = "variable"

# Rename nodes
x_xorig$variable = as.character(x_xorig$variable)
x_xorig$variable = sub("Patient_", "", x_xorig$variable)
x_xorig$variable = sub("UPDRS_", "", x_xorig$variable)
x_xorig$variable = sub("CSF_", "", x_xorig$variable)
x_xorig$variable = sub("RBD_", "", x_xorig$variable)

##for Bar Plot
x_xorig_Bar$variable = as.character(x_xorig_Bar$variable)
x_xorig_Bar$variable = sub("Patient_", "", x_xorig_Bar$variable)
x_xorig_Bar$variable = sub("UPDRS_", "", x_xorig_Bar$variable)
x_xorig_Bar$variable = sub("CSF_", "", x_xorig_Bar$variable)
x_xorig_Bar$variable = sub("RBD_", "", x_xorig_Bar$variable)

png("/Users/Meems/Documents/PhDWork/BayesianNetworkAD/CodeAndResultsforPaper/PaperRevisions/Variable_density_generative.png", width = 12, height = 12, units = 'in', res = 600)
ggplot(data = x_xorig, aes(value,colour=Category)) +  facet_wrap(~ variable, scales = "free") +theme(strip.text.x = element_text(size=4)) + geom_density()+theme_classic()          
dev.off()   
png("/Users/Meems/Documents/PhDWork/BayesianNetworkAD/CodeAndResultsforPaper/PaperRevisions/Variable_density_generative.png", width = 12, height = 12, units = 'in', res = 600)
ggplot(data = x_xorig_Bar, aes(value,fill=Category, colour=Category)) +  facet_wrap(~ variable, scales = "free") +theme(strip.text.x = element_text(size=4)) + geom_histogram(position="dodge",stat="count")+theme_classic() + theme(axis.text = element_text(size=10),
                                                                                                                                                                                                                                  axis.text.x = element_text(size = 10,  angle=90, hjust = 1),
                                                                                                                                                                                                                                  axis.title = element_text(size = 10, face = "bold"),
                                                                                                                                                                                                                                 strip.text = element_text(size = 10))
dev.off()

# =========== NOT Generative method =======================
#load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/3_VirtualPatient_211.RData")
load("~/Documents/PhDWork/BayesianNetworkAD/CodeAndResultsforPaper/PPMI/Workspace/3_VirtualPatient_211.RData")
rm(list = setdiff(ls() , c("dics_data", "dt_bl", "dt_wl","real","real","virtual")))
xorig = real
xsimulated = virtual

# remove aux
remove_row = grep("aux", colnames(xorig) , value = TRUE)
xorig <- xorig[, ! colnames(xorig) %in% remove_row ]

##remove aux
remove_row2 = grep("aux", colnames(xsimulated) , value = TRUE)
xsimulated <- xsimulated[, ! colnames(xsimulated) %in% remove_row2 ]

#print(setdiff(colnames(xorig), colnames(xsimulated)))

# Reshape
xorig = as.data.frame(sapply(xorig, as.numeric))
real_melt = melt(xorig)
real_melt$Category = "Real"

xsimulated = as.data.frame(sapply(xsimulated, as.numeric))
simulated_melt = melt(xsimulated)
simulated_melt$Category = "Simulated"

x_xorig = rbind(real_melt,simulated_melt)

# Rename nodes
x_xorig$variable = as.character(x_xorig$variable)
x_xorig$variable = sub("Patient_", "", x_xorig$variable)
x_xorig$variable = sub("UPDRS_", "", x_xorig$variable)
x_xorig$variable = sub("CSF_", "", x_xorig$variable)
x_xorig$variable = sub("RBD_", "", x_xorig$variable)


png("Variable_density_Not_generative.png", width = 12, height = 12, units = 'in', res = 600)
ggplot(data = x_xorig, aes(value,colour=Category)) +  facet_wrap(~ variable, scales = "free") +theme(strip.text.x = element_text(size=4)) + geom_density()+theme_classic()          
dev.off()  
#===== make side by side plot======
load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/generativeVP.RData")
generative_aucplot = aucplot
generative_bx = bx
generative_mca_plot = mca_plot

load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/not_generativeVP.RData")

ggsave("boxplot_auc", arrangeGrob(aucplot, generative_aucplot))


require(gridExtra)
plot1 <- aucplot
plot2 <- generative_aucplot
grid.arrange(plot1, plot2, ncol=2)

par(mfrow=c(1,2))
generative_bx
bx

##mutual informartion between real and virtual subjects##
library(infotheo)
load("~/Documents/PhDWork/BayesianNetworkAD/CodeAndResultsforPaper/PPMI/Workspace/3_VirtualPatient_120.RData")
rm( list = setdiff(ls() , c("dics_data", "dt_bl", "dt_wl","real", "RealPep", "VirtPep") ))
VirtPep = VirtPep[1:362,]

mi <- sapply(RealPep, function(u){
  sapply(VirtPep, function(v){
    (mutinformation(u, v, method="emp")) # Function with column from x and column from y as inputs
  })
})

library(stringr)
colnames(mi) <- str_c(colnames(mi), '.real')
rownames(mi) <- str_c(rownames(mi), '.sim')

##remove aux columns##
rowAux <- grep("aux", rownames(mi), value = TRUE)
colAux <- grep("aux", colnames(mi), value = TRUE)
mi <- mi[setdiff(rownames(mi),rowAux), setdiff(colnames(mi),colAux)]
write.csv(mi,"~/Documents/PhDWork/BayesianNetworkAD/CodeAndResultsforPaper/PaperRevisions/mutualInfoPPMI.csv")
