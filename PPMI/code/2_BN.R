load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/allData_meta.RData")
source('~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Code/main.R', echo=FALSE)
library(mclust)
library(binst)
library(dplyr)

dics_data = allData_meta
names(dics_data)[names(dics_data) == 'Patient_Demographic_V00'] <- 'Patient_Demographic'
names(dics_data)[names(dics_data) == 'Patient_PDhistory_V00'] <- 'Patient_PDhistory'

factor_columns =  dics_data[,which(sapply(dics_data,class) == 'factor')]
num_columns = dics_data[, !(colnames(dics_data) %in% colnames(factor_columns))]
num_columns = num_columns[, sapply(num_columns, function(col) length(unique(col))) > 1]

nn = num_columns
bin_value = list()
bin_variance = list()
melt_df = setNames(data.frame(matrix(ncol = 3, nrow = 0)) , c("Value", "Node", "Group") )

#pdf(paste("plot", "mclust2",".pdf",sep=""))

for(i in colnames(nn)){
  bics_mm =  mclust::mclustBIC(nn[,i] , modelNames = c("V", "E"))
  fit_mm = mclust::Mclust(nn[,i] , x= bics_mm)
  bin_value[[i]] = as.data.frame(fit_mm$parameters$mean)[,1]
  bin_variance[[i]] = as.data.frame(fit_mm$parameters$variance$scale)
  nn[,i] = fit_mm$classification
  
  meta_df = data.frame("Value"= num_columns[,i],"Node" = i,"Group" = "metaValue" )
  meta_bin =  data.frame("Value"= fit_mm$parameters$mean,"Node" = i,"Group" = "mean" )
  melt_df = rbind(melt_df , rbind(meta_df ,meta_bin ))
  
  
}
#dev.off()
num_columns = as.data.frame(sapply(nn, as.factor))  # removed Patient_ENROLL_AGE, Patient_EDUCYRS 
num_columns = num_columns[, sapply(num_columns, function(col) length(unique(col))) > 1]

# Hard code groups for age

num_columns$Patient_ENROLL_AGE = as.factor(create_bins(dics_data$Patient_ENROLL_AGE, breaks = c(40, 60,80, 100), method = "cuts"))
dics_data = cbind(factor_columns,num_columns )

# ====== blacklist_whitelist=====

rm( list = setdiff( ls() , c("dics_data", "allData_meta","blacklist_whitelist")))

bl_wl = blacklist_whitelist(discPCA2 = dics_data)
dt_bl = bl_wl$blacklist
dt_wl = bl_wl$whitelist
discPCA2 = dics_data 

#Compare algorithms

registerDoMC(cores=3)
cl = makeCluster(3)

cvres1 = bn.cv(discPCA2, "rsmax2", runs=10, fit="bayes", loss="logl",  algorithm.args = list( blacklist=dt_bl, whitelist=dt_wl), cluster=cl) 
cvres2 = bn.cv(discPCA2, "mmhc", runs=10, fit="bayes", loss="logl",  algorithm.args = list(blacklist=dt_bl,  whitelist=dt_wl), cluster=cl) 
cvres3 = bn.cv(discPCA2, "hc", runs=10, fit="bayes", loss="logl", algorithm.args = list(maxp=5, blacklist=dt_bl,  whitelist=dt_wl, restart=10, score="bic"), cluster=cl) 
cvres4 = bn.cv(discPCA2, "tabu", runs=10, fit="bayes", loss="logl", algorithm.args = list(maxp=5, blacklist=dt_bl,  whitelist=dt_wl, restart=10, score="bic"), cluster=cl)
cvres5 = bn.cv(discPCA2, "si.hiton.pc", runs=10, fit="bayes", loss="logl", algorithm.args = list(blacklist=dt_bl,  whitelist=dt_wl, undirected=FALSE), cluster=cl)
cvres6 = bn.cv(discPCA2, "mmpc", runs=10, fit="bayes", loss="logl", algorithm.args = list(blacklist=dt_bl, whitelist=dt_wl, undirected=FALSE), cluster=cl)

pdf(paste("plot", "BNcv",".pdf",sep=""))
plot( cvres1, cvres2, cvres3,cvres4, cvres5 ,cvres6,
      xlab=c( "rsmax2", "mmhc","hc","tabu", "si.hiton.pc", "mmpc"))  # TABU :)
dev.off()
stopCluster(cl)

# Final bayesian network

registerDoMC(cores=3)
cl = makeCluster(3)

boot.stren = boot.strength(discPCA2, algorithm="tabu", R=1000, algorithm.args = list(maxp=5, blacklist=dt_bl, whitelist=dt_wl, restart=50, score="bic"), cluster=cl)
finalBN = tabu(discPCA2, maxp=5, blacklist=dt_bl,whitelist=dt_wl,  score="bic")
stopCluster(cl)

boot_subgraph = boot.stren[boot.stren$strength >0.5 & boot.stren$direction >0.5 , ]
boot_subgraph$strength = round(boot_subgraph$strength, 4 )
boot_subgraph$Group_from = str_extract(boot_subgraph$from, "[^_]+")
boot_subgraph$Group_to = str_extract(boot_subgraph$to, "[^_]+")


# ======== for plottings====
# boot_subgraph$from = sub("UPDRS_", "", boot_subgraph$from)
# boot_subgraph$to = sub("UPDRS_", "", boot_subgraph$to)
# boot_subgraph$from = sub("RBD_", "", boot_subgraph$from)
# boot_subgraph$to = sub("RBD_", "", boot_subgraph$to)
# boot_subgraph$from = sub("CSF_", "", boot_subgraph$from)
# boot_subgraph$to = sub("CSF_", "", boot_subgraph$to)
boot_subgraph$Network = "Bootstrapped"
write.csv(boot_subgraph, file= "BNetwork_GMM_5.csv")

BN_arc = as.data.frame(finalBN$arcs)
BN_arc$Group_from = str_extract(BN_arc$from, "[^_]+")
BN_arc$Group_to = str_extract(BN_arc$to, "[^_]+")
# BN_arc$from = sub("UPDRS_", "", BN_arc$from)
# BN_arc$to = sub("UPDRS_", "", BN_arc$to)
# BN_arc$from = sub("RBD_", "", BN_arc$from)
# BN_arc$to = sub("RBD_", "", BN_arc$to)
# BN_arc$from = sub("CSF_", "", BN_arc$from)
# BN_arc$to = sub("CSF_", "", BN_arc$to)
BN_arc$Network = "FinalBN"
write.csv(BN_arc, file= "FinalBN_arc.csv")

match11 = merge(boot_subgraph,BN_arc, by.x= c("from","to"), by.y= c("from","to"), all = FALSE)
match11 = select(match11 , c("from","to","strength","direction","Group_from.x","Group_to.x"))

Not_match = setdiff(boot_subgraph[,c("from","to")], BN_arc[,c("from","to")])
Not_match$strength = 0 
Not_match$direction = 0
Not_match$Group_from.x = str_extract(Not_match$from, "[^_]+")
Not_match$Group_to.x = str_extract(Not_match$to, "[^_]+")

match11$Category = "Overlap"
Not_match$Category = "NotOverlap"
matched = rbind(match11,Not_match)

write.csv(matched, file= "Overlapping_network.csv")

# ======== for plottings====


# method 1 
boot.stren_df = as.data.frame(boot.stren)
ggplot(boot.stren_df, aes(strength)) + stat_ecdf(geom = "step") +labs(title="Empirical Cumulative Density Function", y = "Threshold", x="Edges")+theme_classic()

# method 2 
library(ggplot2)
threshold.edges.df.new <-  data.frame("Threshold" = c(unique(boot.stren$strength)))
for(i in 1:nrow(threshold.edges.df.new)){
  thresh = threshold.edges.df.new[i,1]
  threshold.edges.df.new[i,"Edges"] <- nrow(boot.stren[boot.stren$strength >= thresh & boot.stren$direction >=thresh, ])
}
ggplot(threshold.edges.df.new, aes(Edges)) + stat_ecdf(geom = "step")+
  labs(title="Empirical Cumulative Density Function",
       y = "Threshold", x="Edges")+
  theme_classic()



save.image("stable_network.RData")

# Validation of network by 10 fold cross-validation

cfm_group = data.frame()
get.pro.average_auto = function(preg_grep , response_grep){
  # preg_grep = "Patient_|SNP"
  # response_grep = "_BL"
  freq_list = list()
  #step1: Subset Predciors 
  subset_discPCA = discPCA2
  
  independent_var = grep( preg_grep , colnames(subset_discPCA) , value = TRUE)
  independent_var = subset_discPCA[ , which(colnames(subset_discPCA) %in% independent_var)]
  to_predict_baseline_names =  grep( response_grep , colnames(subset_discPCA) , value = TRUE)
  to_predict_baseline_names = setdiff(to_predict_baseline_names , grep("aux" , colnames(subset_discPCA), value = TRUE))
  
  # Storage variable
  cfm_table_baseline_overall =  setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("Accuracy", "Kappa", "AccuracyLower" ,"AccuracyUpper" , "AccuracyNull" ,  "AccuracyPValue" , "McnemarPValue","Node"))
  cfm_table_baseline_overAllFold =  setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("Accuracy", "Kappa", "AccuracyLower" ,"AccuracyUpper" , "AccuracyNull" ,"AccuracyPValue" , "McnemarPValue"))
  all_baseline_prediction =  setNames(data.frame(matrix(ncol = 3, nrow = 0)), c( "pred", "oberseved","Node"))
  
  
  
  for(i in to_predict_baseline_names){
    print(i)
    
    to_predict_baseline = subset_discPCA[ , which(colnames(subset_discPCA) %in% i)]
    to_predict_baseline = as.data.frame(to_predict_baseline)
    colnames(to_predict_baseline) = i
    
    #Step3:  Add that baseline feature to predictor dataframe
    inloop_df = as.data.frame(cbind(independent_var , to_predict_baseline))
    
    # Step4: Create blacklist and whitelost for this 
    subset_blacklist_whitelist = blacklist_whitelist(inloop_df)
    
    here_bl = subset_blacklist_whitelist$blacklist
    here_wl = subset_blacklist_whitelist$whitelist
    #============= Predictor column added================
    
    entire_prediction = data.frame()
    cfm_table_baseline_avg =  setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("Accuracy", "Kappa", "AccuracyLower" ,"AccuracyUpper" ,"AccuracyNull" , "AccuracyPValue" , "McnemarPValue"))
    
    fold.index = createMultiFolds(inloop_df[, colnames(to_predict_baseline) ], k=10, times = 1)
    freq_major = data.frame() 
    
    for(j in 1:length(fold.index)){
      f_index = unlist(fold.index[[j]])
      train.set1 = inloop_df[f_index,]
      test.set1  = inloop_df[-f_index,]
      test.set2 <- test.set1[!names(test.set1) %in% colnames(to_predict_baseline)]
      
      
      #frequency of most majority class in test set
      for(a in colnames(test.set1)){
        y  = as.data.frame(table(test.set1[,a]))
        y = y[which.max(y$Freq),]
        y$Freq  =  (y$Freq)/ nrow(test.set1)
        rownames(y) = a
        y$Node = rownames(y)
        freq_major =  rbind(freq_major, y)
      }
      #freq_list[[ as.character(j)]] = list(freq_major) 
      
      
      #learn str
      finalBN_loop = tabu(train.set1, maxp=5, blacklist = here_bl,  score="bic")
      
      # update
      if(dim(here_wl)[1] == 0){
        finalBN_loop = tabu(train.set1, maxp=5, blacklist = here_bl,  score="bic")
      }else{
        finalBN_loop = tabu(train.set1, maxp=5, blacklist = here_bl,  whitelist= here_wl, score="bic")
        print("BN updated")
      }
      
      #learn parameters
      fitted = bn.fit(finalBN_loop, train.set1 , method = "bayes")
      
      #Make prediction
      pred = predict(fitted ,node = colnames(to_predict_baseline),  test.set2)
      
      #Save prediction
      prediction_observation = as.data.frame(cbind("pred"  = pred, "oberseved" = test.set1[, colnames(to_predict_baseline)]))
      
      
      #Confusion matrix
      if(length(unique(prediction_observation$oberseved)) >1 ){
        print("okay")
        prediction_observation = as.data.frame(sapply(prediction_observation, as.factor))
        cfm_table_baseline_avg1 = t(as.data.frame(caret::confusionMatrix(prediction_observation$pred,prediction_observation$oberseved)$overall))
        rownames(cfm_table_baseline_avg1) = i
        #cfm_table_baseline_avg1["Group"] = rownames(cfm_table_baseline_avg1)
        cfm_table_baseline_avg = rbind(cfm_table_baseline_avg , cfm_table_baseline_avg1)
      } else{cfm_table_baseline_avg1 = data.frame("Accuracy" = 0 ,
                                                  "Kappa" = 0 ,
                                                  "AccuracyLower" = 0 ,
                                                  "AccuracyUpper" = 0  ,
                                                  "AccuracyNull" = 0  ,
                                                  "AccuracyPValue" = 0 ,
                                                  "McnemarPValue"= 0)
      print("oops")
      rownames(cfm_table_baseline_avg1) = i
      cfm_table_baseline_avg = rbind(cfm_table_baseline_avg , cfm_table_baseline_avg1)
      }
      
      node = as.data.frame(rep(i, times= dim(prediction_observation)[1]) ) 
      colnames(node) = "Node"
      prediction_observation = cbind(prediction_observation , node)
      entire_prediction = rbind( entire_prediction ,prediction_observation)
      remove(node, prediction_observation)
      print(dim(entire_prediction))
      
      
      # Add group name
      # cfm_table_baseline_avg["Node"] = as.character(i)
    }
    
    #Save all predictions in a dataframe
    #colnames(entire_prediction) = c(paste("pred",i) , i ) 
    all_baseline_prediction = rbind(all_baseline_prediction ,entire_prediction )
    
    #save all internal confusion matrixs
    dummy2 = cfm_table_baseline_avg
    dummy2["Group"] = i
    cfm_group = rbind(cfm_group ,dummy2 )
    
    cfm_table_baseline_overall = rbind(cfm_table_baseline_overall ,cfm_table_baseline_avg )
    
    
  }
  
  
  
  output_list = list(all_baseline_prediction,cfm_table_baseline_overall, freq_major,cfm_group)
  names(output_list) = c("all_baseline_prediction" , "cfm_table_baseline_overall","freq_major","cfm_group")
  return(output_list)
}


# get seq prediction


prediction_bl = get.pro.average_auto( response_grep = "_V00" , preg_grep = "Patient_|SNP_" )  
prediction_01 = get.pro.average_auto( response_grep = "_V01" , preg_grep = "Patient_|SNP_|_V00")  
prediction_02 = get.pro.average_auto( response_grep = "_V02" , preg_grep = "Patient_|SNP_|_V00|_V01") 
prediction_03 = get.pro.average_auto( response_grep = "_V03" , preg_grep = "Patient_|SNP_|_V00|_V01|_V02") 
prediction_04 = get.pro.average_auto( response_grep = "_V04" , preg_grep = "Patient_|SNP_|_V00|_V01|_V02|_V03" )
prediction_05 = get.pro.average_auto( response_grep = "_V05" , preg_grep = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04")
prediction_06 = get.pro.average_auto( response_grep = "_V06" , preg_grep = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05")
prediction_07 = get.pro.average_auto( response_grep = "_V07" , preg_grep = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05|_V06")
prediction_08 = get.pro.average_auto( response_grep = "_V08" , preg_grep = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05|_V06|_V07" )
prediction_09 = get.pro.average_auto( response_grep = "_V09" , preg_grep = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05|_V06|_V07|_V08")


# plot 
library(plotrix)
library(ggplot2)
# Mean and sd.error for plots for majority class
stacked_freq_mean = aggregate( this_fold_fre[, c( "Freq")], list(this_fold_fre$Node), mean)
colnames(stacked_freq_mean) = c("Node", "Accuracy")
stacked_freq_sd = aggregate( this_fold_fre[, c( "Freq")], list(this_fold_fre$Node), std.error)
colnames(stacked_freq_sd) = c("Node", "SD")

stacked_chance = merge(stacked_freq_mean, stacked_freq_sd, by= "Node" )
stacked_chance$Class = "Chance level accuracy"


# Mean and sd.error for plots after prediction
stacked_cal_mean = aggregate( this_fold[, c( "Accuracy")], list(this_fold$Group), mean)
colnames(stacked_cal_mean) = c("Node", "Accuracy")
stacked_cal_sd = aggregate( this_fold[, c( "Accuracy")], list(this_fold$Group), std.error)
colnames(stacked_cal_sd) = c("Node", "SD")

stacked_calculated = merge(stacked_cal_mean, stacked_cal_sd, by= "Node" )
stacked_calculated$Class = "Calculated accuracy"

# Make columns names 
all_stacked = rbind(stacked_chance,stacked_calculated )
all_stacked$Visit = sub("^.*_","", all_stacked$Node )

remove_row = grep("aux", all_stacked$Node , value = TRUE)
all_stacked <- all_stacked[ ! all_stacked$Node %in% remove_row, ]


# plot for UPDRS
subset_updrs = all_stacked[all_stacked$Node %in% grep("UPDRS", all_stacked$Node, value= TRUE), ] 
subset_updrs = subset_updrs[!subset_updrs$Visit == "V09",]    #subset_updrs[!(subset_updrs$TargetNode=="UPDRS_V09 " & subset_updrs$Visit=="V09"),]
subset_updrs$Node = sub("UPDRS_", "", subset_updrs$Node)


png("UPDRS_vp.png", width = 8, height = 8, units = 'in', res = 600)
ggplot(subset_updrs, aes(x= subset_updrs$Node, y= subset_updrs$Accuracy, group= Class,
                         color = Class))+   
  labs(title = "Predicition performance of PPMI model")+
  xlab("UPDRS at each visit") +
  ylab("Accuracy") +
  geom_point() +
  geom_point(size=2, shape=23)+
  geom_errorbar(aes(ymax=subset_updrs$Accuracy + subset_updrs$SD,
                    ymin=subset_updrs$Accuracy - subset_updrs$SD),width=.2,
                position=position_dodge(0.01))  + theme_classic() +
  theme(text=element_text(size=14, face = "bold"),axis.text.x = element_text(angle=90, hjust=1))
dev.off()


# plot for Medical history
subset_Medical = all_stacked[all_stacked$Node %in% grep("MedicalHistory", all_stacked$Node, value= TRUE), ] 
subset_Medical = subset_Medical[!subset_Medical$Visit == "V09",] 

png("Medical1_accuracy.png", width = 8, height = 8, units = 'in', res = 600)
ggplot(subset_Medical, aes(x= subset_Medical$Node, y= subset_Medical$Accuracy, group= Class,
                           color = Class))+   
  labs(title = "Predicition performance of PPMI model")+
  xlab("Medical History at each visit") +
  ylab("Accuracy") +
  geom_point() +
  geom_point(size=2, shape=23)+
  geom_errorbar(aes(ymax=subset_Medical$Accuracy + subset_Medical$SD,
                    ymin=subset_Medical$Accuracy - subset_Medical$SD),width=.2,
                position=position_dodge(0.01))  + theme_classic() +
  theme(text=element_text(size=14, face = "bold"),axis.text.x = element_text(angle=90, hjust=1))
dev.off()


# plot for Biological
subset_Biological = all_stacked[all_stacked$Node %in% grep("Biological", all_stacked$Node, value= TRUE), ] 

png("Biological1_accuracy.png", width = 8, height = 8, units = 'in', res = 600)
ggplot(subset_Biological, aes(x= subset_Biological$Node, y= subset_Biological$Accuracy, group= Class,
                              color = Class))+   
  labs(title = "Predicition performance of PPMI model")+
  xlab("Biological at each visit") +
  ylab("Accuracy") +
  geom_point() +
  geom_point(size=2, shape=23)+
  geom_errorbar(aes(ymax=subset_Biological$Accuracy + subset_Biological$SD,
                    ymin=subset_Biological$Accuracy - subset_Biological$SD),width=.2,
                position=position_dodge(0.01))  + theme_classic() +
  theme(text=element_text(size=14, face = "bold"),axis.text.x = element_text(angle=90, hjust=1))
dev.off()


# plot for Non-motor
subset_NonMotor = all_stacked[all_stacked$Node %in% grep("NonMotor", all_stacked$Node, value= TRUE), ] 

png("NonMotor1_accuracy.png", width = 8, height = 8, units = 'in', res = 600)
ggplot(subset_NonMotor, aes(x= subset_NonMotor$Node, y= subset_NonMotor$Accuracy, group= Class,
                            color = Class))+   
  labs(title = "Predicition performance of PPMI model")+
  xlab("Non-Motor at each visit") +
  ylab("Accuracy") +
  geom_point() +
  geom_point(size=2, shape=23)+
  geom_errorbar(aes(ymax=subset_NonMotor$Accuracy + subset_NonMotor$SD,
                    ymin=subset_NonMotor$Accuracy - subset_NonMotor$SD),width=.2,
                position=position_dodge(0.01))  + theme_classic() +
  theme(text=element_text(size=14, face = "bold"),axis.text.x = element_text(angle=90, hjust=1))

dev.off()


# plot for CSF
subset_CSF = all_stacked[all_stacked$Node %in% grep("CSF_", all_stacked$Node, value= TRUE), ] 

png("CSF_accuracy.png", width = 8, height = 8, units = 'in', res = 600)
ggplot(subset_CSF, aes(x= subset_CSF$Node, y= subset_CSF$Accuracy, group= Class,
                       color = Class))+   
  labs(title = "Predicition performance of PPMI model")+
  xlab("CSF at each visit") +
  ylab("Accuracy") +
  geom_point() +
  geom_point(size=2, shape=23)+
  geom_errorbar(aes(ymax=subset_CSF$Accuracy + subset_CSF$SD,
                    ymin=subset_CSF$Accuracy - subset_CSF$SD),width=.2,
                position=position_dodge(0.01))  + theme_classic() +
  theme(text=element_text(size=14, face = "bold"),axis.text.x = element_text(angle=90, hjust=1))

dev.off()

# plot for RBD
subset_RBD = all_stacked[all_stacked$Node %in% grep("RBD_", all_stacked$Node, value= TRUE), ] 
subset_RBD$Node = sub("RBD_", "", subset_RBD$Node)


png("RBD_accuracy.png", width = 8, height = 8, units = 'in', res = 600)
ggplot(subset_RBD, aes(x= subset_RBD$Node, y= subset_RBD$Accuracy, group= Class,
                       color = Class))+   
  labs(title = "Predicition performance of PPMI model")+
  xlab("RBD at each visit") +
  ylab("Accuracy") +
  geom_point() +
  geom_point(size=2, shape=23)+
  geom_errorbar(aes(ymax=subset_RBD$Accuracy + subset_RBD$SD,
                    ymin=subset_RBD$Accuracy - subset_RBD$SD),width=.2,
                position=position_dodge(0.01))  + theme_classic() +
  theme(text=element_text(size=14, face = "bold"),axis.text.x = element_text(angle=90, hjust=1))

dev.off()




