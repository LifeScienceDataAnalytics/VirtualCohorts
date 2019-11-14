load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/stable_network.RData")
source('~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Code/main.R', echo=TRUE)
load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/Imputeddata_500.RData")
#Change CSF value

CSF_90 = select(all_visit, setdiff(grep("CSF", colnames(all_visit), value = TRUE), grep("aux",grep("CSF", colnames(all_visit), value = TRUE), value = TRUE)) )
CSF_90 = 0.1*CSF_90
colnames(CSF_90) = sub(".V", "_V",colnames(CSF_90) )

for(i in colnames(CSF_90)){
  bics_mm =  mclust::mclustBIC(CSF_90[,i] , modelNames = c("V", "E"))
  fit_mm = mclust::Mclust(CSF_90[,i] , x= bics_mm)
  CSF_90[,i] = fit_mm$classification
}

# Replace old CSF column with one 

dics_alpha = discPCA2[!names(discPCA2) %in% colnames(CSF_90)]
dics_alpha = cbind(dics_alpha,CSF_90 )

# ========================== ========================== ========================== ==========================
# predict using disc_meta_csf as new data
# ========================== ========================== ========================== ==========================

to_predict2 = colnames(dics_alpha)
fitted = bn.fit(finalBN, discPCA2 , method = "bayes")
fitted_grain = as.grain(fitted)
csf_prediction = setNames(data.frame(matrix(ncol = 1, nrow = 362)), "Remove")

for(j in to_predict2){
  
  pred = predict(fitted_grain, response =  j,
                 predictors = setdiff(to_predict2,j),
                 newdata = dics_alpha, 
                 type = "class" , method = "bayes-lw")
  
  pred = as.data.frame( pred$pred)
  csf_prediction = cbind(csf_prediction , pred)
}
csf_prediction$Remove = NULL 

#Melt result into one dataframe
library(reshape)
csf_prediction = as.data.frame(sapply(csf_prediction , as.numeric))
melt_csf_prediction = melt(csf_prediction)
melt_csf_prediction$Category = "Predicted"  

csf_observed = as.data.frame(sapply(dics_alpha, as.numeric))
melt_csf_observed = melt(csf_observed)
melt_csf_observed$Category = "Observed"   

melted_csf = as.data.frame(rbind(melt_csf_prediction, melt_csf_observed))
melted_csf$variable = as.character(melted_csf$variable)
melted_csf$visit = stringr::str_extract(melted_csf$variable ,"V00|V01|V02|V03|V04|V05|V06|V07|V08|V09|V10|V11")
melted_csf$Group = stringr::str_extract(melted_csf$variable ,"UPDRS|MedicalHistory|NonMotor|RBD|CSF|Biological|Imaging|Patient")
melted_csf = melted_csf[ -which(melted_csf$variable %in% grep("aux", melted_csf$variable, value = TRUE)), ]

melted_csf$variable = sub("UPDRS_", "", melted_csf$variable)
melted_csf$variable = sub("RBD_", "", melted_csf$variable)
#===Make plots for UPDRS
plot_list = list()
for(i in unique(melted_csf$Group)){
  print(i)
  Group_data = grep(i,melted_csf$variable, value = TRUE) 
  Group_melted = melted_csf[ melted_csf$variable %in% Group_data, ]
  
  ggp = ggplot() +
    geom_histogram(data = Group_melted, aes(x=value, fill=Category, y=..density..),binwidth=.5, alpha=.5, position="identity") +
    geom_density(data = Group_melted, aes(x=value, fill=Category,colour=Category), alpha=.3)+  facet_wrap(~variable)
   
   plot_list[[i]] = ggp
  # file_name = paste(i,"png",sep = ".")
  # png(file_name, width = 5, height = 5, units = 'in', res = 600)
   #ggplot(Group_melted, aes(value, colour=Category)) + geom_density()  + facet_wrap(~variable) +theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) 
  # dev.off()
  # 
}

# Save plots to tiff. Makes a separate file for each plot.
for (i in 1:length(plot_list)) {
  file_name = paste("csf_",  names(plot_list[i]), ".png", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}


# ========================== ========================== ========================== ==========================
#   Countrafactual age
# ========================== ========================== ========================== ==========================

load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/stable_network.RData")
source('~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Code/main.R', echo=TRUE)
load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/Imputeddata_500.RData")
rm( list = setdiff(ls() , c("disc_meta", "dt_bl", "dt_wl","dics_data" ,"discPCA2","finalBN","all_visit") ))


Age_20 =  all_visit$Patient_ENROLL_AGE-20
Age_20 = as.factor(create_bins(Age_20, breaks = c(40, 60,80, 100), method = "cuts"))


# Replace age

dics_20 = discPCA2
dics_20$Patient_ENROLL_AGE = Age_20

# predict using dics_20 as new data
# ========================== ========================== ========================== ==========================

to_predict2 = colnames(dics_20)
fitted = bn.fit(finalBN, discPCA2 , method = "bayes")
fitted_grain = as.grain(fitted)
age_prediction = setNames(data.frame(matrix(ncol = 1, nrow = 362)), "Remove")

for(j in to_predict2){
  
  pred = predict(fitted_grain, response =  j,
                 predictors = setdiff(to_predict2,j),
                 newdata = dics_20, 
                 type = "class" , method = "bayes-lw")
  
  pred = as.data.frame( pred$pred)
  age_prediction = cbind(age_prediction , pred)
}
age_prediction$Remove = NULL 

#Melt result into one dataframe
library(reshape)
age_prediction = as.data.frame(sapply(age_prediction , as.numeric))
melt_age_prediction = melt(age_prediction)
melt_age_prediction$Category = "Predicted"  

age_observed = as.data.frame(sapply(dics_20, as.numeric))
melt_age_observed = melt(age_observed)
melt_age_observed$Category = "Observed"   

melted_age = as.data.frame(rbind(melt_age_prediction, melt_age_observed))
melted_age$variable = as.character(melted_age$variable)
melted_age$visit = stringr::str_extract(melted_age$variable ,"V00|V01|V02|V03|V04|V05|V06|V07|V08|V09|V10|V11")
melted_age$Group = stringr::str_extract(melted_age$variable ,"UPDRS|MedicalHistory|NonMotor|RBD|CSF|Biological|Imaging|Patient")
melted_age = melted_age[ -which(melted_age$variable %in% grep("aux", melted_age$variable, value = TRUE)), ]

melted_age$variable = sub("UPDRS_", "", melted_age$variable)
melted_age$variable = sub("RBD_", "", melted_age$variable)
#===Make plots for UPDRS
plot_list = list()
for(i in unique(melted_age$Group)){
  print(i)
  Group_data = grep(i,melted_age$variable, value = TRUE) 
  Group_melted = melted_age[ melted_age$variable %in% Group_data, ]
  
  ggp = ggplot() +
    geom_histogram(data = Group_melted, aes(x=value, fill=Category, y=..density..),binwidth=.5, alpha=.5, position="identity") +
    geom_density(data = Group_melted, aes(x=value, fill=Category,colour=Category), alpha=.3)+  facet_wrap(~variable)
  
  plot_list[[i]] = ggp
  # file_name = paste(i,"png",sep = ".")
  # png(file_name, width = 5, height = 5, units = 'in', res = 600)
  #ggplot(Group_melted, aes(value, colour=Category)) + geom_density()  + facet_wrap(~variable) +theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) 
  # dev.off()
  # 
}

# Save plots to tiff. Makes a separate file for each plot.
for (i in 1:length(plot_list)) {
  file_name = paste("Age_",  names(plot_list[i]), ".png", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}

