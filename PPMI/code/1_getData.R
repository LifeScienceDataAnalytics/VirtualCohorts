setwd("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature")
source('~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Code/main.R', echo=FALSE)
library(parallel)
library(dplyr)
library(missForest)
library(parallel)
library(doMC)

varDoc1  = varDoc
varDoc1$Feature=rownames(varDoc)
pd <- getCohort(cohort="PD", enrolledOnly=TRUE)

# baseline features with less than 50% missing value

baseline_feature<- as.data.frame(extractVariables(
  patients= pd,
  variables= setdiff( rownames(varDoc1) , rownames(varDoc1[varDoc1$Group== "Patient" ,])),
  events= "BL"
))
baseline_feature = baseline_feature[ , -which(colMeans(is.na(baseline_feature)) > 0.5)]

#1
updrs_score <- as.data.frame(extractVariables(
  patients= pd,
  variables= intersect(colnames(baseline_feature), c("UPDRS","UPDRS1","UPDRS2", "UPDRS3","UPDRS","UPDRS4")) ,
  events=c("BL","V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V13","V14")
))

#2
Medicalhistory_score <- as.data.frame(extractVariables(
  patients= pd,
  variables= intersect(colnames(baseline_feature), rownames(findVarByGroup("Medical history"))),
  events=c("BL","V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V13","V14")
))

#3
nmotor_score <- as.data.frame(extractVariables(
  patients= pd,
  variables= intersect(colnames(baseline_feature), c("DVT_TOTAL_RECALL","DVS_LNS", "QUIP","SCOPA","STA", "ESS" )),
  events=c("BL","V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V13","V14")
))

RBD_score <- as.data.frame(extractVariables(
  patients= pd,
  variables= intersect(colnames(baseline_feature),c("RBD")),
  events=c("BL","V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V13","V14")
))

#4
CSF_score <- as.data.frame(extractVariables(
  patients= pd,
  variables= intersect(colnames(baseline_feature),c(rownames(varDoc1[varDoc1$`Sub-group`== "Cerebrospinal fluid" ,]), rownames(varDoc1[varDoc1$`Sub-group`== "Cerebrospinal Fluid" ,]))),
  events=c("BL","V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V13","V14")
))

Biological_score <- as.data.frame(extractVariables(
  patients= pd,
  variables= intersect(colnames(baseline_feature),setdiff(rownames(findVarByGroup("Biological")) , c(rownames(varDoc1[varDoc1$`Sub-group`== "Cerebrospinal fluid" ,]), rownames(varDoc1[varDoc1$`Sub-group`== "Cerebrospinal Fluid" ,])))),
  events=c("BL","V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V13","V14")
))

#5
Imaging_score <- as.data.frame(extractVariables(
  patients= pd,
  variables= intersect(colnames(baseline_feature),rownames(findVarByGroup("Imaging"))),
  events=c("BL","V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V13","V14")
))

#6.1 
Patient_Demographic <- as.data.frame(extractVariables(
  patients= pd,
  variables=  setdiff(rownames(varDoc[varDoc$Group == "Patient" & varDoc$`Sub-group` == "Demographic",]) , c("APPRDX","BIRTHDT","GENDER","SimpleGender","ENROLL_AGE" ))
))
colnames(Patient_Demographic) = paste0("Patient_",colnames(Patient_Demographic))

#6.2
Patient_PDhistory <- as.data.frame(extractVariables(
  patients= pd,
  variables = c(rownames(varDoc[varDoc$Group == "Patient" & varDoc$`Sub-group` == "PD history",]), 
                rownames(varDoc[varDoc$Group == "Patient" & varDoc$`Sub-group` == "Socio-economic",]))
))
colnames(Patient_PDhistory) = paste0("Patient_",colnames(Patient_PDhistory))

#6.3
Patient_AgeGender <- as.data.frame(extractVariables(
  patients= pd,
  variables= c("GENDER","SimpleGender","ENROLL_AGE")
))
Patient_AgeGender$GENDER = ifelse(Patient_AgeGender$GENDER == "Female of child bearing potential", 1, 0)
Patient_AgeGender$GENDER = as.factor(Patient_AgeGender$GENDER)
colnames(Patient_AgeGender) = paste0("Patient_",colnames(Patient_AgeGender))

# Add group name to columns
colnames(updrs_score) = paste0("UPDRS_",colnames(updrs_score))
colnames(Medicalhistory_score) = paste0("MedicalHistory_",colnames(Medicalhistory_score))
colnames(nmotor_score) = paste0("NonMotor_",colnames(nmotor_score))
colnames(RBD_score) = paste0("RBD_",colnames(RBD_score))
colnames(CSF_score) = paste0("CSF_",colnames(CSF_score))
colnames(Biological_score) = paste0("Biological_",colnames(Biological_score))
colnames(Imaging_score) = paste0("Imaging_",colnames(Imaging_score))


# Remove columns with less than 50% missing value
updrs_score = updrs_score[ , -which(colMeans(is.na(updrs_score)) > 0.5)]
Medicalhistory_score = Medicalhistory_score[ , -which(colMeans(is.na(Medicalhistory_score)) > 0.5)]
nmotor_score = nmotor_score[ , -which(colMeans(is.na(nmotor_score)) > 0.5)]
RBD_score = RBD_score[ , -which(colMeans(is.na(RBD_score)) > 0.5)]
CSF_score = CSF_score[ , -which(colMeans(is.na(CSF_score)) > 0.5)]
Biological_score = Biological_score[ , -which(colMeans(is.na(Biological_score)) > 0.5)]
Imaging_score = Imaging_score[ , -which(colMeans(is.na(Imaging_score)) > 0.5)]
Imaging_score = as.data.frame(Imaging_score)
colnames(Imaging_score) = "Imaging.BL"

allData = cbind(updrs_score,Medicalhistory_score,nmotor_score,RBD_score,CSF_score,Biological_score,Imaging_score)
remove(updrs_score, Medicalhistory_score,nmotor_score,RBD_score,CSF_score,Biological_score,Imaging_score )

visitData = list("visitbl" = allData[, grep(".BL", colnames(allData), value = TRUE)],
                 "visit1" = allData[, grep(".V01", colnames(allData), value = TRUE)],
                 "visit2" = allData[, grep(".V02", colnames(allData), value = TRUE)],
                 "visit3" = allData[, grep(".V03", colnames(allData), value = TRUE)],
                 "visit4" = allData[, grep(".V04", colnames(allData), value = TRUE)],
                 "visit5" = allData[, grep(".V05", colnames(allData), value = TRUE)],
                 "visit6" = allData[, grep(".V06", colnames(allData), value = TRUE)],
                 "visit7" = allData[, grep(".V07", colnames(allData), value = TRUE)],
                 "visit8" = allData[, grep(".V08", colnames(allData), value = TRUE)],
                 "visit9" = allData[, grep(".V09", colnames(allData), value = TRUE)],
                 "visit10" = allData[, grep(".V10", colnames(allData), value = TRUE)],
                 "visit11" = allData[, grep(".V11", colnames(allData), value = TRUE)],
                 "patientData" = cbind(Patient_Demographic, Patient_PDhistory, Patient_AgeGender))

# Change .BL to V00 
colnames(visitData$visitbl) = sub("BL", "V00", colnames(visitData$visitb)) 


# Replace INF value with NA
visitData = sapply(visitData, remove.inf)
visitData_stat = sapply(visitData,get_stats_on_column_number)

# Create auxillary columns - group-wise and visit-wise
visitData_aux = sapply(visitData, get_aux_all_groups)
visitData_aux$patientData = visitData$patientData

#Impute value (visit wise)
set.seed(123)
visitData_imputed = imputedData = sapply(visitData_aux, function(x)missForest::missForest(x, ntree = 500)[1])
set.seed(123)
visitData_OBB = sapply(visitData_aux, function(x) missForest::missForest(x, ntree = 500)[2])

# Autoencode MedicalHistory, NonMotor and Biological groups at every visit
visitbl = visitData_imputed$visitbl.ximp
visit1 = visitData_imputed$visit1.ximp
visit2 = visitData_imputed$visit2.ximp
visit3 = visitData_imputed$visit3.ximp
visit4 = visitData_imputed$visit4.ximp
visit5 = visitData_imputed$visit5.ximp
visit6 = visitData_imputed$visit6.ximp
visit7 = visitData_imputed$visit7.ximp
visit8 = visitData_imputed$visit8.ximp
visit9 = visitData_imputed$visit9.ximp
visit10 = visitData_imputed$visit10.ximp
visit11 = visitData_imputed$visit11.ximp
patientData = visitData_imputed$patientData.ximp
# Separate aux columns before autoencoding

all_visit = cbind(visitbl,visit1,visit2,visit3,visit4,visit5,visit6,visit7,visit8,visit9,visit10,visit11,patientData)
aux_col = select(all_visit, grep("aux", colnames(all_visit), value = TRUE))

visitbl = visitbl[,-which(names(visitbl) %in% colnames(aux_col))]
visit1 = visit1[,-which(names(visit1) %in% colnames(aux_col))]
visit2 = visit2[,-which(names(visit2) %in% colnames(aux_col))]
visit3 = visit3[,-which(names(visit3) %in% colnames(aux_col))]
visit4 = visit4[,-which(names(visit4) %in% colnames(aux_col))]
visit5 = visit5[,-which(names(visit5) %in% colnames(aux_col))]
visit6 = visit6[,-which(names(visit6) %in% colnames(aux_col))]
visit7 = visit7[,-which(names(visit7) %in% colnames(aux_col))]
visit8 = visit8[,-which(names(visit8) %in% colnames(aux_col))]
visit9 = visit9[,-which(names(visit9) %in% colnames(aux_col))]
visit10 = visit10[,-which(names(visit10) %in% colnames(aux_col))]
visit11 = visit11[,-which(names(visit11) %in% colnames(aux_col))]

# Get autoencoded MedicalHistory, NonMotor and Biological groups at every visit
# visitbl_meta = visit_autoencoded(cohort_data = visitbl)
# visit1_meta = visit_autoencoded(cohort_data = visit1)
# visit2_meta = visit_autoencoded(visit2)
# visit3_meta = visit_autoencoded(visit3)
# visit4_meta = visit_autoencoded(visit4)
# visit5_meta = visit_autoencoded(visit5)
# visit6_meta = visit_autoencoded(visit6)
# visit7_meta = visit_autoencoded(visit7)
# visit8_meta = visit_autoencoded(visit8)
# visit9_meta = visit_autoencoded(visit9)
# visit10_meta = visit_autoencoded(visit10)
# visit11_meta = visit_autoencoded(visit11)

# visitbl_meta = visit_autoencoded(cohort_data = visitbl)

#====== visit visitbl =====
cohort_data = visitbl
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))

Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

visitbl_list = list("Biological_meta" = Biological_meta,
                 "NonMotor_meta" = NonMotor_meta,
                 "MedicalHistory_meta" = MedicalHistory_meta)

visitbl = select(cohort_data,-c(colnames(Biological_features), colnames(NonMotor_features), colnames(MedicalHistory_features)))
visitbl = cbind(visitbl, Biological_meta$meta_feature["Biological_V00"], NonMotor_meta$meta_feature["NonMotor_V00"],MedicalHistory_meta$meta_feature["MedicalHistory_V00"])

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#====== visit visitbl =====

#====== visit visit1 =====
cohort_data = visit1
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
# NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

visit1_list = list( "MedicalHistory_meta" = MedicalHistory_meta)

cohort_data = cohort_data[, -which(colnames(cohort_data) %in% colnames(MedicalHistory_features))]
visit1 = cbind(cohort_data,MedicalHistory_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#====== visit visit1 =====

#====== visit visit2 =====
cohort_data = visit2
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

visit2_list = list("NonMotor_meta" = NonMotor_meta, "MedicalHistory_meta" = MedicalHistory_meta)

cohort_data = cohort_data[, -which(colnames(cohort_data) %in% c(colnames(MedicalHistory_features), colnames(NonMotor_features)))]
visit2 = cbind(cohort_data,MedicalHistory_meta$meta_feature , NonMotor_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#====== visit visit2 =====

#====== visit3 =====
cohort_data = visit3
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
#NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

visit3_list = list( "MedicalHistory_meta" = MedicalHistory_meta)

cohort_data = cohort_data[, -which(colnames(cohort_data) %in% colnames(MedicalHistory_features))]
visit3 = cbind(cohort_data,MedicalHistory_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#======  visit3 =====

#====== visit4 =====
cohort_data = visit4
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

visit4_list = list( "MedicalHistory_meta" = MedicalHistory_meta, "NonMotor_meta" = NonMotor_meta)

cohort_data = cohort_data[, -which(colnames(cohort_data) %in% c(colnames(MedicalHistory_features), colnames(NonMotor_features)))]
visit4 = cbind(cohort_data,MedicalHistory_meta$meta_feature,NonMotor_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#======  visit4 =====

#====== visit5 =====
cohort_data = visit5
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
#NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

visit5_list = list( "MedicalHistory_meta" = MedicalHistory_meta)
cohort_data = cohort_data[, -which(colnames(cohort_data) %in% colnames(MedicalHistory_features))]
visit5 = cbind(cohort_data,MedicalHistory_meta$meta_feature)

#visit5_list = list( "MedicalHistory_meta" = MedicalHistory_meta, "NonMotor_meta" = NonMotor_meta)
# cohort_data = cohort_data[, -which(colnames(cohort_data) %in% c(colnames(MedicalHistory_features), colnames(NonMotor_features)))]
# visit5 = cbind(cohort_data,MedicalHistory_meta$meta_feature,NonMotor_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#======  visit5 =====

#====== visit6 =====
cohort_data = visit6
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

# visit6_list = list( "MedicalHistory_meta" = MedicalHistory_meta)
# cohort_data = cohort_data[, -which(colnames(cohort_data) %in% colnames(MedicalHistory_features))]
# visit6 = cbind(cohort_data,MedicalHistory_meta$meta_feature)

visit6_list = list( "MedicalHistory_meta" = MedicalHistory_meta, "NonMotor_meta" = NonMotor_meta)
cohort_data = cohort_data[, -which(colnames(cohort_data) %in% c(colnames(MedicalHistory_features), colnames(NonMotor_features)))]
visit6 = cbind(cohort_data,MedicalHistory_meta$meta_feature,NonMotor_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#======  visit6 =====

#====== visit7 =====
cohort_data = visit7
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
#NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

visit7_list = list( "MedicalHistory_meta" = MedicalHistory_meta)
cohort_data = cohort_data[, -which(colnames(cohort_data) %in% colnames(MedicalHistory_features))]
visit7 = cbind(cohort_data,MedicalHistory_meta$meta_feature)

# visit7_list = list( "MedicalHistory_meta" = MedicalHistory_meta, "NonMotor_meta" = NonMotor_meta)
# cohort_data = cohort_data[, -which(colnames(cohort_data) %in% c(colnames(MedicalHistory_features), colnames(NonMotor_features)))]
# visit7 = cbind(cohort_data,MedicalHistory_meta$meta_feature,NonMotor_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#======  visit7 =====

#====== visit8 =====
cohort_data = visit8
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

# visit7_list = list( "MedicalHistory_meta" = MedicalHistory_meta)
# cohort_data = cohort_data[, -which(colnames(cohort_data) %in% colnames(MedicalHistory_features))]
# visit7 = cbind(cohort_data,MedicalHistory_meta$meta_feature)

visit8_list = list( "Biological_meta" = Biological_meta, "MedicalHistory_meta" = MedicalHistory_meta, "NonMotor_meta" = NonMotor_meta)
cohort_data = cohort_data[, -which(colnames(cohort_data) %in% c(colnames(MedicalHistory_features), colnames(NonMotor_features), colnames(Biological_features)))]
visit8 = cbind(cohort_data,MedicalHistory_meta$meta_feature,NonMotor_meta$meta_feature, Biological_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#======  visit8 =====

#====== visit9 =====
cohort_data = visit9
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
#NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

visit9_list = list( "MedicalHistory_meta" = MedicalHistory_meta)
cohort_data = cohort_data[, -which(colnames(cohort_data) %in% colnames(MedicalHistory_features))]
visit9 = cbind(cohort_data,MedicalHistory_meta$meta_feature)

# visit7_list = list( "MedicalHistory_meta" = MedicalHistory_meta, "NonMotor_meta" = NonMotor_meta)
# cohort_data = cohort_data[, -which(colnames(cohort_data) %in% c(colnames(MedicalHistory_features), colnames(NonMotor_features)))]
# visit7 = cbind(cohort_data,MedicalHistory_meta$meta_feature,NonMotor_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#======  visit9 =====

#====== visit10 =====
cohort_data = visit10
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))
# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

# visit10_list = list( "MedicalHistory_meta" = MedicalHistory_meta)
# cohort_data = cohort_data[, -which(colnames(cohort_data) %in% colnames(MedicalHistory_features))]
# visit10 = cbind(cohort_data,MedicalHistory_meta$meta_feature)

visit10_list = list( "MedicalHistory_meta" = MedicalHistory_meta, "NonMotor_meta" = NonMotor_meta)
cohort_data = cohort_data[, -which(colnames(cohort_data) %in% c(colnames(MedicalHistory_features), colnames(NonMotor_features)))]
visit10 = cbind(cohort_data,MedicalHistory_meta$meta_feature,NonMotor_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#======  visit10 =====

#====== visit10 =====
cohort_data = visit11
Visit = str_extract(colnames(cohort_data), "V[0-9][0-9]")[1]

Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE )) 

cat(dim(Biological_features))
cat(dim(NonMotor_features))
cat(dim(MedicalHistory_features))

# Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = "Biological")
#NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = "NonMotor")
MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = "MedicalHistory")

visit11_list = list( "MedicalHistory_meta" = MedicalHistory_meta)
cohort_data = cohort_data[, -which(colnames(cohort_data) %in% colnames(MedicalHistory_features))]
visit11 = cbind(cohort_data,MedicalHistory_meta$meta_feature)

# visit11_list = list( "MedicalHistory_meta" = MedicalHistory_meta, "NonMotor_meta" = NonMotor_meta)
# cohort_data = cohort_data[, -which(colnames(cohort_data) %in% c(colnames(MedicalHistory_features), colnames(NonMotor_features)))]
# visit11 = cbind(cohort_data,MedicalHistory_meta$meta_feature,NonMotor_meta$meta_feature)

remove(cohort_data,Visit,Biological_features,NonMotor_features,MedicalHistory_features,Biological_meta,NonMotor_meta,MedicalHistory_meta )
#======  visit10 =====

#====== Patient group =====

Patient_Demographic_meta  = get_meta_feature_autoencoder(group_data = Patient_Demographic , timepoint = "V00", groupname = "Patient_Demographic")
Patient_PDhistory_meta  = get_meta_feature_autoencoder(group_data = Patient_PDhistory , timepoint = "V00", groupname = "Patient_PDhistory")
Patient_data = cbind(Patient_AgeGender,Patient_Demographic_meta$meta_feature["Patient_Demographic_V00"],Patient_PDhistory_meta$meta_feature["Patient_PDhistory_V00"]  )
#====== Patient group =====


#=================================================================================
# Change column suffix 
colnames(visitbl) = gsub(".V00", "_V00", colnames(visitbl))
colnames(visit1) = gsub(".V01", "_V01", colnames(visit1))
colnames(visit2) = gsub(".V02", "_V02", colnames(visit2))
colnames(visit3) = gsub(".V03", "_V03", colnames(visit3))
colnames(visit4) = gsub(".V04", "_V04", colnames(visit4))
colnames(visit5) = gsub(".V05", "_V05", colnames(visit5))
colnames(visit6) = gsub(".V06", "_V06", colnames(visit6))
colnames(visit7) = gsub(".V07", "_V07", colnames(visit7))
colnames(visit8) = gsub(".V08", "_V08", colnames(visit8))
colnames(visit9) = gsub(".V09", "_V09", colnames(visit9))
colnames(visit10) = gsub(".V10", "_V10", colnames(visit10))
colnames(visit11) = gsub(".V11", "_V11", colnames(visit11))

#=================================================================================
aux_col = aux_col[, sapply(aux_col, function(col) length(unique(col))) > 1]
aux_col = as.data.frame(sapply(aux_col, as.factor))
allData_meta = cbind(visitbl,visit1,visit2,visit3,visit4,visit5,visit6,visit7,visit8,visit9,visit10,visit11,Patient_data,aux_col )
save(allData_meta, file = "allData_meta.RData")



