# Replace INF value with NA
remove.inf = function(cohortdata){
  is.na(cohortdata) = sapply(cohortdata , is.nan)
  is.na(cohortdata) = sapply(cohortdata , is.infinite)
  return(cohortdata)
}

# Get stats on group
get_stats_on_column_number = function(cohortdata){
  mydata = cohortdata
  updrs_data = grep("UPDRS_", colnames(mydata), value = TRUE)
  Medicalhistory_data = grep("MedicalHistory_", colnames(mydata), value = TRUE)
  nmotor_data = grep("NonMotor_", colnames(mydata), value = TRUE)
  RBD_data = grep("RBD_", colnames(mydata), value = TRUE)
  CSF_data = grep("CSF_", colnames(mydata), value = TRUE)
  Biological_data = grep("Biological_", colnames(mydata), value = TRUE)
  Imaging_data = grep("Imaging.", colnames(mydata), value = TRUE)
  count_matrix  = setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("UPDRS", "MedicalHistory","NonMotor", "RBD", "CSF","Biological","Imaging"))
  count_matrix = rbind(count_matrix, data.frame("UPDRS" = length(updrs_data),
                                                "MedicalHistory"= length(Medicalhistory_data),
                                                "NonMotor" = length(nmotor_data),
                                                "RBD" = length(RBD_data),
                                                "CSF" =  length(CSF_data),
                                                "Biological" = length(Biological_data),
                                                "Imaging" = length(Imaging_data)))
  
  #rownames(count_matrix) = gsub(".*_.*_", "", colnames(mydata)[1])
  return(count_matrix)
}
# Auxiliary variables keep track of visit-wise and group-wise patient dropout. 
# Measurements of features are marked by value missing not at random (MNAR).
# MNAR  results from a systematic absence of subject data for a measurement type (feature group). 

get_aux_all_groups = function( cohortdata){
  
  mysample = cohortdata
  timepoint = str_extract(colnames(mysample), "V[0-9][0-9]")[1]
  
  UPDRS = select(mysample,grep( "UPDRS",colnames(mysample),value=TRUE)) 
  MedicalHistory = select(mysample,grep( "MedicalHistory",colnames(mysample),value=TRUE)) 
  NonMotor = select(mysample,grep( "NonMotor",colnames(mysample),value=TRUE)) 
  CSF = select(mysample,grep( "CSF",colnames(mysample),value=TRUE)) 
  RBD = select(mysample,grep( "RBD",colnames(mysample),value=TRUE)) 
  Biological = select(mysample,grep( "Biological",colnames(mysample),value=TRUE)) 
  Imaging = select(mysample,grep( "Imaging",colnames(mysample),value=TRUE)) 
  
  output_aux = function(mysubsample){
    
    #return_df = data.frame()
    #groupname = deparse(substitute(a))
    
    if(dim(mysubsample)[2] != 0 ){
            if(dim(mysubsample)[2]== 1){
              mysubsample = mysubsample
            }else{
              # Add a new column for AUX
              new = "new"
              in.loop = mysubsample
              mysubsample[new] <- 0
              
              # Get rownames where all value is NA 
              mysubsample_NA = which(apply(in.loop, 1, function(x) all(is.na(x))))
              mysubsample_pat = names(mysubsample_NA)
                            if(length(mysubsample_pat) !=0 ){
                              mysubsample[which(rownames(mysubsample) %in% mysubsample_pat ),]$new <- 1
                              
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

  UPDRS_aux = output_aux(mysubsample = UPDRS)
  MedicalHistory_aux = output_aux(mysubsample = MedicalHistory)
  NonMotor_aux = output_aux(mysubsample = NonMotor)
  CSF_aux = output_aux(mysubsample = CSF)
  RBD_aux = output_aux(mysubsample = RBD)
  Biological_aux = output_aux(mysubsample = Biological)
  Imaging_aux = output_aux(mysubsample = Imaging)
  
  outputdf <- data.frame(matrix("removelater", ncol = 1, nrow = nrow(mysample)))
  names(outputdf)[1]<- "toremove"
  if(dim(UPDRS)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , UPDRS_aux))
  } else{
    print(paste0("UPDRS data unavailable for visit", timepoint))
  }
  if(dim(MedicalHistory)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , MedicalHistory_aux))
  } else{
    print(paste0("MedicalHistory data unavailable for visit", timepoint))
  }
  if(dim(NonMotor)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , NonMotor_aux))
  } else{
    print(paste0("NonMotor data unavailable for visit", timepoint))
  }
  if(dim(CSF)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , CSF_aux))
  } else{
    print(paste0("CSF data unavailable for visit", timepoint))
  }
  
  if(dim(RBD)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , RBD_aux))
  } else{
    print(paste0("RBD data unavailable for visit", timepoint))
  }
  
  if(dim(Biological)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , Biological_aux))
  } else{
    print(paste0("Biological data unavailable for visit", timepoint))
  }
  
  if(dim(Imaging)[2] != 0 ){
    outputdf = as.data.frame(cbind(outputdf , Imaging_aux))
  } else{
    print(paste0("Imaging data unavailable for visit", timepoint))
  }
  
  outputdf$toremove = NULL
  #print("Aux done")
  
  return(outputdf)
}

# Autoencode
get_meta_feature_autoencoder = function(group_data , timepoint, groupname ){
  #groupdata =  as.data.frame(cbind(group_data, data.frame(y=y_variable) ))
  h2o.init(nthreads = -1)
  myx = as.h2o(group_data)
  n = round(dim(group_data)[2] /2)
  m = round(dim(group_data)[2] /4)
  #n = runif(1, min=20000000, max=99999999)
  r = sample(20:70000000, 1)
  
  hyper_params <- list(activation=c("Rectifier","TanhWithDropout"),  #RectifierWithDropout","TanhWithDropout","MaxoutWithDropout"
                       hidden = list(1, c(n, 1), c(n,m,1), c(m,1)),  # make it dynamic
                       input_dropout_ratio=c(0, 0.05, 0.2, 0.5),
                       #l1=seq(0,1e-4,1e-6),
                       #l2=seq(0,1e-4,1e-6)
                       l2=10^c(-4:4))
  
  
  grid = h2o.grid("deeplearning",
                  grid_id = paste("mygrid", r, sep="_"), 
                  x=colnames(myx), 
                  autoencoder = TRUE,
                  training_frame=myx,
                  seed=1234567, 
                  stopping_metric="MSE", 
                  stopping_rounds=5, 
                  epochs=500,
                  hyper_params = hyper_params, categorical_encoding = "AUTO")
  
  gbm_sorted_grid <- h2o.getGrid(grid_id = paste("mygrid", r, sep="_"), sort_by = "mse")
  print(gbm_sorted_grid)
  fit <- h2o.getModel(gbm_sorted_grid@model_ids[[1]])
  nlayers = length(strsplit(substr(gbm_sorted_grid@summary_table[1,1], 2, nchar(gbm_sorted_grid@summary_table[1,1])-1), ",")[[1]])
  newvar = as.data.frame(h2o.deepfeatures(fit, myx, nlayers))
  colnames(newvar) = paste(groupname , timepoint, sep= "_")
  
  output = list("model" = fit, "meta_feature" = newvar, "nlayers" = nlayers, "modelID" = gbm_sorted_grid@model_ids[[1]] )            
  h2o.saveModel(fit, path = "Auto_Model" ) 
  h2o.shutdown(prompt = FALSE)
  return(output)
}

# Autencoded value for  MedicalHistory, NonMotor and Biological groups at every visit
visit_autoencoded = function(cohort_data){
  
  Biological_features = select(cohort_data, grep("Biological", colnames(cohort_data), value = TRUE ))
  NonMotor_features = select(cohort_data, grep("NonMotor", colnames(cohort_data), value = TRUE ))
  MedicalHistory_features = select(cohort_data, grep("MedicalHistory", colnames(cohort_data), value = TRUE ))  
  
  if( dim(Biological_features)[2] > 1){
    groupName = sub("_.*$", "", colnames(Biological_features)[1])
    Visit = str_extract(colnames(Biological_features), "V[0-9][0-9]")[1]
    Biological_meta  = get_meta_feature_autoencoder(group_data = Biological_features , timepoint = Visit, groupname = groupName)
  }
  
  gc()
  if( dim(NonMotor_features)[2] > 1){
    groupName = sub("_.*$", "", colnames(NonMotor_features)[1])
    Visit = str_extract(colnames(NonMotor_features), "V[0-9][0-9]")[1]
    NonMotor_meta  = get_meta_feature_autoencoder(group_data = NonMotor_features , timepoint = Visit, groupname = groupName)
  }
  
  gc()
  if( dim(MedicalHistory_features)[2] > 1){
    groupName = sub("_.*$", "", colnames(MedicalHistory_features)[1])
    Visit = str_extract(colnames(MedicalHistory_features), "V[0-9][0-9]")[1]
    MedicalHistory_meta  = get_meta_feature_autoencoder(group_data = MedicalHistory_features , timepoint = Visit, groupname = groupName)
  }
  
  gc()
  meta_list = list("Biological_meta" = Biological_meta,
                   "NonMotor_meta" = NonMotor_meta,
                   "MedicalHistory_meta" = MedicalHistory_meta)
  
  return(meta_list)
}


# blacklist-whitelist
blacklist_whitelist = function(discPCA2){
   #discPCA2 = dics_data
  all.columns = colnames(discPCA2)
  
  visit11 = grep(".*_V11",all.columns,value=TRUE)  
  visit10 = grep(".*_V10",all.columns,value=TRUE)
  visit09 = grep(".*_V09",all.columns,value=TRUE)
  visit08 = grep(".*_V08",all.columns,value=TRUE)
  visit07 = grep(".*_V07",all.columns,value=TRUE)
  visit06 = grep(".*_V06",all.columns,value=TRUE)
  visit05 = grep(".*_V05",all.columns,value=TRUE)
  visit04 = grep(".*_V04",all.columns,value=TRUE)
  visit03 = grep(".*_V03",all.columns,value=TRUE)
  visit02 = grep(".*_V02",all.columns,value=TRUE)
  visit01 = grep(".*_V01",all.columns,value=TRUE)
  visitbl = grep(".*_V00",all.columns,value=TRUE) 
  visitSNP = grep("SNP_",all.columns,value=TRUE) 
  visitAUX = grep("aux",all.columns,value=TRUE) 
  visitpat = grep("Patient",all.columns,value=TRUE) 
  visitbl =  setdiff(visitbl, visitpat)
  
  
  # features not covered in thes vists are :
  all = c(visitbl,visit01, visit02,visit03,visit04,visit05, visit06, visit07, visit08, visit09 , visit10 , visit11)
  
  #---- Blacklist timepoint t+1 to t ------
  bl = data.frame()
  # From 11
  from11 = c(visitbl, visit01, visit02, visit03, visit04, visit05, visit06, visit07, visit08, visit09, visit10  )
  for(im in from11){
    bl = rbind.data.frame(bl, data.frame( from=visit11, to=rep(im, each=length(visit11)))) 
  }
  
  from10 = c(visitbl, visit01, visit02, visit03, visit04, visit05, visit06, visit07, visit08, visit09 )
  for(im in from10){
    bl = rbind.data.frame(bl, data.frame( from=visit10, to=rep(im, each=length(visit10)))) 
  }
  
  from09 = c(visitbl, visit01, visit02, visit03, visit04, visit05, visit06, visit07, visit08 )
  for(im in from09){
    bl = rbind.data.frame(bl, data.frame(from=visit09, to=rep(im, each=length(visit09)))) 
  }
  
  from08 = c(visitbl, visit01, visit02, visit03, visit04, visit05, visit06, visit07 )
  for(im in from08){
    bl = rbind.data.frame(bl, data.frame(from=visit08, to=rep(im, each=length(visit08)))) 
  }
  
  from07 =  c(visitbl, visit01, visit02, visit03, visit04, visit05, visit06 )
  for(im in from07){
    bl = rbind.data.frame(bl, data.frame(from=visit07, to=rep(im, each=length(visit07)))) 
  }
  
  from06 = c(visitbl, visit01, visit02, visit03, visit04, visit05)
  for(im in from06){
    bl = rbind.data.frame(bl, data.frame(from=visit06, to=rep(im, each=length(visit06)))) 
  }
  
  from05 = c(visitbl, visit01, visit02, visit03, visit04)
  for(im in from05){
    bl = rbind.data.frame(bl, data.frame(from=visit05, to=rep(im, each=length(visit05)))) 
  }
  
  from04 = c(visitbl, visit01, visit02, visit03)
  for(im in from04){
    bl = rbind.data.frame(bl, data.frame(from=visit04, to=rep(im, each=length(visit04)))) 
  }
  
  from03 = c(visitbl, visit01, visit02)
  for(im in from03){
    bl = rbind.data.frame(bl, data.frame(from=visit03, to=rep(im, each=length(visit03)))) 
  }
  
  from02 = c(visitbl, visit01)
  for(im in from02){
    bl = rbind.data.frame(bl, data.frame(from=visit02, to=rep(im, each=length(visit02)))) 
  }
  
  from01 = c(visitbl)
  for(im in from01){
    bl = rbind.data.frame(bl, data.frame(from=visit01, to=rep(im, each=length(visit01)))) 
  }
  
  
  bl$from <- as.character(bl$from)
  bl$to <- as.character(bl$to)
  #---------------------------------------------------------------------------------------------------------
  # From all longitudinal data to non-longitudinal data----
  
  
  bl2 = data.frame()
  
  nonL.data = setdiff(all.columns , all)
  from.all.Ldata = all
  for(im in nonL.data){
    
    bl2 = rbind.data.frame(bl2, data.frame(from = from.all.Ldata, to = rep(im, each=length(from.all.Ldata)))) 
  }
  bl2$from <- as.character(bl2$from)
  bl2$to <- as.character(bl2$to)
  
  #---------- Blacklist within biomarker group-------
  
  bl4 = data.frame()
  
  # updrs
  from.updrs =  grep("UPDRS",all.columns,value=TRUE)
  to.updrs = grep("Patient|MedicalHistory|NonMotor|Biological|SNP|CSF|RBD",all.columns,value=TRUE)
  
  for(im in to.updrs){
    bl4 = rbind.data.frame(bl4, data.frame(from=from.updrs, to=rep(im, each=length(from.updrs)))) 
  }
  
  #bio  CHECKß
  from.bio =  grep("Biological|CSF",all.columns,value=TRUE)
  to.bio = grep("Patient|SNP",all.columns,value=TRUE)
  
  for(im in to.bio){
    bl4 = rbind.data.frame(bl4, data.frame(from = from.bio, to=rep(im, each=length(from.bio)))) 
  }
  
  #Img   ADD later
  from.img =  grep("Imaging",all.columns,value=TRUE)
  to.img = grep("Biological|CSF|MedicalHistory|NonMotor|RBD|UPDRS|Patient|SNP",all.columns,value=TRUE)
  
  for(im in to.img){
    bl4 = rbind.data.frame(bl4, data.frame(from = from.img, to=rep(im, each=length(from.img))))
  }
  
  
  #Patient 
  from.patient = grep("Patient_", all.columns, value = TRUE)
  to.patient = grep("SNP|Patient", all.columns, value = TRUE)
  
  for(im in to.patient){
    bl4 = rbind.data.frame(bl4, data.frame(from = from.patient, to=rep(im, each=length(from.patient)))) 
  }
  
  #Non motor 
  from.nmotor =  grep("NonMotor|RBD",all.columns,value=TRUE)
  to.nmotor = grep("SNP|Patient|MedicalHistory|Biological|CSF",all.columns,value=TRUE)
  
  for(im in to.nmotor){
    bl4 = rbind.data.frame(bl4, data.frame(from=from.nmotor, to=rep(im, each=length(from.nmotor)))) 
  }
  
  # MedicalHistory 
  from.medicalhist =  grep("MedicalHistory",all.columns,value=TRUE)
  to.medicalhis = grep("SNP|Patient|Biological|CSF",all.columns,value=TRUE)
  
  for(im in to.medicalhis){
    bl4 = rbind.data.frame(bl4, data.frame(from=from.medicalhist, to=rep(im, each=length(from.medicalhist)))) 
  }
  
  #SNP 
  from.SNP = grep("SNP", all.columns, value = TRUE)
  to.SNP = grep("Patient_GENDER|Patient_SimpleGender|Patient_ENROLL_AGE", all.columns, value = TRUE)
  
  for(im in to.SNP){
    bl4 = rbind.data.frame(bl4, data.frame(from = from.SNP, to=rep(im, each=length(from.SNP)))) 
  }
  
  #Patient_GENDER can only influence  Patient_SimpleGender 
  from.gender = grep("Patient_GENDER", all.columns, value = TRUE)
  to.gender = setdiff(all.columns , "Patient_SimpleGender" )
  
  for(im in to.gender){
    bl4 = rbind.data.frame(bl4, data.frame(from = from.gender, to=rep(im, each=length(from.gender)))) 
  }
  
  
  bl4$from <- as.character(bl4$from)
  
  bl4$to <- as.character(bl4$to)

  #======================== From aux to aux----
  bl3 = data.frame()
  
  #1. General aux to aux
  aux_columns = grep("aux" ,all.columns, value= TRUE )
  all_aux_again = grep("aux" ,all.columns, value= TRUE )
  
  non_aux = setdiff(all.columns, aux_columns)
  
  #non_aux = setdiff(grep("UPDRS", all.columns, value = TRUE), aux_columns)
  for(im in aux_columns){
    
    bl3 = rbind.data.frame(bl3, data.frame(from = all_aux_again, to = rep(im, each=length(all_aux_again)))) 
  }
  
  
  for(imm in all_aux_again){
    
    bl3 = rbind.data.frame(bl3, data.frame(from = non_aux, to = rep(imm, each=length(non_aux)))) 
  }
  
  for(immm in non_aux){
    
    bl3 = rbind.data.frame(bl3, data.frame(from = all_aux_again, to = rep(immm, each=length(all_aux_again)))) 
  }
  
  bl3$from <- as.character(bl3$from)
  bl3$to <- as.character(bl3$to)
  
  # remove whitelist pair from blacklist pairs
  
  
  
  #================================================================================
  dt_bl = rbind(bl,bl2,bl4,bl3)
  
  #Whitelist aux column to those that created them. 
  #====== whitelist=======
  aux_columns = grep("aux" ,all.columns, value= TRUE )
  non_aux = setdiff(all.columns, aux_columns)
  dt_wl = data.frame()
  
  for(i in aux_columns){
    groupname = stringr::str_extract(i, "CSF|Biological|UPDRS|MedicalHistory|NonMotor")
    timepoint_loop = stringr::str_extract(i, "V00|V01|V02|V03|V04|V05|V06|V07|V08|V09|V10|V11")
    get_node = grep(timepoint_loop , grep(groupname , non_aux, value = TRUE), value = TRUE)
    dt_wl = rbind.data.frame(dt_wl, data.frame(from = i, to = get_node) )
    
  }
  
  output = list("blacklist" = dt_bl, "whitelist" = dt_wl)
  
  return(output)
}




