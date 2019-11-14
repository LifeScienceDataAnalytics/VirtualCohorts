load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Workspace/stable_network.RData")
source('~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/PPMI_selectedFeature/Code/main.R', echo=TRUE)

train1 =  dics_data[-4,]
test1 = dics_data[4,]
BN_fold  = finalBN
# learn parameters
fitted = bn.fit(BN_fold, train1 , method = "bayes")
fitted_grain = as.grain(fitted)

# For saving prediction results
prediction_fold = setNames(data.frame(matrix(ncol = 3, nrow = 0)),  c( "Feature","Values", "Class" ))

sequence_prediction_one = function(grep_to_predict, grep_predictors , fitted_grain ,test1 ){
  # Prediction - save visit wise and later combine them
  # grep_to_predict = "_BL" # node to predict
  # grep_predictors = "Patient_|SNP_"  # predictor
  
  
  nodes_predict= setdiff(grep(grep_to_predict, colnames(train1) , value = TRUE) , grep("aux_", colnames(train1) , value = TRUE))
  #nodes_predict = setdiff(nodes_predict, "Patient_RABLACK")
  predictors_var = grep(grep_predictors, colnames(train1) , value = TRUE)
  
  
  df_loop = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Node", "prob", "level"))
  saveped = list()
  for(b in nodes_predict){
    pred = predict(fitted_grain, response = b,  predictors = predictors_var ,newdata = test1,  method = "bayes-lw", type = "distribution")
    df_loop1 = as.data.frame( t(as.data.frame(pred$pred)))
    df_loop1$level = sub( paste(b,".", sep = ""), "", rownames(df_loop1))
    df_loop1 = data.frame("Node" = rownames(df_loop1), "prob" = df_loop1$V1 , "level" =  df_loop1$level)
    df_loop = rbind(df_loop, df_loop1)
    
  }
  
  #save evidence 
  df_loop_evidence = setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Node", "prob"))
  saveped = list()
  for(b in nodes_predict){
    pred = predict(fitted_grain, response = b,  predictors = predictors_var ,newdata = test1,  method = "bayes-lw", type = "distribution")
    df_loop2  = pred$pEvidence
    df_loop2 = data.frame("Node" = b, "prob" = df_loop2)
    df_loop_evidence = rbind(df_loop_evidence , df_loop2)
  }
  
  #df_loop$level = sub('.*\\.', '', df_loop$Node)
  df_loop_melt = melt(df_loop)
  df_loop_melt$group = gsub("(.+?)(\\_.*)", "\\1", df_loop_melt$Node)
  
  predicted_node  = df_loop  # as.data.frame(t(as.data.frame(pred$pred)))
  #colnames(predicted_node) = "Values"
  predicted_node$Class = "Predicted"
  predicted_node$Feature = rownames(predicted_node)
  
  observed_node = as.data.frame(test1[, nodes_predict])
  observed_node = as.data.frame(t(as.data.frame(observed_node)))
  colnames(observed_node) = "level"
  observed_node$Node = rownames(observed_node)
  observed_node$Class = "Observed"
  observed_node$Feature = rownames(observed_node)
  observed_node$prob = 1
  
  prediction_fold_loop = rbind(predicted_node ,observed_node)
  prediction_fold_loop$group = gsub("(.+?)(\\_.*)", "\\1", prediction_fold_loop$Node)
  prediction_fold_loop$Visit = str_extract(prediction_fold_loop$Node, "V00|V01|V02|V03|V04|V05|V06|V07|V08|V09|V10|V11")
  prediction_fold = rbind(prediction_fold ,prediction_fold_loop )
  
  #output_for_one_visit = list("prediction_fold" = prediction_fold )
  
  return(prediction_fold)
}

prediction_bl = sequence_prediction_one( grep_to_predict = "_V00" , grep_predictors = "Patient_|SNP_", fitted_grain=fitted_grain ,test1=test1 )  
prediction_01 = sequence_prediction_one( grep_to_predict = "_V01" , grep_predictors = "Patient_|SNP_|_V00", fitted_grain=fitted_grain ,test1=test1 )  
prediction_02 = sequence_prediction_one( grep_to_predict = "_V02" , grep_predictors = "Patient_|SNP_|_V00|_V01", fitted_grain=fitted_grain ,test1=test1 ) 
prediction_03 = sequence_prediction_one( grep_to_predict = "_V03" , grep_predictors = "Patient_|SNP_|_V00|_V01|_V02", fitted_grain=fitted_grain ,test1=test1 ) 
prediction_04 = sequence_prediction_one( grep_to_predict = "_V04" , grep_predictors = "Patient_|SNP_|_V00|_V01|_V02|_V03", fitted_grain=fitted_grain ,test1=test1 )
prediction_05 = sequence_prediction_one( grep_to_predict = "_V05" , grep_predictors = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04", fitted_grain=fitted_grain ,test1=test1 )
prediction_06 = sequence_prediction_one( grep_to_predict = "_V06" , grep_predictors = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05", fitted_grain=fitted_grain ,test1=test1 )
prediction_07 = sequence_prediction_one( grep_to_predict = "_V07" , grep_predictors = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05|_V06", fitted_grain=fitted_grain ,test1=test1 )
prediction_08 = sequence_prediction_one( grep_to_predict = "_V08" , grep_predictors = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05|_V06|_V07", fitted_grain=fitted_grain ,test1=test1 )
prediction_09 = sequence_prediction_one( grep_to_predict = "_V09" , grep_predictors = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05|_V06|_V07|_V08", fitted_grain=fitted_grain ,test1=test1 )
prediction_10 = sequence_prediction_one( grep_to_predict = "_V10" , grep_predictors = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05|_V06|_V07|_V08|_V09", fitted_grain=fitted_grain ,test1=test1 )
prediction_11 = sequence_prediction_one( grep_to_predict = "_V11" , grep_predictors = "Patient_|SNP_|_V00|_V01|_V02|_V03|_V04|_V05|_V06|_V07|_V08|_V09|_V10", fitted_grain=fitted_grain ,test1=test1 )

this_fold_one = rbind(prediction_bl,
                      prediction_01,
                      prediction_02,
                      prediction_03,
                      prediction_04,
                      prediction_05,
                      prediction_06,
                      prediction_07,
                      prediction_08,
                      prediction_09,
                      prediction_10,
                      prediction_11)

rownames(this_fold_one) <- 1:nrow(this_fold_one)
this_fold_one$Feature = NULL

# plot 
this_fold_one$Node = as.character(this_fold_one$Node)

png("oneP.png", width = 20, height = 12, units = 'in', res = 600)
ggplot(this_fold_one, aes(Node, prob)) + 
  geom_point() +facet_grid(vars(group), vars(Visit)) + 
  geom_point(aes(colour = factor(Class)))+ 
  geom_bar(stat="identity" , color = this_fold_one$level,position = "dodge2") + 
  geom_text(aes(label=level), vjust=0.5, check_overlap = TRUE, position = position_nudge(x = 8,y = 0))
dev.off()



#  save.image("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/bitcluster/GMM_bn/generativeBN_185.RData")
#  load("~/Documents/Masters_thesis/Markdown_pages/CompareAllMethod/bitcluster/GMM_bn/generativeBN_185.RData")
rm( list = setdiff(ls() , c("disc_meta", "dt_bl", "dt_wl" , "real", "realCopy", "final_VP_BN" , "boot.stren","this_fold_one") ))