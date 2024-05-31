#' Separates the patients into two groups based on the score of interest
#'
#' This function separates the patients of a dataset into two groups: one which lowly express the desired metric and the other one which highly express the desired metric.
#'
#' @param applyFileUsed A character value indicating which `*.apply` file to use for patient separation. If set to `classification`, then the Classification Signature `*.apply` file located in the `output > apply` folder is used. If set to `custom`, then each `*.apply` file from the `output > apply > customSignatures > fullTables` folder are iteratively used.
#'
#' @param metricToUse A character value indicating which metric to use for patients separation. If set to `CS`, then the Clinical Scores for each patient will be used. If set to `QS`, then the Quantitative Scores for each patient will be used for patients separation.
#'
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @return Generated output files are exported to the `output > class` folder. Their precise amounts and organization varies according to the `applyFileUsed` argument (see the online tutorial for more details).
#'
#' @importFrom foreach %do%
#'
#' @export

separatePatients = function(applyFileUsed = "classification", metricToUse = "CS", verbose = TRUE)
{
  if(applyFileUsed == "classification")
  {
    signatureResultsToUseList = list.files(file.path("output", "apply"), pattern = "(.*)classificationSignature(.*).apply")
  } else if(applyFileUsed == "custom")
  {
    signatureResultsToUseList = list.files(file.path("output", "class", "input"), pattern = "*.apply")
  }
  
  a = NULL
  
  foreach::foreach(a = 1:length(signatureResultsToUseList)) %do%
    {
      currentSignatureResultsFile = signatureResultsToUseList[a]
      currentSignatureName = gsub("(.+)Sign = (.+)\\)(.+)", "\\2", currentSignatureResultsFile)
      
      if(verbose == TRUE)
      {	
        print(paste("## Currently working on the class : ", currentSignatureName, " (", a, "/", length(signatureResultsToUseList), ") ##", sep = ""))
      }		
      
      if(dir.exists(file.path("output", "class", currentSignatureName)) == TRUE)
      {
        unlink(file.path("output", "class", currentSignatureName, "*"), recursive = TRUE, force = TRUE)
      } else
      {
        dir.create(file.path("output", "class", currentSignatureName), showWarnings = TRUE, recursive = FALSE)
      }
      
      if(applyFileUsed == "classification")
      {
        currentSignatureResultsPath = file.path("output", "apply", currentSignatureResultsFile)
      } else if(applyFileUsed == "custom")
      {
        currentSignatureResultsPath = file.path("output", "class", "input", currentSignatureResultsFile)
      }
      
      currentSignatureResults = data.table::fread(input = currentSignatureResultsPath, sep = "\t", header = TRUE, data.table = FALSE)
      
      currentSignatureResults$vitalStatus = gsub("Alive", "0", as.character(currentSignatureResults$vitalStatus))
      currentSignatureResults$vitalStatus = as.numeric(gsub("Dead", "1", currentSignatureResults$vitalStatus))
      
      if("score" %in% colnames(currentSignatureResults))
      {
        colIDtoChange = which(colnames(currentSignatureResults) == "score")
        colnames(currentSignatureResults)[colIDtoChange] = "clinicalScore"
      }
      
      if(metricToUse == "QS")
      {
        columnIDToUse = which(colnames(currentSignatureResults) == "quantitativeScore")
      } else if (metricToUse == "CS")
      {
        columnIDToUse = which(colnames(currentSignatureResults) == "clinicalScore")
      } else
      {
        columnIDToUse = which(colnames(currentSignatureResults) == "clinicalScore")
      }
      
      # if(verbose == TRUE)
      # {
      #   print("## Random sampling of patients ##")
      # }
      
      # currentSignatureResults = currentSignatureResults[sample(nrow(currentSignatureResults)),]
      
      options(scipen=999)
      
      totalCases = length(rownames(currentSignatureResults))
      
      # if(verbose == TRUE)
      # {
      #   print("## Defining of training and validation datasets ##")
      # }
      
      # casesForTraining = as.integer(c(1:(round(totalCases/2))))
      # casesForValidation = as.integer(c((round(totalCases/2)+1):totalCases))
      # 
      # dataTraining = data.frame(currentSignatureResults[casesForTraining,])
      # 
      patientsScoreOutput = currentSignatureResults
      patientsScoreOutput = patientsScoreOutput[order(patientsScoreOutput[,columnIDToUse]),]
      patientsScoreOutput$group = ""
      
      if(verbose == TRUE)
      {
        print("## Creating ROC curves and determining the best cutoff value possible ##")
      }			
      
      cutoffResults = data.frame()
      
    
        if(metricToUse == "QS")
        {
          scoresSequence = sort(patientsScoreOutput[, columnIDToUse])
        } else if(metricToUse == "CS")
        {
          scoresSequence = sort(patientsScoreOutput[, columnIDToUse])
        } else
        {
          scoresSequence = sort(patientsScoreOutput[, columnIDToUse])
        }
      
      
      # o = NULL
      # 
      # foreach::foreach(o = seq_len(length(cutoffSequence))) %do%	
      #   {
      #     cutoff = cutoffSequence[o]
      #     
      #     patientsScoreOutputTraining[patientsScoreOutputTraining[, columnIDToUse] >= cutoff, "group"] = 1
      #     patientsScoreOutputTraining[patientsScoreOutputTraining[, columnIDToUse] < cutoff, "group"] = 0
      #     
      #     if(length(unique(patientsScoreOutputTraining$group)) > 1)
      #     {
      #       
            # The positive reference here is the LTS group (>= threshold), where group == 1
            # The negative reference here is the STS group (< threshold), where group == 0
            
            # currentThresholdTruePositives = nrow(patientsScoreOutputTraining[patientsScoreOutputTraining$group == 1 & patientsScoreOutputTraining$vitalStatus == 0,])
            # currentThresholdTrueNegatives = nrow(patientsScoreOutputTraining[patientsScoreOutputTraining$group == 0 & patientsScoreOutputTraining$vitalStatus == 1,])
            # 
            # currentThresholdFalsePositives = nrow(patientsScoreOutputTraining[patientsScoreOutputTraining$group == 1 & patientsScoreOutputTraining$vitalStatus == 1,])
            # currentThresholdFalseNegatives = nrow(patientsScoreOutputTraining[patientsScoreOutputTraining$group == 0 & patientsScoreOutputTraining$vitalStatus == 0,])
            # 
            # currentThresholdTruePositiveRate = currentThresholdTruePositives / (currentThresholdTruePositives + currentThresholdFalseNegatives)
            # currentThresholdFalsePositiveRate = currentThresholdFalsePositives / (currentThresholdFalsePositives + currentThresholdTrueNegatives)
            #   
              
            # pathToSave = file.path("output", "class", currentSignatureName, "cutoffTrainingCurves", paste("survivalCurve_cutoff=", cutoff , ".pdf", sep = ""))
            # 
            # pdf(pathToSave)
            # 
            # patientsScoreOutputTraining$group = as.numeric(patientsScoreOutputTraining$group)
            # patientsScoreOutputTraining$survivedDays = as.numeric(patientsScoreOutputTraining$survivedDays)
            # patientsScoreOutputTraining$vitalStatus = as.numeric(patientsScoreOutputTraining$vitalStatus)
            # 
            # survivalObjectGlobal = survival::Surv(time = patientsScoreOutputTraining$survivedDays, event = patientsScoreOutputTraining$vitalStatus)
            # survivalCurves = survival::survfit(survivalObjectGlobal ~ group, data = patientsScoreOutputTraining, conf.type = "log-log")
            # survivalTests = survival::survdiff(survivalObjectGlobal ~ group, data = patientsScoreOutputTraining, rho = 0)
            # 
            # lowAndHigh_results = survRM2::rmst2(patientsScoreOutputTraining$survivedDays, patientsScoreOutputTraining$vitalStatus, patientsScoreOutputTraining$group)
            # 
            # RMST_low = as.numeric(lowAndHigh_results$RMST.arm0$rmst[1])
            # RMST_high = as.numeric(lowAndHigh_results$RMST.arm1$rmst[1])
            # 
            # finalMetric = abs(RMST_low - RMST_high)
            # 
            # n1 = as.numeric(nrow(patientsScoreOutputTraining[patientsScoreOutputTraining$group == 0,]))
            # n2 = as.numeric(nrow(patientsScoreOutputTraining[patientsScoreOutputTraining$group == 1,]))
            # 
            # alivePatientsOfGoodPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 0 & patientsScoreOutputTraining$group == 1,]))
            # alivePatientsOfBadPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 0 & patientsScoreOutputTraining$group == 0,]))
            # deadPatientsOfGoodPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 1 & patientsScoreOutputTraining$group == 1,]))
            # deadPatientsOfBadPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 1 & patientsScoreOutputTraining$group == 0,]))
            # realAlivePatientsOfGoodPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 0,]))
            # 
            # legendText = paste("(blue = ", n1, " et red = ", n2, ")", sep = "")
            # legendKeys = c("<= to cutoff", "> to cutoff")
            # 
            # cutoffResults[o, "cutoff"] = cutoff
            # cutoffResults[o, "lessOrEqualThanThreshold"] = n1
            # cutoffResults[o, "moreThanThreshold"] = n2
            # cutoffResults[o, "n1DivByn2"] = n1/n2
            # cutoffResults[o, "truePositives"] = currentThresholdTruePositives
            # cutoffResults[o, "trueNegatives"] = currentThresholdTrueNegatives
            # cutoffResults[o, "falsePositives"] = currentThresholdFalsePositives
            # cutoffResults[o, "falseNegatives"] = currentThresholdFalseNegatives
            # cutoffResults[o, "truePositiveRate"] = currentThresholdTruePositiveRate
            # cutoffResults[o, "falsePositiveRate"] = currentThresholdFalsePositiveRate
            # 
            # currentPlot = graphics::plot(
            #   survivalCurves,
            #   col = c("blue", "red"),
            #   xlab = "Time (days)",
            #   ylab = "Overall Survival",
            #   lty=1,
            #   lwd=2,
            #   log = FALSE,
            #   main = paste("Cutoff : ", cutoff, sep = ""),
            # )
            # graphics::legend(100, 0.2, legendKeys, lty = 1:1, lwd = 2:2, col = c("blue", "red"))
            # graphics::mtext(legendText, side=3, adj=0.5, cex=0.8)
            # 
            # grDevices::dev.off()					
        #   }
        # }
        # 
      # cutoffResults = cutoffResults[!is.na(cutoffResults$cutoff),]
      # oldCutoffResults = cutoffResults
      # 
      

      ROC_object = ROCit::rocit(score = scoresSequence, class = patientsScoreOutput$vitalStatus, negref = 1)
      ROC_CI = ROCit::ciROC(ROC_object, level = 0.95)
      ROC_AUC_CI = ROCit::ciAUC(ROC_object, level = 0.95)
      ROC_dataframe = data.frame(Cutoff = ROC_object$Cutoff, TPR = ROC_object$TPR, FPR = ROC_object$FPR)
      
      pathToSave = file.path("output", "class", currentSignatureName, paste("fullDataset_ROC-Curve.pdf", sep = ""))
      
      grDevices::pdf(pathToSave)
      
      ROC_plot = plot(ROC_object, values = TRUE)
      
      graphics::lines(ROC_CI$UpperTPR~ROC_CI$FPR, col = 2, lty = 2)
      graphics::lines(ROC_CI$LowerTPR~ROC_CI$FPR, col = 2, lty = 2)
      
      # plot(ROC_CI, col = 1, legend = FALSE)
      
      # lines(ROC_CI$TPR~ROC_CI$FPR, col = 2, lwd = 2)
      # lines(ROC_CI$LowerTPR~ROC_CI$FPR, col = 2, lty = 2)
      # 
      # lines(ROC_CI$UpperTPR~ROC_CI$FPR, col = 2, lty = 2)
      # 
      grDevices::dev.off()
      
      ROC_values = data.frame()
      ROC_values[1, "AUC"] = as.numeric(ROC_plot$AUC)
      ROC_values[1, "AUC_CI95-Lower"] = as.numeric(ROC_AUC_CI$lower)
      ROC_values[1, "AUC_CI95-Upper"] = as.numeric(ROC_AUC_CI$upper)
      ROC_values[1, "OptimalYoudenIndexPoint_TPR"] = as.numeric(ROC_plot$"optimal Youden Index point"["TPR"])
      ROC_values[1, "OptimalYoudenIndexPoint_FPR"] = as.numeric(ROC_plot$"optimal Youden Index point"["FPR"])
      ROC_values[1, "OptimalYoudenIndexPoint_cutoff"] = as.numeric(ROC_plot$"optimal Youden Index point"["cutoff"])
     
      utils::write.table(ROC_values, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "fullDataset_ROC-GeneralResults.txt"), sep = "\t")
      
      utils::write.table(ROC_dataframe, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "fullDataset_ROC-FullResults.txt"), sep = "\t")
      # 
      # dataTrainingAlivePatientsNumber = length(rownames(dataTraining[dataTraining$vitalStatus == 0,]))
      # dataTrainingDeadPatientsNumber = length(rownames(dataTraining[dataTraining$vitalStatus == 1,]))
      # 
      # predictedAliveToDeadPatientsRatio = dataTrainingAlivePatientsNumber/dataTrainingDeadPatientsNumber
      # predictedDeadToAlivePatientsRatio = dataTrainingDeadPatientsNumber/dataTrainingAlivePatientsNumber
      # 
      # cutoffResults = cutoffResults[cutoffResults$n1DivByn2 <= max(predictedAliveToDeadPatientsRatio, predictedDeadToAlivePatientsRatio) & cutoffResults$n1DivByn2 >= min(predictedAliveToDeadPatientsRatio, predictedDeadToAlivePatientsRatio),]
      # 
      # if(length(rownames(cutoffResults)) == 0)
      # {
      #   cutoffResults = oldCutoffResults[oldCutoffResults$Metric == max(oldCutoffResults$Metric), "cutoff"]
      #   bestCutoff = min(cutoffResults)
      # } else
      # {
      #   cutoffResults = cutoffResults[cutoffResults$Metric == max(cutoffResults$Metric), "cutoff"]
      #   bestCutoff = min(cutoffResults)
      # }
      
      
      bestCutoff = as.numeric(ROC_plot$"optimal Youden Index point"["cutoff"])
      
      if(verbose == TRUE)
      {
        print(paste("## The best cutoff value is : ", bestCutoff, " ##", sep = ""))				
      }	
      
      pathToSave = file.path("output", "class", currentSignatureName, paste("fullDataset_SurvivalCurve_Cutoff=", bestCutoff , ".pdf", sep = ""))

      grDevices::pdf(pathToSave)
      
      patientsScoreOutput[patientsScoreOutput[, columnIDToUse] >= bestCutoff, "group"] = 1
      patientsScoreOutput[patientsScoreOutput[, columnIDToUse] < bestCutoff, "group"] = 0
      

      patientsScoreOutput$group = as.numeric(patientsScoreOutput$group)
      patientsScoreOutput$survivedDays = as.numeric(patientsScoreOutput$survivedDays)
      patientsScoreOutput$vitalStatus = as.numeric(patientsScoreOutput$vitalStatus)

      survivalObjectGlobal = survival::Surv(time = patientsScoreOutput$survivedDays, event = patientsScoreOutput$vitalStatus)
      survivalCurves = survival::survfit(survivalObjectGlobal ~ group, data = patientsScoreOutput, conf.type = "log-log")
      survivalTests = survival::survdiff(survivalObjectGlobal ~ group, data = patientsScoreOutput, rho = 0)
      survivalTests_khi = survivalTests$chisq
      
      survivalTests_pValue = 1 - stats::pchisq(survivalTests$chisq, length(survivalTests$n) - 1)
      
      
      lowAndHigh_results = survRM2::rmst2(patientsScoreOutput$survivedDays, patientsScoreOutput$vitalStatus, patientsScoreOutput$group)

      RMST_low = as.numeric(lowAndHigh_results$RMST.arm0$rmst[1])
      RMST_high = as.numeric(lowAndHigh_results$RMST.arm1$rmst[1])

      finalMetric = abs(RMST_low - RMST_high)

      
      
      
      n1 = as.numeric(nrow(patientsScoreOutput[patientsScoreOutput$group == 0,]))
      n2 = as.numeric(nrow(patientsScoreOutput[patientsScoreOutput$group == 1,]))

      legendKeys = c("< to cutoff", ">= to cutoff")

     
      currentPlot = graphics::plot(
        survivalCurves,
        col = c("blue", "red"),
        xlab = "Time (days)",
        ylab = "Overall Survival",
        lty=1,
        lwd=2,
        log = FALSE,
        main = list(paste("Cutoff = ", bestCutoff, " (blue = ", n1, " et red = ", n2, ")\n",
                          "True Positive Rate (%) = ", round(as.numeric(ROC_plot$"optimal Youden Index point"["TPR"])*100,2), "\n",
                          "False Positive Rate (%) = ", round(as.numeric(ROC_plot$"optimal Youden Index point"["FPR"])*100,2), "\n",
                          
                          "RMST abs difference (days) : ", round(RMST_high, 2), " (red) - ", round(RMST_low, 2), " (blue) = ", round(finalMetric, 2), "\n",
                          "Khi2 value between curves : ", survivalTests_khi, " (p-value = ", survivalTests_pValue, ")",
                          sep = ""), cex = 0.7)
      )
      
      graphics::legend(100, 0.2, legendKeys, lty = 1:1, lwd = 2:2, col = c("blue", "red"))
     

      grDevices::dev.off()
      
      if(verbose == TRUE)
      {	
        print("## Exporting other results ##")
      }		
      
      patients_STS = patientsScoreOutput[patientsScoreOutput$group == 0, ]
      patients_LTS = patientsScoreOutput[patientsScoreOutput$group == 1, ]
      
      utils::write.table(patients_LTS, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "LTSpatients_scores.txt"), sep = "\t")
      utils::write.table(patients_STS, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "STSpatients_scores.txt"), sep = "\t")
      
      
      patientsGroupsCaseUUID = data.frame(matrix(0, nrow = max(c(nrow(patients_STS), nrow(patients_LTS))), ncol = 2), stringsAsFactors = FALSE)
      colnames(patientsGroupsCaseUUID) = c("groupSTS_CaseUUID", "groupLTS_CaseUUID")
      
      patientsGroupsCaseUUID$groupSTS_CaseUUID = c(patients_STS$CaseUUID, rep("", nrow(patientsGroupsCaseUUID) - length(patients_STS$CaseUUID)))
      patientsGroupsCaseUUID$groupLTS_CaseUUID = c(patients_LTS$CaseUUID, rep("", nrow(patientsGroupsCaseUUID) - length(patients_LTS$CaseUUID)))
      
      utils::write.table(patientsGroupsCaseUUID, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "STS-LTS_groups_caseUUID.txt"), sep = "\t")
      
      
      patientsGroupsLegend = data.frame(matrix(0, nrow = 2, ncol = 2), stringsAsFactors = FALSE)
      colnames(patientsGroupsLegend) = c("groupID", "groupLegend")
      patientsGroupsLegend[1, ] = c(0, "Short-term survivors (STS)")
      patientsGroupsLegend[2, ] = c(1, "Long-term survivors (LTS)")
      
      
      utils::write.table(patientsGroupsLegend, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "STS-LTS_groups_legend.txt"), sep = "\t")
      
      utils::write.table(patientsScoreOutput, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "fullOutputTable_allPatients.class"), sep = "\t")
      
   
      
    }
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
}




# separatePatients = function(applyFileUsed = "classification", metricToUse = "CS", verbose = TRUE)
# {
#   if(applyFileUsed == "classification")
#   {
#     signatureResultsToUseList = list.files(file.path("output", "apply"), pattern = "(.*)classificationSignature(.*).apply")
#   } else if(applyFileUsed == "custom")
#   {
#     signatureResultsToUseList = list.files(file.path("output", "class", "input"), pattern = "*.apply")
#   }
#   
#   a = NULL
#   
#   foreach::foreach(a = 1:length(signatureResultsToUseList)) %do%
#     {
#       currentSignatureResultsFile = signatureResultsToUseList[a]
#       currentSignatureName = gsub("(.+)Sign = (.+)\\)(.+)", "\\2", currentSignatureResultsFile)
#       
#       if(verbose == TRUE)
#       {	
#         print(paste("## Currently working on the class : ", currentSignatureName, " (", a, "/", length(signatureResultsToUseList), ") ##", sep = ""))
#       }		
#       
#       if(dir.exists(file.path("output", "class", currentSignatureName)) == TRUE)
#       {
#         unlink(file.path("output", "class", currentSignatureName, "*"), recursive = TRUE, force = TRUE)
#         dir.create(file.path("output", "class", currentSignatureName, "cutoffTrainingCurves"), showWarnings = TRUE, recursive = FALSE)
#       } else
#       {
#         dir.create(file.path("output", "class", currentSignatureName), showWarnings = TRUE, recursive = FALSE)
#         dir.create(file.path("output", "class", currentSignatureName, "cutoffTrainingCurves"), showWarnings = TRUE, recursive = FALSE)
#       }
#       
#       if(applyFileUsed == "classification")
#       {
#         currentSignatureResultsPath = file.path("output", "apply", currentSignatureResultsFile)
#       } else if(applyFileUsed == "custom")
#       {
#         currentSignatureResultsPath = file.path("output", "class", "input", currentSignatureResultsFile)
#       }
#       
#       currentSignatureResults = data.table::fread(input = currentSignatureResultsPath, sep = "\t", header = TRUE, data.table = FALSE)
#       
#       currentSignatureResults$vitalStatus = gsub("Alive", "0", as.character(currentSignatureResults$vitalStatus))
#       currentSignatureResults$vitalStatus = as.numeric(gsub("Dead", "1", currentSignatureResults$vitalStatus))
#       
#       if("score" %in% colnames(currentSignatureResults))
#       {
#         colIDtoChange = which(colnames(currentSignatureResults) == "score")
#         colnames(currentSignatureResults)[colIDtoChange] = "clinicalScore"
#       }
#       
#       if(metricToUse == "QS")
#       {
#         columnIDToUse = which(colnames(currentSignatureResults) == "quantitativeScore")
#       } else if (metricToUse == "CS")
#       {
#         columnIDToUse = which(colnames(currentSignatureResults) == "clinicalScore")
#       } else
#       {
#         columnIDToUse = which(colnames(currentSignatureResults) == "clinicalScore")
#       }
#       
#       if(verbose == TRUE)
#       {
#         print("## Random sampling of patients ##")
#       }
#       
#       currentSignatureResults = currentSignatureResults[sample(nrow(currentSignatureResults)),]
#       
#       options(scipen=999)
#       
#       totalCases = length(rownames(currentSignatureResults))
#       
#       if(verbose == TRUE)
#       {
#         print("## Defining of training and validation datasets ##")
#       }
#       
#       casesForTraining = as.integer(c(1:(round(totalCases/2))))
#       casesForValidation = as.integer(c((round(totalCases/2)+1):totalCases))
#       
#       dataTraining = data.frame(currentSignatureResults[casesForTraining,])
#       
#       patientsScoreOutputTraining = dataTraining
#       patientsScoreOutputTraining = patientsScoreOutputTraining[order(patientsScoreOutputTraining[,columnIDToUse]),]
#       patientsScoreOutputTraining$group = ""
#       
#       if(verbose == TRUE)
#       {
#         print("## Determining the best cutoff value possible on the training dataset ##")
#       }			
#       
#       cutoffResults = data.frame()
#       
#       if("score" %in% colnames(patientsScoreOutputTraining))
#       {
#         cutoffSequence = sort(patientsScoreOutputTraining$score)
#       } else
#       {
#         if(metricToUse == "QS")
#         {
#           cutoffSequence = sort(patientsScoreOutputTraining$quantitativeScore)
#         } else if(metricToUse == "CS")
#         {
#           cutoffSequence = sort(patientsScoreOutputTraining$clinicalScore)
#         } else
#         {
#           cutoffSequence = sort(patientsScoreOutputTraining$clinicalScore)
#         }
#       }
#       
#       o = NULL
#       
#       foreach::foreach(o = seq_len(length(cutoffSequence))) %do%	
#         {
#           cutoff = cutoffSequence[o]
#           
#           patientsScoreOutputTraining[patientsScoreOutputTraining[, columnIDToUse] >= cutoff, "group"] = 1
#           patientsScoreOutputTraining[patientsScoreOutputTraining[, columnIDToUse] < cutoff, "group"] = 0
#           
#           if(length(unique(patientsScoreOutputTraining$group)) > 1)
#           {
#             pathToSave = file.path("output", "class", currentSignatureName, "cutoffTrainingCurves", paste("survivalCurve_cutoff=", cutoff , ".pdf", sep = ""))
#             
#             pdf(pathToSave)
#             
#             patientsScoreOutputTraining$group = as.numeric(patientsScoreOutputTraining$group)
#             patientsScoreOutputTraining$survivedDays = as.numeric(patientsScoreOutputTraining$survivedDays)
#             patientsScoreOutputTraining$vitalStatus = as.numeric(patientsScoreOutputTraining$vitalStatus)
#             
#             survivalObjectGlobal = survival::Surv(time = patientsScoreOutputTraining$survivedDays, event = patientsScoreOutputTraining$vitalStatus)
#             survivalCurves = survival::survfit(survivalObjectGlobal ~ group, data = patientsScoreOutputTraining, conf.type = "log-log")
#             survivalTests = survival::survdiff(survivalObjectGlobal ~ group, data = patientsScoreOutputTraining, rho = 0)
#             
#             lowAndHigh_results = survRM2::rmst2(patientsScoreOutputTraining$survivedDays, patientsScoreOutputTraining$vitalStatus, patientsScoreOutputTraining$group)
#             
#             RMST_low = as.numeric(lowAndHigh_results$RMST.arm0$rmst[1])
#             RMST_high = as.numeric(lowAndHigh_results$RMST.arm1$rmst[1])
#             
#             finalMetric = abs(RMST_low - RMST_high)
#             
#             n1 = as.numeric(survivalTests$n[1])
#             n2 = as.numeric(survivalTests$n[2])
#             
#             alivePatientsOfGoodPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 0 & patientsScoreOutputTraining$group == 1,]))
#             alivePatientsOfBadPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 0 & patientsScoreOutputTraining$group == 0,]))
#             deadPatientsOfGoodPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 1 & patientsScoreOutputTraining$group == 1,]))
#             deadPatientsOfBadPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 1 & patientsScoreOutputTraining$group == 0,]))
#             realAlivePatientsOfGoodPrognosis = length(rownames(patientsScoreOutputTraining[patientsScoreOutputTraining$vitalStatus == 0,]))
#             
#             legendText = paste("(blue = ", n1, " et red = ", n2, ")", sep = "")
#             legendKeys = c("<= to cutoff", "> to cutoff")
#             
#             cutoffResults[o, "cutoff"] = cutoff
#             cutoffResults[o, "n1"] = n1
#             cutoffResults[o, "n2"] = n2
#             cutoffResults[o, "n1DivByn2"] = n1/n2
#             cutoffResults[o, "RMST-InfCutoff"] = RMST_low
#             cutoffResults[o, "RMST-SupCutoff"] = RMST_high
#             cutoffResults[o, "Metric"] = finalMetric
#             cutoffResults[o, "PrecisionOfPrediction"] = ((alivePatientsOfGoodPrognosis + deadPatientsOfBadPrognosis)/(n1 + n2))*100
#             cutoffResults[o, "SpecificityOfPrediction"] = (alivePatientsOfGoodPrognosis/realAlivePatientsOfGoodPrognosis)*100
#             cutoffResults[o, "PPVOfPrediction"] = (deadPatientsOfBadPrognosis/(alivePatientsOfBadPrognosis + deadPatientsOfBadPrognosis))*100
#             
#             currentPlot = graphics::plot(
#               survivalCurves,
#               col = c("blue", "red"),
#               xlab = "Time (days)",
#               ylab = "Overall Survival",
#               lty=1,
#               lwd=2,
#               log = FALSE,
#               main = paste("Cutoff : ", cutoff, sep = ""),
#             )
#             graphics::legend(100, 0.2, legendKeys, lty = 1:1, lwd = 2:2, col = c("blue", "red"))
#             graphics::mtext(legendText, side=3, adj=0.5, cex=0.8)
#             
#             grDevices::dev.off()					
#           }
#         }
#       
#       cutoffResults = cutoffResults[!is.na(cutoffResults$cutoff),]
#       oldCutoffResults = cutoffResults
#       
#       utils::write.table(cutoffResults, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "fullCutoffResults.txt"), sep = "\t")
#       
#       dataTrainingAlivePatientsNumber = length(rownames(dataTraining[dataTraining$vitalStatus == 0,]))
#       dataTrainingDeadPatientsNumber = length(rownames(dataTraining[dataTraining$vitalStatus == 1,]))
#       
#       predictedAliveToDeadPatientsRatio = dataTrainingAlivePatientsNumber/dataTrainingDeadPatientsNumber
#       predictedDeadToAlivePatientsRatio = dataTrainingDeadPatientsNumber/dataTrainingAlivePatientsNumber
#       
#       cutoffResults = cutoffResults[cutoffResults$n1DivByn2 <= max(predictedAliveToDeadPatientsRatio, predictedDeadToAlivePatientsRatio) & cutoffResults$n1DivByn2 >= min(predictedAliveToDeadPatientsRatio, predictedDeadToAlivePatientsRatio),]
#       
#       if(length(rownames(cutoffResults)) == 0)
#       {
#         cutoffResults = oldCutoffResults[oldCutoffResults$Metric == max(oldCutoffResults$Metric), "cutoff"]
#         bestCutoff = min(cutoffResults)
#       } else
#       {
#         cutoffResults = cutoffResults[cutoffResults$Metric == max(cutoffResults$Metric), "cutoff"]
#         bestCutoff = min(cutoffResults)
#       }
#       
#       if(verbose == TRUE)
#       {
#         print(paste("## The best cutoff value is : ", bestCutoff, " ##", sep = ""))				
#         print("## Applying the best cutoff value on the validation dataset ##")
#       }		
#       
#       patientsScoreOutputValidation = data.frame(currentSignatureResults[casesForValidation,])			
#       patientsScoreOutputValidation = patientsScoreOutputValidation[order(patientsScoreOutputValidation[, columnIDToUse]),]
#       patientsScoreOutputValidation$group = ""
#       
#       patientsScoreOutputValidation[patientsScoreOutputValidation[, columnIDToUse] >= bestCutoff, "group"] = 1
#       patientsScoreOutputValidation[patientsScoreOutputValidation[, columnIDToUse] < bestCutoff, "group"] = 0
#       
#       patientsScoreOutputValidation$group = as.numeric(patientsScoreOutputValidation$group)
#       patientsScoreOutputValidation$survivedDays = as.numeric(patientsScoreOutputValidation$survivedDays)
#       patientsScoreOutputValidation$vitalStatus = as.numeric(patientsScoreOutputValidation$vitalStatus)
#       
#       survivalObjectGlobal = survival::Surv(time = patientsScoreOutputValidation$survivedDays, event = patientsScoreOutputValidation$vitalStatus)
#       survivalCurves = survival::survfit(survivalObjectGlobal ~ group, data = patientsScoreOutputValidation, conf.type = "log-log")
#       
#       n1 = nrow(patientsScoreOutputValidation[patientsScoreOutputValidation$group == 0,])
#       n2 = nrow(patientsScoreOutputValidation[patientsScoreOutputValidation$group == 1,])
#       
#       legendText = paste("(blue = ", n1, " et red = ", n2, ")", sep = "")
#       legendKeys = c("<= to cutoff", "> to cutoff")
#       
#       correlationValidation = NULL
#       correlationValidation = cbind(correlationValidation, patientsScoreOutputValidation[, columnIDToUse], patientsScoreOutputValidation$survivedDays)
#       colnames(correlationValidation) = c("clinicalScore", "survivedDays")
#       
#       utils::write.table(correlationValidation, row.names = FALSE, col.names = FALSE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "correlationValidationPatients.txt"), sep = "\t")
#       
#       if(verbose == TRUE)
#       {		
#         print("## Exporting results ##")
#       }		
#       
#       bestPatientsDays = patientsScoreOutputValidation[patientsScoreOutputValidation[, columnIDToUse] >= bestCutoff, "survivedDays"]
#       worsePatientsDays = patientsScoreOutputValidation[patientsScoreOutputValidation[, columnIDToUse] < bestCutoff, "survivedDays"]
#       utils::write.table(bestPatientsDays, row.names = FALSE, col.names = FALSE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "bestPatientsDays_validationPatients.txt"), sep = "\t")
#       utils::write.table(worsePatientsDays, row.names = FALSE, col.names = FALSE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "worsePatientsDays_validationPatients.txt"), sep = "\t")
#       
#       pathToSave = file.path("output", "class", currentSignatureName, paste("survivalCurveValidation_cutoff=", bestCutoff , ".pdf", sep = ""))
#       
#       pdf(pathToSave)
#       
#       graphics::plot(
#         survivalCurves,
#         col = c("blue", "red"),
#         xlab = "Time (days)",
#         ylab = "Overall Survival",
#         lty=1,
#         lwd=2,
#         log = FALSE,
#         main = paste("Cutoff : ", bestCutoff, sep = ""),
#       )
#       
#       graphics::legend(100, 0.2, legendKeys, lty = 1:1, lwd = 2:2, col = c("blue", "red"))
#       graphics::mtext(legendText, side=3, adj=0.5, cex=0.8)
#       
#       grDevices::dev.off()
#       
#       if(verbose == TRUE)
#       {	
#         print("## Applying the best cutoff value on all available patients ##")
#       }		
#       
#       patientsScoreOutputAllPatients = currentSignatureResults
#       patientsScoreOutputAllPatients = patientsScoreOutputAllPatients[order(patientsScoreOutputAllPatients[, columnIDToUse]),]
#       patientsScoreOutputAllPatients$group = ""
#       
#       patientsScoreOutputAllPatients[patientsScoreOutputAllPatients[, columnIDToUse] >= bestCutoff, "group"] = 1
#       patientsScoreOutputAllPatients[patientsScoreOutputAllPatients[, columnIDToUse] < bestCutoff, "group"] = 0
#       
#       patientsScoreOutputAllPatients$group = as.numeric(patientsScoreOutputAllPatients$group)
#       patientsScoreOutputAllPatients$survivedDays = as.numeric(patientsScoreOutputAllPatients$survivedDays)
#       patientsScoreOutputAllPatients$vitalStatus = as.numeric(patientsScoreOutputAllPatients$vitalStatus)
#       
#       alivePatientsOfGoodPrognosis = length(rownames(patientsScoreOutputAllPatients[patientsScoreOutputAllPatients$vitalStatus == 0 & patientsScoreOutputAllPatients$group == 1,]))
#       alivePatientsOfBadPrognosis = length(rownames(patientsScoreOutputAllPatients[patientsScoreOutputAllPatients$vitalStatus == 0 & patientsScoreOutputAllPatients$group == 0,]))
#       deadPatientsOfGoodPrognosis = length(rownames(patientsScoreOutputAllPatients[patientsScoreOutputAllPatients$vitalStatus == 1 & patientsScoreOutputAllPatients$group == 1,]))
#       deadPatientsOfBadPrognosis = length(rownames(patientsScoreOutputAllPatients[patientsScoreOutputAllPatients$vitalStatus == 1 & patientsScoreOutputAllPatients$group == 0,]))
#       realAlivePatientsOfGoodPrognosis = length(rownames(patientsScoreOutputAllPatients[patientsScoreOutputAllPatients$vitalStatus == 0,]))		
#       
#       survivalObjectGlobal = survival::Surv(time = patientsScoreOutputAllPatients$survivedDays, event = patientsScoreOutputAllPatients$vitalStatus)
#       survivalCurves = survival::survfit(survivalObjectGlobal ~ group, data = patientsScoreOutputAllPatients, conf.type = "log-log")
#       
#       KhiTwoBetweenCurves = survival::survdiff(survival::Surv(survivedDays, vitalStatus) ~ group, data = patientsScoreOutputAllPatients, rho = 0)
#       KhiTwoBetweenCurves = KhiTwoBetweenCurves$chisq
#       
#       n1 = nrow(patientsScoreOutputAllPatients[patientsScoreOutputAllPatients$group == 0,])
#       n2 = nrow(patientsScoreOutputAllPatients[patientsScoreOutputAllPatients$group == 1,])
#       
#       lowAndHigh_results = survRM2::rmst2(patientsScoreOutputAllPatients$survivedDays, patientsScoreOutputAllPatients$vitalStatus, patientsScoreOutputAllPatients$group)
#       
#       RMST_low = as.numeric(lowAndHigh_results$RMST.arm0$rmst[1])
#       RMST_high = as.numeric(lowAndHigh_results$RMST.arm1$rmst[1])
#       
#       finalMetric = abs(RMST_low - RMST_high)
#       
#       legendText = paste("(blue = ", n1, " et red = ", n2, ")", sep = "")
#       legendKeys = c("<= to cutoff", "> to cutoff")
#       
#       correlationAllPatients = NULL
#       correlationAllPatients = cbind(correlationAllPatients, patientsScoreOutputAllPatients[, columnIDToUse], patientsScoreOutputAllPatients$survivedDays)
#       colnames(correlationAllPatients) = c("clinicalScore", "survivedDays")
#       
#       utils::write.table(correlationAllPatients, row.names = FALSE, col.names = FALSE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "correlationAllPatients.txt"), sep = "\t")
#       
#       if(verbose == TRUE)
#       {	
#         print("## Exporting results ##")
#       }		
#       
#       bestPatientsDays = patientsScoreOutputAllPatients[patientsScoreOutputAllPatients[, columnIDToUse] >= bestCutoff, "survivedDays"]
#       worsePatientsDays = patientsScoreOutputAllPatients[patientsScoreOutputAllPatients[, columnIDToUse] < bestCutoff, "survivedDays"]
#       utils::write.table(bestPatientsDays, row.names = FALSE, col.names = FALSE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "bestPatientsDays_allPatients.txt"), sep = "\t")
#       utils::write.table(worsePatientsDays, row.names = FALSE, col.names = FALSE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "worsePatientsDays_allPatients.txt"), sep = "\t")
#       
#       pathToSave = file.path("output", "class", currentSignatureName, paste("survivalCurveAllPatients_cutoff=", bestCutoff , ".pdf", sep = ""))
#       
#       pdf(pathToSave)
#       
#       graphics::par(mar=c(4,5,8,2))
#       graphics::plot(
#         survivalCurves,
#         col = c("blue", "red"),
#         xlab = "Time (days)",
#         ylab = "Overall Survival",
#         lty=1,
#         lwd=2,
#         log = FALSE,
#         main = list(paste("Cutoff = ", bestCutoff , "\n", 
#                           "(blue = ", n1, " et red = ", n2, ")\n",
#                           "Blue : <= cutoff & red : > cutoff\n",
#                           "Precision of prediction (%) = ", round(((alivePatientsOfGoodPrognosis + deadPatientsOfBadPrognosis)/(length(rownames(patientsScoreOutputAllPatients))))*100,2), "\n",
#                           "Specificity of prediction (%) = ", round((alivePatientsOfGoodPrognosis/realAlivePatientsOfGoodPrognosis)*100,2), "\n",
#                           "PPV of prediction (%) = ", round((deadPatientsOfBadPrognosis/(alivePatientsOfBadPrognosis + deadPatientsOfBadPrognosis))*100, 2), "\n",
#                           "RMST abs difference (days) : ", round(RMST_high, 2), " (red) - ", round(RMST_low, 2), " (blue) = ", round(finalMetric, 2), "\n",
#                           "Khi2 value between curves : ", KhiTwoBetweenCurves,
#                           sep = ""), cex = 0.8)
#       )
#       
#       grDevices::dev.off()
#       
#       utils::write.table(patientsScoreOutputAllPatients, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "fullOutputTable_allPatients.class"), sep = "\t")
#       
#       allPatientsStatistics = data.frame()
#       allPatientsStatistics[1, "PrecisionOfPrediction"] = ((alivePatientsOfGoodPrognosis + deadPatientsOfBadPrognosis)/(n1 + n2))*100
#       allPatientsStatistics[1, "SpecificityOfPrediction"] = (alivePatientsOfGoodPrognosis/realAlivePatientsOfGoodPrognosis)*100
#       allPatientsStatistics[1, "PPVOfPrediction"] = (deadPatientsOfBadPrognosis/(alivePatientsOfBadPrognosis + deadPatientsOfBadPrognosis))*100
#       allPatientsStatistics[1, "RMST-InfToCutoff"] = RMST_low
#       allPatientsStatistics[1, "RMST-SupToCutoff"] = RMST_high
#       
#       utils::write.table(allPatientsStatistics, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "class", currentSignatureName, "finalCurvesStatistics_allPatients.txt"), sep = "\t")
#       
#       if(verbose == TRUE)
#       {		
#         print("## Zipping the cutoffTrainingCurves folder ##")
#       }
#       
#       previouswd = getwd()
#       setwd(file.path("output", "class", currentSignatureName))
#       filesToZip = list.files(path = "cutoffTrainingCurves", pattern = "*", recursive = TRUE)
#       zip::zipr(zipfile = "cutoffTrainingCurves.zip", files = file.path("cutoffTrainingCurves", filesToZip))
#       
#       unlink(file.path("cutoffTrainingCurves", filesToZip), force = TRUE, recursive = TRUE)
#       directoriesToDelete = gsub("\\./", "", list.dirs()[-1])
#       unlink(directoriesToDelete, force = TRUE, recursive = TRUE)
#       setwd(previouswd)
#       
#     }
#   
#   sink(file.path("silentOutput.log"))
#   gc()
#   sink()
# }