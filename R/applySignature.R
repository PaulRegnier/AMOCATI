#' Compute the Quantitative Scores and the Clinical Scores for each patient and for one or several gene signature(s)
#'
#' This function computes the Quantitative Scores and the Clinical Scores individually for each patient of a whole dataset, either for the Classification Signature or for one or several custom gene signatures.
#'
#' @param signatureUsed A character vector which tells rather to use the Classification Signature (if set to `classification`) or one or several custom gene signatures (if set to `custom`). The Classification Signature is read from the `output > signatures > classificationSignature.sign` text file, whereas the custom signatures are read from the `output > signatures > customSignatures.sign` text file. 
#'
#' @param distinguishingGenes A boolean determining if one or several gene(s) should be used (if set to `TRUE`) to further discriminate the patients according to their relative low, intermediate or high expression of this/these gene(s).
#'
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @return Generated output files are exported to the `output > apply` folder. Their precise amounts and organization varies according to the `signatureUsed` and `distinguishingGenes` arguments (see the online tutorial for more details).
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#'
#' @export

applySignature = function(signatureUsed = "classification", distinguishingGenes = FALSE, verbose = TRUE)
{
  basisFolder = getwd()
  
  lm_eqn = function(df)
  {
    m = stats::lm(y ~ x, df)
    a = format(stats::coef(m)[1], digits = 2)
    b = format(stats::coef(m)[2], digits = 2)
    r2 = format(summary(m)$r.squared, digits = 4)
    
    pearsonCorrelation = stats::cor.test(df$x, df$y, method = "pearson")
    pearsonScore = round(as.numeric(pearsonCorrelation$estimate), 3)
    pearsonPvalue = format(as.numeric(pearsonCorrelation$p.value), scientific = TRUE, digits = 2)
    spearmanCorrelation = stats::cor.test(df$x, df$y, method = "spearman")
    spearmanScore = round(as.numeric(spearmanCorrelation$estimate), 3)
    spearmanPvalue = format(as.numeric(spearmanCorrelation$p.value), scientific = TRUE, digits = 2)
    
    return(list(b = as.numeric(b), a = as.numeric(a), r2 = as.numeric(r2), pearsonScore = as.numeric(pearsonScore), pearsonPvalue = as.numeric(pearsonPvalue), spearmanScore = as.numeric(spearmanScore), spearmanPvalue = as.numeric(spearmanPvalue)))
  }
  
  if(verbose == TRUE)
  {
    print("## Opening metaResults file ##")
  }
  
  fileToLoad = list.files(file.path("output", "metaResults"), pattern = "*.meta")[1]
  metaResults_allGenesPooledPath = file.path("output", "metaResults", fileToLoad)
  metaResults_allGenesPooled = data.table::fread(input = metaResults_allGenesPooledPath, sep = "\t", header = TRUE, data.table = FALSE)
  metaResults_filtered = metaResults_allGenesPooled[order(metaResults_allGenesPooled$GeneScore, decreasing = TRUE), ]
  rownames(metaResults_filtered) = NULL
  genesRanksTotal = length(which((as.character(metaResults_filtered$AssociatedSetting) != "N/A") == TRUE))
  
  if(signatureUsed == "classification")
  {
    distinguishingGenes = FALSE
  }
  
  if(distinguishingGenes == FALSE)
  {
    distinguishingGenesMode = FALSE
    if(signatureUsed == "classification")
    {
      distinguishingGenesDf = data.frame(HGNC_GeneSymbol = "classification", GeneSymbol = "classification", stringsAsFactors = FALSE)
    } else if(signatureUsed == "custom")
    {
      distinguishingGenesDf = data.frame(HGNC_GeneSymbol = "noDG", GeneSymbol = "noDG", stringsAsFactors = FALSE)
    }
  } else
  {	
    distinguishingGenesMode = TRUE
    distinguishingGenesFileToOpen = list.files(file.path("output", "apply", "input"), pattern = "*.txt")[1]
    distinguishingGenesDf = data.table::fread(input = file.path("output", "apply", "input", distinguishingGenesFileToOpen), sep="\t", header=TRUE, data.table = FALSE)
  }
  
  distinguishingGenesList = as.vector(distinguishingGenesDf$HGNC_GeneSymbol)
  
  if(verbose == TRUE)
  {
    print("## Opening fullData file ##")
  }
  
  fileToLoad = list.files(file.path("output", "data"), pattern = "*.data")[1]
  fullGenesDataPath = file.path("output", "data", fileToLoad)
  fullGenesData = data.table::fread(input = fullGenesDataPath, sep = "\t", header = TRUE, data.table = FALSE)
  
  genesNumber = length(colnames(fullGenesData))-3
  totalCases = length(rownames(fullGenesData))
  
  if(signatureUsed == "classification")
  {
    if(verbose == TRUE)
    {
      print("## Importing classification signature ##")
    }
    
    signaturesToMatchFile = file.path("output", "signatures", "classificationSignature.sign")
    signaturesToMatch = data.table::fread(input = signaturesToMatchFile, sep="\t", header=TRUE, data.table = FALSE)	
  } else if(signatureUsed == "custom")
  {
    if(verbose == TRUE)
    {
      print("## Importing custom genes signatures ##")
    }
    
    signaturesToMatchFile = file.path("output", "signatures", "customSignatures.sign")
    signaturesToMatch = data.table::fread(input = signaturesToMatchFile, sep="\t", header=TRUE, data.table = FALSE)
    
    if(verbose == TRUE)
    {
      print("## Emptying output folder ##")
    }
    
    if(distinguishingGenes == TRUE)
    {
      unlink(file.path("output", "apply", "customSignatures", "CS"), recursive = TRUE, force = TRUE)
      unlink(file.path("output", "apply", "customSignatures", "QS"), recursive = TRUE, force = TRUE)
      
      zipFilesToErase = list.files(file.path("output", "apply", "customSignatures", "fullTables"), pattern = "*.zip")
      unlink(file.path("output", "apply", "customSignatures", "fullTables", zipFilesToErase), force = TRUE)
      
      unlink(file.path("output", "apply", "customSignatures", "CS VS QS"), recursive = TRUE, force = TRUE)
      
      synthesisFolderToErase = list.files(file.path("output", "apply", "customSignatures", "synthesis"))
      fileToKeep = grep("\\.txt", synthesisFolderToErase)
      synthesisFolderToErase = synthesisFolderToErase[-fileToKeep]
      
      y = NULL
      
      foreach::foreach(y = 1:length(synthesisFolderToErase)) %do%
        {
          unlink(file.path("output", "apply", "customSignatures", "synthesis", synthesisFolderToErase[y]), recursive = TRUE, force = TRUE)
        }
      
      dir.create(file.path("output", "apply", "customSignatures", "CS"))
      dir.create(file.path("output", "apply", "customSignatures", "QS"))
      dir.create(file.path("output", "apply", "customSignatures", "CS VS QS"))
      # dir.create(file.path("output", "apply", "customSignatures", "synthesis"))
    } else
    {
      unlink(file.path("output", "apply", "customSignatures", "fullTables", "noDG"), recursive = TRUE, force = TRUE)
      dir.create(file.path("output", "apply", "customSignatures", "fullTables", "noDG"))
    }
  }
  
  g = NULL
  
  foreach::foreach(g = 1:length(distinguishingGenesList)) %do%
    {
      distinguishingGene = distinguishingGenesList[g]
      
      if(signatureUsed == "custom" & distinguishingGenesMode == TRUE)
      {
        dir.create(file.path("output", "apply", "customSignatures", "CS", distinguishingGene))
        dir.create(file.path("output", "apply", "customSignatures", "CS", distinguishingGene, "graphs"))
        dir.create(file.path("output", "apply", "customSignatures", "CS", distinguishingGene, "stats"))
        
        dir.create(file.path("output", "apply", "customSignatures", "QS", distinguishingGene))
        dir.create(file.path("output", "apply", "customSignatures", "QS", distinguishingGene, "graphs"))
        dir.create(file.path("output", "apply", "customSignatures", "QS", distinguishingGene, "stats"))
        
        dir.create(file.path("output", "apply", "customSignatures", "CS VS QS", distinguishingGene))
        dir.create(file.path("output", "apply", "customSignatures", "CS VS QS", distinguishingGene, "graphs"))
        dir.create(file.path("output", "apply", "customSignatures", "CS VS QS", distinguishingGene, "stats"))
        
        dir.create(file.path("output", "apply", "customSignatures", "synthesis", distinguishingGene))
        dir.create(file.path("output", "apply", "customSignatures", "fullTables", distinguishingGene))
      }
      
      currentGenePooledResults_rownames = c("CS_LowMean", "CS_InterMean", "CS_HighMean", "QS_LowMean", "QS_InterMean", "QS_HighMean", "LowVSInter_CS", "LowVSInter_QS", "InterVSHigh_CS", "InterVSHigh_QS", "LowVSHigh_CS", "LowVSHigh_QS", "CS-VS-QS_Pearson", "CS-VS-QS_PearsonPvalue", "CS-VS-QS_Spearman", "CS-VS-QS_SpearmanPvalue", "CS-VS-QS_a", "CS-VS-QS_b", "CS-VS-QS_r2")
      
      if(verbose == TRUE)
      {
        if(distinguishingGenesMode == FALSE)
        {
          print("## Analysis ran without distinguishing gene ##")
          
        } else
        {
          print(paste("## Distinguishing gene (", distinguishingGene, ") : ", g, "/", length(distinguishingGenesList), " ##", sep = ""))
          
        }
      }
      
      coresNumber = parallel::detectCores()
      cl = parallel::makeCluster(coresNumber)
      doSNOW::registerDoSNOW(cl)
      
      ntasks = length(colnames(signaturesToMatch))
      pb = tcltk::tkProgressBar(max=ntasks)
      progress = function(n) tcltk::setTkProgressBar(pb, n, label = paste(round(n/ntasks*100, 0), "% done"))
      opts = list(progress=progress)
      
      z = NULL
      
      currentGenePooledResults = foreach::foreach(z = seq_len(length(colnames(signaturesToMatch))), .packages = c("foreach", "survival", "tcltk", "data.table", "ggplot2"), .combine=cbind, .options.snow = opts) %dopar%
        {
          currentGeneAllResults = data.frame(matrix("", nrow = length(currentGenePooledResults_rownames), ncol = 1), stringsAsFactors = FALSE)
          signatureToMatch = as.vector(signaturesToMatch[,z])
          signatureToMatch = signatureToMatch[signatureToMatch != ""]
          currentSignatureName = as.vector(colnames(signaturesToMatch)[z])
          
          colnames(currentGeneAllResults)[1] = currentSignatureName
          rownames(currentGeneAllResults) = currentGenePooledResults_rownames
          
          resultingMetaResults_signatureGenes = metaResults_filtered[metaResults_filtered$HGNC_GeneSymbol %in% signatureToMatch,]
          resultingGenesData_signatureGenes = cbind(fullGenesData[,1:3], fullGenesData[,colnames(fullGenesData) %in% signatureToMatch])
          rownames(resultingGenesData_signatureGenes) = NULL
          
          if(distinguishingGenesMode == FALSE)
          {
            resultingGenesData_distinguishGene = fullGenesData[,1:3]
            resultingGenesData_distinguishGene$NO_DG = ""
          } else
          {
            resultingMetaResults_distinguishGene = metaResults_filtered[metaResults_filtered$HGNC_GeneSymbol == distinguishingGene,]
            tempDataFrame = NULL
            
            resultingGenesData_distinguishGene = cbind(fullGenesData[,1:3], fullGenesData[,colnames(fullGenesData) == distinguishingGene])
            
            rownames(resultingGenesData_distinguishGene) = NULL
            colnames(resultingGenesData_distinguishGene)[4] = distinguishingGene
            
            resultingGenesData_distinguishGene$group = "Inter"
            
            currentGeneSetting = as.vector(resultingMetaResults_distinguishGene$AssociatedSetting)
            currentCutoffToUse = as.numeric(resultingMetaResults_distinguishGene$Cutoff)
            
            if(length(currentCutoffToUse) == 0)
            {
              currentCutoffToUse = 15
            }
            
            currentGeneData = resultingGenesData_distinguishGene[,colnames(resultingGenesData_distinguishGene) == distinguishingGene]		
            currentQuantiles = stats::quantile(as.numeric(currentGeneData), probs = c(currentCutoffToUse/100, 1-currentCutoffToUse/100), type = 7, na.rm = TRUE)
            
            lowThreshold = as.numeric(currentQuantiles[1])
            highThreshold = as.numeric(currentQuantiles[2])
            
            rowsPatientsHIGH = which(currentGeneData >= highThreshold)
            rowsPatientsLOW = which(currentGeneData <= lowThreshold)
            
            resultingGenesData_distinguishGene$group[rowsPatientsHIGH] = "High"
            resultingGenesData_distinguishGene$group[rowsPatientsLOW] = "Low"
          }
          
          currentSignatureRealLength = 0
          currentSignatureTotalClinicalScore = 0
          currentSignatureTotalQuantitativeScore = 0
          
          b = NULL
          
          foreach::foreach(b = seq_len(length(signatureToMatch))) %do%
            {
              currentGeneOfSignature = signatureToMatch[b]
              currentGeneMetaResults = resultingMetaResults_signatureGenes[resultingMetaResults_signatureGenes$HGNC_GeneSymbol == currentGeneOfSignature,]
              
              currentGeneSetting = as.vector(currentGeneMetaResults$AssociatedSetting)
              currentGeneData = resultingGenesData_signatureGenes[,colnames(resultingGenesData_signatureGenes) == currentGeneOfSignature]
              
              if(length(currentGeneSetting) > 0 & length(currentGeneData) > 0)
              {
                currentCutoffToUse = as.numeric(currentGeneMetaResults$Cutoff)
                
                if(is.na(currentCutoffToUse) == TRUE)
                {
                  currentCutoffToUse = 15
                }
                
                currentQuantiles = stats::quantile(as.numeric(currentGeneData), probs = c(currentCutoffToUse/100, 1-currentCutoffToUse/100), type = 7, na.rm = TRUE)
                
                lowThreshold = as.numeric(currentQuantiles[1])
                highThreshold = as.numeric(currentQuantiles[2])
                
                currentSignatureRealLength = currentSignatureRealLength + 1
                
                if(currentGeneSetting == "correlation")
                {
                  highSurvivalPatients = which(currentGeneData >= highThreshold)
                  interSurvivalPatients = which(currentGeneData < highThreshold & currentGeneData > lowThreshold)
                  lowSurvivalPatients = which(currentGeneData <= lowThreshold)
                } else if(currentGeneSetting == "inverseCorrelation")
                {
                  lowSurvivalPatients = which(currentGeneData >= highThreshold)
                  interSurvivalPatients = which(currentGeneData < highThreshold & currentGeneData > lowThreshold)
                  highSurvivalPatients = which(currentGeneData <= lowThreshold)
                } else if(currentGeneSetting == "paradox")
                {
                  highSurvivalPatients = which(currentGeneData >= highThreshold | currentGeneData <= lowThreshold)
                  interSurvivalPatients = NULL
                  lowSurvivalPatients = which(currentGeneData < highThreshold & currentGeneData > lowThreshold)
                } else if(currentGeneSetting == "inverseParadox")
                {
                  lowSurvivalPatients = which(currentGeneData >= highThreshold | currentGeneData <= lowThreshold)
                  interSurvivalPatients = NULL
                  highSurvivalPatients = which(currentGeneData < highThreshold & currentGeneData > lowThreshold)
                } else
                {
                  randomParameter = sample(1:4, 1)
                  
                  if(randomParameter == 1)
                  {
                    highSurvivalPatients = which(currentGeneData >= highThreshold)
                    interSurvivalPatients = which(currentGeneData < highThreshold & currentGeneData > lowThreshold)
                    lowSurvivalPatients = which(currentGeneData <= lowThreshold)
                  } else if(randomParameter == 2)
                  {
                    lowSurvivalPatients = which(currentGeneData >= highThreshold)
                    interSurvivalPatients = which(currentGeneData < highThreshold & currentGeneData > lowThreshold)
                    highSurvivalPatients = which(currentGeneData <= lowThreshold)
                  } else if(randomParameter == 3)
                  {
                    highSurvivalPatients = which(currentGeneData >= highThreshold | currentGeneData <= lowThreshold)
                    interSurvivalPatients = NULL
                    lowSurvivalPatients = which(currentGeneData < highThreshold & currentGeneData > lowThreshold)
                  } else if(randomParameter == 4)
                  {
                    lowSurvivalPatients = which(currentGeneData >= highThreshold | currentGeneData <= lowThreshold)
                    interSurvivalPatients = NULL
                    highSurvivalPatients = which(currentGeneData < highThreshold & currentGeneData > lowThreshold)
                  }
                }
                
                currentSignatureGeneQuantitativeScore = currentGeneData
                
                currentSignatureGeneClinicalScore = 1*currentGeneData
                
                currentSignatureGeneClinicalScore[highSurvivalPatients] = 1*currentSignatureGeneClinicalScore[highSurvivalPatients]
                currentSignatureGeneClinicalScore[lowSurvivalPatients] = -1*currentSignatureGeneClinicalScore[lowSurvivalPatients]
                currentSignatureGeneClinicalScore[interSurvivalPatients] = 0
                
                currentSignatureTotalClinicalScore = currentSignatureTotalClinicalScore + currentSignatureGeneClinicalScore
                currentSignatureTotalQuantitativeScore = currentSignatureTotalQuantitativeScore + currentSignatureGeneQuantitativeScore
              }
            }
          
          clinicalScore = NULL
          quantitativeScore = NULL
          group = NULL
          
          resultingGenesData_distinguishGene[,"clinicalScore"] = currentSignatureTotalClinicalScore/currentSignatureRealLength
          resultingGenesData_distinguishGene[,"quantitativeScore"] = currentSignatureTotalQuantitativeScore/currentSignatureRealLength
          
          resultingGenesData_distinguishGene = resultingGenesData_distinguishGene[order(resultingGenesData_distinguishGene$CaseUUID),]
          resultingGenesData_distinguishGene = resultingGenesData_distinguishGene[!is.na(resultingGenesData_distinguishGene$quantitativeScore),]
          resultingGenesData_distinguishGene = resultingGenesData_distinguishGene[!is.na(resultingGenesData_distinguishGene$clinicalScore),]
          
          if(signatureUsed == "custom" & distinguishingGenesMode == FALSE)
          {				
            currentRawDataPlotting = as.data.frame(cbind(x = resultingGenesData_distinguishGene$clinicalScore, y = resultingGenesData_distinguishGene$quantitativeScore))
            currentGeomTextResult = lm_eqn(currentRawDataPlotting)
            
            resultsToExport = t(as.data.frame(unlist(currentGeomTextResult)))
            rownames(resultsToExport) = NULL
            resultsToExport = as.data.frame(resultsToExport[,colnames(resultsToExport) != "label"])
            colnames(resultsToExport) = "results"
            rownames(resultsToExport) = c("a", "b", "r2", "pearsonScore", "pearsonPvalue", "spearmanScore", "spearmanPvalue")
            
            currentGeneAllResults[c("CS-VS-QS_Pearson", "CS-VS-QS_PearsonPvalue", "CS-VS-QS_Spearman", "CS-VS-QS_SpearmanPvalue", "CS-VS-QS_a", "CS-VS-QS_b", "CS-VS-QS_r2"), 1] = c(currentGeomTextResult$pearsonScore, currentGeomTextResult$pearsonPvalue, currentGeomTextResult$spearmanScore, currentGeomTextResult$spearmanPvalue, currentGeomTextResult$a, currentGeomTextResult$b, currentGeomTextResult$r2)
            
            currentGeneAllResults[c("CS_LowMean", "CS_InterMean", "CS_HighMean", "LowVSInter_CS", "InterVSHigh_CS", "LowVSHigh_CS"), 1] = c("", "", "", "", "", "")
          } else if(signatureUsed == "custom" & distinguishingGenesMode == TRUE)
          {
            currentRawDataPlotting = as.data.frame(cbind(x = resultingGenesData_distinguishGene$clinicalScore, y = resultingGenesData_distinguishGene$quantitativeScore))
            currentGeomTextResult = lm_eqn(currentRawDataPlotting)
            
            resultsToExport = t(as.data.frame(unlist(currentGeomTextResult)))
            rownames(resultsToExport) = NULL
            resultsToExport = as.data.frame(resultsToExport[,colnames(resultsToExport) != "label"])
            colnames(resultsToExport) = "results"
            rownames(resultsToExport) = c("a", "b", "r2", "pearsonScore", "pearsonPvalue", "spearmanScore", "spearmanPvalue")
            
            currentGeneAllResults[c("CS-VS-QS_Pearson", "CS-VS-QS_PearsonPvalue", "CS-VS-QS_Spearman", "CS-VS-QS_SpearmanPvalue", "CS-VS-QS_a", "CS-VS-QS_b", "CS-VS-QS_r2"), 1] = c(currentGeomTextResult$pearsonScore, currentGeomTextResult$pearsonPvalue, currentGeomTextResult$spearmanScore, currentGeomTextResult$spearmanPvalue, currentGeomTextResult$a, currentGeomTextResult$b, currentGeomTextResult$r2)
            
            utils::write.table(resultsToExport, row.names = TRUE, col.names = TRUE, quote = FALSE, file = file.path(basisFolder, "output", "apply", "customSignatures", "CS VS QS", distinguishingGene, "stats", paste(distinguishingGene, "_plot-QS-VS-CS_stats_(Sign = ", currentSignatureName, ").txt", sep = "")), sep = "\t")
            
            currentPlot = ggplot2::ggplot(resultingGenesData_distinguishGene, ggplot2::aes(clinicalScore, quantitativeScore)) +
              ggplot2::geom_point(ggplot2::aes(colour = factor(group)), size = 2) +
              ggplot2::scale_color_manual(values=c("Low"="green", "Inter"="black","High"="red")) +
              ggplot2::theme(axis.line = ggplot2::element_line(colour = "black", size = 1, linetype = "solid"), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white"), legend.position="none", aspect.ratio=1, axis.text = ggplot2::element_text(face="bold", color="black", size=10), axis.ticks = ggplot2::element_line(size = 1)) +
              ggplot2::geom_vline(xintercept = 0, na.rm = FALSE, show.legend = NA, color = "#6768AB", size = 1) +
              ggplot2::xlab("Clinical Score") +
              ggplot2::ylab("Quantitative Score") +
              ggplot2::geom_smooth(method="lm", color = "#6768AB", fill = "#E8E8FD")
            
            ggplot2::ggsave(plot = currentPlot, path = file.path(basisFolder, "output", "apply", "customSignatures", "CS VS QS", distinguishingGene, "graphs"), height = 12, width = 12, device = "pdf", filename = paste(distinguishingGene, "_plot-QS-VS-CS_(Sign = ", currentSignatureName, ").pdf", sep = ""))
            
            resultingGenesData_distinguishGene$group = factor(resultingGenesData_distinguishGene$group, levels = c("Low", "Inter", "High"))
            
            settingsToExport = c("clinicalScore", "quantitativeScore")
            
            l = NULL
            
            foreach::foreach(l = 1:length(settingsToExport)) %do%
              {
                currentSetting = settingsToExport[l]
                
                if(currentSetting == "clinicalScore")
                {
                  exportShortcut = "CS"
                } else
                {
                  exportShortcut = "QS"
                }
                
                scoresLow = resultingGenesData_distinguishGene[resultingGenesData_distinguishGene$group == "Low", currentSetting]
                scoresInter = resultingGenesData_distinguishGene[resultingGenesData_distinguishGene$group == "Inter", currentSetting]
                scoresHigh = resultingGenesData_distinguishGene[resultingGenesData_distinguishGene$group == "High", currentSetting]
                shapiroPValue = stats::shapiro.test(c(scoresLow, scoresInter, scoresHigh))$p.value
                
                scoresLow_mean = mean(scoresLow)
                scoresInter_mean = mean(scoresInter)
                scoresHigh_mean = mean(scoresHigh)
                
                scoresLow_median = stats::median(scoresLow)
                scoresInter_median = stats::median(scoresInter)
                scoresHigh_median = stats::median(scoresHigh)
                
                scoresLow_sem = stats::sd(scoresLow)/sqrt(length(scoresLow))
                scoresInter_sem = stats::sd(scoresInter)/sqrt(length(scoresInter))
                scoresHigh_sem = stats::sd(scoresHigh)/sqrt(length(scoresHigh))
                
                if(shapiroPValue > 0.05)
                {
                  pValue_LowVSInter = stats::t.test(scoresLow, scoresInter)$p.value
                  pValue_InterVSHigh = stats::t.test(scoresInter, scoresHigh)$p.value
                  pValue_LowVSHigh = stats::t.test(scoresLow, scoresHigh)$p.value
                  testUsed = "t.test"
                } else
                {
                  pValue_LowVSInter = stats::wilcox.test(scoresLow, scoresInter)$p.value
                  pValue_InterVSHigh = stats::wilcox.test(scoresInter, scoresHigh)$p.value
                  pValue_LowVSHigh = stats::wilcox.test(scoresLow, scoresHigh)$p.value
                  testUsed = "wilcox.test"
                }
                
                if(length(scoresLow) >= 3 & length(scoresInter) >= 3)
                {
                  if(shapiroPValue > 0.05)
                  {
                    pValue_LowVSInter = stats::t.test(scoresLow, scoresInter)$p.value
                    testUsed = "t.test"
                  } else
                  {
                    pValue_LowVSInter = stats::wilcox.test(scoresLow, scoresInter)$p.value
                    testUsed = "wilcox.test"
                  }
                } else
                {
                  pValue_LowVSInter = ""
                  testUsed = "none"
                }
                
                if(length(scoresInter) >= 3 & length(scoresHigh) >= 3)
                {
                  if(shapiroPValue > 0.05)
                  {
                    pValue_InterVSHigh = stats::t.test(scoresInter, scoresHigh)$p.value
                    testUsed = "t.test"
                  } else
                  {
                    pValue_InterVSHigh = stats::wilcox.test(scoresInter, scoresHigh)$p.value
                    testUsed = "wilcox.test"
                  }
                } else
                {
                  pValue_InterVSHigh = ""
                  testUsed = "none"
                }
                
                if(length(scoresLow) >= 3 & length(scoresHigh) >= 3)
                {
                  if(shapiroPValue > 0.05)
                  {
                    pValue_LowVSHigh = stats::t.test(scoresLow, scoresHigh)$p.value
                    testUsed = "t.test"
                  } else
                  {
                    pValue_LowVSHigh = stats::wilcox.test(scoresLow, scoresHigh)$p.value
                    testUsed = "wilcox.test"
                  }
                } else
                {
                  pValue_LowVSHigh = ""
                  testUsed = "none"
                }
                
                if(length(scoresLow) >= 1)
                {
                  scoresLow_mean = mean(scoresLow)
                  scoresLow_median = stats::median(scoresLow)
                } else
                {
                  scoresLow_mean = ""
                  scoresLow_median = ""
                }
                
                if(length(scoresInter) >= 1)
                {
                  scoresInter_mean = mean(scoresInter)
                  scoresInter_median = stats::median(scoresInter)
                } else
                {
                  scoresInter_mean = ""
                  scoresInter_median = ""
                }
                
                if(length(scoresHigh) >= 1)
                {
                  scoresHigh_mean = mean(scoresHigh)
                  scoresHigh_median = stats::median(scoresHigh)
                } else
                {
                  scoresHigh_mean = ""
                  scoresHigh_median = ""
                }
                
                data_export = data.frame(t(data.frame(c("Comparison =>", "Low VS Inter", "Inter VS High", "Low VS High"), stringsAsFactors = FALSE)), stringsAsFactors = FALSE)
                rownames(data_export) = NULL
                colnames(data_export) = NULL
                data_export = rbind(data_export, c("P-values =>", pValue_LowVSInter, pValue_InterVSHigh, pValue_LowVSHigh))
                data_export = rbind(data_export, c("Group =>", "Low", "Inter", "High"))
                data_export = rbind(data_export, c("Mean =>", scoresLow_mean, scoresInter_mean, scoresHigh_mean))
                data_export = rbind(data_export, c("Median =>", scoresLow_median, scoresInter_median, scoresHigh_median))
                data_export = rbind(data_export, c("Test used =>", testUsed, testUsed, testUsed))
                
                currentGeneAllResults[c("QS_LowMean", "QS_InterMean", "QS_HighMean", "LowVSInter_QS", "InterVSHigh_QS", "LowVSHigh_QS"), 1] = c(scoresLow_mean, scoresInter_mean, scoresHigh_mean, pValue_LowVSInter, pValue_InterVSHigh, pValue_LowVSHigh)
                
                utils::write.table(data_export, row.names = FALSE, col.names = FALSE, quote = FALSE, file = file.path(basisFolder, "output", "apply", "customSignatures", exportShortcut, distinguishingGene, "stats", paste(distinguishingGene, "_bars-", exportShortcut, "_stats_(Sign = ", currentSignatureName, ").txt", sep = "")), sep = "\t")
                
                summary_group = rbind("Low", "Inter", "High")
                summary_mean = rbind(scoresLow_mean, scoresInter_mean, scoresHigh_mean)
                summary_sem = rbind(scoresLow_sem, scoresInter_sem, scoresHigh_sem)
                
                scores_summary = data.frame(cbind(summary_group, summary_mean, summary_sem))
                
                rownames(scores_summary) = NULL
                SEM = NULL
                colnames(scores_summary) = c("group", "mean", "SEM")
                scores_summary$group = as.factor(as.character(scores_summary$group))
                
                scores_summary$mean = as.numeric(as.character(scores_summary$mean))
                scores_summary$SEM = as.numeric(as.character(scores_summary$SEM))
                
                scores_summary$group = factor(scores_summary$group,as.character(scores_summary$group))
                
                currentPlot = ggplot2::ggplot(scores_summary, ggplot2::aes(y=mean, x=group)) +
                  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - SEM, ymax = mean + SEM), width=0.2, size = 1, color = c("Low"="green", "Inter"="black","High"="red")) +
                  ggplot2::geom_point(ggplot2::aes(color = group), stat="identity", size=4) +
                  ggplot2::scale_color_manual(values=c("Low"="green", "Inter"="black","High"="red")) +
                  ggplot2::xlab("Group") +
                  ggplot2::ylab(currentSetting) +
                  ggplot2::theme(axis.line = ggplot2::element_line(colour = "black", size = 1, linetype = "solid"), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white"), legend.position="none", aspect.ratio=2, axis.text = ggplot2::element_text(face="bold", color="black", size=10), axis.ticks = ggplot2::element_line(size = 1)) +
                  if(currentSetting == "clinicalScore")
                  {
                    ggplot2::geom_hline(yintercept = 0, na.rm = FALSE, show.legend = NA, color = "#6768AB", size = 1)
                  }
                
                ggplot2::ggsave(plot = currentPlot, path = file.path(basisFolder, "output", "apply", "customSignatures", exportShortcut, distinguishingGene, "graphs"), height = 12, width = 12, device = "pdf", filename = paste(distinguishingGene, "_bars-", exportShortcut, "_(Sign = ", currentSignatureName, ").pdf", sep = ""))
                
                if(l == 1)
                {
                  currentGeneAllResults[c("CS_LowMean", "CS_InterMean", "CS_HighMean", "LowVSInter_CS", "InterVSHigh_CS", "LowVSHigh_CS"), 1] = c(scoresLow_mean, scoresInter_mean, scoresHigh_mean, pValue_LowVSInter, pValue_InterVSHigh, pValue_LowVSHigh)
                } else
                {
                  currentGeneAllResults[c("QS_LowMean", "QS_InterMean", "QS_HighMean", "LowVSInter_QS", "InterVSHigh_QS", "LowVSHigh_QS"), 1] = c(scoresLow_mean, scoresInter_mean, scoresHigh_mean, pValue_LowVSInter, pValue_InterVSHigh, pValue_LowVSHigh)
                }
                
                utils::write.table(data_export, row.names = FALSE, col.names = FALSE, quote = FALSE, file = file.path(basisFolder, "output", "apply", "customSignatures", exportShortcut, distinguishingGene, "stats", paste(distinguishingGene, "_bars-", exportShortcut, "_stats_(Sign = ", currentSignatureName, ").txt", sep = "")), sep = "\t")
              }
            
            currentGeneAllResults = c(colnames(currentGeneAllResults)[1], currentGeneAllResults[, 1])
          }
          
          if(signatureUsed == "custom")
          {
            utils::write.table(resultingGenesData_distinguishGene, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path(basisFolder, "output", "apply", "customSignatures", "fullTables", distinguishingGene, paste(distinguishingGene, "_fullOutputData_(Sign = ", currentSignatureName, ").apply", sep = "")), sep = "\t")
          } else if(signatureUsed == "classification")
          {
            utils::write.table(resultingGenesData_distinguishGene, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path(basisFolder, "output", "apply", paste(distinguishingGene, "_fullOutputData_(Sign = ", currentSignatureName, ").apply", sep = "")), sep = "\t")
          }
          
          return(currentGeneAllResults)
        }
      
      close(pb)
      parallel::stopCluster(cl)
      
      if(signatureUsed == "custom")
      {
        currentGenePooledResults = data.frame(currentGenePooledResults, stringsAsFactors = FALSE)
        
        if(distinguishingGenesMode == FALSE)
        {
          currentGenePooledResults = rbind("", currentGenePooledResults)
          currentGenePooledResults[1,] = colnames(currentGenePooledResults)
          colnames(currentGenePooledResults) = NULL
          rownames(currentGenePooledResults)[1] = "Signatures"
        } else
        {
          rownames(currentGenePooledResults) = c("Signatures", currentGenePooledResults_rownames)
        }
        
        if(g == 1)
        {
          allSignaturesPooledResults = currentGenePooledResults[-c(2:13),]
          utils::write.table(allSignaturesPooledResults, row.names = TRUE, col.names = FALSE, quote = FALSE, file = file.path(basisFolder, "output", "apply", "customSignatures", "synthesis", "correlations_CS-VS-QS_allSignaturesSynthesis.txt"), sep = "\t")
        }
        
        if(distinguishingGenesMode == TRUE)
        {
          currentGenePooledResults = currentGenePooledResults[c(1:13),]
          
          utils::write.table(currentGenePooledResults, row.names = TRUE, col.names = FALSE, quote = FALSE, file = file.path(basisFolder, "output", "apply", "customSignatures", "synthesis", distinguishingGene, paste(distinguishingGene, "_allSignaturesSynthesis.txt", sep = "")), sep = "\t")
          
          if(verbose == TRUE)
          {		
            print(paste("## Zipping the ", distinguishingGene, " folders ##", sep = ""))
          }
          
          foldersToZip = c("CS", "QS", "CS VS QS", "fullTables")
          foldersToZipShortName = c("CS", "QS", "CS-VS-QS", "fullTables")
          
          foreach::foreach(z = 1:length(foldersToZip)) %do%
            {
              currentFolderToZip = foldersToZip[z]
              currentFolderToZipShortName = foldersToZipShortName[z]
              
              previouswd = getwd()
              setwd(file.path("output", "apply", "customSignatures", currentFolderToZip))
              
              filesToZip = list.files(path = file.path(distinguishingGene), pattern = "*", recursive = TRUE)
              zip::zipr(zipfile = paste(currentFolderToZipShortName, "_", distinguishingGene, ".zip", sep = ""), files = file.path(distinguishingGene, filesToZip))
              
              unlink(file.path(distinguishingGene, filesToZip), force = TRUE, recursive = TRUE)
              directoriesToDelete = gsub("\\./", "", list.dirs()[-1])
              keepNoDGFolder = grep("noDG", directoriesToDelete)
              if(length(keepNoDGFolder) > 0)
              {
                directoriesToDelete = directoriesToDelete[-keepNoDGFolder]
                
              }
              unlink(directoriesToDelete, force = TRUE, recursive = TRUE)
              
              setwd(previouswd)
            }
        }
      }
    }
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
}
