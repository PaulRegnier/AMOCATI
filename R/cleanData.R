#' Cleans the dataset once the related data (RNA-Seq + clinical) are downloaded, matched and processed
#'
#' This function cleans the dataset after its preprocessing. It is automatically called by the related pooling functions and notably remove columns that are not genes, removes duplicated/incorrect patients, formats gene names to the latest HGNC format, merges duplicated genes, normalizes and transforms the dataset, exports the final dataset and zips the original RNA-Seq and clinical data.
#'
#' @param fullData A single table with patients in rows and genes in columns.
#'
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @return Generated output files are exported to the `output > data` folder.
#'
#' @importFrom foreach %do%
#'
#' @export

cleanData = function(fullData = fullData, verbose = TRUE)
{
  if(verbose == TRUE)
  {		
    print("## Removing columns that are not genes ##")
  }
  
  totalGenesNames = colnames(fullData)
  genesToKeepID = grep("ENSG", totalGenesNames)
  fullData = fullData[, c(1:3, genesToKeepID)]
  colnames(fullData) = gsub("\\.[0-9]+", "", colnames(fullData))
  
  if(verbose == TRUE)
  {		
    print("## Removing duplicated patients ##")
  }
  
  duplicatedPatientsID = which(duplicated(fullData$CaseUUID))
  duplicatedPatientsCaseUUID = unique(fullData$CaseUUID[duplicatedPatientsID])
  associatedRowsToDelete = which(fullData$CaseUUID %in% duplicatedPatientsCaseUUID)
  fullData$survivedDays = as.numeric(fullData$survivedDays)
  fullData = data.frame(fullData, stringsAsFactors = FALSE)
  
  fullData[, -c(1:3)] = apply(fullData[, -c(1:3)], 2, as.vector)
  fullData[, -c(1:3)] = apply(fullData[, -c(1:3)], 2, as.numeric)

  if(length(associatedRowsToDelete) > 0)
  {
    d = NULL
    
    foreach::foreach(d = 1:length(duplicatedPatientsCaseUUID)) %do%
      {
        currentDuplicatedCaseUUID = duplicatedPatientsCaseUUID[d]
        currentDuplicatedData = fullData[fullData$CaseUUID == currentDuplicatedCaseUUID,]
        currentDuplicatedData = currentDuplicatedData[order(currentDuplicatedData$survivedDays, decreasing = TRUE),]
        
        if(length(unique(currentDuplicatedData$vitalStatus)) == 1 & length(unique(currentDuplicatedData$survivedDays)) == 1)
        {
          currentRemainingData = cbind(as.data.frame(currentDuplicatedData[1, c(1:3)], stringsAsFactors = FALSE), t(as.data.frame(round(colMeans(currentDuplicatedData[, -c(1:3)])), stringsAsFactors = FALSE)))
        } else
        {
          currentRemainingData = currentDuplicatedData[1, ]
        }
        
        if(d == 1)
        {
          totalDuplicatedRemainingData = currentRemainingData
        } else
        {
          totalDuplicatedRemainingData = rbind(totalDuplicatedRemainingData, currentRemainingData)
        }
      }
    
    fullData = fullData[-associatedRowsToDelete, ]
    fullData = rbind(fullData, totalDuplicatedRemainingData)
    fullData = fullData[order(fullData$CaseUUID),]
    rownames(fullData) = NULL
  }
  
  if(verbose == TRUE)
  {		
    print("## Removing incorrect patients ##")
  }
  
  survivedDaysNAPatients = which(is.na(fullData$survivedDays))
  
  if(length(survivedDaysNAPatients) > 0)
  {
    fullData = fullData[-survivedDaysNAPatients, ]
  }
  
  survivedDaysNullPatients = which(is.null(fullData$survivedDays))
  
  if(length(survivedDaysNullPatients) > 0)
  {
    fullData = fullData[-survivedDaysNullPatients, ]
  }
  
  survivedDaysZeroPatients = which(fullData$survivedDays <= 0)
  
  if(length(survivedDaysZeroPatients) > 0)
  {
    fullData = fullData[-survivedDaysZeroPatients, ]
  }
  
  vitalStatusCorrectPatients = which(fullData$vitalStatus %in% c("Dead", "Alive"))
  
  if(length(vitalStatusCorrectPatients) > 0)
  {
    fullData = fullData[vitalStatusCorrectPatients, ]
  }
  
  if(verbose == TRUE)
  {
    print("## Getting annotations from BiomaRt ##")
  }
  
  martAnnotations = biomaRt::useDataset("hsapiens_gene_ensembl", biomaRt::useMart("ensembl"))
  genesNamesCorrespondence = biomaRt::getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"), values = colnames(fullData)[-c(1:3)], mart = martAnnotations)
  
  if(verbose == TRUE)
  {
    print("## Removing ENSG genes not matching with an official Gene Symbol ##")
  }
  
  duplicatedENSGGenes = which(duplicated(genesNamesCorrespondence[, "ensembl_gene_id"]))
  
  if(length(duplicatedENSGGenes) > 0)
  {
    rowsToDelete = NULL
    f = NULL
    
    foreach::foreach(f = 1:length(duplicatedENSGGenes)) %do%
      {
        currentDuplicatedGeneID = duplicatedENSGGenes[f]
        currentDuplicatedGeneENSG = genesNamesCorrespondence[currentDuplicatedGeneID, "ensembl_gene_id"]
        currentDuplicatedGeneRows = which(genesNamesCorrespondence$ensembl_gene_id %in% currentDuplicatedGeneENSG)
        rowsToDelete = c(rowsToDelete, currentDuplicatedGeneRows[-1])
      }
    
    genesNamesCorrespondence = genesNamesCorrespondence[-rowsToDelete, ]
  }
  
  columnsToKeep = as.vector(c(colnames(fullData)[1:3], genesNamesCorrespondence[genesNamesCorrespondence$hgnc_symbol != "", "ensembl_gene_id"]))
  
  fullData = fullData[, columnsToKeep]
  
  if(verbose == TRUE)
  {
    print("## Converting remaining ENSG gene names to official Gene Symbol ##")
  }
  
  e = NULL
  
  foreach::foreach(e = 4:ncol(fullData)) %do%
    {
      currentGeneToRename = colnames(fullData)[e]
      currentGene_matchingLine = which(genesNamesCorrespondence$ensembl_gene_id == currentGeneToRename)
      currentGene_symbol = genesNamesCorrespondence[currentGene_matchingLine, "hgnc_symbol"]
      colnames(fullData)[e] = currentGene_symbol
    }
  
  if(verbose == TRUE)
  {
    print("## Merging counts of duplicated official Gene Symbol ##")
  }
  
  fullData[, -c(1:3)] = apply(fullData[, -c(1:3)], 2, as.numeric)
  
  duplicatedGeneSymbolGenes = which(duplicated(colnames(fullData)))
  
  if(length(duplicatedGeneSymbolGenes) > 0)
  {
    columnsToDelete = NULL
    columnsToAdd = NULL
    g = NULL
    
    foreach::foreach(g = 1:length(duplicatedGeneSymbolGenes)) %do%
      {
        currentDuplicatedGeneID = duplicatedGeneSymbolGenes[g]
        currentDuplicatedGeneSymbol = colnames(fullData)[currentDuplicatedGeneID]
        currentDuplicatedGeneColumnsID = which(colnames(fullData) %in% currentDuplicatedGeneSymbol)
        currentDuplicatedGeneColumns = fullData[, currentDuplicatedGeneColumnsID]
        columnsToDelete = c(columnsToDelete, currentDuplicatedGeneColumnsID)
        columnsToAdd = cbind(columnsToAdd, rowSums(currentDuplicatedGeneColumns))
        colnames(columnsToAdd)[g] = currentDuplicatedGeneSymbol
      }
    
    fullData = fullData[, -columnsToDelete]
    fullData = cbind(fullData, columnsToAdd)
  }
  
  if(verbose == TRUE)
  {		
    print("## Normalizing and transforming data ##")
  }
  
  currentRNASeqData_data = sapply(data.table::transpose(fullData[,-c(1:3)]), as.numeric)
  condition = factor(rep(c("A"), length(colnames(currentRNASeqData_data))))
  currentDESeqDataSet = DESeq2::DESeqDataSetFromMatrix(currentRNASeqData_data, S4Vectors::DataFrame(condition), design = ~ 1)
  
  currentDESeqDataSet = DESeq2::varianceStabilizingTransformation(currentDESeqDataSet, blind = TRUE)	
  currentDESeqDataSet = SummarizedExperiment::assay(currentDESeqDataSet)
  
  currentDESeqDataSet = data.frame(t(currentDESeqDataSet), stringsAsFactors = FALSE)
  colnames(currentDESeqDataSet) = colnames(fullData[,-c(1:3)])
  
  currentDESeqDataSet_untouched = cbind(fullData[,c(1:3)], currentDESeqDataSet)
  
  if(verbose == TRUE)
  {		
    print("## Exporting results (may be long) ##")
  }			
  
  utils::write.table(currentDESeqDataSet_untouched, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "fullData_untouched.data"), sep = "\t")
  
  nonGenomicData = fullData[,c(1:3)]
  
  rowsToChangeDead = as.numeric(rownames(nonGenomicData[as.numeric(nonGenomicData$survivedDays) > 1825 & nonGenomicData$vitalStatus == "Dead", ]))
  rowsToChangeAlive = as.numeric(rownames(nonGenomicData[as.numeric(nonGenomicData$survivedDays) > 1825 & nonGenomicData$vitalStatus == "Alive", ]))
  
  rowsToChange = c(rowsToChangeDead, rowsToChangeAlive)
  rowsToChange = rowsToChange[!is.na(rowsToChange)]
  nonGenomicData[rownames(nonGenomicData) %in% rowsToChange, "survivedDays"] = 1825
  nonGenomicData[rownames(nonGenomicData) %in% rowsToChange, "vitalStatus"] = "Alive"
  
  currentDESeqDataSet_fittedTo5YearsSurvival = cbind(nonGenomicData, currentDESeqDataSet)
  
  utils::write.table(currentDESeqDataSet_fittedTo5YearsSurvival, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "fullData.data"), sep = "\t")
  
  currentDESeqDataSet_fittedTo5YearsSurvival_simplified = currentDESeqDataSet_fittedTo5YearsSurvival
  
  currentDESeqDataSet_fittedTo5YearsSurvival_clinicalData = currentDESeqDataSet_fittedTo5YearsSurvival_simplified[, c("CaseUUID", "vitalStatus", "survivedDays")]
  columnsToDelete = which(colnames(currentDESeqDataSet_fittedTo5YearsSurvival_simplified) %in% c("vitalStatus", "survivedDays"))
  currentDESeqDataSet_fittedTo5YearsSurvival_simplified = currentDESeqDataSet_fittedTo5YearsSurvival_simplified[, -columnsToDelete]
  rownames(currentDESeqDataSet_fittedTo5YearsSurvival_simplified) = currentDESeqDataSet_fittedTo5YearsSurvival_simplified$CaseUUID
  currentDESeqDataSet_fittedTo5YearsSurvival_simplified = currentDESeqDataSet_fittedTo5YearsSurvival_simplified[, colnames(currentDESeqDataSet_fittedTo5YearsSurvival_simplified) != "CaseUUID"]
  currentDESeqDataSet_fittedTo5YearsSurvival_simplified = t(currentDESeqDataSet_fittedTo5YearsSurvival_simplified)
  currentDESeqDataSet_fittedTo5YearsSurvival_simplified = rbind(colnames(currentDESeqDataSet_fittedTo5YearsSurvival_simplified), currentDESeqDataSet_fittedTo5YearsSurvival_simplified)
  currentDESeqDataSet_fittedTo5YearsSurvival_simplified = cbind(rownames(currentDESeqDataSet_fittedTo5YearsSurvival_simplified), currentDESeqDataSet_fittedTo5YearsSurvival_simplified)
  currentDESeqDataSet_fittedTo5YearsSurvival_simplified[1,1] = "Symbol"
  
  utils::write.table(currentDESeqDataSet_fittedTo5YearsSurvival_simplified, row.names = FALSE, col.names = FALSE, quote = FALSE, file = file.path("output", "data", "fullData_simplifiedReversed.data"), sep = "\t")
  
  if(verbose == TRUE)
  {		
    print("## Zipping the data (Clinical + RNASeq) folder ##")
  }
  
  previouswd = getwd()
  setwd(file.path("data"))
  dirsToZip = list.dirs(recursive = TRUE)
  dirsToZip = dirsToZip[dirsToZip != "."]
  zip::zipr(zipfile = "data.zip", files = dirsToZip)
  
  unlink(file.path(dirsToZip, "*.*"), force = TRUE)
  setwd(previouswd)
  
  # Empty temp folder aftewards
  
  unlink(file.path("temp"), force = TRUE, recursive = TRUE)
  dir.create(file.path("temp"))
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
}