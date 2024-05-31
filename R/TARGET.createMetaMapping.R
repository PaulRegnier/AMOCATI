#' Create the metamapping for TARGET projects
#'
#' This function creates the metamapping related to TARGET projects from GDC.
#'
#' @param projectID A character value determining which project to use.
#' 
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @return Generated output file is exported to the `output > data` folder.
#'
#' @importFrom foreach %do%
#'
#' @export

TARGET.createMetaMapping = function(projectID = NULL, verbose = TRUE)
{
  if(verbose == TRUE)
  {
    print("## Opening clinical data file ##")
  }
  
  targetXlsxData = as.data.frame(readxl::read_excel(file.path("data", "Clinical", "TARGETClinicalData.xlsx")))	
  columnsToKeep = c("TARGET USI", "Vital Status", "Overall Survival Time in Days")
  targetXlsxData = targetXlsxData[columnsToKeep]
  colnames(targetXlsxData) = c("TARGET-ID", "vitalStatus", "survivedDays")
  
  duplicatedPatientsID = which(duplicated(targetXlsxData$"TARGET-ID"))


  if(length(duplicatedPatientsID) > 0)
  {
    rowsToDelete = NULL
    a = NULL
    
    foreach::foreach(a = 1:length(duplicatedPatientsID)) %do%
      {
        currentDuplicatedPatientID = duplicatedPatientsID[a]
        currentDuplicatedPatient = targetXlsxData[currentDuplicatedPatientID, "TARGET-ID"]
        currentDuplicatedPatientData = targetXlsxData[targetXlsxData$"TARGET-ID" == currentDuplicatedPatient, ]
        
        rowsToDelete = c(rowsToDelete, rownames(currentDuplicatedPatientData)[-1])
      }
    
    rowsToDelete = unique(as.numeric(rowsToDelete))
    
    targetXlsxData = targetXlsxData[-rowsToDelete, ]
    
  }
  
 
  newXlsxData = data.frame()
  
  if(verbose == TRUE)
  {
    print("## Data matching on case UUID (may be long) ##")
  }
  
  progressBar = tcltk::tkProgressBar(title = "Work in progress, please wait for completion", min = 0, max = length(rownames(targetXlsxData)))
  
  i = NULL
  
  foreach::foreach(i = 1:nrow(targetXlsxData)) %do%
    {
      tcltk::setTkProgressBar(progressBar, i, label = paste(round(i/length(rownames(targetXlsxData))*100, 0), "% done"))
      
      currentRowTargetID = as.vector(targetXlsxData[i,]$"TARGET-ID")
      
      metaMappingTargetRaw = GenomicDataCommons::files() %>% GenomicDataCommons::filter(~ cases.submitter_id == currentRowTargetID & data_type == "Gene Expression Quantification" & cases.submitter_id == currentRowTargetID) %>% GenomicDataCommons::select(c("file_id", "cases.case_id", "data_type", "file_name", "file_size")) %>% GenomicDataCommons::response_all()
      
      if(length(metaMappingTargetRaw$results) > 0)
      {
        metaMappingTargetRaw = metaMappingTargetRaw$results
        newXlsxData[i, "TARGET-ID"] = currentRowTargetID
        newXlsxData[i, "vitalStatus"] = as.vector(targetXlsxData[targetXlsxData$"TARGET-ID" == currentRowTargetID,]$"vitalStatus")
        newXlsxData[i, "survivedDays"] = as.vector(targetXlsxData[targetXlsxData$"TARGET-ID" == currentRowTargetID,]$"survivedDays")
        newXlsxData[i, "RNASeq_FileName"] = gsub("\\.gz", "", as.vector(metaMappingTargetRaw$file_name)[1])
        newXlsxData[i, "RNASeq_Case_UUID"] = as.vector(unlist(metaMappingTargetRaw$cases))[1]
        newXlsxData[i, "RNASeq_File_UUID"] = as.vector(metaMappingTargetRaw$file_id)[1]
        newXlsxData[i, "RNASeq_DataType"] = as.vector(metaMappingTargetRaw$data_type)[1]
        newXlsxData[i, "RNASeq_FileSize"] = as.vector(metaMappingTargetRaw$file_size)[1]
      }
    }
  
  close(progressBar)
  
  if(verbose == TRUE)
  {
    print("## Eliminating non-matched patients ##")
  }
  
  newXlsxData = newXlsxData[!is.na(newXlsxData$"RNASeq_FileName"),]
  
  if(verbose == TRUE)
  {
    print("## Exporting final metamapping ##")
  }
  
  utils::write.table(newXlsxData, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "finalMetaMapping.txt"), sep = "\t")
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
}