#' Create the metamapping for TCGA or CGCI projects
#'
#' This function creates the metamapping related to TCGA or CGCI projects from GDC.
#' 
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @return Generated output file is exported to the `output > data` folder.
#'
#' @importFrom foreach %do%
#'
#' @export

TCGA_CGCI.createMetaMapping = function(verbose = TRUE)
{
  metaMappingRNASeqRawResults = utils::read.table(file.path("output", "data", "metaMappingRNASeqRaw.txt"), sep = "\t", header = TRUE)
  metaMappingClinicalRawResults = utils::read.table(file.path("output", "data", "metaMappingClinicalRaw.txt"), sep = "\t", header = TRUE)
  
  if(verbose == TRUE)
  {
    print("## Renaming columns ##")
  }
  
  colnames(metaMappingRNASeqRawResults) = gsub("file_name", "RNASeq_FileName", colnames(metaMappingRNASeqRawResults))
  colnames(metaMappingRNASeqRawResults) = gsub("cases", "RNASeq_Case_UUID", colnames(metaMappingRNASeqRawResults))
  colnames(metaMappingRNASeqRawResults) = gsub("file_id", "RNASeq_File_UUID", colnames(metaMappingRNASeqRawResults))
  colnames(metaMappingRNASeqRawResults) = gsub("data_type", "RNASeq_DataType", colnames(metaMappingRNASeqRawResults))
  colnames(metaMappingRNASeqRawResults) = gsub("file_size", "RNASeq_FileSize", colnames(metaMappingRNASeqRawResults))
  
  colnames(metaMappingClinicalRawResults) = gsub("file_name", "Clinical_FileName", colnames(metaMappingClinicalRawResults))
  colnames(metaMappingClinicalRawResults) = gsub("cases", "Clinical_Case_UUID", colnames(metaMappingClinicalRawResults))
  colnames(metaMappingClinicalRawResults) = gsub("file_id", "Clinical_File_UUID", colnames(metaMappingClinicalRawResults))
  colnames(metaMappingClinicalRawResults) = gsub("data_type", "Clinical_DataType", colnames(metaMappingClinicalRawResults))
  colnames(metaMappingClinicalRawResults) = gsub("file_size", "Clinical_FileSize", colnames(metaMappingClinicalRawResults))
  
  metaMappingClinicalRawResults = metaMappingClinicalRawResults[grep("(nationwidechildrens.org_clinical)|(genome.wustl.edu_clinical)", metaMappingClinicalRawResults$Clinical_FileName), ]
  
  if(verbose == TRUE)
  {
    print("## Data matching on case UUID ##")
  }	
  
  a = NULL
  
  foreach::foreach(a = 1:nrow(metaMappingRNASeqRawResults)) %do%
    {
      currentRowCase = as.vector(metaMappingRNASeqRawResults[a,]$"RNASeq_Case_UUID")
      
      if(currentRowCase %in% metaMappingClinicalRawResults$"Clinical_Case_UUID")
      {
        metaMappingRNASeqRawResults[metaMappingRNASeqRawResults$"RNASeq_Case_UUID" == currentRowCase, "Clinical_Case_UUID"] = metaMappingClinicalRawResults[metaMappingClinicalRawResults$"Clinical_Case_UUID" == currentRowCase,]$"Clinical_Case_UUID"[1]
        
        metaMappingRNASeqRawResults[metaMappingRNASeqRawResults$"RNASeq_Case_UUID" == currentRowCase, "Clinical_FileName"] = gsub("\\.gz", "", metaMappingClinicalRawResults[metaMappingClinicalRawResults$"Clinical_Case_UUID" == currentRowCase,]$"Clinical_FileName"[1])
        
        metaMappingRNASeqRawResults[metaMappingRNASeqRawResults$"RNASeq_Case_UUID" == currentRowCase, "Clinical_File_UUID"] = metaMappingClinicalRawResults[metaMappingClinicalRawResults$"Clinical_Case_UUID" == currentRowCase,]$"Clinical_File_UUID"[1]
        
        metaMappingRNASeqRawResults[metaMappingRNASeqRawResults$"RNASeq_Case_UUID" == currentRowCase, "Clinical_DataType"] = metaMappingClinicalRawResults[metaMappingClinicalRawResults$"Clinical_Case_UUID" == currentRowCase,]$"Clinical_DataType"[1]
        
        metaMappingRNASeqRawResults[metaMappingRNASeqRawResults$"RNASeq_Case_UUID" == currentRowCase, "Clinical_FileSize"] = metaMappingClinicalRawResults[metaMappingClinicalRawResults$"Clinical_Case_UUID" == currentRowCase,]$"Clinical_FileSize"[1]
      } else
      {
        metaMappingRNASeqRawResults = metaMappingRNASeqRawResults[-a,]
      }
    }
  
  metaMappingRNASeqRawResults$"RNASeq_FileName" = gsub("\\.gz", "", metaMappingRNASeqRawResults$"RNASeq_FileName")
  
  if(verbose == TRUE)
  {
    print("## Eliminating non-matched patients ##")
  }
  
  metaMappingRNASeqRawResults = metaMappingRNASeqRawResults[is.na(metaMappingRNASeqRawResults$"Clinical_Case_UUID") == FALSE,]
  
  if(verbose == TRUE)
  {
    print("## Exporting final metamapping ##")
  }
  
  utils::write.table(metaMappingRNASeqRawResults, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "finalMetaMapping.txt"), sep = "\t")
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
}