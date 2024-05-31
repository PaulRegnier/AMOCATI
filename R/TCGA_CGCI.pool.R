#' Pools the RNA-Seq and clinical data from the desired TCGA or CGCI project
#'
#' This function opens and pools the RNA-Seq and clinical data from the desired GDC TCGA or CGCI project.
#'
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @importFrom foreach %dopar%
#'
#' @export

TCGA_CGCI.pool = function(verbose = TRUE)
{
  filteredData = data.frame(NULL)
  
  if(verbose == TRUE)
  {
    print("## Opening of metamapping file ##")
  }
  
  metaMapping = utils::read.table(file.path("output", "data", "finalMetaMapping.txt"), sep = "\t", header = TRUE)
  totalCases = length(rownames(metaMapping))
  
  if(verbose == TRUE)
  {
    print("## Beginning of data filtering and pooling ##")
  }
  
  basisFolder = getwd()
  
  coresNumber = parallel::detectCores()
  cl = parallel::makeCluster(coresNumber)
  doSNOW::registerDoSNOW(cl)
  
  ntasks = totalCases
  pb = tcltk::tkProgressBar(max = ntasks)
  progress = function(n) tcltk::setTkProgressBar(pb, n, label = paste(round(n/ntasks*100, 0), "% done"))
  opts = list(progress = progress)
  
  b = NULL
  
  fullData = foreach::foreach(b = seq_len(totalCases), .packages = c("XML", "data.table"), .combine = rbind, .options.snow = opts) %dopar%
    {
      currentCase = as.vector(metaMapping[b, "Clinical_Case_UUID"])
      
      currentFilenameRNAseq = as.vector(metaMapping[b, "RNASeq_FileName"])
      
      currentRNASeqDataPath = file.path(basisFolder, "data", "RNASeq", currentFilenameRNAseq)
      currentRNASeqData = data.table::fread(input = currentRNASeqDataPath, sep = "\t", header = FALSE, data.table = FALSE)
      
      currentRNASeqData = currentRNASeqData[, currentRNASeqData[1, ] %in% c("gene_id", "unstranded")]
      
      colnames(currentRNASeqData) = c("ENSG_GeneName", "GeneExpressionValue")
      currentRNASeqData$"ENSG_GeneName" = gsub("\\..+", "", currentRNASeqData$"ENSG_GeneName")
      
      currentRNASeqData = as.data.frame(t(currentRNASeqData))
      colnames(currentRNASeqData) = as.vector(t(currentRNASeqData["ENSG_GeneName",]))
      currentRNASeqData = currentRNASeqData[rownames(currentRNASeqData) != "ENSG_GeneName",]
      
      currentFilenameClinical = as.vector(metaMapping[b,"Clinical_FileName"])
      
      if(length(grep("(.+)\\.txt", currentFilenameClinical)) == 0)
      {
        currentClinicalData = XML::xmlParseDoc(file.path(basisFolder, "data", "Clinical", currentFilenameClinical))
        
        if(length(XML::getNodeSet(currentClinicalData, "//clin_shared:vital_status")) > 0)
        {
          currentClinical_survivedDaysIfDead_matchingNodes = as.numeric(XML::xmlValue(XML::getNodeSet(currentClinicalData, "//clin_shared:days_to_death")))
          currentClinical_survivedDaysIfAlive_matchingNodes = as.numeric(XML::xmlValue(XML::getNodeSet(currentClinicalData, "//clin_shared:days_to_last_followup")))
                                                                                                                                                    
          if(all(is.na(currentClinical_survivedDaysIfDead_matchingNodes)) == TRUE)
          {
            currentClinical_mostRecentOccurenceID = which.max(currentClinical_survivedDaysIfAlive_matchingNodes)
          } else if(all(is.na(currentClinical_survivedDaysIfAlive_matchingNodes)) == TRUE)
          {
            currentClinical_mostRecentOccurenceID = which.max(currentClinical_survivedDaysIfDead_matchingNodes)
          } else
          {
            currentClinical_mostRecentOccurenceID = which.max(currentClinical_survivedDaysIfDead_matchingNodes)
          }
          
          if(length(currentClinical_mostRecentOccurenceID) > 0)
          {
            currentClinical_vitalStatus = XML::xmlValue(XML::getNodeSet(currentClinicalData, "//clin_shared:vital_status")[[currentClinical_mostRecentOccurenceID]])
            
            if(currentClinical_vitalStatus != "Alive")
            {
              currentClinical_survivedDays = as.numeric(XML::xmlValue(XML::getNodeSet(currentClinicalData, "//clin_shared:days_to_death")[[currentClinical_mostRecentOccurenceID]]))
            } else
            {
              currentClinical_survivedDays = as.numeric(XML::xmlValue(XML::getNodeSet(currentClinicalData, "//clin_shared:days_to_last_followup")[[currentClinical_mostRecentOccurenceID]]))
            }
            
            filteredData[b, "CaseUUID"] = currentCase
            filteredData[b, "vitalStatus"] = currentClinical_vitalStatus
            filteredData[b, "survivedDays"] = currentClinical_survivedDays
            
            temporaryBinding = cbind(filteredData[b,], as.data.frame(currentRNASeqData))
            return(temporaryBinding)
          }
        }
      }		
    }
  
  close(pb)
  parallel::stopCluster(cl)
  
  cleanData(verbose = verbose, fullData = fullData)
}