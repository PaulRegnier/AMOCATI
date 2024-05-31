#' Pools the RNA-Seq and clinical data from the desired TARGET project
#'
#' This function opens and pools the RNA-Seq and clinical data from the desired GDC TARGET project.
#'
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @importFrom foreach %dopar%
#'
#' @export

TARGET.pool = function(verbose = TRUE)
{
  filteredData = data.frame(NULL)
  tempDataFrame = NULL
  
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
  
  j = NULL
  
  fullData = foreach::foreach(j = seq_len(totalCases), .packages = c("data.table"), .combine = rbind, .options.snow = opts) %dopar%
    {
      currentRowCase = metaMapping[j,]
      currentFilenameRNAseq = as.vector(currentRowCase$"RNASeq_FileName")
      
      if(file.exists(file.path(basisFolder, "data", "RNASeq", currentFilenameRNAseq)) == TRUE)
      {
        currentRNASeqDataPath = file.path(basisFolder, "data", "RNASeq", currentFilenameRNAseq)
        currentRNASeqData = data.table::fread(input = currentRNASeqDataPath, sep = "\t", header = FALSE, data.table = FALSE)
        
        currentRNASeqData = currentRNASeqData[, currentRNASeqData[1, ] %in% c("gene_id", "unstranded")]
        
        colnames(currentRNASeqData) = c("ENSG_GeneName", "GeneExpressionValue")
        currentRNASeqData$"ENSG_GeneName" = gsub("\\..+", "", currentRNASeqData$"ENSG_GeneName")
        
        currentRNASeqData = as.data.frame(t(currentRNASeqData))
        colnames(currentRNASeqData) = as.vector(t(currentRNASeqData["ENSG_GeneName",]))
        currentRNASeqData = currentRNASeqData[rownames(currentRNASeqData) != "ENSG_GeneName",]
        
        filteredData[j, "CaseUUID"] = as.vector(currentRowCase$"RNASeq_Case_UUID")
        
        if(is.na(currentRowCase$vitalStatus) == TRUE)
        {
          filteredData[j, "vitalStatus"] = "N/A"
          
        } else
        {
          
          if(length(grep("^[Aa][Ll][Ii][Vv][Ee]$", as.vector(currentRowCase$"vitalStatus"))) > 0)
          {
            filteredData[j, "vitalStatus"] = "Alive"
          } else if(length(grep("^[Dd][Ee][Aa][Dd]$", as.vector(currentRowCase$"vitalStatus"))) > 0)
          {
            filteredData[j, "vitalStatus"] = "Dead"
          } else
          {
            filteredData[j, "vitalStatus"] = "N/A"
            
          }
          
          
        }
        
        filteredData[j, "survivedDays"] = as.vector(currentRowCase$"survivedDays")
        
        temporaryBinding = cbind(filteredData[j,], as.data.frame(currentRNASeqData))
        return(rbind(tempDataFrame,temporaryBinding))
      }
    }
  
  close(pb)	
  parallel::stopCluster(cl)
  
  cleanData(verbose = verbose, fullData = fullData)
}