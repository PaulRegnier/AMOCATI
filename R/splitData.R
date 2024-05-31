#' Split the dataset into two new subdataset
#'
#' This function splits the dataset into two new subdatasets (`A` and `B`), using any proportion of interest. This is useful in datasets with a lot of samples to construct intra-cohort training/validation datasets.
#'
#' @param datasetAProportion A numeric value indicating the proportion of patients to put into the new subdataset `A`. The remaining patients will be put in the new subdataset `B`.
#' 
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @return Generated output files are exported to the `output > data` folder.
#'
#' @export

splitData = function(datasetAProportion = 0.5, verbose = TRUE)
{
  if(verbose == TRUE)
  {
    print("## fullData file opening ##")
  }
  
  fileToLoad = list.files(file.path("output", "data"), pattern = "fullData.data")
  fullGenesDataPath = file.path("output", "data", fileToLoad)
  fullGenesData = data.table::fread(input = fullGenesDataPath, sep = "\t", header = TRUE, data.table = FALSE)
  
  if(verbose == TRUE)
  {
    print("## Randomly sampling dataset's rows several times ##")
  }
  
  fullGenesData = fullGenesData[sample(nrow(fullGenesData)),]
  fullGenesData = fullGenesData[sample(nrow(fullGenesData)),]
  fullGenesData = fullGenesData[sample(nrow(fullGenesData)),]
  fullGenesData = fullGenesData[sample(nrow(fullGenesData)),]
  rownames(fullGenesData) = NULL
  
  if(verbose == TRUE)
  {
    print("## Splitting the dataset into 2 parts ##")
  }
  
  patientsNumberA = round(datasetAProportion*nrow(fullGenesData))
  datasetA = fullGenesData[1:patientsNumberA,]
  datasetB_from = patientsNumberA+1
  datasetB_to = nrow(fullGenesData)
  datasetB = fullGenesData[datasetB_from:datasetB_to,]
  
  if(verbose == TRUE)
  {
    print("## Exporting splitted data files ##")
  }
  
  utils::write.table(datasetA, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "fullData_datasetA.data"), sep = "\t")
  utils::write.table(datasetB, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "data", "fullData_datasetB.data"), sep = "\t")
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
}