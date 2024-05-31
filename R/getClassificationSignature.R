#' Extract the Classification Signature from the metaResults table
#'
#' This function helps to visually extract the most relevant genes that will be part of the Classification Signature.
#'
#' @param GeneScoreThreshold A numeric value indicating the threshold above the one genes are kept relative to their GeneScore.
#'
#' @param Gene_SNR_ExpressionThreshold A numeric value indicating the threshold above the one genes are kept relative to their Gene_SNR_Expression.
#'
#' @param exportSignature A boolean indicating if the function should actually export the Classification Signature (if set to `TRUE`) or not (if set to `FALSE`).
#'
#' @param verbose A boolean determining if the function should regularly output its progression in the console (if set to `TRUE`).
#'
#' @return Generated output files are exported to the `output > signatures` folder.
#'
#' @export

getClassificationSignature = function(GeneScoreThreshold = NULL, Gene_SNR_ExpressionThreshold = NULL, exportSignature = TRUE, verbose = TRUE)
{
  GeneScore = NULL
  Gene_SNR_Expression = NULL
  
  if(verbose == TRUE)
  {
    print("## Opening metaResults file ##")
  }
  
  fileToLoad = list.files(file.path("output", "metaResults"), pattern = "*.meta")[1]
  metaResults_allGenesPooledPath = file.path("output", "metaResults", fileToLoad)
  metaResults_allGenesPooled = data.table::fread(input = metaResults_allGenesPooledPath, sep = "\t", header = TRUE, data.table = FALSE)
  
  metaResults_allGenesPooled$GeneScore = log(as.numeric(metaResults_allGenesPooled$GeneScore))
  metaResults_allGenesPooled$Gene_SNR_Expression = log(as.numeric(metaResults_allGenesPooled$Gene_SNR_Expression))
  
  rowsToDelete = which(abs(metaResults_allGenesPooled$GeneScore) == "Inf")
  rowsToDelete = c(rowsToDelete, which(abs(metaResults_allGenesPooled$Gene_SNR_Expression) == "Inf"))
  
  if(length(rowsToDelete) > 0)
  {
    metaResults_allGenesPooled = metaResults_allGenesPooled[-rowsToDelete, ]
  }
  
  if(is.null(GeneScoreThreshold) == TRUE | is.null(Gene_SNR_ExpressionThreshold) == TRUE)
  {
    if(verbose == TRUE)
    {
      print("## One or more thresholds missing : only showing plot ##")
    }
    
    rawPlot = ggplot2::ggplot(metaResults_allGenesPooled, ggplot2::aes(x = GeneScore, y = Gene_SNR_Expression) ) +
      ggplot2::geom_bin2d(bins = 100) +
      ggplot2::scale_fill_continuous(type = "viridis") +
      ggplot2::xlab("log(GeneScore)") +
      ggplot2::ylab("log(Gene SNR Expression)") +
      ggplot2::theme_bw()
    
    plot(rawPlot)
  } else
  {
    
    bestSignature = as.data.frame(metaResults_allGenesPooled[metaResults_allGenesPooled$GeneScore >= GeneScoreThreshold & metaResults_allGenesPooled$Gene_SNR_Expression >= Gene_SNR_ExpressionThreshold, "HGNC_GeneSymbol"])		
    
    colnames(bestSignature) = "classificationSignature"
    
    genesNumber = nrow(bestSignature)
    
    classificationSignaturePlot = ggplot2::ggplot(metaResults_allGenesPooled, ggplot2::aes(x = GeneScore, y = Gene_SNR_Expression) ) +
      ggplot2::geom_bin2d(bins = 100) +
      ggplot2::scale_fill_continuous(type = "viridis") +
      ggplot2::xlab("log(GeneScore)") +
      ggplot2::ylab("log(Gene SNR Expression)") +
      ggplot2::geom_hline(yintercept = Gene_SNR_ExpressionThreshold, linetype = "dashed", color = "red", size = 1) +
      ggplot2::geom_vline(xintercept = GeneScoreThreshold, linetype = "dashed", color = "red", size = 1) +
      ggplot2::ggtitle(paste("Gene Score threshold : log = ", GeneScoreThreshold, ", original = ", round(exp(1)^GeneScoreThreshold, digits = 3), "\n", "Gene SNR Expression threshold : log =  ", Gene_SNR_ExpressionThreshold, ", original = ", round(exp(1)^Gene_SNR_ExpressionThreshold, digits = 3), "\n", "Number of genes that passed both thresholds = ", genesNumber, sep = "")) +
      ggplot2::theme_bw()
    
    classificationSignaturePlot = classificationSignaturePlot + ggplot2::theme(plot.title = ggplot2::element_text(color = "black", size = 12, hjust = 0.5))
    
    plot(classificationSignaturePlot)
    
    if(verbose == TRUE)
    {
      print(paste("## ", genesNumber, " genes passed both thresholds ##", sep = ""))
    }
    
    if(exportSignature == TRUE)
    {
      if(verbose == TRUE)
      {
        print("## Exporting classification signature ##")
      }
      
      utils::write.table(bestSignature, row.names = FALSE, col.names = TRUE, quote = FALSE, file = file.path("output", "signatures", "classificationSignature.sign"), sep = "\t")
      
      if(verbose == TRUE)
      {
        print("## Saving plot with thresholds ##")
      }
      
      ggplot2::ggsave(filename = "classificationSignature_plot.pdf", plot = classificationSignaturePlot, device = "pdf", path = file.path("output", "signatures"))
      
    }
  }
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
}