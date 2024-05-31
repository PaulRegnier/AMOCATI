#' Lists the projects that are downloadable through AMOCATI
#'
#' This function lists the GDC projects that are currently suitable for downloading and processing through AMOCATI.
#'
#' @return Outputs its results directly in the console.
#' 
#' @importFrom dplyr %>%
#'
#' @export

listProjectsAttributes = function()
{
  # projectsList = data.frame(matrix("", ncol = 5, nrow = 0), stringsAsFactors = FALSE)
  # colnames(projectsList) = c("Project", "ProjectID", "SamplesType", "ProjectName", "RNASeqWorkflows")
  # 
  # projectsList[1, ] = c("TCGA", "TCGA-BRCA", "Primary Tumor", "Breast Invasive Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[2, ] = c("TCGA", "TCGA-COAD", "Primary Tumor", "Colon Adenocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[3, ] = c("TCGA", "TCGA-ACC", "Primary Tumor", "Adrenocortical Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[4, ] = c("TCGA", "TCGA-SKCM", "Primary Tumor", "Skin Cutaneous Melanoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[5, ] = c("TCGA", "TCGA-SKCM", "Metastatic", "Skin Cutaneous Melanoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[6, ] = c("TCGA", "TCGA-UCEC", "Primary Tumor", "Uterine Corpus Endometrial Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[7, ] = c("TCGA", "TCGA-UCS", "Primary Tumor", "Uterine Carcinosarcoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[8, ] = c("TCGA", "TCGA-UVM", "Primary Tumor", "Uveal Melanoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[9, ] = c("TCGA", "TCGA-THYM", "Primary Tumor", "Thymoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[10, ] = c("TCGA", "TCGA-THCA", "Primary Tumor", "Thyroid Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[11, ] = c("TCGA", "TCGA-TGCT", "Primary Tumor", "Testicular Germ Cell Tumors", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[12, ] = c("TCGA", "TCGA-STAD", "Primary Tumor", "Stomach Adenocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[13, ] = c("TCGA", "TCGA-SARC", "Primary Tumor", "Sarcoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[14, ] = c("TCGA", "TCGA-READ", "Primary Tumor", "Rectum Adenocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[15, ] = c("TCGA", "TCGA-PRAD", "Primary Tumor", "Prostate Adenocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[16, ] = c("TCGA", "TCGA-PCPG", "Primary Tumor", "Pheochromocytoma and Paraganglioma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[17, ] = c("TCGA", "TCGA-PAAD", "Primary Tumor", "Pancreatic Adenocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[18, ] = c("TCGA", "TCGA-OV", "Primary Tumor", "Ovarian Serous Cystadenocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[19, ] = c("TCGA", "TCGA-PAAD", "Primary Tumor", "Pancreatic Adenocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[20, ] = c("TCGA", "TCGA-MESO", "Primary Tumor", "Mesothelioma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[21, ] = c("TCGA", "TCGA-LUAD", "Primary Tumor", "Lung Adenocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[22, ] = c("TCGA", "TCGA-LUSC", "Primary Tumor", "Lung Squamous Cell Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[23, ] = c("TCGA", "TCGA-LIHC", "Primary Tumor", "Liver Hepatocellular Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[24, ] = c("TCGA", "TCGA-LGG", "Primary Tumor", "Brain Lower Grade Glioma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[25, ] = c("TCGA", "TCGA-KIRP", "Primary Tumor", "Kidney Renal Papillary Cell Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[26, ] = c("TCGA", "TCGA-KIRC", "Primary Tumor", "Kidney Renal Clear Cell Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[27, ] = c("TCGA", "TCGA-KICH", "Primary Tumor", "Kidney Chromophobe", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[28, ] = c("TCGA", "TCGA-HNSC", "Primary Tumor", "Head and Neck Squamous Cell Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[29, ] = c("TCGA", "TCGA-GBM", "Primary Tumor", "Glioblastoma Multiforme", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[30, ] = c("TCGA", "TCGA-ESCA", "Primary Tumor", "Esophageal Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[31, ] = c("TCGA", "TCGA-DLBC", "Primary Tumor", "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[32, ] = c("TCGA", "TCGA-CHOL", "Primary Tumor", "Cholangiocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[33, ] = c("TCGA", "TCGA-CESC", "Primary Tumor", "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[34, ] = c("TCGA", "TCGA-BLCA", "Primary Tumor", "Bladder Urothelial Carcinoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[35, ] = c("TARGET", "TARGET-NBL", "Primary Tumor", "Neuroblastoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[36, ] = c("TARGET", "TARGET-WT", "Primary Tumor", "High-Risk Wilms Tumor", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[37, ] = c("TARGET", "TARGET-OS", "Primary Tumor", "Osteosarcoma", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[38, ] = c("TARGET", "TARGET-RT", "Primary Tumor", "Rhabdoid Tumor", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[39, ] = c("CGCI", "CGCI-HTMCP-CC", "Primary Tumor", "HIV+ Tumor Molecular Characterization Project - Cervical Cancer", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[40, ] = c("TCGA", "TCGA-LAML", "Primary Blood Derived Cancer - Peripheral Blood", "Acute Myeloid Leukemia", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[41, ] = c("TARGET", "TARGET-AML", "Primary Blood Derived Cancer - Peripheral Blood, Recurrent Blood Derived Cancer - Peripheral Blood", "Acute Myeloid Leukemia (Peripheral Blood)", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[42, ] = c("TARGET", "TARGET-AML", "Primary Blood Derived Cancer - Bone Marrow, Recurrent Blood Derived Cancer - Bone Marrow", "Acute Myeloid Leukemia (Bone Marrow)", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ")
  # projectsList[43, ] = c("TARGET", "TARGET-ALL-P2", "Primary Blood Derived Cancer - Peripheral Blood, Recurrent Blood Derived Cancer - Peripheral Blood", "Acute Lymphoblastic Leukemia - Phase II (Peripheral Blood)", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ") 
  # projectsList[44, ] = c("TARGET", "TARGET-ALL-P2", "Primary Blood Derived Cancer - Bone Marrow, Recurrent Blood Derived Cancer - Bone Marrow", "Acute Lymphoblastic Leukemia - Phase II (Bone Marrow)",  "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ") 
  # projectsList[45, ] = c("TARGET", "TARGET-ALL-P3", "Primary Blood Derived Cancer - Bone Marrow, Recurrent Blood Derived Cancer - Bone Marrow", "Acute Lymphoblastic Leukemia - Phase III (Bone Marrow)", "HTSeq - Counts,  HTSeq - FPKM, HTSeq - FPKM-UQ") 
  
  projectsListOnline = GenomicDataCommons::projects() %>% GenomicDataCommons::response_all()
  projectsListOnline = projectsListOnline$results
  projectsListOnline = projectsListOnline[projectsListOnline$released == TRUE & projectsListOnline$state == "open", ]
  
  projectsToKeep_ID = grep("TCGA|TARGET|CGCI", projectsListOnline$id)
  
  projectsListOnline = projectsListOnline[projectsToKeep_ID, ]
  
  projectsList = data.frame(ProjectID = projectsListOnline$project_id, ProjectName = projectsListOnline$name)
  
  projectsList = projectsList[order(projectsList$ProjectID), ]
  rownames(projectsList) = NULL
  print(projectsList)
  
  sink(file.path("silentOutput.log"))
  gc()
  sink()
}




