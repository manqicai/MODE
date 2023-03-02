source("D:/Manqi/Meeting/20221027/utilities.R")
source("D:/Manqi/Meeting/20221023/utilities.R")
source("D:/Manqi/Package/MODE/R/EnsDeconv_wrap.R")
load("D:/Manqi/Meeting/20221023/minfi_DLPFC_ref.RData")
source("D:/Manqi/Meeting/20221027/utilities.R")
source("D:/Manqi/Meeting/20221023/utilities.R")
source("D:/Manqi/Package/MODE/R/Multi_mrk.R")
source("D:/Manqi/Package/MODE/R/Multi_mrk2.R")
source("D:/Manqi/EnsDeconv/R/CTS_EnsDeconv_wrapper.R")
source("D:/Manqi/EnsDeconv/R/CTS_EnsDeconv_LS.R")
source("D:/Manqi/deconv.ensemble/R/process_all.R")
source("D:/Manqi/deconv.ensemble/R/analysis_wrap.R")
source("D:/Manqi/deconv.ensemble/R/ensemble_prep.R")
source("D:/Manqi/deconv.ensemble/R/analyze_dset.R")
source("D:/Manqi/deconv.ensemble/R/get_metric_wrap.R")
source("D:/Manqi/Package/MODE/R/Multi_mrk_wrapper.R")
source("D:/Manqi/Package/MODE/R/Multi_mrk3.R")
load("D:/Manqi/data/DNAm/minfi_Gasparoni_ref.RData")
library(readxl)
GSE108066_knnWkernel_1000_2 <- readRDS("D:/Manqi/Meeting/20221103_07/GSE108066_knnWkernel_1000_2.rds")

load("D:/Manqi/Meeting/20221103_07/knn_GSE108066.RData")

meta_GSE108066 <- read_excel("D:/Manqi/Meeting/20221103_07/13059_2019_1747_MOESM1_ESM.xlsx",
                             sheet = "Table S2 - Covariates")
load("D:/Manqi/data/ROS/ROSDNA_bulk_meta_trueProp.RData")
ROS_true_4ct <- ROS_true/rowSums(ROS_true)
ROS_true_4ct[,2] <- ROS_true_4ct[,2]+ROS_true_4ct[,3]
ROS_true_4ct <- ROS_true_4ct[,-3]
colnames(ROS_true_4ct)[c(2,4)] = c("Micro_Endo","Oligo_MBP")
#' Title
#'
#' @param sig
#'
#' @return
#' @export
#'
#' @examples
sig_benchmarking_2 <- function(sig,include_minfi = F){

  B <- ros_dna_mtx_sub
  fractrue  <-  ROS_true_4ct

  cat("ROS")
  meta_sig <- data.frame(CellType = colnames(sig),Sample_Name =colnames(sig) )
  rownames(meta_sig) = colnames(sig)

  res_sig_ros <- EnsDeconv_wrap(B = B,frac_true = fractrue,
                                 sig = sig[intersect(rownames(sig),
                                                                          rownames(B)),],
                                include_minfi = include_minfi)
  rm(B)
  gc()
  DLPFC_minfi_truefrac <- matrix(0,58,2)
  colnames(DLPFC_minfi_truefrac) = unique(DLPFC_minfi_ref_meta$CellType)
  for (i in 1:nrow(DLPFC_minfi_truefrac)) {
    if(DLPFC_minfi_ref_meta$CellType[i] == "NeuN_pos"){
      DLPFC_minfi_truefrac[i,1] = 1
    }else{
      DLPFC_minfi_truefrac[i,2] = 1
    }
  }

  rownames(DLPFC_minfi_truefrac) = colnames(DLPFC_minfi_beta)

  Gasparoni_minfi_ref_meta_sub <- Gasparoni_minfi_ref_meta %>% filter(cell != "Bulk") %>%
    mutate(CellType = ifelse(cell == "Neuron","NeuN_pos","NeuN_neg"))

  Gasparoni_truefrac <- matrix(0,62,2)
  for (i in 1:nrow(Gasparoni_truefrac)) {
    if(Gasparoni_minfi_ref_meta_sub$CellType[i] == "NeuN_pos"){
      Gasparoni_truefrac[i,1] = 1
    }else{
      Gasparoni_truefrac[i,2] = 1
    }
  }
  colnames(Gasparoni_truefrac) = c("NeuN_pos","NeuN_neg" )
  rownames(Gasparoni_truefrac) = Gasparoni_minfi_ref_meta_sub$V1

  B <- DLPFC_minfi_beta[rowSums(is.na(DLPFC_minfi_beta))==0,]
  metasub <- DLPFC_minfi_ref_meta %>% dplyr::select(CellType,sex,ethnicity)
  fractrue <- DLPFC_minfi_truefrac
  cat("Guitivano")
  meta_sig <- data.frame(CellType = colnames(sig),Sample_Name =colnames(sig) )
  rownames(meta_sig) = colnames(sig)

  res_sig_gui <- EnsDeconv_wrap(B = B,frac_true = fractrue,
                                 sig = sig[intersect(rownames(sig),rownames(B)),],
                                include_minfi = include_minfi)

  cat("Gasparoni")
  B <- Gasparoni_minfi_beta[rowSums(is.na(Gasparoni_minfi_beta))==0, Gasparoni_minfi_ref_meta_sub$V1]
  metasub <- Gasparoni_minfi_ref_meta_sub %>% dplyr::select(CellType,AD)
  fractrue <- Gasparoni_truefrac

  res_sig_gas <- EnsDeconv_wrap(B = B,frac_true = fractrue,
                                 sig = sig[intersect(rownames(sig),rownames(B)),],
                                include_minfi = include_minfi)

  cat("GSE108066")

  # Deconvolve GSE108066

  nam <- str_split(files_name,"_")
  colnames(GSE108066_knnWkernel_1000_2) = paste0(sapply(nam, '[[',2),"_",sapply(nam, '[[',3))
  rownames(GSE108066_knnWkernel_1000_2) = meta450k_sub$TargetID
  meta_GSE108066 <- meta_GSE108066[match(colnames(GSE108066_knnWkernel_1000_2),meta_GSE108066$Sample),] %>% filter(Diagnosis == "CTL") %>% mutate(CellType = ifelse(Cell == "NeuN","NeuN_pos","NeuN_neg"))

  GSE108066_knnWkernel_1000_2_ctl <- GSE108066_knnWkernel_1000_2[,meta_GSE108066$Sample]

  GSE108066_truefrac <- matrix(0,45,2)
  colnames(GSE108066_truefrac) = unique(DLPFC_minfi_ref_meta$CellType)
  for (i in 1:nrow(GSE108066_truefrac)) {
    if(meta_GSE108066$CellType[i] == "NeuN_pos"){
      GSE108066_truefrac[i,1] = 1
    }else{
      GSE108066_truefrac[i,2] = 1
    }
  }



  B <- GSE108066_knnWkernel_1000_2_ctl[rowSums(is.na(GSE108066_knnWkernel_1000_2_ctl))< 25,]

  for (j in 1:nrow(B)) {
    for (i in which(is.na(B[j, ]))) {
      B[j, i] <- mean(B[j,which(meta_GSE108066$CellType == meta_GSE108066$CellType[i])],  na.rm = TRUE)
    }
  }

  rownames(GSE108066_truefrac) = colnames(B)
  true_frac <- GSE108066_truefrac

  res_sig_gse <- EnsDeconv_wrap(B = B,frac_true = GSE108066_truefrac,
                                 sig = sig[intersect(rownames(sig),rownames(B)),],
                                include_minfi = include_minfi)

  return(list(res_sig_ros = res_sig_ros, res_sig_gse= res_sig_gse,res_sig_gui = res_sig_gui,res_sig_gas = res_sig_gas))
}
