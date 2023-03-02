#' Multi_mrk is used to generate dataframe that saved the multi-markers
#'
#' @param B
#' @param frac_true
#' @param include_minfi
#' @param dmet_list
#'
#' @return
#' @importFrom ACAT ACAT
#' @importFrom MIND est_frac
#' @export
#'
#' @examples
EnsDeconv_wrap <- function(B = ros_dna_mtx_sub,frac_true = ROS_true_4ct,
                      include_minfi = TRUE, sig ,dmet_list = c("CIBERSORT","EPIC")){
  # source("D:/Manqi/Meeting/20221027/utilities.R")
  # source("D:/Manqi/Meeting/20221023/utilities.R")
  # source("D:/Manqi/Package/MODE/R/Multi_mrk.R")
  # source("D:/Manqi/Package/MODE/R/Multi_mrk2.R")
  # source("D:/Manqi/Package/MODE/R/Multi_mrk3.R")
  # source("D:/Manqi/EnsDeconv/R/CTS_EnsDeconv_wrapper.R")
  # source("D:/Manqi/EnsDeconv/R/CTS_EnsDeconv_LS.R")
  # source("D:/Manqi/deconv.ensemble/R/process_all.R")
  # source("D:/Manqi/deconv.ensemble/R/analysis_wrap.R")
  # source("D:/Manqi/deconv.ensemble/R/ensemble_prep.R")
  # source("D:/Manqi/deconv.ensemble/R/analyze_dset.R")
  # source("D:/Manqi/deconv.ensemble/R/get_metric_wrap.R")

  #message(paste0("Number of markers: ", nrow(sig)))
  ref_list = list(list(ref_matrix = sig, meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                       data_name = "sig"))
  names(ref_list)[1] = "sig"
  suppressWarnings(res_ROS <- EnsDeconv(count_bulk = B,
                                        ref_list  = ref_list,
                                        params = get_params(data_type = "singlecell-rna", data_name = "sig",
                                                            n_markers = 50, Marker.Method = "none",
                                                            TNormalization = "none", CNormalization = "none",
                                                            Scale = "linear", dmeths = dmet_list),
                                        outpath = "D:/Manqi/tmp/"))
  frac_est = est_frac(sig, B)

  new_list <-res_ROS[["allgene_res"]][[1]]

  p <- res_ROS[["allgene_res"]][[1]]$p
  p$dmeths = "NNLS"
  a <- res_ROS[["allgene_res"]][[1]][["a"]]
  a[["p_hat"]][[1]][[1]] = frac_est
  names(a[["p_hat"]][[1]]) = "NNLS"
  new_list <- list(a = a, p =p)

  res_ROS[["allgene_res"]] <- append(res_ROS[["allgene_res"]],list(new_list))

  if(include_minfi){
    minfi_res <- DNAm_minfi_deconv_general(dat = sig,metaref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                                           bulk =B,true_frac = NULL,findMrks = F)

    p <- res_ROS[["allgene_res"]][[1]]$p
    p$dmeths = "minfi"
    p$data_name = "sig"
    a <- res_ROS[["allgene_res"]][[1]][["a"]]
    a[["p_hat"]][[1]][[1]] = minfi_res[[2]]$counts
    names(a[["p_hat"]][[1]]) = "minfi"
    new_list <- list(a = a, p =p)

    res_ROS[["allgene_res"]] <- append(res_ROS[["allgene_res"]],list(new_list))

    p <- res_ROS[["allgene_res"]][[1]]$p
    p$dmeths = "minfi_M"
    p$data_name = "sig"
    a <- res_ROS[["allgene_res"]][[1]][["a"]]
    a[["p_hat"]][[1]][[1]] = minfi_res[[1]]$counts
    names(a[["p_hat"]][[1]]) = "minfi_M"
    new_list <- list(a = a, p =p)


    res_ROS[["allgene_res"]] <- append(res_ROS[["allgene_res"]],list(new_list))
  }

  EnsDeconv_p = CTS_EnsDeconv_wrapper(res_ROS[["allgene_res"]])
  phat_list <- lapply(res_ROS[["allgene_res"]], function(i){
    as.matrix(i[["a"]][["p_hat"]][[1]][[1]])
  })

  phat_list <- append(phat_list,list(EnsDeconv_p$ensemble_p))

  res_sig <- sapply(phat_list, Get_resTbl,frac_true,inputList = F)
  cat("Result for sig to deconvolve")
  print(res_sig)

  #print(Get_resTbl(EnsDeconv_p$ensemble_p,true_frac = frac_true,inputList = F))

  #message(paste0("Mean cor: ",mean(diag(cor(res_ROS$EnsDeconv$ensemble_p, frac_true, method = 's')))))
  return(res_sig)
}
