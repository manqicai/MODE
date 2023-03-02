#' Multi_mrk_wrapper  Wrapper function for multi-markers deconvolution with EnsDeconv
#'
#' @param com_params
#' @param my_DM_df
#' @param my_eQTM_df
#' @param my_DE_df
#' @param my_DM_df_onesided
#' @param my_DE_df_onesided
#' @param ct_ind
#' @param beta_mtx
#' @param B
#' @param frac_true
#' @param ncore
#' @param include_minfi
#' @param dmet_list
#' @param probeInfo
#'
#' @return
#' @export
#'
#' @examples
Multi_mrk_wrapper <- function(com_params,my_DM_df,my_eQTM_df,my_DE_df,my_DM_df_onesided,
                              my_DE_df_onesided,ct_ind = c("Astro","Micro_Endo", "Neuro","Oligo_MBP"),
                              beta_mtx = sig_all,
                              B = ros_dna_mtx_sub,frac_true = ROS_true_4ct,
                              ncore = 20,dmet_list = c("CIBERSORT","EPIC","FARDEEP"),probeInfo = "ROSMAP",padjust_method = "bonferroni"){

  pb <- progress_bar$new(
    format = "Current : :current [:bar] :elapsed | percent: :percent",
    total = nrow(com_params),
    clear = FALSE,
    force = TRUE,
    width = 60)

  progress_letter <- rep(1:10, 10)  # token reported in progress bar

  # allowing progress bar to be used in foreach -----------------------------
  progress <- function(n){
    pb$tick(tokens = list(letter = progress_letter[n]))
  }


  opts <- list(progress = progress)

  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  res <- foreach(rownum = 1:nrow(com_params), .export=c("Multi_mrk","Multi_mrk2","Multi_mrk3",
                                                        "DNAm_minfi_deconv_general","MyestimateCellCounts",
                                                        "pickCompProbes","projectCellType","data.table"),.packages=c("dplyr","reshape2","tidyverse","EnsDeconv","ACAT","sparseMatrixStats","tidyverse","MIND","minfi","bumphunter","quadprog","EpiSCORE"),.options.snow=opts, .errorhandling = 'pass') %dopar% {
    sig_name = paste0(unlist(paste0(colnames(com_params),"_", com_params[rownum,])), collapse = "")
    source("D:/Manqi/Meeting/20221027/utilities.R")
    source("D:/Manqi/Meeting/20221023/utilities.R")
    source("D:/Manqi/Package/MODE/R/Multi_mrk.R")
    source("D:/Manqi/Package/MODE/R/Multi_mrk2.R")
    source("D:/Manqi/Package/MODE/R/Multi_mrk3.R")
    source("D:/Manqi/EnsDeconv/R/CTS_EnsDeconv_wrapper.R")
    source("D:/Manqi/EnsDeconv/R/CTS_EnsDeconv_LS.R")
    source("D:/Manqi/deconv.ensemble/R/process_all.R")
    source("D:/Manqi/deconv.ensemble/R/analysis_wrap.R")
    source("D:/Manqi/deconv.ensemble/R/ensemble_prep.R")
    source("D:/Manqi/deconv.ensemble/R/analyze_dset.R")
    source("D:/Manqi/deconv.ensemble/R/get_metric_wrap.R")
    output <- Multi_mrk(DM_df = my_DM_df,eQTM_df = my_eQTM_df,DE_df = my_DE_df,padjust_method = padjust_method,
                        ct_ind = ct_ind,pval_threshold = com_params$pval_threshold[rownum],
                        filter_type = com_params$filter_type[rownum],cutoff_val = com_params$cutoff_val[rownum],DE_type = com_params$DE_type[rownum],cutoff_DE = com_params$cutoff_DE[rownum],beta_mtx = beta_mtx,
                        sig_name = sig_name, Gene_group =com_params$Gene_group[rownum],B = B,frac_true = frac_true,
                        include_minfi = com_params$include_minfi[rownum], dmet_list = dmet_list,probeInfo = probeInfo)
    output_onesided <- Multi_mrk(DM_df = my_DM_df_onesided,eQTM_df = my_eQTM_df,DE_df = my_DE_df_onesided,padjust_method = padjust_method,
                                 ct_ind = ct_ind,pval_threshold = com_params$pval_threshold[rownum],
                                 filter_type = com_params$filter_type[rownum],cutoff_val = com_params$cutoff_val[rownum],DE_type = com_params$DE_type[rownum],cutoff_DE = com_params$cutoff_DE[rownum],beta_mtx = beta_mtx,
                                 sig_name = sig_name, Gene_group =com_params$Gene_group[rownum],B = B,frac_true = frac_true,
                                 include_minfi = com_params$include_minfi[rownum], dmet_list = dmet_list,probeInfo = probeInfo )

    output2 <- Multi_mrk2(DM_df = my_DM_df,eQTM_df = my_eQTM_df,DE_df = my_DE_df,padjust_method = padjust_method,
                          ct_ind = ct_ind,pval_threshold = com_params$pval_threshold[rownum],
                          filter_type = com_params$filter_type[rownum],cutoff_val = com_params$cutoff_val[rownum],DE_type = com_params$DE_type[rownum],cutoff_DE = com_params$cutoff_DE[rownum],beta_mtx = beta_mtx,
                          sig_name = sig_name, Gene_group =com_params$Gene_group[rownum],B = B,frac_true = frac_true,
                          include_minfi =  com_params$include_minfi[rownum] , dmet_list = dmet_list,probeInfo = probeInfo )
    output_onesided2 <- Multi_mrk2(DM_df = my_DM_df_onesided,eQTM_df = my_eQTM_df,DE_df = my_DE_df_onesided,padjust_method = padjust_method,
                                   ct_ind = ct_ind,pval_threshold = com_params$pval_threshold[rownum],
                                   filter_type = com_params$filter_type[rownum],cutoff_val = com_params$cutoff_val[rownum],DE_type = com_params$DE_type[rownum],cutoff_DE = com_params$cutoff_DE[rownum],beta_mtx = beta_mtx,
                                   sig_name = sig_name, Gene_group =com_params$Gene_group[rownum],B = B,frac_true = frac_true,
                                   include_minfi = com_params$include_minfi[rownum], dmet_list = dmet_list,probeInfo = probeInfo )
    output3 <- Multi_mrk3(DM_df = my_DM_df,eQTM_df = my_eQTM_df,DE_df = my_DE_df,padjust_method = padjust_method,
                          ct_ind = ct_ind,pval_threshold = com_params$pval_threshold[rownum],
                          filter_type = com_params$filter_type[rownum],cutoff_val = com_params$cutoff_val[rownum],DE_type = com_params$DE_type[rownum],cutoff_DE = com_params$cutoff_DE[rownum],beta_mtx = beta_mtx,
                          sig_name = sig_name, Gene_group =com_params$Gene_group[rownum],B = B,frac_true = frac_true,
                          include_minfi = com_params$include_minfi[rownum], dmet_list = dmet_list,probeInfo = probeInfo )
    output_onesided3 <- Multi_mrk3(DM_df = my_DM_df_onesided,eQTM_df = my_eQTM_df,DE_df = my_DE_df_onesided,padjust_method = padjust_method,
                                   ct_ind = ct_ind,pval_threshold = com_params$pval_threshold[rownum],
                                   filter_type = com_params$filter_type[rownum],cutoff_val = com_params$cutoff_val[rownum],DE_type = com_params$DE_type[rownum],cutoff_DE = com_params$cutoff_DE[rownum],beta_mtx = beta_mtx,
                                   sig_name = sig_name, Gene_group =com_params$Gene_group[rownum],B = B,frac_true = frac_true,
                                   include_minfi = com_params$include_minfi[rownum], dmet_list = dmet_list,probeInfo = probeInfo )
    gc()
    return(list(output = output,output_onesided = output_onesided,
                output2 = output2,output_onesided2 = output_onesided2,
                output3 = output3,output_onesided3 = output_onesided3))
  }
  stopCluster(cl)
  gc()
  complete_res = which(sapply(res,class) == 'list')
  com_params_sub = com_params[complete_res,]
  res = res[complete_res]
  print_text <- c("1st alg 2-sided","1st alg 1-sided","2nd alg 2-sided","2nd alg 1-sided","3rd alg 2-sided","3rd alg 1-sided")
  #p1list <- list()
  plist <- list()
  sub_res2list <- list()
  sub_reslist <- list()
  for (j in 1:6) {
    sub_reslist[[j]] <- sapply(res, function(i){
      i[[j]][[3]]
    })

    sub_res2list[[j]]<- lapply(res, function(i){
      i[[j]][[2]][[2]]
    })

    sub_res2list[[j]] <- unlist(sub_res2list[[j]],recursive = F)
    ind = sapply(sub_res2list[[j]], function(x){
      length(x[["a"]][["p_hat"]][[1]])
    })
    sub_res2list[[j]] = sub_res2list[[j]][which(ind == 1)]
    EnsDeconv_p = CTS_EnsDeconv_wrapper(sub_res2list[[j]])
    print(paste0(print_text[j],range(sub_reslist[[j]])))
    # p1list[[j]] <- plot(sub_reslist[[j]])

    new_list <-sub_res2list[[j]][[1]]

    p <- sub_res2list[[j]][[1]]$p
    p$dmeths = "EnsDeconv"
    a <- sub_res2list[[j]][[1]][["a"]]
    a[["p_hat"]][[1]][[1]] = EnsDeconv_p$ensemble_p
    names(a[["p_hat"]][[1]]) = "EnsDeconv"
    new_list <- list(a = a, p =p)

    sub_res2list[[j]] <- append(sub_res2list[[j]],list(new_list))

    # ggerrorplot
    cl <- makeCluster(20)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    res_ROS_out_all <- foreach(i = 1:length(sub_res2list[[j]]), .export=c("get_metric_wrap", "ExportRes_all_new","getOutput_ensemble_cv","ord_name","sum_to_one"),.packages=c("dplyr","reshape2","tidyverse")) %dopar% {
      ress <- get_metric_wrap(sub_res2list[[j]][[i]], true = frac_true,trueMet ="ROS" )
      ress <- ExportRes_all_new(ress)
      return(ress)
    }
    stopCluster(cl)

    res_ROS_out_all <- bind_rows(res_ROS_out_all)
    celltype_col = pal_lancet("lanonc")(4)
    new = res_ROS_out_all %>% group_by(CellType,Method) %>%
      mutate(Spearman_mean = mean(Spearman,na.rm = T),
             Method = ifelse(Method == "EnsDeconv","EnsDeconv",Method))

    checkkkk =new %>% group_by(CellType,Method) %>%
      mutate(p_mean = mean(Spearman_mean,na.rm = T),
             Method = ifelse(Method == "EnsDeconv","EnsDeconv",Method)) %>%
      dplyr::filter(!str_detect(Method,"fold|allgene|mean|Mean"))

    plist[[j]] <- new %>% ggerrorplot("Method", "Spearman_mean",order = levels(reorder(checkkkk$Method, checkkkk$p_mean, median, order = TRUE)),
                        desc_stat = "median_mad", error.plot = "errorbar")+
      geom_point(aes(color =CellType),size = 4)+      scale_color_manual(
        values=celltype_col)+ theme(legend.title = element_blank(),
                                    legend.text=element_text(size=12),plot.title = element_text(hjust = 0.5))+
      xlab("Deconvolution Method") + ylab("")+
      stat_summary(
        geom = "point",
        shape = 3,
        size = 5,
        fun = "median",show.legend = F)+
      coord_flip()+ggtitle(print_text[j])

    }
  return(list(res = res, plist = plist,sub_reslist = sub_reslist,sub_res2list = sub_res2list,com_params_sub = com_params_sub ))
}
