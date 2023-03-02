#' Multi_mrk is used to generate dataframe that saved the multi-markers
#'
#' @param DM_df dataframes,where the first column is TargetID (cpg), the second column is Refgene, and
#' the rest are DM pvalues for each cell type
#' @param eQTM_df dataframes,where the first column is TargetID (cpg), the second column is Refgene, and
#' the 3rd column is eQTM pvalues
#' @param DE_df dataframes,where the first column is Refgene, and the rest are DE pvalues for each cell type
#' @param padjust_method correction method, a character string. Can be abbreviated. See \cite{p.adjust}
#' @param ct_ind character vectors for cell type
#' @param pval_threshold numeric. pvalue threshold for count significant DM within gene
#' @param filter_type Character count or proportion
#' @param cutoff_val Must be consistent with filter_type
#' @param DE_type Character top or fdr
#' @param n_mrk Must be consistent with DE_type
#' @param beta_mtx
#' @param sig_name
#' @param probeInfo ROSMAP or 450k or 850k
#' @param Gene_group
#' @param B bulk data
#' @param frac_true
#' @param include_minfi
#' @param dmet_list
#' @param sc_signal pre-processed DM signals of all single cell and only kept 867 gene
#' @param fin_ncpg number of cpg per gene
#' @param cpg_pergene Logical indicator. Get cpg per gene or directly top cpgs
#'
#' @return
#' @importFrom ACAT ACAT
#' @importFrom MIND est_frac
#' @export
#'
#' @examples
Multi_mrk_new <- function(DM_df = DM_df,eQTM_df = eQTM_df,DE_df = DE_df,padjust_method = "bonferroni",
                      ct_ind = ct_ind,pval_threshold = 0.05,
                      filter_type = "prop",cutoff_val = 0.9,DE_type = "top",n_mrk = 100,
                      beta_mtx = sig_all,sig_name,probeInfo = "ROSMAP",Gene_group = "1stExon",
                      B = ros_dna_mtx_sub,frac_true = ROS_true_4ct,
                      include_minfi = TRUE,dmet_list = c("CIBERSORT","EPIC","FARDEEP"),
                      sc_signal = NULL,cpg_pergene = TRUE,fin_ncpg = 1, cutoff_val2 = 0.3){
  if(probeInfo == "ROSMAP"){
    probeInfo <- as.data.frame(data.table::fread('E:\\ROSMAP\\ROSMAP_arrayMethylation_metaData.tsv')) %>%
      dplyr::select(TargetID,GeneGroup=Gene_Feature)
  }else if(probeInfo == "450k"){
    probeInfo <- as.data.frame(EpiSCORE::probeInfo450k.lv)  %>%
      mutate(GeneGroup = recode(GeneGroup,
                                "1" = "TSS1500",
                                "2" = "TSS200",
                                "3" = "5UTR",
                                "4" = "1stExon",
                                "5" = "Body",
                                "6" = "3UTR")) %>% dplyr::select(TargetID = probeID,GeneGroup)
  }else{
    probeInfo <- as.data.frame(EpiSCORE::probeInfo850k.lv)  %>%
      mutate(GeneGroup = recode(GeneGroup,
                                "1" = "TSS1500",
                                "2" = "TSS200",
                                "3" = "5UTR",
                                "4" = "1stExon",
                                "5" = "Body",
                                "6" = "3UTR")) %>% dplyr::select(TargetID = probeID,GeneGroup)
  }

  # First combine DM and eQTM
  suppressMessages(com1 <- left_join(DM_df,eQTM_df) %>% left_join(probeInfo) %>%
                     filter(complete.cases(.)))
  eqval <- com1$eQTM_pval
  com1 <-com1 %>% ungroup() %>% mutate_at({{ct_ind}},function(p) {
    p[p>1]=1
    return(p)}) %>%
    mutate_at({{ct_ind}},function(p) {
      if(sum(p==1) >0 & sum(p==0) >0){
        p[p==0] = 1e-130
      }
      ACAT(t(matrix(c(p,eqval),ncol = 2)))})

  # Adjust pvalue
  com1 <-com1 %>%   mutate_at({{ct_ind}},function(p) {p.adjust(p,method = padjust_method)})
  if(is.null(sc_signal)){
    # Count significant DM within gene
    if(filter_type == "prop"){
      com2 <- com1 %>% group_by(RefGene) %>% mutate(Num_cpg = n()) %>%
        mutate(across(ct_ind, function(i){
          sum(i < pval_threshold)/Num_cpg
        }, .names = "{.col}_DMR"))
    }else{
      com2 <- com1 %>% group_by(RefGene) %>%
        mutate(across(ct_ind, function(i){
          sum(i < pval_threshold)
        }, .names = "{.col}_DMR"))
    }
  }else{
    if(filter_type == "prop"){
      sub_signal <- sc_signal %>% group_by(RefGene) %>%mutate(Num_cpg = n()) %>%
        mutate(across(ct_ind, function(i){
          i <- p.adjust(i,padjust_method)
          sum(i < pval_threshold)/Num_cpg
        }, .names = "{.col}_DMR"))%>% dplyr::select(-ct_ind) %>% unique()
    }else if(filter_type == "count"){
      sub_signal <- sc_signal %>% group_by(RefGene) %>%
        mutate(across(ct_ind, function(i){
          i <- p.adjust(i,padjust_method)
          sum(i < pval_threshold)
        }, .names = "{.col}_DMR"))%>% dplyr::select(-ct_ind) %>% unique()
    }else{
      sub_signal <- sc_signal %>% group_by(RefGene)%>% mutate(Num_cpg = n()) %>%
        mutate(across(ct_ind, function(i){
          i <- p.adjust(i,padjust_method)
          sum(i < pval_threshold)
        }, .names = "{.col}_DMR")) %>%
        mutate(across(ct_ind, function(i){
          i <- p.adjust(i,padjust_method)
          sum(i < pval_threshold)/Num_cpg
        }, .names = "{.col}_DMR2"))%>% dplyr::select(-ct_ind) %>% unique()
    }

    com2 <- right_join(com1,sub_signal)
  }


  output_mrk <- lapply(ct_ind, function(i){

    sub = paste0(i,"_DMR")
    sub2 = paste0(i,"_DMR2")
    sub_de = paste0(i,"_DE")

    # First filter DE:
    DE_df_sub <- DE_df %>% arrange(.data[[sub_de]]) %>%
      dplyr::slice(1:n_mrk)
    com3 <- com2 %>% ungroup() %>% left_join(DE_df_sub)%>%
      filter(complete.cases(.))

    #filter(.data[[sub]]>cutoff_val & .data[[i]] <pval_threshold)%>%

    if(cpg_pergene){

      if(filter_type %in% c("count","prop")){
      com4 <- com3 %>% dplyr::select(TargetID,RefGene, starts_with(i),GeneGroup) %>% ungroup() %>%
        filter(.data[[sub]]>cutoff_val& .data[[i]] <pval_threshold)%>%
        filter( str_detect(GeneGroup,as.character(Gene_group))) %>%
        group_by(RefGene) %>%
        slice_min(.data[[i]], n = fin_ncpg) %>% ungroup() %>% arrange(.data[[i]]) %>%
        dplyr::slice(1:n_mrk)
      }else{
        com4 <- com3 %>% dplyr::select(TargetID,RefGene, starts_with(i),GeneGroup) %>% ungroup() %>%
          filter(.data[[sub]]>cutoff_val&.data[[sub2]]>cutoff_val2& .data[[i]] <pval_threshold)%>%
          filter( str_detect(GeneGroup,as.character(Gene_group))) %>%
          group_by(RefGene) %>%
          slice_min(.data[[i]], n = fin_ncpg) %>% ungroup() %>% arrange(.data[[i]]) %>%
          dplyr::slice(1:n_mrk)
      }
    }else{
      if(filter_type %in% c("count","prop")){
      com4 <- com3 %>% dplyr::select(TargetID,RefGene, starts_with(i),GeneGroup) %>% ungroup() %>%
        filter(.data[[sub]]>cutoff_val& .data[[i]] <pval_threshold)%>%
        filter( str_detect(GeneGroup,as.character(Gene_group))) %>% arrange(RefGene) %>%
        dplyr::slice(1:n_mrk)
      }else{
        com4 <- com3 %>% dplyr::select(TargetID,RefGene, starts_with(i),GeneGroup) %>% ungroup() %>%
          filter(.data[[sub]]>cutoff_val&.data[[sub2]]>cutoff_val2& .data[[i]] <pval_threshold)%>%
          filter( str_detect(GeneGroup,as.character(Gene_group))) %>% arrange(RefGene) %>%
          dplyr::slice(1:n_mrk)
      }
    }


    return(com4$TargetID)
  })
  cat("New first alg: Number of markers per ct: ", sapply(output_mrk,length))
  output <- unique(unlist(output_mrk))
  #output <- output[!duplicated(output)]

  sig <- beta_mtx[intersect(output,rownames(beta_mtx)),]
  return(list(sig = sig,output_mrk = output_mrk))
  # if(nrow(sig)<4){
  #   return(list(sig = sig, res_ROS = NULL,
  #               mean_cor = NA,output_mrk,
  #               bench = NULL))
  # }else{
  # #stopifnot(nrow(sig) > 50)
  # cat("New first alg: Number of markers: ", nrow(sig))
  # ref_list = list(list(ref_matrix = sig, meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
  #                 data_name = sig_name))
  # names(ref_list)[1] = sig_name
  # suppressWarnings(res_ROS <- EnsDeconv(count_bulk = B,
  #                     ref_list  = ref_list,
  #                     params = get_params(data_type = "singlecell-rna", data_name = sig_name,
  #                                         n_markers = 50, Marker.Method = "none",
  #                                         TNormalization = "none", CNormalization = "none",
  #                                         Scale = "linear", dmeths = dmet_list),
  #                     outpath = "D:/Manqi/tmp1/"))
  # frac_est = est_frac(sig, B)
  #
  # new_list <-res_ROS[["allgene_res"]][[1]]
  #
  # p <- res_ROS[["allgene_res"]][[1]]$p
  # p$dmeths = "NNLS"
  # a <- res_ROS[["allgene_res"]][[1]][["a"]]
  # a[["p_hat"]][[1]][[1]] = frac_est
  # names(a[["p_hat"]][[1]]) = "NNLS"
  # new_list <- list(a = a, p =p)
  #
  # res_ROS[["allgene_res"]] <- append(res_ROS[["allgene_res"]],list(new_list))
  #
  # if(include_minfi){
  #   minfi_res <- DNAm_minfi_deconv_general(dat = sig,metaref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
  #                                          bulk =B,true_frac = NULL,findMrks = F)
  #
  #   p <- res_ROS[["allgene_res"]][[1]]$p
  #   p$dmeths = "minfi"
  #   p$data_name = sig_name
  #   a <- res_ROS[["allgene_res"]][[1]][["a"]]
  #   a[["p_hat"]][[1]][[1]] = minfi_res[[2]]$counts
  #   names(a[["p_hat"]][[1]]) = "minfi"
  #   new_list <- list(a = a, p =p)
  #
  #   res_ROS[["allgene_res"]] <- append(res_ROS[["allgene_res"]],list(new_list))
  #
  #   p <- res_ROS[["allgene_res"]][[1]]$p
  #   p$dmeths = "minfi_M"
  #   p$data_name = sig_name
  #   a <- res_ROS[["allgene_res"]][[1]][["a"]]
  #   a[["p_hat"]][[1]][[1]] = minfi_res[[1]]$counts
  #   names(a[["p_hat"]][[1]]) = "minfi_M"
  #   new_list <- list(a = a, p =p)
  #
  #
  #   res_ROS[["allgene_res"]] <- append(res_ROS[["allgene_res"]],list(new_list))
  # }
  #
  # bench <- sig_benchmarking(sig)
  #
  #
  # #message(paste0("Mean cor: ",mean(diag(cor(res_ROS$EnsDeconv$ensemble_p, frac_true, method = 's')))))
  # return(list(sig = sig, res_ROS = res_ROS,
  #             mean_cor = mean(diag(cor(res_ROS$EnsDeconv$ensemble_p, frac_true, method = 's'))),output_mrk,
  #             bench = bench))
  # }
}


