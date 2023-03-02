#' Multi_mrk2 is used to generate dataframe that saved the multi-markers, it filter DMR first and then
#' combine DM, DE and eQTM
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
#' @param cutoff_DE Must be consistent with DE_type
#' @param beta_mtx
#' @param sig_name
#' @param probeInfo
#' @param Gene_group
#' @param B
#' @param frac_true
#' @param include_minfi
#'
#' @return
#' @importFrom ACAT ACAT
#' @importFrom MIND est_frac
#' @export
#'
#' @examples
Multi_mrk3 <- function(DM_df = DM_df,eQTM_df = eQTM_df,DE_df = DE_df,padjust_method = "bonferroni",
                       ct_ind = ct_ind,pval_threshold = 0.05,
                       filter_type = "prop",cutoff_val = 0.9,DE_type = "top",cutoff_DE = 100,
                       beta_mtx = sig_all,sig_name,probeInfo = "ROSMAP",Gene_group = "1stExon",
                       B = ros_dna_mtx_sub,frac_true = ROS_true_4ct,
                       include_minfi = TRUE,dmet_list = c("CIBERSORT","EPIC","FARDEEP"),cutoff_prop = 0){
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

  suppressMessages(com1 <- left_join(DM_df,eQTM_df)%>% left_join(DE_df)%>% left_join(probeInfo)%>%
                     filter(complete.cases(.)))


  # First count DMR using DM
  if(filter_type == "prop"){
    com1 <- com1 %>% group_by(RefGene) %>% mutate(Num_cpg = n()) %>%
      mutate(across(ct_ind, function(i){
        sum(p.adjust(i,method = 'fdr') < pval_threshold)/Num_cpg
      }, .names = "{.col}_DMR"))
  }else if(filter_type == "count"){
    com1 <- com1 %>% group_by(RefGene) %>%
      mutate(across(ct_ind, function(i){
        sum(p.adjust(i,method = 'fdr')  < pval_threshold)
      }, .names = "{.col}_DMR"))
  }else{
    com1 <- com1 %>% group_by(RefGene)%>% mutate(Num_cpg = n())  %>%
      mutate(eQTM_DMR = sum(eQTM_pval < pval_threshold),
             eQTM_DMR_prop = sum(eQTM_pval < pval_threshold)/Num_cpg)%>%
      filter(eQTM_DMR>cutoff_val|eQTM_DMR_prop > cutoff_prop )
  }

  # Second combine DM, DE and eQTM

  com2 <-com1 %>% ungroup() %>% mutate_at({{ct_ind}},function(p) {
    p[p>1]=1
    return(p)})
  for (i in ct_ind) {
    com2 <- com2 %>% filter(!((.data[[paste0(i,"_DE")]]==1 & .data[[i]] ==0)|
                                (.data[[paste0(i,"_DE")]]==0 & .data[[i]] ==1)))
    eqval <- com2$eQTM_pval
    com2[[i]]=ACAT(t(matrix(c(com2[[i]],eqval,com2[[paste0(i,"_DE")]]),ncol = 3)))
  }


  # Third Adjust pvaluec & filter DMR and select top markers per cell type
  com3 <-com2 %>%   mutate_at({{ct_ind}},function(p) {p.adjust(p,method = padjust_method)})


  output_mrk <- lapply(ct_ind, function(i){
    sub = paste0(i,"_DMR")
    sub_de = paste0(i,"_DE")

    fin <- com3 %>% dplyr::select(TargetID,RefGene, starts_with(i),GeneGroup) %>% ungroup() %>%
      filter(.data[[sub]]>cutoff_val& .data[[i]] <pval_threshold)%>%
      filter( str_detect(GeneGroup,as.character(Gene_group)))  %>% arrange(.data[[i]]) %>%
      dplyr::slice(1:cutoff_DE)
    return(fin$TargetID)
  })
  output <- unlist(output_mrk)
  output <- output[!duplicated(output)]

  sig <- beta_mtx[intersect(output,rownames(beta_mtx)),]
  #message(paste0("Number of markers: ", nrow(sig)))
  ref_list = list(list(ref_matrix = sig, meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                       data_name = sig_name))
  names(ref_list)[1] = sig_name
  suppressWarnings(res_ROS <- EnsDeconv(count_bulk = B,
                                        ref_list  = ref_list,
                                        params = get_params(data_type = "singlecell-rna",
                                                            data_name = sig_name, n_markers = 50,
                                                            Marker.Method = "none",
                                                            TNormalization = "none",
                                                            CNormalization = "none",
                                                            Scale = "linear",
                                                            dmeths = dmet_list ),
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
    p$data_name = sig_name
    a <- res_ROS[["allgene_res"]][[1]][["a"]]
    a[["p_hat"]][[1]][[1]] = minfi_res[[2]]$counts
    names(a[["p_hat"]][[1]]) = "minfi"
    new_list <- list(a = a, p =p)

    res_ROS[["allgene_res"]] <- append(res_ROS[["allgene_res"]],list(new_list))

    p <- res_ROS[["allgene_res"]][[1]]$p
    p$dmeths = "minfi_M"
    p$data_name = sig_name
    a <- res_ROS[["allgene_res"]][[1]][["a"]]
    a[["p_hat"]][[1]][[1]] = minfi_res[[1]]$counts
    names(a[["p_hat"]][[1]]) = "minfi_M"
    new_list <- list(a = a, p =p)


    res_ROS[["allgene_res"]] <- append(res_ROS[["allgene_res"]],list(new_list))
  }
  #message(paste0("Mean cor: ",mean(diag(cor(res_ROS$EnsDeconv$ensemble_p, frac_true, method = 's')))))
  return(list(sig = sig, res_ROS = res_ROS,
              mean_cor = mean(diag(cor(res_ROS$EnsDeconv$ensemble_p, frac_true, method = 's'))),output_mrk))
}

