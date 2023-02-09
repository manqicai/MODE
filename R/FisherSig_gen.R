#' FisherSig_gen
#'
#' @param cpgpath
#'
#' @return
#' @export
#'
#' @examples
FisherSig_gen <- function(cpgpath ="E:/Luo2022/cpg_snm3C.rds"){
  celltype_all = openxlsx::read.xlsx('https://ars.els-cdn.com/content/image/1-s2.0-S2666979X22000271-mmc7.xlsx', startRow = 3)
  rownames(celltype_all) = celltype_all$Cell.ID
  celltype_all$Cell.ID = NULL
  celltype_all$cellType = celltype_all$MajorType
  celltype_all$cellType[grep('Exc', celltype_all$cellType)] = 'Exc'
  celltype_all$cellType[grep('Inh', celltype_all$cellType)] = 'Inh'
  celltype_all$cellType[grep('Oligo', celltype_all$cellType)] = 'Oligo'
  celltype_all$cellType[grep('Astro', celltype_all$cellType)] = 'Astro'
  celltype_all$cellType[grep('Micro', celltype_all$cellType)] = 'Micro-Endo'
  celltype_all$cellType[grep('Outlier', celltype_all$cellType)] = NA # n = 88
  celltype_all = celltype_all[!is.na(celltype_all$cellType),]
  celltype_all$tech = sapply(rownames(celltype_all), function(x) unlist(strsplit(x, '_'))[2])
  celltype_all$tech[celltype_all$tech == 'NDARKD326LNK'] = 'snmC-seq'
  celltype_all$tech[celltype_all$tech == 'NDARKJ183CYT'] = 'snmC-seq2'
  #celltype_all$tech[celltype_all$tech == 'UMB4540'] = ' snmC-seq'
  celltype_all$tech[celltype_all$tech == ' snmC-seq'] = 'snmC-seq'

  celltype_all <- celltype_all %>% filter(tech != "UMB4540")
  cpg_list <- readRDS(cpgpath)
  ncpg = sapply(cpg_list, nrow)

  celltype_all = celltype_all[names(cpg_list),]
  # hist(log10(ncpg))
  # boxplot(ncpg ~ celltype_all$cellType)
  celltype_all$Missing = ncpg
  # 2886 cells have no 450k coverage
  table(ncpg!=0, celltype_all$tech)


  cpg_list = cpg_list[ncpg != 0]
  celltype_all = celltype_all[names(cpg_list),]
  ncpg = ncpg[names(cpg_list)]
  meta = as.data.frame(fread('E:\\ROSMAP\\ROSMAP_arrayMethylation_metaData.tsv'))
  meta$CHR = as.character(meta$CHR)
  meta <- as.data.frame(meta)
  rownames(meta) = paste0(meta$CHR, ":", meta$MAPINFO)
  meth = coverage = matrix(NA, nrow(meta), length(cpg_list))
  rownames(meth) = rownames(coverage) = paste0(meta$CHR, ":", meta$MAPINFO)
  colnames(meth) = colnames(coverage) = names(cpg_list)

  for(i in 1:length(cpg_list)) {
    if(i %% 100 == 0) print(i)
    id = paste0(cpg_list[[i]]$CHR, ":", cpg_list[[i]]$MAPINFO)
    meth[id, i] = cpg_list[[i]]$V5
    coverage[id, i] = cpg_list[[i]]$V6
  }

  id = which(rowSums(!is.na(meth)) != 0)
  coverage = coverage[id,]
  meth = meth[id,]

  beta = meth/coverage
  mean(beta != 0 & beta != 1, na.rm = T)

  subid <- rownames(meth)
  meta = meta[subid,]

  meta$RefGene = sapply(meta$RefGene, function(x) paste0(unique(unlist(strsplit(x, ';'))), collapse = ';'))
  celltype_all$celltype = celltype_all$cellType
  celltype_all$celltype[celltype_all$celltype%in% c('Exc', 'Inh')] = 'Neuro'
  meth4 = cbind(rowSums(meth[,celltype_all$celltype=='Astro'], na.rm = T),
                rowSums(meth[,celltype_all$celltype=='Micro-Endo'], na.rm = T),
                rowSums(meth[,celltype_all$celltype=='Neuro'], na.rm = T),
                rowSums(meth[,celltype_all$celltype=='Oligo'], na.rm = T))

  coverage4 = cbind(rowSums(coverage[,celltype_all$celltype=='Astro'], na.rm = T),
                    rowSums(coverage[,celltype_all$celltype=='Micro-Endo'], na.rm = T),
                    rowSums(coverage[,celltype_all$celltype=='Neuro'], na.rm = T),
                    rowSums(coverage[,celltype_all$celltype=='Oligo'], na.rm = T))

  colnames(coverage4) = colnames(meth4) = c('Astro', 'Micro_Endo', 'Neuro', 'Oligo_MBP')

  Mrk <- lapply(1:4, function(j){
    dnam_pval <- sapply(1:nrow(meth4), function(i) {
      m = c(meth4[i,j], sum(meth4[i,-j]))
      n = c(coverage4[i,j], sum(coverage4[i,-j]))
      if(all(n > 0)) return(fisher.test(cbind(m, n - m), alternative = 'less')$p.value) else return(NA)
    })
    names(dnam_pval) = rownames(meth4)
    return(sort(dnam_pval))
  })

  Mrk_twosided <- lapply(1:4, function(j){
    dnam_pval <- sapply(1:nrow(meth4), function(i) {
      m = c(meth4[i,j], sum(meth4[i,-j]))
      n = c(coverage4[i,j], sum(coverage4[i,-j]))
      if(all(n > 0)) return(fisher.test(cbind(m, n - m), alternative = 'two.sided')$p.value) else return(NA)
    })
    names(dnam_pval) = rownames(meth4)
    return(sort(dnam_pval))
  })
  return(list(Mrk = Mrk, Mrk_twosided= Mrk_twosided))
}
