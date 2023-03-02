gen_cpg <- function(cell_dir,meta = meta_850k){

  file = dir(cell_dir, '.gz', recursive = T)
  
  cell = sapply(file, function(x) unlist(strsplit(x, 'allc_'))[2])
  cell = gsub('.tsv.gz', '', cell)
  
  # table(cell %in% rownames(celltype_all)[!is.na(celltype_all$cellType)])
  # rownames(celltype_all)[!rownames(celltype_all) %in% cell & celltype_all$tech == 'mCTseq']
  # cell[!cell %in% rownames(celltype_all)]
  rm(file)
  
  cell = cell[cell %in% rownames(celltype_all)]
  # celltype_all = celltype_all[cell,]
  
  library(data.table)
  library(doParallel)
  library(foreach)
  library(doSNOW)
  library(progress)
  # rm(celltype_all); gc()
  
  
  pb <- progress_bar$new(
    format = "Current : :current [:bar] :elapsed | percent: :percent",
    total = length(cell),
    clear = FALSE,
    force = TRUE,
    width = 60)
  
  progress_letter <- rep(1:10, 10)  # token reported in progress bar
  
  # allowing progress bar to be used in foreach -----------------------------
  progress <- function(n){
    pb$tick(tokens = list(letter = progress_letter[n]))
  }
  
  
  opts <- list(progress = progress)
  cl <- makeCluster(80)
  registerDoSNOW(cl)
  
  system.time(cpg <- foreach(i = names(cell), .packages = c("dplyr","data.table","RANN"),.options.snow=opts, .errorhandling = 'pass') %dopar% {
    setwd(cell_dir)
    c1 = fread(i, select = c(1,2,5,6))
    
    colnames(c1)[1:2] = c("CHR","MAPINFO")
    c2 = merge(c1, meta)
    rm(c1)
    gc()
    return(c2)
  })
  names(cpg) = cell
  stopCluster(cl)
  return(cpg)
}
