setwd('...') # Set your working dir

path.risk.components <- 'data/all_participants/VRS_pca.csv'
path.risk.binary <- 'data/all_participants/risk_binary.csv'
path.risk.continuous <- 'data/all_participants/risk_continuous.csv'
path.risk.aggregate <- 'data/all_participants/aggregate VRS.csv'
path.risk.fsrs <- '/Users/lairaichin/Desktop/PRS_VRS_UKB/data/all_participants/risk_FSRS.csv'
path.risk.all <- 'data/all_participants/Vascular factors.csv'
path.brain.idp <-  'data/all_participants/brain_idp.csv'
path.cognition.factors <- 'data/all_participants/cognition_factors.csv'
path.cognition.raw <- 'data/all_participants/cognition_raw.csv'
path.demo <- 'data/all_participants/demo.csv'
path.global.gretna <- 'data/gretna_bn/global.csv'
path.global.gretna.wide <- 'data/gretna_bn/global_wide.csv'
path.nodal.gretna <- 'data/gretna_bn/nodal.csv'
path.modular.gretna <- 'data/gretna_bn/modular_yeo_gretna.csv'
path.description.bna <- 'bna_246_description_BNV.csv'


library(tidyverse)
library(lubridate)
library(psych)
library(reshape)
library(magrittr)
library(corrplot)

# Data reading and writing ---------
create_when_absent <- function(path) {
  if(!dir.exists(path)) dir.create(path)
}


read_connectome_lairai <- function(f) {
  if (endsWith(f, '.txt') | endsWith(f, '.txt*') | endsWith(f, '.edge') | endsWith(f, '.edge*')) {
    return(as.matrix(read.table(f)))
  } 
  if (endsWith(f, '.csv') | endsWith(f, '.csv*')) {
    return(as.matrix(read_table(f)))
  } 
  stop(str_c('Invalid connectome for ', f))
  # call. = F 会不显示报错的调用
}


get_file_name <- function(x) {
  strsplit(x, '\\.')[[1]][1]
}


read_connectome_path <- function(connectome.path, n = 246) {
  files <- list.files(connectome.path, pattern = '*.txt')
  data_3d <-
    array(100, dim = c(length(files), n, n)) # A 3D array consisting all the data, indexing by file, row and then col
  
  # Initialising the array
  print("Initialising the array...")
  progress <- txtProgressBar(min = 0, max = length(files), initial = 0, char = "*", width = NA, style = 3)
  for (f in 1:length(files)) {
    mat_csv <- read_connectome_lairai(str_c(connectome.path, '/', files[f]))
    #print(files[f])
    data_3d[f,,] = mat_csv
    setTxtProgressBar(progress, f)
  }
  close(progress)
  return(data_3d)
}


read_connectome_path_NBR <- function(connectome.path, n = 246) {
  #Returns an array with N*N*n(sample)
  files <- list.files(connectome.path, pattern = '*.txt')
  data_3d <-
    array(100, dim = c(n, n, length(files))) # A 3D array consisting all the data, indexing by file, row and then col
  
  # Initialising the array
  print("Initialising the array...")
  progress <- txtProgressBar(min = 0, max = length(files), initial = 0, char = "*", width = NA, style = 3)
  for (f in 1:length(files)) {
    mat_csv <- read_connectome_lairai(str_c(connectome.path, '/', files[f]))
    #print(files[f])
    data_3d[,,f] = mat_csv
    setTxtProgressBar(progress, f)
  }
  close(progress)
  return(data_3d)
}


sum_nodes_bnv <- function(df, atlas.description, save.size=F, save.color=F, save.label=T, size.col, color.col, arrange.col='region', path.node.bnv, path.node.csv) {
  if (nrow(df) != nrow(atlas.description)) stop("Inconsistent row number!")
  # A df must contains a column called 'region' with the same size as atlas.description and arranged ascendingly by region number
  result <- atlas.description  %>% 
    dplyr::select(
      region, 
      acronym,
      `full name`,
      X, Y, Z
    ) %>% 
    mutate(color = 1, size = 1, label = acronym)
  
  if (save.size == T) result$size = df[[size.col]]
  if (save.color == T) result$color = df[[color.col]]
  if (save.label == F) result$label = '-'
  
  
  write.table(
    result %>% dplyr::select(X, Y, Z, color, size, label),
    path.node.bnv,
    row.names = F, col.names = F, quote = F)
  
  
  full_join(result, df, by = 'region') %>% 
    write_csv(path.node.csv)
}

write_connetome_lairai <- function(network, path.save, .bnv = T) {
  if(.bnv) path.save <- str_replace(path.save, '\\.txt', '\\.edge')
  write.table(network, path.save, row.names = F, col.names = F, quote = F)
}

write_connetome_lairai <- function(network, path.save, .bnv = T) {
  if(.bnv) path.save <- str_replace(path.save, '\\.txt', '\\.edge')
  write.table(network, path.save, row.names = F, col.names = F, quote = F)
}

walk_connectome <- function(path.connectome, .f, ...) {
  walk(list.files(path.connectome), function(connectome) {
    .f(connectome, ...)
  })
}
# Vector tools --------------------------------------------------
scale.lairai <- function(x) {
  return( (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
}

# Tools of two vectors ----------------------------------------------------
fit_Y_from_X <- function(Y, X) {
  mod.temp <- lm(Y~X)
  p.val <-  broom::tidy(mod.temp)$p.value[2]
  if (p.val > 0.05) {
    #print('Insignificant effect, hence did not adjust')
    return(Y)
  }
  estimate <- broom::tidy(mod.temp)$estimate[2]
  mu <- mean(X, na.rm = T)
  return(Y + estimate*(mu-X))
}


