setwd("/home/jalil/Documents/Fusion_multimodal/")
options(digits = 5)
addTaskCallback(function(...) {set.seed(42);TRUE})
pacman::p_load(tidyverse,Hmisc,rstatix,vtreat,easystats,ggpubr,magrittr,jmuOutlier,cowplot,grid,gridExtra,raveio,xlsx)

# Obtaining correlation matrices
path = "xcpengine/woGSR/"
atlas = "desikanKilliany"
atlas = "HarvardOxford"
regions_labels = read_table(paste0(getwd(),"/",path,atlas,"/../../",atlas,"NodeNames.txt"),col_names = F)

TS_files <- list.files(paste0(getwd(),"/",path,atlas), pattern = "1D", full.names = T)
TS_names <- TS_files %>% map(~ .x %>% basename %>% str_remove("_ts.1D"))

#only for DK
TS_corr <- TS_files %>% map(~ .x %>% read.csv(sep = " ",header = F) %>% 
                              select(-c(1,2,3,4,10,11,12,15,18,19,20,21,22,23,
                                        33,34,35,36,37,38,39,40,41,42,43,78)) %>% 
                              select(c(19:86,1:18)) %>%
                                as.matrix() %>% rcorr())
# For HOA
TS_corr <- TS_files %>% map(~ .x %>% read.csv(sep = " ",header = F) %>% 
           as.matrix() %>% rcorr() )

# Obtain the r and p values
TS_corr_r <- TS_corr %>% map(~ .x[["r"]] %>% as.data.frame() %>% 
             round(digits = 5) %>% set_colnames(regions_labels$X1)) %>% set_names(TS_names) 
TS_corr_pvals <- TS_corr %>% map(~ .x[["P"]] %>% as.data.frame() %>% 
                 round(digits = 5) %>% set_colnames(regions_labels$X1)) %>% 
                 map(~replace(.x, is.na(.x), 0)) %>% set_names(TS_names)

TS_corr_r %>% iwalk(~write_csv(.x, paste0(path,"/",atlas,"_TScorr/",.y, "_r.csv"),col_names = T))
TS_corr_pvals %>% iwalk(~write_csv(.x, paste0(path,"/",atlas,"_TScorr/",.y, "_pval.csv"),col_names = T))
