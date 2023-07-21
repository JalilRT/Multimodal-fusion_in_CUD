setwd("/home/jalil/Documents/Psilantro/Fusion_multimodal/")
options(digits = 4)
addTaskCallback(function(...) {set.seed(42);TRUE})
pacman::p_load(tidyverse,readxl,rstatix,vtreat,easystats,ggpubr,magrittr,jmuOutlier,cowplot,grid,gridExtra,raveio,xlsx)

path <- paste0(getwd(),"/Final_analysis")

atlas <- "HOA"
if (atlas=="HOA") {
  regions_labels <- read_table(paste0(path,"/atlas_regions/HarvardOxfordNodeNames.txt"),col_names = F)
} else if(atlas=="DK") {
  regions_labels <- read_table(paste0(path,"/atlas_regions/desikanKillianyNodeNames.txt"),col_names = F)
} else {
  print("error in definition, must be HOA or DK atlas")
}

# Demographic data --------------------------------------------------------
Demographic <-  read_xlsx(paste0(path,"/Demographics/subj_dwi_rsfmri_all.xlsx"), sheet = 1,na = "na") %>% 
  filter(RID != 25 & RID != 43 & RID != 116 & RID != 117) %>% mutate_if(is.character,as.numeric)
Demographic_CUD <- Demographic %>% filter(Group=="2")
Demogr_names <- Demographic[17:35] %>% colnames()

# mcca + jica -------------------------------------------------------------

modal <- "woGSR" #could be woGSR or GSR
reg <- "regout" #could be regressed out (regoutgm) or raw data (rawgm)
struc_names <- str_replace(regions_labels$X1,c("-"),c("_"))

mcca_files <- list.files(paste0(path,"/GT_results/",atlas,"/mcca_jica"), 
                         pattern = paste0(tolower(atlas),modal,"_efmm_",reg),full.names = T) %>% 
              map(read_mat)

mcca_files_IQ <- list.files(paste0(path,"/GT_results/",atlas,"/mcca_jica"), 
                            pattern = paste0(tolower(atlas),modal,"_mccajica_",reg),full.names = T) %>% 
                 map(read_mat)

mccajica_mm_fmri <- names(mcca_files[[1]]) %>% 
  str_detect('fmri') %>%
  keep(mcca_files[[1]], .)
mccajica_mm_smri <- names(mcca_files[[1]]) %>% 
  str_detect('smri') %>%
  keep(mcca_files[[1]], .)

mccajica_IQ <- names(mcca_files_IQ[[1]]) %>% 
  str_detect('Iq') %>%
  keep(mcca_files_IQ[[1]], .)
mccajica_IQ_column <- mccajica_IQ %>% map(~which(.x>0.8))

# CC1 (fmri) --------------------------------------------------------------

mm_fmri_group <- mccajica_mm_fmri %>% names() %>% map(~mccajica_mm_fmri[[.x]] %>% as_tibble() %>% 
                                                        add_column(Group=c(rep(2,63),rep(1,42)))) %>% set_names(mccajica_mm_fmri %>% names())
mm_fmri_byIQ <- 1:length(mm_fmri_group) %>% 
  map(~mm_fmri_group[[.x]] %>% select(c(mccajica_IQ_column[[.x]],"Group"))) %>% 
  set_names(mm_fmri_group %>% names()) 
mm_fmri_byIQ <- mm_fmri_byIQ[lapply(mm_fmri_byIQ,length)>1]

mm_fmri_byIQ_long <- names(mm_fmri_byIQ) %>% map(~ pivot_longer(mm_fmri_byIQ[[.x]], -Group,names_to = "Variables", values_to = "Values")) %>% set_names(names(mm_fmri_byIQ))

mm_fmri_byIQ_mean <- mm_fmri_byIQ_long %>% map(~ .x %>% group_by(Group,Variables) %>% summarise(Mean=mean(Values),Median=median(Values),SD=sd(Values)))

# mm_shapiro <- mm_fmri_byIQ %>% map(~ .x %>% shapiro_test(colnames(.x)[!colnames(.x) %in% c("Group")])) %>% 
#   reduce(rbind) %>% add_column(Regions=names(mm_fmri_byIQ))

mm_fmri.Wtest_Residuals <- mm_fmri_byIQ_long %>% 
  map(~ .x %>% group_by(Variables) %>% 
        wilcox_test(Values~Group) %>% # it is still being non-normal distribution, then we'll use Wilcox
        adjust_pvalue(method = "fdr") %>% add_significance())

mm_fmri.Ttest_Residuals <- mm_fmri_byIQ_long %>% 
  map(~ .x %>% group_by(Variables) %>% 
        t_test(Values~Group) %>% 
        adjust_pvalue(method = "fdr") %>% add_significance())

mm_fmri_byIQ.Permtest_Residuals <- names(mm_fmri_byIQ) %>% 
  map(~ mm_fmri_byIQ[[.x]] %>% select(-"Group") %>% colnames %>% 
        map(function(y) coin::independence_test(formula(paste0(y,'~Group')),
                                                mm_fmri_byIQ[[.x]],teststat = "maximum",
                                                distribution = coin::approximate(nresample = 10000)))) %>% 
  set_names(names(mm_fmri_byIQ))

mm_fmri_byIQ.Permtest <- names(mm_fmri_byIQ.Permtest_Residuals) %>% map(~
                                                                          tibble(z_stat=mm_fmri_byIQ.Permtest_Residuals[[.x]] %>% map(coin::statistic) %>% 
                                                                                   reduce(rbind) %>% c()) %>% 
                                                                          add_column(
                                                                            p_value=mm_fmri_byIQ.Permtest_Residuals[[.x]] %>% map(coin::pvalue) %>% 
                                                                              reduce(rbind) %>% c()) %>% 
                                                                          add_column(
                                                                            p_fdr=p.adjust(mm_fmri_byIQ.Permtest_Residuals[[.x]] %>% map(coin::pvalue) %>% 
                                                                                             reduce(rbind) %>% c(),method = "fdr"))) %>% set_names(names(mm_fmri_byIQ))
mm_fmri_p.adjust <- names(mm_fmri.Wtest_Residuals) %>% 
  map(~ mm_fmri.Wtest_Residuals[[.x]] %>% 
        add_column(modal=names(mm_fmri.Wtest_Residuals[.x]),.before=1)) %>% 
  reduce(rbind) %>% adjust_pvalue(method = "fdr")

dir.create(paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/Group_dif/"),recursive = T)
mm_fmri.Wtest_Residuals %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/","/Group_dif/",.y, "_Wilcox.csv")))
mm_fmri_byIQ.Permtest %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/","/Group_dif/",.y, "_Perm.csv")))

mm_fmri_byIQ_mean %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/",.y, "_mean.csv")))
write.csv(mm_fmri_p.adjust, paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/statistic_adjusted.csv"))

## Correlations
mm_fmri_byIQ_CUD <- mm_fmri_byIQ %>% map(~ .x %>% filter(Group==2) %>% select(-Group)) 
mm_fmri_byIQ_CUD_rho <- mm_fmri_byIQ_CUD %>% names() %>% 
  map(function(met) mm_fmri_byIQ_CUD[[met]] %>% names() %>% 
        map(function(y) Demogr_names %>% 
              map(function(name) cor.test(mm_fmri_byIQ_CUD[[met]][[y]], 
                                          Demographic_CUD[[name]],method = "spearman")) %>% 
              map(extract("estimate")) %>% reduce(rbind) %>% 
              set_rownames(Demogr_names)) %>% reduce(cbind) %>% 
        set_colnames(letters[1:length(mm_fmri_byIQ_CUD[[met]] %>% names())])) %>% 
  set_names(mm_fmri_byIQ_CUD %>% names()) %>% 
  map(~ as_tibble(.x %>% as_tibble()) %>% 
        add_column(Clinicals=rownames(.x),.before = 1))

mm_fmri_byIQ_CUD_rhop <- mm_fmri_byIQ_CUD %>% names() %>% 
  map(function(met) mm_fmri_byIQ_CUD[[met]] %>% names() %>% 
        map(function(y) Demogr_names %>% 
              map(function(name) cor.test(mm_fmri_byIQ_CUD[[met]][[y]], 
                                          Demographic_CUD[[name]],method = "spearman")) %>% 
              map(extract("p.value")) %>% reduce(rbind) %>% 
              set_rownames(Demogr_names)) %>% reduce(cbind) %>% 
        set_colnames(letters[1:length(mm_fmri_byIQ_CUD[[met]] %>% names())])) %>% 
  set_names(mm_fmri_byIQ_CUD %>% names()) %>% 
  map(~ as_tibble(.x %>% as_tibble()) %>% 
        add_column(Clinicals=rownames(.x),.before = 1))

mm_fmri_byIQ_CUD_pval_adj <- mm_fmri_byIQ_CUD_rhop %>% names() %>% map(function(met)
  mm_fmri_byIQ_CUD_rhop[[met]][-1] %>% as.matrix %>% as.vector %>% 
    p.adjust(method='fdr') %>% 
    matrix(ncol=ncol(mm_fmri_byIQ_CUD_rhop[[met]][-1])) %>% as_tibble() %>% 
    add_column(Clinicals=mm_fmri_byIQ_CUD_rhop[[met]][1] %>% unlist(), .before = 1) %>% 
    set_colnames(mm_fmri_byIQ_CUD_rhop[[met]] %>% names())) %>% 
  set_names(mm_fmri_byIQ_CUD %>% names())

# corr permutation

mm_fmri_byIQ_CUD_Permpval <- mm_fmri_byIQ_CUD %>% names() %>% 
  map(function(met) mm_fmri_byIQ_CUD[[met]] %>% names() %>% 
        map(function(y) Demogr_names %>% 
              map(function(name) perm.cor.test(mm_fmri_byIQ_CUD[[met]][[y]], 
                                               Demographic_CUD[[name]],method = "spearman", num.sim = 1000)) %>% 
              map(extract("p.value")) %>% reduce(rbind) %>% 
              set_rownames(Demogr_names)) %>% reduce(cbind) %>% 
        set_colnames(letters[1:length(mm_fmri_byIQ_CUD[[met]] %>% names())])) %>% 
  set_names(mm_fmri_byIQ_CUD %>% names()) %>% 
  map(~ as_tibble(.x %>% as_tibble()) %>% 
        add_column(Clinicals=rownames(.x),.before = 1))

mm_fmri_byIQ_CUD_Permpval_adj <- mm_fmri_byIQ_CUD_Permpval %>% names() %>% map(function(met)
  mm_fmri_byIQ_CUD_Permpval[[met]][-1] %>% as.matrix %>% as.vector %>% 
    p.adjust(method='fdr') %>% 
    matrix(ncol=ncol(mm_fmri_byIQ_CUD_Permpval[[met]][-1])) %>% as_tibble() %>% 
    add_column(Clinicals=mm_fmri_byIQ_CUD_Permpval[[met]][1] %>% unlist(), .before = 1) %>% 
    set_colnames(mm_fmri_byIQ_CUD_Permpval[[met]] %>% names())) %>% 
  set_names(mm_fmri_byIQ_CUD %>% names())

dir.create(paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/Corr/"),recursive = T)
mm_fmri_byIQ_CUD_rho %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/Corr/",.y, "_rho.csv")))
mm_fmri_byIQ_CUD_rhop %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/Corr/",.y, "_pval.csv")))
mm_fmri_byIQ_CUD_pval_adj %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/Corr/",.y, "_padj.csv")))

dir.create(paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/PermCorr"),recursive = T)
mm_fmri_byIQ_CUD_Permpval %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/PermCorr/",.y, "_p.csv")))
mm_fmri_byIQ_CUD_Permpval_adj %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/",modal,"/PermCorr/",.y, "_padj.csv")))

# CC2 (smri) --------------------------------------------------------------

mm_smri_group <- mccajica_mm_smri %>% names() %>% map(~mccajica_mm_smri[[.x]] %>% as_tibble() %>% 
                                                        add_column(Group=c(rep(2,63),rep(1,42)))) %>% set_names(mccajica_mm_smri %>% names())
mm_smri_byIQ <- 1:length(mm_smri_group) %>% 
  map(~mm_smri_group[[.x]] %>% select(c(mccajica_IQ_column[[.x]],"Group"))) %>% 
  set_names(mm_smri_group %>% names()) 
mm_smri_byIQ <- mm_smri_byIQ[lapply(mm_smri_byIQ,length)>1]

mm_smri_byIQ_long <- names(mm_smri_byIQ) %>% map(~ pivot_longer(mm_smri_byIQ[[.x]], -Group,names_to = "Variables", values_to = "Values")) %>% set_names(names(mm_smri_byIQ))

mm_smri_byIQ_mean <- mm_smri_byIQ_long %>% map(~ .x %>% group_by(Group,Variables) %>% summarise(Mean=mean(Values),Median=median(Values),SD=sd(Values)))

mm_smri.Wtest_Residuals <- mm_smri_byIQ_long %>% 
  map(~ .x %>% group_by(Variables) %>% 
        wilcox_test(Values~Group) %>% # it is still being non-normal distribution, then we'll use Wilcox
        adjust_pvalue(method = "fdr") %>% add_significance())

mm_smri_byIQ.Permtest_Residuals <- names(mm_smri_byIQ) %>% 
  map(~ mm_smri_byIQ[[.x]] %>% select(-"Group") %>% colnames %>% 
        map(function(y) coin::independence_test(formula(paste0(y,'~Group')),mm_smri_byIQ[[.x]],
                                                teststat = "maximum",
                                                distribution = coin::approximate(nresample = 10000)))) %>% set_names(names(mm_smri_byIQ))

mm_smri_byIQ.Permtest <- names(mm_smri_byIQ.Permtest_Residuals) %>% 
  map(~tibble(z_stat=mm_smri_byIQ.Permtest_Residuals[[.x]] %>% map(coin::statistic) %>% 
                reduce(rbind) %>% c()) %>% 
        add_column(p_value=mm_smri_byIQ.Permtest_Residuals[[.x]] %>% map(coin::pvalue) %>% 
                     reduce(rbind) %>% c()) %>% 
        add_column(p_fdr=p.adjust(mm_smri_byIQ.Permtest_Residuals[[.x]] %>% map(coin::pvalue) %>% 
                                    reduce(rbind) %>% c(),method = "fdr"))) %>% set_names(names(mm_smri_byIQ))

mm_smri_p.adjust <- names(mm_smri.Wtest_Residuals) %>% 
  map(~ mm_smri.Wtest_Residuals[[.x]] %>% 
        add_column(modal=names(mm_smri.Wtest_Residuals[.x]),.before=1)) %>% 
  reduce(rbind) %>% adjust_pvalue(method = "fdr")

dir.create(paste0(path,"/statistics/",atlas,"/Multimodal/smri/Group_dif"),recursive = T)
mm_smri.Wtest_Residuals %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/smri/Group_dif/",.y, "_Wilcox.csv")))
mm_smri_byIQ.Permtest %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/smri/Group_dif/",.y, "_Perm.csv")))

mm_smri_byIQ_mean %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/smri/",.y, "_mean.csv")))
write.csv(mm_smri_p.adjust, paste0(path,"/statistics/",atlas,"/Multimodal/","smri/statistic_adjusted.csv"))

## Correlations
mm_smri_byIQ_CUD <- mm_smri_byIQ %>% map(~ .x %>% filter(Group==2) %>% select(-Group)) 
mm_smri_byIQ_CUD_rho <- mm_smri_byIQ_CUD %>% names() %>% 
  map(function(met) mm_smri_byIQ_CUD[[met]] %>% names() %>% 
        map(function(y) Demogr_names %>% 
              map(function(name) cor.test(mm_smri_byIQ_CUD[[met]][[y]], 
                                          Demographic_CUD[[name]],method = "spearman")) %>% 
              map(extract("estimate")) %>% reduce(rbind) %>% 
              set_rownames(Demogr_names)) %>% reduce(cbind) %>% 
        set_colnames(letters[1:length(mm_smri_byIQ_CUD[[met]] %>% names())])) %>% 
  set_names(mm_smri_byIQ_CUD %>% names()) %>% 
  map(~ as_tibble(.x %>% as_tibble()) %>% 
        add_column(Clinicals=rownames(.x),.before = 1))

mm_smri_byIQ_CUD_rhop <- mm_smri_byIQ_CUD %>% names() %>% 
  map(function(met) mm_smri_byIQ_CUD[[met]] %>% names() %>% 
        map(function(y) Demogr_names %>% 
              map(function(name) cor.test(mm_smri_byIQ_CUD[[met]][[y]], 
                                          Demographic_CUD[[name]],method = "spearman")) %>% 
              map(extract("p.value")) %>% reduce(rbind) %>% 
              set_rownames(Demogr_names)) %>% reduce(cbind) %>% 
        set_colnames(letters[1:length(mm_smri_byIQ_CUD[[met]] %>% names())])) %>% 
  set_names(mm_smri_byIQ_CUD %>% names()) %>% 
  map(~ as_tibble(.x %>% as_tibble()) %>% 
        add_column(Clinicals=rownames(.x),.before = 1))

mm_smri_byIQ_CUD_pval_adj <- mm_smri_byIQ_CUD_rhop %>% names() %>% map(function(met)
  mm_smri_byIQ_CUD_rhop[[met]][-1] %>% as.matrix %>% as.vector %>% 
    p.adjust(method='fdr') %>% 
    matrix(ncol=ncol(mm_smri_byIQ_CUD_rhop[[met]][-1])) %>% as_tibble() %>% 
    add_column(Clinicals=mm_smri_byIQ_CUD_rhop[[met]][1] %>% unlist(), .before = 1) %>% 
    set_colnames(mm_smri_byIQ_CUD_rhop[[met]] %>% names())) %>% 
  set_names(mm_smri_byIQ_CUD %>% names())

# corr permutation
mm_smri_byIQ_CUD_Permpval <- mm_smri_byIQ_CUD %>% names() %>% 
  map(function(met) mm_smri_byIQ_CUD[[met]] %>% names() %>% 
        map(function(y) Demogr_names %>% 
              map(function(name) perm.cor.test(mm_smri_byIQ_CUD[[met]][[y]], 
                                               Demographic_CUD[[name]],method = "spearman", num.sim = 1000)) %>% 
              map(extract("p.value")) %>% reduce(rbind) %>% 
              set_rownames(Demogr_names)) %>% reduce(cbind) %>% 
        set_colnames(letters[1:length(mm_smri_byIQ_CUD[[met]] %>% names())])) %>% 
  set_names(mm_smri_byIQ_CUD %>% names()) %>% 
  map(~ as_tibble(.x %>% as_tibble()) %>% 
        add_column(Clinicals=rownames(.x),.before = 1))

mm_smri_byIQ_CUD_Permpval_adj <- mm_smri_byIQ_CUD_Permpval %>% names() %>% map(function(met)
  mm_smri_byIQ_CUD_Permpval[[met]][-1] %>% as.matrix %>% as.vector %>% 
    p.adjust(method='fdr') %>% 
    matrix(ncol=ncol(mm_smri_byIQ_CUD_Permpval[[met]][-1])) %>% as_tibble() %>% 
    add_column(Clinicals=mm_smri_byIQ_CUD_Permpval[[met]][1] %>% unlist(), .before = 1) %>% 
    set_colnames(mm_smri_byIQ_CUD_Permpval[[met]] %>% names())) %>% 
  set_names(mm_smri_byIQ_CUD %>% names())

dir.create(paste0(path,"/statistics/",atlas,"/Multimodal/smri/Corr"),recursive = T)
mm_smri_byIQ_CUD_rho %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/smri/Corr/",.y, "_rho.csv")))
mm_smri_byIQ_CUD_rhop %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/smri/Corr/",.y, "_pval.csv")))
mm_smri_byIQ_CUD_pval_adj %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/smri/Corr/",.y, "_padj.csv")))

dir.create(paste0(path,"/statistics/",atlas,"/Multimodal/smri/PermCorr"),recursive = T)
mm_smri_byIQ_CUD_Permpval %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/smri/PermCorr/",.y, "_p.csv")))
mm_smri_byIQ_CUD_Permpval_adj %>% iwalk(~write_csv(.x, paste0(path,"/statistics/",atlas,"/Multimodal/smri/PermCorr/",.y, "_padj.csv")))


# ICs plotting ------------------------------------------------------------

mcca_files_IC <- list.files(paste0(path,"/GT_results/",atlas,"/mcca_jica"), 
                            pattern = paste0(tolower(atlas),modal,"_mccajica_",reg),full.names = T) %>% 
  map(read_mat)

mccajica_IC <- names(mcca_files_IC[[1]]) %>% 
  str_detect('IC') %>%
  keep(mcca_files_IQ[[1]], .) %>% map(~ .x %>% t() %>% as_tibble())

mm_mri_ICbyIQ <- 1:length(mccajica_IC) %>% 
  map(~mccajica_IC[[.x]] %>% select(c(mccajica_IQ_column[[.x]]))) %>% 
  set_names(mccajica_IC %>% names()) 
mm_mri_ICbyIQ <- mm_mri_ICbyIQ[lapply(mm_mri_ICbyIQ,length)>0]

# Slice by modality (first 112 = fmri, last 112 = smri)

mm_fmri_ICbyIQ <- mm_mri_ICbyIQ %>% map(~ .x %>% slice(1:(nrow(.x)/2)) %>% 
                                          scale(center = T,scale = T) %>% 
                                          as_tibble)
mm_smri_ICbyIQ <- mm_mri_ICbyIQ %>% map(~ .x %>% slice(((nrow(.x)/2)+1):nrow(.x)) %>% 
                                          scale(center = T,scale = T) %>% 
                                          as_tibble)

# Correlation between ICs -------------------------------------------------

# Correlation between IC modalities 
cor.test(mm_fmri_byIQ$BC_efmm_fmri$V2,mm_smri_byIQ$BC_efmm_smri$V2,method = "spearman") #Yes
cor.test(mm_fmri_byIQ$PC_efmm_fmri$V1,mm_smri_byIQ$PC_efmm_smri$V1,method = "spearman") #No
cor.test(mm_fmri_byIQ$PC_efmm_fmri$V4,mm_smri_byIQ$PC_efmm_smri$V4,method = "spearman") #Yes

pacman::p_load(ggside)
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(3, 1, 7, 2)]

BC_IC.plot <- tibble(mm_fmri_byIQ$BC_efmm_fmri$V2,mm_smri_byIQ$BC_efmm_smri,
       .name_repair = "minimal") %>% 
  set_colnames(c("fmri","smri","Group")) %>% 
  mutate(Group = case_when(Group == 1 ~ "HC",
                           Group == 2 ~ "CUD"),
         Group = factor(Group,levels=c("HC","CUD"))) %>% 
  ggscatter(x = "smri", y = "fmri", color = "Group") +
  geom_xsideboxplot(aes(y = smri,group = Group, color = Group), 
                    orientation = "y") +
  geom_ysideboxplot(aes(x = fmri, group = Group, color = Group), 
                    orientation = "x") +
  labs(x = "tc-dMRI", y = "fc-fMRI",
       subtitle = "Effective mixing matrices of Betweenness Centrality") +  
  theme_ggside_void()+
  theme(text = element_text(size = 21,family = "times new roman")) +
  scale_color_manual(values = my_pal) +
  geom_smooth(aes(color=Group),method='lm', formula= y~x,se = F)

PC_IC.plot <- tibble(mm_fmri_byIQ$PC_efmm_fmri$V4,mm_smri_byIQ$PC_efmm_smri[2:3],
                     .name_repair = "minimal") %>% 
  set_colnames(c("fmri","smri","Group")) %>% 
  mutate(Group = case_when(Group == 1 ~ "HC",
                           Group == 2 ~ "CUD"),
         Group = factor(Group,levels=c("HC","CUD"))) %>% 
  ggscatter(x = "smri", y = "fmri", color = "Group") +
  geom_xsideboxplot(aes(y = smri,group = Group, color = Group), 
                    orientation = "y") +
  geom_ysideboxplot(aes(x = fmri, group = Group, color = Group), 
                    orientation = "x") +
  labs(x = "tc-dMRI", y = "fc-fMRI",
       subtitle = "Effective mixing matrices of Participation Coefficient") +  
  theme_ggside_void()+
  theme(text = element_text(size = 21,family = "times new roman")) +
  scale_color_manual(values = my_pal) +
  geom_smooth(aes(color=Group),method='lm', formula= y~x,se = F)

IC.plot <- ggarrange(BC_IC.plot,NULL,PC_IC.plot,ncol = 1,
                     labels = c("A","","B"),
          common.legend = T,heights = c(1, 0.05, 1),
          font.label = list(family="times new roman",size=26))

ggsave(units="in", width=10, height=12, dpi=300,
       plot = IC.plot,bg = "white",
       paste0(path,"/images/",atlas,"/Multimodal/","Plot_ICs.png"))

# Creating list of nodes --------------------------------------------------

if (atlas=="HOA") {
HOA_coords <- read_csv(paste0(path,"/atlas_regions/",atlas,"_atlas/","HO_atlas_coords.node"),col_names = F) %>% 
  set_names(c("X","Y","Z","Color","Scale","Region"))
HOA_coords <- HOA_coords[c(rbind(seq(2, nrow(HOA_coords), 2), seq(1, nrow(HOA_coords), 2))),]

BC_fmri_IC <- mm_fmri_ICbyIQ$Bc_IC %>% 
  map(~ HOA_coords %>% mutate(Color=.x)) %>% 
  map(~ .x %>% mutate(Scale=case_when(Color >= 2.3 | Color <= -2.3  ~ 1,Color<= 2.3 & Color >= -2.3 ~ 0))) 
  
PC_fmri_IC <- mm_fmri_ICbyIQ$Pc_IC %>% 
  map(~ HOA_coords %>% mutate(Color=.x)) %>% 
  map(~ .x %>% mutate(Scale=case_when(Color >= 2.3 | Color <= -2.3  ~ 1,Color<= 2.3 & Color >= -2.3 ~ 0))) 

BC_smri_IC <- mm_smri_ICbyIQ$Bc_IC %>% 
  map(~ HOA_coords %>% mutate(Color=.x)) %>% 
  map(~ .x %>% mutate(Scale=case_when(Color >= 2.3 | Color <= -2.3  ~ 1,Color<= 2.3 & Color >= -2.3 ~ 0)))  
  
PC_smri_IC <- mm_smri_ICbyIQ$Pc_IC %>% 
  map(~ HOA_coords %>% mutate(Color=.x)) %>% 
  map(~ .x %>% mutate(Scale=case_when(Color >= 2.3 | Color <= -2.3  ~ 1,Color<= 2.3 & Color >= -2.3 ~ 0))) 
  
dir.create(paste0(path,"/atlas_regions/",atlas,"_atlas/coords_node/mcca_jica/"),recursive = T)
BC_fmri_IC %>% iwalk(~write_delim(.x, paste0(path,"/atlas_regions/",atlas,"_atlas/coords_node/mcca_jica/BC_fmri_coords_",.y,".node"),col_names = F,delim = "\t"))
PC_fmri_IC %>% iwalk(~write_delim(.x, paste0(path,"/atlas_regions/",atlas,"_atlas/coords_node/mcca_jica/PC_fmri_coords_",.y,".node"),col_names = F,delim = "\t"))
BC_smri_IC %>% iwalk(~write_delim(.x, paste0(path,"/atlas_regions/",atlas,"_atlas/coords_node/mcca_jica/BC_smri_coords_",.y,".node"),col_names = F,delim = "\t"))
PC_smri_IC %>% iwalk(~write_delim(.x, paste0(path,"/atlas_regions/",atlas,"_atlas/coords_node/mcca_jica/PC_smri_coords_",.y,".node"),col_names = F,delim = "\t"))

# Figure ------------------------------------------------------------------

BC_fmri_V2 <- ggdraw() + draw_image(paste0(path,"/images/",atlas,"/Multimodal/BC/","BC_fmri_V2.png"))
BC_smri_V2 <- ggdraw() + draw_image(paste0(path,"/images/",atlas,"/Multimodal/BC/","BC_smri_V2.png"))
PC_fmri_V4 <- ggdraw() + draw_image(paste0(path,"/images/",atlas,"/Multimodal/PC/","PC_fmri_V4.png"))
PC_smri_V4 <- ggdraw() + draw_image(paste0(path,"/images/",atlas,"/Multimodal/PC/","PC_smri_V4.png"))

bar <- ggdraw() + draw_image(paste0(path,"/images/",atlas,"/Multimodal/","bar.tif"))

multimodal <- ggarrange(BC_fmri_V2,BC_smri_V2,
                        PC_fmri_V4,PC_smri_V4,
                        ncol = 2,nrow = 2,
                        labels = c("Betweenness Centrality - fMRI",
                                   "Betweenness Centrality - sMRI",
                                   "Participation Coefficient - fMRI",
                                   "Participation Coefficient - sMRI"),
                        hjust = -0.2,font.label = list(size=15,face='bold'))

#source("~/Documents/Psilantro/Fusion_multimodal/Final_analysis/Unimodal_analysis.R")

PC <- ggdraw() +
  draw_plot(plots_met$Pc + labs(y="left Caudate") + 
              theme(text = element_text(size = 25))) +
  draw_image(paste0(path,"/atlas_regions/",atlas,"_atlas/","PC_lh_s.png"),  x = 0.38, y = 0.37, scale = .18)

Figure3 <- ggarrange(PC,ggplot() + theme_void(),multimodal,
                     ggplot() + theme_void(),ggplot() + theme_void(),bar,
                     ncol = 3,nrow=2,labels = c("A","B"),
                     widths = c(0.8,0.06,1.2),heights = c(1,0.042),
                     font.label = list(size=22,face='bold'))

ggsave(units="in", width=17, height=9, dpi=300,
       plot = Figure3,bg = "white",
       paste0(path,"/images/",atlas,"/Figure3-R.png"))

library(svglite)
svglite(filename = paste0(path,"/images/",atlas,"/Figure3.svg"),
        width = 17,height = 9,bg = "white")
Figure3
dev.off()

} else if(atlas=="DK") {
  
} else {
  print("error in definition, must be HOA or DK atlas")
}

