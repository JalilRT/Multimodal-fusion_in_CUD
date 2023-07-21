setwd("/home/jalil/Documents/Psilantro/Fusion_multimodal/")
options(digits = 4)
addTaskCallback(function(...) {set.seed(42);TRUE})
pacman::p_load(tidyverse,readxl,rstatix,vtreat,easystats,ggpubr,magrittr,jmuOutlier,cowplot,grid,gridExtra,raveio,xlsx)

path <- paste0(getwd(),"/Final_analysis")

atlas <- "DK"

if (atlas=="HOA") {
  regions_labels <- read_table(paste0(path,"/atlas_regions/HarvardOxfordNodeNames.txt"),col_names = F)
} else if(atlas=="DK") {
  regions_labels <- read_table(paste0(path,"/atlas_regions/desikanKillianyNodeNames.txt"),col_names = F)
} else {
  print("error in definition, must be HOA atlas or DK atlas")
}

struc_names <- str_replace(regions_labels$X1,c("-"),c("_"))

# Demographic data --------------------------------------------------------
Demographic <-  read_xlsx(paste0(path,"/Demographics/subj_dwi_rsfmri_all.xlsx"), sheet = 1,na = "na") %>% 
  filter(RID != 25 & RID != 43 & RID != 116 & RID != 117) %>% mutate_if(is.character,as.numeric)
Demographic_CUD <- Demographic %>% filter(Group=="2")
Demogr_names <- Demographic[17:28] %>% colnames()

# Demographic stats
Demographic %>% group_by(Group) %>% select(Group,age,educ,coc.age.onset,years.begin,week.dose) %>% 
  summarise(median_age = median(age), n = n(), range = range(age),
            median_educ = median(as.numeric(educ) %>% na.omit()),
            coc_onset =  median(as.numeric(coc.age.onset) %>% na.omit()),
            n_onset = length(as.numeric(coc.age.onset) %>% na.omit()),
            range_onset = range(as.numeric(coc.age.onset) %>% na.omit()),
            coc_years =  median(as.numeric(years.begin) %>% na.omit()),
            n_years = length(as.numeric(years.begin) %>% na.omit()),
            range_years = range(as.numeric(years.begin) %>% na.omit()),
            mean_consumption = length(as.numeric(week.dose) %>% na.omit()))
Demographic %>% group_by(Group) %>% select(handeness, Group) %>% table()
table(Demographic %>% filter(Group==2) %>% select(week.dose))
Demographic %>% shapiro_test(age)
t.test(age~Group,Demographic)
chisq.test(as.numeric(Demographic$educ) %>% na.omit(), Demographic$Group[-c(82,88)])

# GT-Local metrics --------------------------------------------------------

modal <- c("fmri_wogsr","smri")
Net_val <- c("Bc","Dc","Ncc","Nle","Pc")

for (mri in modal){
  for (Net in Net_val) {
    pathAtlas = paste0(path,"/GT_results/",atlas,"/",mri)
    Net_files <- list.files(pathAtlas, pattern = Net, full.names = T)
    Net_lmResiduals <- Net_files %>% map(~ read_xlsx(.x,col_names = F)) %>% reduce(rbind) %>% set_colnames(struc_names) %>% add_column(Group=c(rep(2,63),rep(1,42))) 
    Net_Dtest <- tibble(Net_lmResiduals,Demographic %>% select(age,sex,educ,tobc.totyears),MeanNet=Net_lmResiduals[-c(1,ncol(Net_lmResiduals))] %>% rowMeans())
    
    Net_shapiro <- struc_names %>% map(~ Net_lmResiduals[[.x]] %>% shapiro_test()) %>% reduce(rbind) %>% add_column(Regions=struc_names)
    
    Net_valResiduals_long <- Net_lmResiduals %>% pivot_longer(-Group,names_to = "Variables", values_to = "Values")
    
    Net_stat.Wtest_Residuals <- Net_valResiduals_long %>% 
      group_by(Variables) %>% 
      wilcox_test(Values~Group) %>% # it is still being non-normal distribution, then we'll use Wilcox
      adjust_pvalue(method = "fdr") %>% add_significance()
    
    Net_stat.Ttest_Residuals <- Net_valResiduals_long %>% 
      group_by(Variables) %>% 
      t_test(Values~Group) %>% # it is still being non-normal distribution, then we'll use Wilcox
      adjust_pvalue(method = "fdr") %>% add_significance()
    
    Net_stat.Permtest_Residuals <- struc_names %>% 
      map(function(x) coin::independence_test(formula(paste0(x,'~Group')),Net_lmResiduals,
                                              teststat = "maximum",
                                              distribution = coin::approximate(nresample = 10000))) %>% set_names(struc_names)
    
    Net_stat.Permtest <- tibble(struc_names) %>% add_column(
      z_stat=Net_stat.Permtest_Residuals %>% map(coin::statistic) %>% 
        reduce(rbind) %>% c()) %>% 
      add_column(
        p_value=Net_stat.Permtest_Residuals %>% map(coin::pvalue) %>% 
          reduce(rbind) %>% c()) %>% 
      add_column(
        p_fdr=p.adjust(Net_stat.Permtest_Residuals %>% map(coin::pvalue) %>% 
                         reduce(rbind) %>% c(),method = "fdr"))
    
    dir.create(paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Group_dif/",Net,"/"),recursive = T)
    write_csv(Net_stat.Ttest_Residuals,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Group_dif/",Net,"/",Net, "_T-test.csv"))
    write_csv(Net_stat.Wtest_Residuals,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Group_dif/",Net,"/",Net, "_Wilcox.csv"))
    write_csv(Net_stat.Permtest,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Group_dif/",Net,"/",Net, "_Perm.csv"))
    write_csv(Net_shapiro,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Group_dif/",Net,"/",Net, "_Shapiro.csv"))
    
    # Correlation by spearman with covs 
    
    Net_valResiduals_CUD <- Net_lmResiduals %>% filter(Group=="2")
    
    Net_corr_r <- struc_names %>% map(function(roi) Demogr_names %>% map(function(name) 
      cor.test(Demographic_CUD[name] %>% unlist(),Net_valResiduals_CUD[roi] %>% 
                 unlist(),method = "spearman")) %>% 
        map(extract("estimate")) %>% reduce(rbind) %>% set_rownames(Demogr_names)) %>% 
      reduce(cbind) %>% set_colnames(struc_names)
    Net_corr_r_tibble <- as_tibble(rownames(Net_corr_r)) %>% add_column(Net_corr_r %>% as_tibble())
     
    Net_corr_p <- struc_names %>% map(function(roi) Demogr_names %>% map(function(name) 
      cor.test(Demographic_CUD[name] %>% unlist(),Net_valResiduals_CUD[roi] %>% 
                 unlist(),method = "spearman")) %>% 
        map(extract("p.value")) %>% reduce(rbind) %>% set_rownames(Demogr_names)) %>% 
      reduce(cbind) %>% set_colnames(struc_names) 
    Net_corr_psig <- apply(Net_corr_p,1,function(x) x[which(x <= 0.05)])
    
    Net_corr_pAdj <- Net_corr_p %>% as.matrix %>% as.vector %>% 
      p.adjust(method='fdr') %>% matrix(ncol=ncol(Net_corr_p)) %>% 
      set_colnames(struc_names) %>% set_rownames(Demogr_names)
    Net_corr_psigAdj <- apply(Net_corr_pAdj,1,function(x) x[which(x <= 0.05)])
    Net_corr_pAdj_tibble <- as_tibble(rownames(Net_corr_pAdj)) %>% add_column(Net_corr_pAdj %>% as_tibble())
    
        
    # Correlation by permutation with covs 
    Netperm_corr_p <- struc_names %>% map(function(roi) Demogr_names %>% map(function(name) 
      perm.cor.test(Demographic_CUD[name] %>% unlist(),Net_valResiduals_CUD[roi] %>% 
                      unlist(),method = "spearman",num.sim = 10)) %>% 
        map(extract("p.value")) %>% reduce(rbind) %>% set_rownames(Demogr_names)) %>% 
      reduce(cbind) %>% set_colnames(struc_names)
    Netperm_corr_psig <- apply(Netperm_corr_p,1,function(x) x[which(x <= 0.05)])
    Netperm_corr_p_tibble <- as_tibble(rownames(Netperm_corr_p)) %>% add_column(Netperm_corr_p %>% as_tibble())
    
    Netperm_corr_pAdj <- Netperm_corr_p %>% as.matrix %>% as.vector %>% 
      p.adjust(method='fdr') %>% matrix(ncol=ncol(Netperm_corr_p)) %>% 
      set_colnames(struc_names) %>% set_rownames(Demogr_names)
    Netperm_corr_psigAdj <- apply(Netperm_corr_pAdj,1,function(x) x[which(x <= 0.05)])
    Netperm_corr_pAdj_tibble <- as_tibble(rownames(Netperm_corr_pAdj)) %>% add_column(Netperm_corr_pAdj %>% as_tibble())
    
    # writing

    dir.create(paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Corr/",Net,"/"),recursive = T)
    
    write_csv(Net_corr_r_tibble,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Corr/",Net,"/", Net,"_rho.csv"))
    write_lines(Net_corr_psigAdj %>% as.character(),paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Corr/",Net,"/", Net,"_Corr_pvalueAdjSig.txt"))
    write_csv(Net_corr_pAdj_tibble,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Corr/",Net,"/", Net,"_Corr_pvalueAdj.csv"))
    
    dir.create(paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/PermCorr/",Net,"/"),recursive = T)
    write_csv(Netperm_corr_p_tibble,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/PermCorr/", Net,"/",Net,"_PermCorr_pvalue.csv"))
    write_lines(Netperm_corr_psigAdj %>% as.character(),paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/PermCorr/", Net,"/",Net,"_PermCorr_pvalueAdjSig.txt"))
    write_csv(Netperm_corr_pAdj_tibble,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/PermCorr/", Net,"/",Net,"_PermCorr_pvalueAdj.csv"))
    
  }}

# GT-Global metrics -------------------------------------------------------

for (mri in modal) {
  Global_NetMetric <- list.files(paste0(path,"/GT_results/",atlas,"/",mri), 
                                 pattern = "global", full.names = T) %>% 
    map(read_xlsx) %>% reduce(rbind) %>% 
    add_column(Group=c(rep(2,63),rep(1,42)))
  
  Global_shapiro <- colnames(Global_NetMetric)[-7] %>% 
    map(~ Global_NetMetric[[.x]] %>% shapiro_test()) %>% reduce(rbind) %>% 
    add_column(Metric=colnames(Global_NetMetric)[-7])
  
  Global_NetMetric_long <- Global_NetMetric %>% 
    pivot_longer(-Group,names_to = "Variables", values_to = "Values")
  
  Global_Net.Wtest <- Global_NetMetric_long %>% 
    group_by(Variables) %>% 
    wilcox_test(Values~Group) %>% # it is still being non-normal distribution, then we'll use Wilcox
    adjust_pvalue(method = "fdr") %>% add_significance()
  
  Global_Net.Ttest <- Global_NetMetric_long %>% 
    group_by(Variables) %>% 
    t_test(Values~Group) %>% # it is still being non-normal distribution, then we'll use Wilcox
    adjust_pvalue(method = "fdr") %>% add_significance()
  
  Global_stat.Permtest_Residuals <- colnames(Global_NetMetric)[-7] %>% 
    map(function(x) coin::independence_test(formula(paste0(x,'~Group')),Global_NetMetric,
                                            teststat = "maximum",
                                            distribution = coin::approximate(nresample = 10000)))
  
  Global_stat.Permtest <- tibble(Metric=colnames(Global_NetMetric)[-7]) %>% 
    add_column(
    z_stat=Global_stat.Permtest_Residuals %>% map(coin::statistic) %>% 
      reduce(rbind) %>% c()) %>% 
    add_column(
      p_value=Global_stat.Permtest_Residuals %>% map(coin::pvalue) %>% 
        reduce(rbind) %>% c()) %>% 
    add_column(
      p_fdr=p.adjust(Global_stat.Permtest_Residuals %>% map(coin::pvalue) %>% 
                       reduce(rbind) %>% c(),method = "fdr"))
  
  write_csv(Global_stat.Permtest,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Group_dif/","Perm_Global_stats.csv"))
  write_csv(Global_stat.Permtest,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Group_dif/","Wilcox_Global_stats.csv"))
  write_csv(Global_stat.Permtest,paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Group_dif/","Ttest_Global_stats.csv"))
}

# Visualization -----------------------------------------------------------

adj="p" # corrected ("p.adj") or uncorrected ("p")
mri="fmri_wogsr" # modality: "fmri_wogsr" or "smri"
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(3, 1, 7, 2)]

# Selecting significants with correction regions by metric
Metrics.list <- list.files(paste0(path,"/statistics/",atlas,"/Unimodal/",mri,"/Group_dif"), pattern = "*Wilcox.csv",recursive = T, full.names = T)
Metrics.csv <- Metrics.list %>% map(read_csv) 
Metrics.names <- 1:length(Metrics.list) %>% 
  map(function(x) str_remove(str_split_fixed(Metrics.list[x],"/",13),"_Wilcox.csv")) %>% 
  map(last) %>% map(basename) %>% unlist()
Metric.sig.regions <- 1:length(Metrics.list) %>% 
  map(function(x) Metrics.csv[[x]]$Variables[Metrics.csv[[x]][[adj]] < 0.05]) %>% set_names(Metrics.names)
  
NetMetric.values <- Metrics.names %>% map(function(x) 
  list.files(paste0(path,"/GT_results/",atlas,"/",mri),pattern = x, full.names = T)) %>% 
    map(~ read_xlsx(.x,col_names=F) %>% set_colnames(struc_names) %>% 
    add_column(Group=c(rep("CUD",63),rep("HC",42)))) %>% set_names(Metrics.names)
  
NetMetric.filter <- Metrics.names %>% map(function(x) NetMetric.values[[x]] %>% 
                                              select(Metric.sig.regions[[x]],Group)) %>% 
  set_names(Metrics.names)

NetMetric.filter <- NetMetric.filter[lapply(NetMetric.filter,length)>1]

NetMetric_longer <- NetMetric.filter %>% 
  map(~ .x %>% 
        pivot_longer(-Group,names_to = "Regions", values_to = "Values") %>% 
        group_by(Group,Regions) %>% 
        summarise(Values=mean(Values))) 

#ploting results
plot_fun <- function(met){
  y_metric <- grep('R|L',names(NetMetric.filter[[met]]),value = T)
  
  plot_metric <- NetMetric.filter[[met]] %>% 
    mutate(Group=factor(Group,levels=c("HC","CUD"))) %>% 
    ggplot(aes(x = Group, y = get(y_metric), 
               color = Group, fill = Group)) +
    scale_color_manual(values = my_pal, guide = "none") +
    scale_fill_manual(values = my_pal, guide = "none") +
    ggdist::stat_halfeye(
      adjust = .33, width = .7, fill = "grey85",
      interval_colour = NA, 
      position = position_nudge(x = .01),
      aes(thickness = stat(f*n))) +
    geom_boxplot(width = .1) +
    gghalves::geom_half_point(side = "l", 
                              range_scale = .3, alpha = .4, size = 2) +
    stat_summary(geom = "point",
                 fun = median, shape = 23, fill="white",
                 color = "black",size = 2) + 
    theme_classic(base_size = 15) + theme(axis.title.x=element_blank())+
    labs(y=y_metric)
  
}

if (adj=="p.adj") {
  
  metric.name <- c()
  for (metric in Metrics.names) {
    if (identical(grep('R|L',names(NetMetric.filter[[metric]]),value = T),character(0))==FALSE) {
      metric.name[metric] <- metric}}
  metric.name <- metric.name %>% as.character()
  
  plots_met <- metric.name %>% map(~ plot_fun(.x)) %>% set_names(metric.name)
  # Save images
  dir.create(paste0(path,"/images/",atlas,"/Unimodal/",mri,"/corrected/"),recursive = T)
  names(plots_met) %>% 
    map(function(Met) ggsave(units="in", width=8, height=9, dpi=300,
                             plot = plots_met[[Met]],
                             paste0(path,"/images/",atlas,"/Unimodal/",mri,"/corrected/","Plot_",Met,".png")))
  
} else if(adj=="p") {

  Plots_Metrics <- NetMetric_longer %>% map(function(Met)
    ggdotchart(Met, x = "Regions", y = "Values",
               color = "Group", size = 3,
               add = "segment",rotate = TRUE,
               add.params = list(color = "lightgray", size = 1.5),
               position = position_dodge(0),
               palette = "npg",
               ggtheme = theme_pubr(legend = "bottom")))
  # Save images
  dir.create(paste0(path,"/images/",atlas,"/Unimodal/",mri,"/uncorrected"),recursive = T)
  names(Plots_Metrics) %>% 
    map(function(Met) ggsave(units="in", width=10, height=5, dpi=300,
                             plot = Plots_Metrics[[Met]],
                             paste0(path,"/images/",atlas,"/Unimodal/",mri,"/uncorrected/","Plot_",Met,".png")))
  
} else {
  print("error in adj definition, must be uncorrected (p) or corrected (p.adj) ")
}

# Adding brain to Vis -----------------------------------------------------

if (atlas=="HOA") {
HOA_coords <- read_csv(paste0(path,"/atlas_regions/",atlas,"_atlas/","HOA_coords.node"),col_names = F) %>% 
  set_names(c("X","Y","Z","Scale","Color","Region"))

regions_sig <- plots_met %>% map(~ .x[["data"]] %>% select(-"Group") %>% colnames)
which_line <- regions_sig %>% map(~ struc_names %in% .x %>% which)

# PC
HOA_coords$Scale[which_line$Pc-1] <- 1 #is -1 because R-L are inverted in atlas
NetMetric_longer$Pc %>% filter(Regions==regions_sig[["Pc"]] & Group=="CUD") %>% as.vector()

write_csv(HOA_coords,paste0(path,"/atlas_regions/",atlas,"_atlas/","HO_atlas_PC.node"))

# Plotting
PC <- ggdraw() +
  draw_plot(plots_met$Pc) +
  draw_image(paste0(path,"/atlas_regions/",atlas,"_atlas/","PC_lh_s.png"),  x = 0.38, y = 0.37, scale = .18)

ggsave(units="in", width=10, height=9, dpi=300,
       plot = PC,
       paste0(path,"/images/Unimodal/",atlas,"/Unimodal/",mri,"/corrected/","Plot_PC.png"))
} else if(atlas=="DK") {
  print("No significant results after correction")
} else {
  print("error in definition, must be HOA atlas or DK atlas")
}
