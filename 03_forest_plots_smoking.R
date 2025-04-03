#------------------------------------------------------------------------
#This script produces key figures of results 
#------------------------------------------------------------------------------

source("00_init.R")
library(ggpubr)
library(RColorBrewer)

#select result 
dat <- readRDS("output/smoke_simulation_final.RDS")
dat <- readRDS("output/smoke_simulation_final_burnin.RDS")


dat_mean="pstar_mean"
dat_sd="pstar_sd"

cbp2 <- c( "#000000",brewer.pal(5, "Set2"))
shape_list <- c(16, 16, 17, 17, 17)

dodge <- position_dodge(width=-0.75)  
beta_names <- unique(dat$beta_setting)

background <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1),  
  strip.text = element_text(size = 12, face = "bold", color = "gray40"),  
  strip.background = element_rect(fill = "gray85", color = "grey85"), 
  panel.spacing = unit(1, "lines"), 
  legend.position = "right",
  legend.title = element_text(face = "bold"),  
  legend.key = element_rect(fill = "white", color = "black"), 
  legend.box.spacing = unit(0.5, "lines"), 
  legend.direction = "vertical", 
  panel.grid.major.x = element_line(color = "gray80", linewidth = 0.5, linetype="dashed"), 
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),     
  panel.border =  element_rect(color = "gray85", linewidth = 3),
)

for (scenario in beta_names){
  
  #Plot for subset of results---------------------------------------------------
  plot_dat <- dat %>% filter(method %in% c( c("Permuted Block", "BRAR", "Neyman Allocation", "Gittins Index", "Gittins Index truncated",
                                              "Randomised Gittins Index")))
  
  
  p <- ggplot(data = plot_dat %>% filter(beta_setting == scenario)) + 
    facet_wrap(~method, ncol = 3) +
    geom_linerange(aes(ymin=get(dat_mean)-1.96*get(dat_sd)/sqrt(N), 
                       ymax=get(dat_mean)+1.96*get(dat_sd)/sqrt(N), 
                       x=alpha_setting, color=missing_approach), position=dodge) + 
    geom_rect(data = plot_dat %>% filter(beta_setting == scenario), 
              aes(xmin = as.numeric(factor(alpha_setting)) - 0.5, 
                  xmax = as.numeric(factor(alpha_setting)) + 0.5, 
                  ymin = -Inf, ymax = Inf), 
              color = "grey95", alpha=0, linetype="dashed") +  
    geom_point(aes(y = get(dat_mean), x = alpha_setting, col = missing_approach, 
                   shape = missing_approach), 
               position = position_dodge(width = -0.7), size = 2.5) + 
    coord_flip() +  
    theme_light() +  
    xlab("Missing Data Mechanism") + 
    ylab(expression(p^"*"~(1622))) + 
    guides(fill = "none") +  # Remove alpha guide
    scale_color_manual(name = "Missing Data Approach", values = cbp2) + 
    scale_shape_manual(name = "Missing Data Approach", values = shape_list) +  
    scale_x_discrete(limits=rev)+background
    
  
  p
  
  
  #plot for all methods---------------------------------------------------------
  
  dat$method <- factor(dat$method, 
                          levels=c("Fixed Randomisation", "Permuted Block", "BRAR", "Neyman Allocation", 
                                   "Current Belief", "Current Belief truncated", "Gittins Index",  "Gittins Index truncated",           
                                   "Randomised Belief Index", "Randomised Gittins Index"))
  
  # Updated plot
  p_all <- ggplot(data = dat %>% filter(beta_setting == scenario)) + 
    facet_wrap(~method, ncol = 4) +  
    geom_point(aes(y = get(dat_mean), x = alpha_setting, col = missing_approach, 
                   shape = missing_approach), 
               position = position_dodge(width = -0.7), size = 3) + 
    geom_linerange(aes(ymin=get(dat_mean)-1.96*get(dat_sd)/sqrt(N), 
                       ymax=get(dat_mean)+1.96*get(dat_sd)/sqrt(N), 
                       x=alpha_setting, color=missing_approach), position=dodge) + 
    coord_flip() + 
    theme_light() +  
    xlab("Missing Data Mechanism") +  
    ylab(expression(p^"*"~(1622))) + 
    guides(fill = "none") +  # Remove alpha guide
    scale_color_manual(name = "Missing Data Approach", values = cbp2) +  
    scale_shape_manual(name = "Missing Data Approach", values = shape_list) +  
    scale_x_discrete(limits=rev)+
    geom_rect(data = plot_dat %>% filter(beta_setting == scenario), 
              aes(xmin = as.numeric(factor(alpha_setting)) - 0.5, 
                  xmax = as.numeric(factor(alpha_setting)) + 0.5, 
                  ymin = -Inf, ymax = Inf), 
              color = "grey95", alpha=0, linetype="dashed") +  background
  # Print the plot
  print(p)
  
  #retain correct path 
  png(paste0("output/pstar_forest_smoke/pstar_", scenario, "_smoke.png"), height=21, width=25, units="cm", res=1200)
  #png(paste0("output/pstar_forest_smoke/pstar_", scenario, "_smoke_burnin.png"), height=21, width=25, units="cm", res=1200)
  
  print(p)
  dev.off()
  
  #retain correct path 
  png(paste0("output/pstar_forest_smoke/pstar_", scenario, "_all_smoke.png"), height=24, width=31, units="cm", res=1200)
  #png(paste0("output/pstar_forest_smoke/pstar_", scenario, "_all_smoke_burnin.png"), height=24, width=31, units="cm", res=1200)
  
  print(p_all)
  dev.off()
  
}


#plots for bias----------------------------------------------------------------------------------------------------

mle_dat <- dat %>% select(method, beta_setting, alpha_setting, missing_approach, starts_with("p0"), 
                          starts_with("p1") )%>%
  filter(method %in% c("Permuted Block", "BRAR"))
mle_dat$method <- factor(mle_dat$method)
levels(mle_dat$method) <- c("Permuted Block - MLE", "BRAR - MLE")
IPW_dat <- dat %>% select(method, beta_setting, alpha_setting, missing_approach, starts_with("IPW") ) %>%
  filter(method =="BRAR")
IPW_dat$method <- "BRAR - IPW"


order <- names(mle_dat)
order
names(IPW_dat) = c("method", "beta_setting", "alpha_setting", "missing_approach", 
                   "p0_cc_mean", "p0_cc_sd", "p0_cc_na",
                   "p1_cc_mean", "p1_cc_sd", "p1_cc_na",
                   "p0_single_imp_mean", "p0_single_imp_sd", "p0_single_imp_na",
                   "p1_single_imp_mean", "p1_single_imp_sd", "p1_single_imp_na")  #make sure col names match with mle_dat
IPW_dat <- IPW_dat[,order]
subset_dat <- rbind(mle_dat, IPW_dat)

dodge <- position_dodge(width=-0.9)  


for (scenario in beta_names){
  
  #use this for the permuted block/BRAR cases 
  plot_dat <- subset_dat %>% filter(beta_setting==scenario) 
  cbp2 <- c( "#000000",brewer.pal(4, "Set2"))[-2] #remove color for Complete cases nt
  
  
  #use this for all other plots
  #plot_dat <- dat %>% filter(beta_setting==scenario) 
  #cbp2 <- c( "#000000",brewer.pal(4, "Set2")) 
  
  
  plot0 <- ggplot(data=plot_dat) + 
    geom_point(aes(y=p0_cc_mean, x=alpha_setting, col=missing_approach), position=dodge) +
    geom_rect(data = plot_dat %>% filter(beta_setting == scenario), 
              aes(xmin = as.numeric(factor(alpha_setting)) - 0.5, 
                  xmax = as.numeric(factor(alpha_setting)) + 0.5, 
                  ymin = -Inf, ymax = Inf), 
              color = "grey95", alpha=0, linetype="dashed") + 
    geom_hline(aes(yintercept=true_success[true_success$beta==scenario,"p0"]), col="black")+
    geom_linerange(aes(ymin=p0_cc_mean-1.96*p0_cc_sd/sqrt(N), 
                       ymax=p0_cc_mean+1.96*p0_cc_sd/sqrt(N), x=alpha_setting, color=missing_approach), position=dodge) + 
    
    geom_point(aes(y=p0_single_imp_mean, x=alpha_setting, col=missing_approach), position=dodge, shape=8) +
    geom_linerange(aes(ymin=p0_single_imp_mean-1.96*p0_single_imp_sd/sqrt(N), 
                       ymax=p0_single_imp_mean+1.96*p0_single_imp_sd/sqrt(N), x=alpha_setting, color=missing_approach), position=dodge) + 
    
    facet_wrap(~method, ncol=1) + coord_flip()+ 
    
    xlab("Missing Data Mechanisms") + ylab(expression(hat(p)[0])) + labs(col = "Missing data approach") +
    theme_light() +  # Clean background with light theme
    guides(alpha="none") +theme(axis.text.x = element_text(angle = 45, vjust=0.5)) +
    scale_color_manual(name="Missing data approach", values=cbp2) + scale_x_discrete(limits=rev)+
    background #+ ylim(c(0.08, 0.32))
  
  
  
  plot0
  
  
  
  plot1 <- ggplot(data=plot_dat) + 
    geom_point(aes(y=p1_cc_mean, x=alpha_setting, col=missing_approach), position=dodge) +
    geom_rect(data = plot_dat %>% filter(beta_setting == scenario), 
              aes(xmin = as.numeric(factor(alpha_setting)) - 0.5, 
                  xmax = as.numeric(factor(alpha_setting)) + 0.5, 
                  ymin = -Inf, ymax = Inf), 
              color = "grey95", alpha=0, linetype="dashed") + 
    geom_hline(aes(yintercept=true_success[true_success$beta==scenario,"p1"]), col="black")+
    
    geom_linerange(aes(ymin=p1_cc_mean-1.96*p1_cc_sd/sqrt(N), 
                       ymax=p1_cc_mean+1.96*p1_cc_sd/sqrt(N), x=alpha_setting, color=missing_approach), position=dodge) + 
    
    geom_point(aes(y=p1_single_imp_mean, x=alpha_setting, col=missing_approach), position=dodge, shape=8) +
    geom_linerange(aes(ymin=p1_single_imp_mean-1.96*p1_single_imp_sd/sqrt(N), 
                       ymax=p1_single_imp_mean+1.96*p1_single_imp_sd/sqrt(N), x=alpha_setting, color=missing_approach), position=dodge) + 
    
    facet_wrap(~method, ncol=1) + coord_flip()+ 
    
    xlab("Missing Data Mechanisms") + ylab(expression(hat(p)[1])) + labs(col = "Missing data approach") +
    theme_light() +  # Clean background with light theme
    guides(alpha="none") +theme(axis.text.x = element_text(angle = 45, vjust=0.5)) +
    scale_color_manual(name="Missing data approach", values=cbp2) + scale_x_discrete(limits=rev)+
    background 
  
  plot1
  
  
  p = ggarrange(plot0, plot1, ncol=2, common.legend = T, legend="right")
  
  
  #retain correct path
  png(paste0("output/bias_forest_smoking/bias_", scenario, "_smoking.png"), height=25, width=23, units="cm", res=1200)
  #png(paste0("output/bias_forest_smoking/bias_", scenario, "_smoking_burnin.png"), height=25, width=23, units="cm", res=1200)
  #png(paste0("output/bias_forest_smoking/bias_", scenario, "_all_smoking_burnin.png"), height=35, width=23, units="cm", res=1200)
  #png(paste0("output/bias_forest_smoking/bias_", scenario, "_all_smoking.png"), height=35, width=23, units="cm", res=1200)
  
  print(p)
  dev.off()
  
}





