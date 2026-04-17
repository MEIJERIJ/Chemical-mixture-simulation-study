################################################################################
# Simulation study 
#
# Visualising results
#
# 19/11/2025 - Jonas Meijerink
################################################################################


library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)
library(purrr)
library(forcats)
library(ggdist)


### Load previous results
#-------------------------
out <- checkpoint_results$results |> 
  flatten()
dfs <- purrr::keep(out, is.data.frame)
output <- dplyr::bind_rows(dfs)
out_bkmr <- checkpoint_results_bkmr$results |> 
  flatten()
dfs_bkmr <- purrr::keep(out_bkmr, is.data.frame)
output_bkmr <- dplyr::bind_rows(dfs_bkmr)


output <- rbind(output, output_bkmr)


### Define method order
#----------------------
method_order <- c(
  "Single pollutant linear model",
  "Multiple pollutant linear model",
  "Ridge regression",
  "Lasso regression",
  "Elastic Net regression",
  "Horseshoe regression",
  "Regularised horseshoe regression",
  "BART",
  "DART",
  "SoftBART",
  "SoftBART (non-sparse)",
  "SoftBART semi-parametric",
  "BKMR",
  "Bayesian model averaging",
  "WQS regression (deciles)",
  "WQS regression (quartiles)",
  "WQS regression (standardised continuous)"
)


### Define group order
#---------------------
group_order <- c(
  "OLS methods",
  "Shrinkage methods",
  "Model averaging",
  "Weighted index",
  "Non- or semi-parametric"
)


### Define parameter order
#-------------------------
parameter_order <- c("Joint",
                     "bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda", "pfna", 
                     "pfos")



### Generate a color vector of the same length
#---------------------------------------------
method_colors <- c(
  "Single pollutant linear model"                 = "#000000FF",  
  "Multiple pollutant linear model"               = "#4F4D4ABB",  
  "Ridge regression"                              = "#0E00A8FF",  
  "Lasso regression"                              = "#0C4EC2FF",
  "Elastic Net regression"                        = "#7C5EC2FF",
  "Horseshoe regression"                          = "#9511A1FF",
  "Regularised horseshoe regression"              = "#A62395FF",
  "BART"                                          = "#BC3488FF",
  "DART"                                          = "#CC4678FF",
  "SoftBART"                                      = "#DA596AFF",
  "SoftBART (non-sparse)"                         = "#E56B5DFF",
  "SoftBART semi-parametric"                      = "#F07F4FFF",
  "BKMR"                                          = "#F89441FF",
  "Bayesian model averaging"                      = "#FDAB33FF",
  "WQS regression (deciles)"                      = "#FDC928FF",
  "WQS regression (quartiles)"                    = "#F9ED65FF",
  "WQS regression (standardised continuous)"      = "#F0F921FF"
)
names(method_colors) <- method_order


### Assign groups
#----------------
output <- output %>%
  mutate(Method = factor(Method, levels = method_order),  # Correct method order
         Group = case_when(
           Method %in% c("Multiple pollutant linear model","Single pollutant linear model") ~ "OLS methods",
           Method %in% c("Elastic Net regression",
                         "Lasso regression", "Ridge regression",
                         "Horseshoe regression", "Regularised horseshoe regression") ~ "Shrinkage methods",
           Method %in% c("BKMR", "BART",
                         "DART","SoftBART","SoftBART (non-sparse)",
                         "SoftBART semi-parametric") ~ "Non- or semi-parametric",
           Method == "Bayesian model averaging" ~ "Model averaging",
           TRUE ~ "Weighted index"
         ),
         Group = factor(Group, levels = group_order)) %>%  # Correct group order
  mutate(parameter = factor(parameter, levels = parameter_order))


################################################################################
#
# Joint effect overview Bias, CI width
#
################################################################################


### Filter joint effect and leaf out some methods
#------------------------------------------------
joint <- subset(output, parameter == "Joint" 
                & Method != "Regularised horseshoe regression"
                & Method != "SoftBART (non-sparse)")


### Plot bias
#------------
p_joint_bias <- joint %>%
  filter(  Method != "Bayesian model averaging" 
         & Method != "WQS regression (deciles)" 
         & Method != "WQS regression (quartiles)" 
         & Method != "WQS regression (standardised continuous)" ) %>%
  ggplot( aes(x = Method, y = Estimate + 1.19922, fill = Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  ggh4x::facet_nested(
    cols = vars(round(SNR,2)),
    scales = "free_x",
    space ="free"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  labs(fill = 'Method:') +
  ggtitle("Bias") +
  ylab("Bias") +
  geom_hline(yintercept = 0, col = "red") 


### Plot the CI width
#--------------------
p_joint_width <- joint %>%
  filter(  Method != "Bayesian model averaging" 
         & Method != "WQS regression (deciles)" 
         & Method != "WQS regression (quartiles)" 
         & Method != "WQS regression (standardised continuous)" ) %>%
  ggplot(aes(x=Method, y=CI_width, fill=Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)+
  ggh4x::facet_nested(
    cols = vars(round(SNR,2)),
    scales = "free_x",
    space = "free"
  ) +
  theme_bw() +
  theme(
    legend.position   = "right",
    strip.background  = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )  +
  guides(fill = guide_legend(nrow = 2)) +  
  labs(fill = 'Method:')+
  ggtitle("Width of 95% confidence/credible interval") +
  xlab("") + ylab("Interval width")


### Combine plots with a common legend
#-------------------------------------
final_plot <- ggarrange(
  p_joint_bias + theme(legend.position = "none") ,
  p_joint_width + theme(legend.position = "none"),
  ncol = 1, nrow = 2, heights = c(0.5, 0.5),
  common.legend = TRUE,     
  legend = "bottom")
final_plot <- annotate_figure(
  final_plot,
  top = text_grob(
    "Joint effect by signal-to-noise ratio", 
    face = "bold", size = 16
  )
)
final_plot


################################################################################
#
# PFOS effect overview Bias, CI width
#
################################################################################


### Filter joint effect and leaf out some methods
#------------------------------------------------
pfos <- subset(output, parameter == "pfos")


### Plot bias
#------------
p_pfos_bias <- pfos %>%
  filter(Method != "Bayesian model averaging" & Group != "Non- or semi-parametric" &
           Group != "Weighted index") %>%
  ggplot(aes(x = Method, y = Estimate + 0.4, fill = Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  ggh4x::facet_nested(
    cols = vars( round(SNR,2)),
    scales = "free_x",
    space ="free"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  labs(fill = 'Method:') +
  ggtitle("Bias") +
  ylab("Bias") +
  geom_hline(yintercept = 0, col = "red") 


### Plot the CI width
#--------------------
p_pfos_width <- pfos %>%
  filter(Method != "Bayesian model averaging" & Group != "Non- or semi-parametric" &
           Group != "Weighted index") %>%
  ggplot(aes(x=Method, y=CI_width, fill=Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)+
  ggh4x::facet_nested(
    cols = vars( round(SNR,2)),
    scales = "free_x",
    space = "free"
  ) +
  theme_bw() +
  theme(
    legend.position   = "right",
    strip.background  = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )  +
  guides(fill = guide_legend(nrow = 2)) +  
  labs(fill = 'Method:')+
  ggtitle("Width of 95% confidence/credible interval") +
  xlab("") + ylab("Interval width")


### Combine plots with a common legend
#-------------------------------------
pfos_plot <- ggarrange(
  p_pfos_bias + theme(legend.position = "none") ,
  p_pfos_width + theme(legend.position = "none"),
  ncol = 1, nrow = 2, heights = c(0.5, 0.5),
  common.legend = TRUE,     
  legend = "bottom")
pfos_plot <- annotate_figure(
  pfos_plot,
  top = text_grob(
    "PFOS effect by signal-to-noise ratio", 
    face = "bold", size = 16
  )
)
pfos_plot


################################################################################
#
# PFOA effect overview Bias, CI width
#
################################################################################


### Filter lbpfoa and leaf out some methods
#------------------------------------------
pfoa <- subset(output, parameter == "lbpfoa")


### Plot bias
#------------
p_pfoa_bias <- pfoa %>%
  filter(Method != "Bayesian model averaging" & Group != "Non- or semi-parametric" &
           Group != "Weighted index") %>%
  ggplot(aes(x = Method, y = Estimate + 0.1, fill = Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  ggh4x::facet_nested(
    cols = vars( round(SNR,2)),
    scales = "free_x",
    space ="free"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  labs(fill = 'Method:') +
  ggtitle("Bias") +
  ylab("Bias") +
  geom_hline(yintercept = 0, col = "red") 


### Plot the CI width
#--------------------
p_pfoa_width <- pfoa %>%
  filter(Method != "Bayesian model averaging" & Group != "Non- or semi-parametric" &
           Group != "Weighted index") %>%
  ggplot(aes(x=Method, y=CI_width, fill=Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)+
  ggh4x::facet_nested(
    cols = vars( round(SNR,2)),
    scales = "free_x",
    space = "free"
  ) +
  theme_bw() +
  theme(
    legend.position   = "right",
    strip.background  = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )  +
  guides(fill = guide_legend(nrow = 2)) +  
  labs(fill = 'Method:')+
  ggtitle("Width of 95% confidence/credible interval") +
  xlab("") + ylab("Interval width")


### Combine plots with a common legend
#-------------------------------------
pfoa_plot <- ggarrange(
  p_pfoa_bias + theme(legend.position = "none") ,
  p_pfoa_width + theme(legend.position = "none"),
  ncol = 1, nrow = 2, heights = c(0.5, 0.5),
  common.legend = TRUE,     
  legend = "bottom")
pfoa_plot <- annotate_figure(
  pfoa_plot,
  top = text_grob(
    "PFOA effect by signal-to-noise ratio", 
    face = "bold", size = 16
  )
)
pfoa_plot


################################################################################
#
# PFOS (branched) effect overview Bias, CI width
#
################################################################################


### Filter joint effect and leaf out some methods
#------------------------------------------------
bpfos <- subset(output, parameter == "bpfos")


### Plot bias
#------------
p_bpfos_bias <- bpfos %>%
  filter(Method != "Bayesian model averaging" & Group != "Non- or semi-parametric" &
           Group != "Weighted index") %>%
  ggplot(aes(x = Method, y = Estimate + 0.1, fill = Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  ggh4x::facet_nested(
    cols = vars( round(SNR,2)),
    scales = "free_x",
    space ="free"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  labs(fill = 'Method:') +
  ggtitle("Bias") +
  ylab("Bias") +
  geom_hline(yintercept = 0, col = "red") 


### Plot the CI width
#--------------------
p_bpfos_width <- bpfos %>%
  filter(Method != "Bayesian model averaging" & Group != "Non- or semi-parametric" &
           Group != "Weighted index") %>%
  ggplot(aes(x=Method, y=CI_width, fill=Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)+
  ggh4x::facet_nested(
    cols = vars( round(SNR,2)),
    scales = "free_x",
    space = "free"
  ) +
  theme_bw() +
  theme(
    legend.position   = "right",
    strip.background  = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )  +
  guides(fill = guide_legend(nrow = 2)) +  
  labs(fill = 'Method:')+
  ggtitle("Width of 95% confidence/credible interval") +
  xlab("") + ylab("Interval width")


### Combine plots with a common legend
#-------------------------------------
bpfos_plot <- ggarrange(
  p_bpfos_bias + theme(legend.position = "none") ,
  p_bpfos_width + theme(legend.position = "none"),
  ncol = 1, nrow = 2, heights = c(0.5, 0.5),
  common.legend = TRUE,     
  legend = "bottom")
bpfos_plot <- annotate_figure(
  bpfos_plot,
  top = text_grob(
    "PFOS (branched) effect by signal-to-noise ratio", 
    face = "bold", size = 16
  )
)
bpfos_plot


################################################################################
#
# PFBA effect overview Bias, CI width
#
################################################################################


### Filter joint effect and leaf out some methods
#------------------------------------------------
pfba <- subset(output, parameter == "pfba")


### Plot bias
#------------
p_pfba_bias <- pfba %>%
  filter(Method != "Bayesian model averaging" & Group != "Non- or semi-parametric" &
           Group != "Weighted index") %>%
  ggplot(aes(x = Method, y = Estimate + 0.1, fill = Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  ggh4x::facet_nested(
    cols = vars( round(SNR,2)),
    scales = "free_x",
    space ="free"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  labs(fill = 'Method:') +
  ggtitle("Bias") +
  ylab("Bias") +
  geom_hline(yintercept = 0, col = "red") 


### Plot the CI width
#--------------------
p_pfba_width <- pfba %>%
  filter(Method != "Bayesian model averaging" & Group != "Non- or semi-parametric" &
           Group != "Weighted index") %>%
  ggplot(aes(x=Method, y=CI_width, fill=Method)) +
  geom_jitter(aes(color = Method), width = 0.1, alpha = 0.05, 
              show.legend = FALSE) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)+
  ggh4x::facet_nested(
    cols = vars( round(SNR,2)),
    scales = "free_x",
    space = "free"
  ) +
  theme_bw() +
  theme(
    legend.position   = "right",
    strip.background  = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )  +
  guides(fill = guide_legend(nrow = 2)) +  
  labs(fill = 'Method:')+
  ggtitle("Width of 95% confidence/credible interval") +
  xlab("") + ylab("Interval width")


### Combine plots with a common legend
#-------------------------------------
pfba_plot <- ggarrange(
  p_pfba_bias + theme(legend.position = "none") ,
  p_pfba_width + theme(legend.position = "none"),
  ncol = 1, nrow = 2, heights = c(0.5, 0.5),
  common.legend = TRUE,     
  legend = "bottom")
pfba_plot <- annotate_figure(
  pfba_plot,
  top = text_grob(
    "PFBA effect by signal-to-noise ratio", 
    face = "bold", size = 16
  )
)
pfba_plot


################################################################################

##############################  Power overview  ################################

################################################################################


### Summarize empirical power using exact binomial confidence intervals
#----------------------------------------------------------------------
power_summary <- output %>%
  filter(Method != "Bayesian model averaging") %>%
  filter(
    parameter == "Joint" |
      (parameter != "Joint" &
         !Group %in% c("Weighted index", "Non- or semi-parametric"))
  ) %>%
  
  group_by(Method, Size,Group, SNR, parameter) %>%
  summarise(
    Positives = sum(Power == 1, na.rm = TRUE), 
    n = n(),                                   
    Power = Positives / n,                     
    CI_lower = binom.test(Positives, n, conf.level = 0.95)$conf.int[1],
    CI_upper = binom.test(Positives, n, conf.level = 0.95)$conf.int[2],
    .groups = 'drop'
  ) %>%
  mutate(parameter = dplyr::recode(parameter,
                                 "bpfos"   = "PFOS (branched)",
                                 "lbpfhxs" = "PFHXS (total)",
                                 "lbpfoa"  = "PFOA (total)",
                                 "pfba"    = "PFBA",
                                 "pfda"    = "PFDA",
                                 "pfna"    = "PFNA",
                                 "pfos"    = "PFOS",
                                 "Joint"   = "Joint"
  )) 


### Plot all parameter for OLS and shrinkage methods
#--------------------------------------------------
shape_vals <- c(16, 17, 15, 7, 8)

power_summary %>%
  filter(Group == "OLS methods" | Group == "Shrinkage methods") %>%
  filter(Method != "Regularised horseshoe regression" & Method != "Lasso regression") %>%
  ggplot(aes(x = factor(round(SNR, 2)), y = Power,
             color = Method, group = Method, shape = Method)) +
  geom_point(size = 2, alpha = 1, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.5), width = 0.5,
                alpha = 1, show.legend = FALSE, linewidth = 0.7) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = shape_vals) +
  facet_wrap(~parameter, ncol = 3) +
  theme_bw(base_size = 11) +
  ylim(0, 1.01) +
  ggtitle("Empirical power") +
  ylab("Power") + xlab("Signal-to-noise ratio") +
  guides(
    color = guide_legend(ncol = 1, override.aes = list(linewidth = 1.2, size = 3)),
    shape = guide_legend(ncol = 1)
  ) +
  labs(color = NULL, shape = NULL) +       
  theme(
    strip.background   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(size = 16),
    axis.title.y       = element_text(size = 14),
    axis.title.x       = element_text(size = 14),
    legend.position    = "right",           e
    legend.text        = element_text(size = 12),
    legend.key.width   = unit(1.2, "cm")
  )




### Plot all parameter for non- or semi-parametric
#-------------------------------------------------
shape_vals <- c(18, 0, 1, 2, 4, 5, 6, 9, 10, 11, 12)
power_summary %>%
  filter(Group != "OLS methods" & Group != "Shrinkage methods") %>%
  ggplot(aes(x=factor(round(SNR,2)), y=Power, color=Method, group=Method, shape = Method)) +
  geom_point(size = 2, alpha = 1, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.5), width = 0.5, alpha=1, show.legend = FALSE, size = 0.7) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~parameter, ncol = 3)+
  theme_bw(base_size = 11) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = shape_vals) +
  ylim(0,1) + 
  guides(color = guide_legend(ncol = 3)) + 
  ggtitle("Empirical power") +
  ylab("Power") + xlab("Signal-to-noise ratio")  +
  labs(fill = 'Method:') + 
  theme(
    strip.background  = element_blank(),           
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.position = "bottom",          
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 11.5)
  ) 


### Plot all parameter all methods
#---------------------------------
shape_vals <- c(16, 17, 15, 4, 7, 8, 18, 0, 1, 2, 4, 5, 6, 9, 10, 11, 12)
power_summary %>%
  ggplot(aes(x=factor(round(SNR,2)), y=Power, color=Method, group=Method, shape = Method)) +
  geom_point(size = 2, alpha = 1, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.5), width = 0.5, alpha=1, show.legend = FALSE, size = 0.7) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~parameter, ncol = 2)+
  theme_bw(base_size = 11) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = shape_vals) +
  ylim(0,1) + 
  guides(color = guide_legend(ncol = 3)) + 
  ggtitle("Empirical power") +
  ylab("Power") + xlab("Signal-to-noise ratio")  +
  labs(fill = 'Method:') + 
  theme(
    strip.background  = element_blank(),            
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.position = "bottom",           
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 11.5)
  ) 



################################################################################

############################  Coverage overview  ###############################

################################################################################


### Summarize empirical power using exact binomial confidence intervals
#----------------------------------------------------------------------
coverage_summary <- output %>%
  filter(Method != "Bayesian model averaging" & Group != "Weighted index") %>%
  filter(
    parameter == "Joint" |
      (parameter != "Joint" &
         !Group %in% c("Weighted index", "Non- or semi-parametric"))
  ) %>%
  
  group_by(Method, Size,Group, SNR, parameter) %>%
  summarise(
    Positives = sum(CI_coverage == 1, na.rm = TRUE),  
    n = n(),                                    
    Coverage = Positives / n,                      # Empirical power
    CI_lower = binom.test(Positives, n, conf.level = 0.95)$conf.int[1],
    CI_upper = binom.test(Positives, n, conf.level = 0.95)$conf.int[2],
    .groups = 'drop'
  ) %>%
  mutate(parameter = dplyr::recode(parameter,
                                   "bpfos"   = "PFOS (branched)",
                                   "lbpfhxs" = "PFHXS (total)",
                                   "lbpfoa"  = "PFOA (total)",
                                   "pfba"    = "PFBA",
                                   "pfda"    = "PFDA",
                                   "pfna"    = "PFNA",
                                   "pfos"    = "PFOS",
                                   "Joint"   = "Joint"
  )) 


### Plot all parameter all methods
#---------------------------------
shape_vals <- c(16, 17, 15, 4, 7, 8, 18, 0, 1, 2, 4, 5, 6, 9, 10, 11, 12)
coverage_summary %>%
  ggplot(aes(x=factor(round(SNR,2)), y=Coverage, color=Method, group=Method, shape = Method)) +
  geom_point(size = 2, alpha = 1, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(0.5), width = 0.5, alpha=1, show.legend = FALSE, size = 0.7) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~parameter, ncol = 3)+
  theme_bw(base_size = 11) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = shape_vals) +
  ylim(0,1) + 
  guides(color = guide_legend(ncol = 3)) +  # Legend with 4 rows
  ggtitle("Empirical 95% credible/confidence interval coverage") +
  ylab("Coverage") + xlab("Signal-to-noise ratio")  +
  labs(fill = 'Method:') + 
  theme(
    strip.background  = element_blank(),           
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.position = "bottom",           
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 11.5)
  ) +
  geom_hline(yintercept = 0.95, col = "red") 




################################################################################

########################  Variable importance  #################################

################################################################################



### Filter on PIP
#----------------
pip <- output %>% filter(!is.na(PIP))


### Scale for BART and DART
#--------------------------
pip2 <- pip %>%
  group_by(Method) %>%
  mutate(
    PIP_plot = if_else(
      Method %in% c("BART", "DART"),
      (PIP - min(PIP, na.rm = TRUE)) /
        (max(PIP, na.rm = TRUE) - min(PIP, na.rm = TRUE)),
      PIP
    )
  ) %>%
  ungroup() %>%
  mutate(parameter = dplyr::recode(parameter,
                                   "bpfos"   = "PFOS (branched)",
                                   "lbpfhxs" = "PFHXS (total)",
                                   "lbpfoa"  = "PFOA (total)",
                                   "pfba"    = "PFBA",
                                   "pfda"    = "PFDA",
                                   "pfna"    = "PFNA",
                                   "pfos"    = "PFOS",
                                   "Joint"   = "Joint"
  )) 


### Plot 
#-------
ggplot(pip2, aes(y = parameter, x = PIP_plot, color = Method)) +
  stat_halfeye()+
  geom_point(size = 0.7, alpha = 0.15, show.legend = FALSE) +
  facet_grid(Method ~ round(SNR,2), drop = TRUE, scales = "free",space = "free") + 
  theme_bw(base_size = 11) +
  scale_fill_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  theme(
    legend.position   = "right",
    strip.background  = element_blank(),           
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    #panel.border = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 5)) + 
  ggtitle("Variable importance") +
  ylab("") + xlab("")  +
  labs(fill = 'Method:')

