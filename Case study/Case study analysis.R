################################################################################
# 
# Case study: Association between PFAS mixture and BMI z-score
# 
# Data used from the 3M hotspot study
# 
# 16/03/2026
################################################################################


rm(list = ls())


### Load libraries
#-----------------
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(gWQS)
library(anthroplus)


### Load data
#------------
df_original <- read_excel("~/Dataset_HBM_3M.xlsx")


### Set working directory
#------------------------
setwd("~/Method")


### Load functions
#-----------------
source("~/Method/slr.R")
source("~/Method/mlr.R")
source("~/Method/ridge.R")
source("~/Method/lasso.R")
source("~/Method/enet.R")
source("~/Method/horseshoe.R")
source("~/Method/reghorseshoe.R")
source("~/Method/softbart.R")
source("~/Method/bkmr.R")
source("~/Method/bas.R")
source("~/Method/rhwqs.R")


################################################################################
#                                                                              #
#---------------------------- Data manipulations ------------------------------#
#                                                                              #
################################################################################


### Calculate the BMI z-score
#----------------------------
df_original$zbmi <- anthroplus_zscores(
  sex = df_original$GESL,
  age_in_months = df_original$LEEFTIJD*12,
  height_in_cm = df_original$OB_Lengte,
  weight_in_kg = df_original$OB_Gewicht)[["zbfa"]]


### Filter on exclusion parameters
#---------------------------------
df <- df_original %>%
  filter(RB_D1==0 & GENEESMIDDEL_SCHILDKLIERAAND==0 & GENEESMIDDEL_DIABETES==0 & 
           GENEESMIDDEL_NIERZIEKTE==0)
# CORTICOSTEROID had only missing values after excluding previous ones


### Make gender a dummy variable (0 male/1 female)
#-------------------------------------------------
df$GESL <- (df$GESL-1)


### Create dummy variables for age (1=[12.5, 14.5] age;2=]14.5, 15.5] age;3=> 15.5 age)
#--------------------------------------------------------------------------------------
df$Age2 <- ifelse(df$AGE==2,1,0)
df$Age3 <- ifelse(df$AGE==3,1,0)


### Create dummy variables for income 
# 1=heel erg moeilijk tot moeilijk;
# 2=lukt om rond te komen;
# 3=lukt om comfortabel te leven
#----------------------------------------------------------------------
df$Income2 <- ifelse(df$RONDKOMEN2==2,1,0)
df$Income3 <- ifelse(df$RONDKOMEN2==3,1,0)


### Create dataset
#-----------------
obs <- data.frame(
  "Y" = df$zbmi,
  "bpfos" = scale(log(df$bpfos_imp)) , 
  "lbpfhxs" = scale(log(df$lbpfhxs_imp)), 
  "lbpfoa"=scale(log(df$lbpfoa_imp)),
  "pfba" = scale(log(df$pfba_imp)), 
  "pfda" = scale(log(df$pfda_imp)) , 
  "pfna" = scale(log(df$pfna_imp)) , 
  "pfos" = scale(log(df$pfos_imp)),
  "age2"=df$Age2, 
  "age3"=df$Age3,
  "gender"=df$GESL, 
  "income2"=df$Income2,
  "income3"=df$Income3
)
obs<-na.omit(obs)


### Select original exposures
#----------------------------
exposures <- data.frame("bpfos" = (df$bpfos_imp), "lbpfhxs" = (df$lbpfhxs_imp), 
                        "lbpfoa"=(df$lbpfoa_imp), "pfba" = (df$pfba_imp), 
                        "pfda" = (df$pfda_imp), "pfna" = (df$pfna_imp), 
                        "pfos" = (df$pfos_imp))
exposures <- na.omit(exposures)


### Extract SDs used for scaling
#-------------------------------
scaled_data <- scale(log(exposures[1:7]))
sds <- attr(scaled_data, "scaled:scale")
mu <- attr(scaled_data, "scaled:center")


### Construct the log transformed scaled quantiles
#-------------------------------------------------
q_low <- (log(apply(exposures[, 1:7], 2, quantile, probs = 0.25))-mu)/sds
q_high <- (log(apply(exposures[, 1:7], 2, quantile, probs = 0.75))-mu)/sds



### Remove all unnecessary variables/data
#----------------------------------------
rm(list=setdiff(ls(), c("obs", "exposures","q_low", "q_high", 
                        "slr", "mlr", "ridge", "lasso", "enet", 
                        "horseshoe", "reghorseshoe", "softbart",
                        "bkmr", "bas", "rhwqs"
                        )))


################################################################################
#                                                                              #
#----------------------- Apply the different methods --------------------------#
#                                                                              #
################################################################################


exposure_names <- c("bpfos", "lbpfhxs", "lbpfoa", "pfba", "pfda",  "pfna", "pfos")


### Define seed and run methods
#------------------------------
set.seed(20032026)
coef_slr <- do.call(rbind, lapply(exposure_names, slr, data = obs))
set.seed(20032026)
coef_mlr <- mlr(df = obs, qlow = q_low, qhigh = q_high)
set.seed(20032026)
coef_ridge <- ridge(obs, qlow = q_low, qhigh = q_high)
set.seed(20032026)
coef_lasso <- lasso(obs, qlow = q_low, qhigh = q_high)
set.seed(20032026)
coef_enet <- enet(obs, qlow = q_low, qhigh = q_high)
set.seed(20032026)
coef_horseshoe <- horseshoe(df = obs, qlow = q_low, qhigh = q_high)
set.seed(20032026)
coef_reghorseshoe <- reghorseshoe(df = obs, qlow = q_low, qhigh = q_high)
set.seed(20032026)
coef_softbart <- softbart(df = obs, q_low = q_low, q_high = q_high)
set.seed(20032026)
coef_bkmr <- bkmr(df = obs, qlow = q_low, qhigh = q_high)
set.seed(20032026)
coef_bas <- bas(obs = obs)
set.seed(20032026)
coef_rhwqs <- rhwqs(obs = obs)


################################################################################
#                                                                              #
#----------------------------- Visualise results ------------------------------#
#                                                                              #
################################################################################


### Bind all results
#-------------------
out <- rbind(coef_slr, coef_mlr, coef_ridge, coef_lasso, coef_enet, coef_horseshoe,
             coef_reghorseshoe, coef_softbart, coef_bkmr, coef_rhwqs, coef_bas)


### Add stars as these methods have slightly different interpretation
#--------------------------------------------------------------------
out <- out |>
  mutate(Method = case_match(Method,
                             "Single pollutant linear model" ~ "Single pollutant linear model*",
                             "WQS regression (quartiles)"    ~ "WQS regression (quartiles)**",
                             "SoftBART (sparse)"             ~ "SoftBART (sparse)***",
                             .default = Method
  ))


### Fix ordering of the exposures
#--------------------------------
exposure_order_fixed <- c("PFOS (branched)", "PFHxS (total)", "PFOA (total)", "PFBA", "PFDA", "PFNA", "PFOS")


### Fix ordering of the methods
#------------------------------
method <- rev(c(
  "Single pollutant linear model*",
  "Multiple pollutant linear model",
  "Ridge regression",
  "Lasso regression",
  "Elastic Net regression",
  "Horseshoe regression",
  "Regularised horseshoe regression",
  "SoftBART (sparse)***",
  "BKMR",
  "WQS regression (quantiles)**"
))


### Apply fixed order
#--------------------
ind_data <- out |>
  filter(parameter != "Joint") |>
  filter(!is.na(Estimate)) |>
  mutate(Method = factor(Method,
                        levels = rev(method))) |>
  mutate(parameter = dplyr::recode(parameter,
                                 "bpfos"   = "PFOS (branched)",
                                 "lbpfhxs" = "PFHxS (total)",
                                 "lbpfoa"  = "PFOA (total)",
                                 "pfba"    = "PFBA",
                                 "pfda"    = "PFDA",
                                 "pfna"    = "PFNA",
                                 "pfos"    = "PFOS",
                                 "Joint"   = "Joint"
  )) |>
  mutate(
    parameter = factor(parameter,
                     levels = rev(exposure_order_fixed))
  ) 

### Color palette — one per method
#---------------------------------
method_colors <- c(
  "Single pollutant linear model*"     = "#264653",
  "Multiple pollutant linear model"    = "#2A9D8F",
  "Ridge regression"                   = "#57CC99",
  "Lasso regression"                   = "#E9C46A",
  "Elastic Net regression"             = "#F4A261",
  "Horseshoe regression"               = "#E76F51",
  "Regularised horseshoe regression"   = "#9B2226",
  "SoftBART (sparse)***"               = "#C77DFF",
  "BKMR"                               = "#457B9D",
  "WQS regression (quantiles)**"       = "#E63946"
)


### A: first 4 exposures
#-----------------------
p_ind_left <- ind_data |>
  filter(parameter %in% c("PFOA (total)", "PFOS (branched)", "PFHXS (total)",
                           "PFNA")) |>
  ggplot(aes(x     = Estimate,
             y     = parameter,
             color = Method)) +
  geom_vline(xintercept = 0,
             linetype   = "dashed",
             color      = "grey50",
             linewidth  = 0.7) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    position  = position_dodge(width = 0.7),
    width     = 0.3,
    linewidth = 0.8
  ) +
  xlim(
    -0.32,0.33
  )+
  geom_point(
    size     = 3,
    position = position_dodge(width = 0.7)
  ) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(
    x     = "BMI z-score per 1 SD increase in log(PFAS)",
    y     = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 10),
    axis.text.y        = element_text(size = 11, face = "bold"),
    axis.text.x        = element_text(size = 9),
    legend.position    = "none"
  )


### B: last 3 exposures
#----------------------
p_ind_right <- ind_data |>
  filter(parameter %in% c("PFOS", "PFBA", "PFDA")) |>
  ggplot(aes(x     = Estimate,
             y     = parameter,
             color = Method)) +
  geom_vline(xintercept = 0,
             linetype   = "dashed",
             color      = "grey50",
             linewidth  = 0.7) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    position  = position_dodge(width = 0.7),
    width     = 0.3,
    linewidth = 0.8
  ) +
  xlim(
    -0.32,0.33
  )+
  geom_point(
    size     = 3,
    position = position_dodge(width = 0.7)
  ) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(
    x     = "BMI z-score per 1 SD increase in log(PFAS)",
    y     = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 10),
    axis.text.y        = element_text(size = 11, face = "bold"),
    axis.text.x        = element_text(size = 9),
    legend.position    = "none"
  )


### Collect shared legend — 2 rows below A and B
#-----------------------------------------------
p_top_horizontal <- (p_ind_left | p_ind_right) +
  plot_layout(guides = "collect") &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(0.45, "cm"),
    plot.margin      = margin(t = 0, r = 0, b = 0, l = 0)
  ) &
  guides(color = guide_legend(nrow         = 2,
                              override.aes = list(size = 5)))


### Variable importance plot
#---------------------------
pip_ranked <- out |>
  filter(parameter != "Joint") |>
  filter(!is.na(PIP)) |>
  filter(Method != "Ridge regression") |>
  mutate(parameter = dplyr::recode(parameter,
                                   "bpfos"   = "PFOS (branched)",
                                   "lbpfhxs" = "PFHxS (total)",
                                   "lbpfoa"  = "PFOA (total)",
                                   "pfba"    = "PFBA",
                                   "pfda"    = "PFDA",
                                   "pfna"    = "PFNA",
                                   "pfos"    = "PFOS"
  )) |>
  group_by(Method) |>
  mutate(rank = rank(-PIP, ties.method = "min")) |>
  ungroup() |>
  mutate(
    Method = factor(Method,
                    levels = out |>
                      filter(parameter != "Joint", !is.na(PIP), Method != "Ridge regression") |>
                      group_by(Method) |>
                      summarise(m = mean(PIP), .groups = "drop") |>
                      arrange(desc(m)) |>
                      pull(Method)
    ),
    parameter = factor(parameter,                  # ← uses already-recoded values
                       levels = group_by(pick(parameter, PIP), parameter) |>
                         summarise(m = mean(PIP), .groups = "drop") |>
                         arrange(m) |>
                         pull(parameter)
    ),
    text_color = ifelse(rank <= 2, "white", "grey20")
  )
pip_ranked <- pip_ranked |>
  mutate(Method = stringr::str_wrap(as.character(Method), width = 12))
p_pip <- ggplot(pip_ranked,
                aes(x = Method, y = parameter, fill = rank)) +
  geom_tile(color = "white", linewidth = 0) +
  geom_text(
    aes(label = round(PIP, 2), color = text_color),
    size     = 3
  ) +
  scale_fill_gradient(
    low    = "#C0392B",
    high   = "#FDDBC7",
    limits = c(1, 7),
    breaks = 1:7,
    labels = paste0("#", 1:7),
    name   = "Rank\n(1=most\nimportant)"
  ) +
  scale_color_identity() +
  labs(
    title = "Method specific variable importance",
    x     = "",
    y     = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_text(size = 12.5),
    axis.text.x     = element_text(angle     = 0,
                                   hjust     = 0.5,
                                   vjust     = 1,
                                   size      = 12.5,
                                   lineheight = 0.85),
    plot.title      = element_text(face = "bold", size = 12),
    legend.position   = "right",
    legend.key.height = unit(1.0, "cm")
  )



### Joint effect plot
#--------------------
p_joint <- out |>
  filter(parameter == "Joint") |>
  filter(!is.na(Estimate)) |>
  ggplot(aes(x     = Estimate,
             y     = reorder(Method, Estimate),
             color = Method)) +
  geom_vline(xintercept = 0,
             linetype   = "dashed",
             color      = "grey50",
             linewidth  = 0.7) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width     = 0.3,
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  scale_color_manual(values = method_colors, guide = "none") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Joint mixture effect (75th vs 25th percentile)",
    x     = "Difference in BMI z-score (95% CI)",
    y     = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 12),
    axis.text          = element_text(size = 11)
  )
