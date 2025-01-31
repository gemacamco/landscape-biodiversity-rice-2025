# Clean complete code: loop

# install and load packages
install.packages("metafor")
install.packages("dplyr")
install.packages("snow")
install.packages("jtools")
install.packages("ggplot2")

library(metafor)
library(dplyr)
library(snow)
library(jtools)
library(ggplot2)

################################################################################

# load dataset
setwd("C:/R_files/Meta-analisis/data/meta")
data <- read.csv("meta_ajust.csv", sep = ";") 

################################################################################

# calculate effect size
data <- escalc(measure = "ZCOR", ri = r_dir, ni = n_total, data = data)

# function to build the model
build_model <- function(data, mods = NULL) {
  # Only include the 'mods' argument if it's not NULL
  if (is.null(mods)) {
    # Call rma.mv without the mods argument
    rma.mv(yi, 
           vi, 
           random = ~ 1 | ref / id_r, 
           data = data, 
           method = "REML", 
           level = 95)
  } else {
    # Call rma.mv with the mods argument if it's provided
    rma.mv(yi, 
           vi, 
           random = ~ 1 | ref / id_r, 
           data = data, 
           mods = mods, 
           method = "REML", 
           level = 95)
  }
}

# meta-analysis options

# 1.1. Landscape complexity effect on biodiversity -> "general"                   # Faltaría mod8 y mod9
# 1.2. Landscape complexity effect on vertebrates -> "land_vert"
# 1.3. Landscape complexity effect on invertebrates -> "land_invert"

# 2.1. Landscape dimensions effect on biodiversity -> "dim_gen"
# 2.2. Landscape dimensions effect on vertebrates -> "dim_vert"
# 2.3. Landscape dimensions effect on invertebrates -> "dim_invert"

# 3.1. Landscape dimensions effect on invertebrate's orders -> "dim_order"

# Create an empty list to store models for each analysis
model_list <- list()

# Select the analysis from options above
LB <- "dim_order"

# Use switch to select the appropriate data and model structure
invisible(
  switch(LB,
         "general" = {
           data_subset <- data
           mod <- build_model(data_subset, mods = NULL)
           model_list[["general"]] <- list(mod = mod, data = data_subset)  # Store the model and data
         },
         "land_vert" = {
           data_subset <- filter(data, tax == "vertebrate")
           mod <- build_model(data_subset, mods = NULL)
           model_list[["land_vert"]] <- list(mod = mod, data = data_subset)  # Store the model and data
         },
         "land_invert" = {
           data_subset <- filter(data, tax == "invertebrate")
           mod <- build_model(data_subset, mods = NULL)
           model_list[["land_invert"]] <- list(mod = mod, data = data_subset)  # Store the model and data
         },
         "dim_gen" = {
           data_subset <- data
           mod <- build_model(data_subset, mods = ~g_indicator)
           model_list[["dim_gen"]] <- list(mod = mod, data = data_subset)  # Store the model and data
         },
         "dim_vert" = {
           data_subset <- filter(data, tax == "vertebrate")
           mod <- build_model(data_subset, mods = ~g_indicator)
           model_list[["dim_vert"]] <- list(mod = mod, data = data_subset)  # Store the model and data
         },
         "dim_invert" = {
           data_subset <- filter(data, tax == "invertebrate")
           mod <- build_model(data_subset, mods = ~g_indicator)
           model_list[["dim_invert"]] <- list(mod = mod, data = data_subset)  # Store the model and data
         },
         "dim_order" = { 
           ordenes <- c("Araneae", "Coleoptera", "Hemiptera", "Lepidoptera", "Hymenoptera", "Odonata") 
           
           for (nivel in ordenes) {
             data_subset <- subset(data, orden == nivel)
             
             if (length(unique(data_subset$g_indicator)) > 1) {
               mod <- build_model(data_subset, mods = ~g_indicator)
               model_list[[nivel]] <- list(mod = mod, data = data_subset)  # Store the model and data
             } else {
               cat("Skipping order:", nivel, "due to insufficient levels of 'g_indicator'.\n")
             }
           }
         }
  )
)

# Now process the models stored in the model_list
for (order_name in names(model_list)) {
  mod_data <- model_list[[order_name]]
  mod <- mod_data$mod
  data_subset <- mod_data$data
  
  # 1. Model summary
  cat("Summary for", order_name, ":\n")
  print(summary(mod))
  
  # 2. Extract coefficients and back-transform them
  coef <- data.frame(
    Term = rownames(mod$b),  
    Estimate = as.numeric(mod$b),
    CI.Lower = as.numeric(mod$ci.lb),
    CI.Upper = as.numeric(mod$ci.ub),
    pval = as.numeric(mod$pval)
  )
  
  coef_bt <- data.frame(
    Term = coef$Term,
    Estimate = transf.ztor(coef$Estimate),
    CI.Lower = transf.ztor(coef$CI.Lower),
    CI.Upper = transf.ztor(coef$CI.Upper),
    pval = coef$pval
  )
  
  # 3. Profile likelihood plots
  par(mfrow=c(1,1))
  profile(mod, sigma2=1)  # Reference
  profile(mod, sigma2=2)  # id_r
  
  # 4. Forest plot
  forest(mod,
         addpred=TRUE, 
         header=TRUE, 
         atransf=transf.ztor, 
         annotate=TRUE, 
         order="obs", 
         xlim=c(-1.5, 1.5), 
         addfit=TRUE, 
         slab=data_subset$id_r)  # Use the subsetted data
  
  # 5. I2 for multilevel and multivariate models
  W <- diag(1/mod$vi)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  i2 <- 100 * sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P)))
  cat("I2 for", order_name, ":", i2, "\n")
  
  # 6. Graphical display
  plot <- ggplot(coef_bt, aes(x = Estimate, y = Term)) +
    geom_point(size = 3, color = "purple") +
    geom_errorbarh(aes(xmin = CI.Lower, xmax = CI.Upper), height = 0.1, color = "purple", linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    theme_minimal() +
    labs(title = paste("Effect of Landscape Complexity on", order_name), x = "Pearson's r", y = NULL)
  
  print(plot)
}

