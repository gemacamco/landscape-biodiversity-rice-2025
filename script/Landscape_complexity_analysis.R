
## Article title: Effects of landscape complexity on biodiversity of rice agroecosystems: a meta-analysis 
## Authors: Gema Cambero-Conejero, Carles Alcaraz, Néstor Pérez-Méndez
## Contact: Gema Cambero Conejero; gema.cambero@irta.cat
##          Néstor Pérez Méndez; nestor.perez@irta.cat


################################################################################

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

############################ RUN THIS ##########################################

  # load dataset
  setwd("C:/R_files/git_open/meta/data")
  data <- read.csv("dataset.csv", sep = ";") 
  
  # calculate effect size
    data <- escalc(measure = "ZCOR", ri = r_dir, ni = n_total, data = data)
    
  # function to build the model
      build_model <- function(data, mods = NULL) {
        if(length(unique(data$ref)) <= 2) {  # if there are only two studies, don't include random effects
          random_effect <- ~ 1
        } else {
          random_effect <- ~ 1 | ref / id_r 
        }
        
        if (is.null(mods)) {   # only include the 'mods' argument if it's not NULL
          rma.mv(yi,           # call rma.mv without the mods argument
                 vi, 
                 random = ~ 1 | ref / id_r, 
                 data = data, 
                 method = "REML", 
                 level = 95)
        } else {               # Call rma.mv with the mods argument if it's provided
          rma.mv(yi,           
                 vi, 
                 random = ~ 1 | ref / id_r, 
                 data = data, 
                 mods = mods, 
                 method = "REML", 
                 level = 95)
        }
      }

################################################################################
      
  # META-ANALYSIS LABELS - select the labels of the analysis you want to run 
  
    # 1.1. Landscape complexity effect on biodiversity  -> "general"                   
    # 1.2. Landscape complexity effect on vertebrates   -> "land_vert"
    # 1.3. Landscape complexity effect on invertebrates -> "land_invert"

    # 2.1. Landscape dimensions effect on biodiversity  -> "dim_gen"
    # 2.2. Landscape dimensions effect on vertebrates   -> "dim_vert"
    # 2.3. Landscape dimensions effect on invertebrates -> "dim_invert"

    # 3.1. Landscape dimensions effect on Araneae     -> "dim_Araneae"
    # 3.2. Landscape dimensions effect on Coleoptera  -> "dim_Coleoptera"
    # 3.3. Landscape dimensions effect on Hemiptera   -> "dim_Hemiptera"
    # 3.4. Landscape dimensions effect on Lepidoptera -> "dim_Lepidoptera"
    # 3.5. Landscape dimensions effect on Hymenoptera -> "dim_Hymenoptera"
    # 3.6. Landscape dimensions effect on Odonata     -> "dim_Odonata"
      
    # 4. Landscape complexity effect on biodiversity metrics -> "land_desbio"
      
    # 5. Landscape complexity effect on aquatic and terrestrial biodiversity -> "land_habitat"  

  
    # Insert the choosen label
    LB <- "general"  
    
    # Run this to select the appropriate data and model structure
      switch(LB,
             "general" = {
               data_model <- data
               mods <- NULL
               mod <- build_model(data_model, mods = mods)
             },
             "land_vert" = {
               data_model <- filter(data, tax == "vertebrate")
               mods <- NULL
               mod <- build_model(data_model, mods = mods)
             },
             "land_invert" = {
               data_model <- filter(data, tax == "invertebrate")
               mods <- NULL
               mod <- build_model(data_model, mods = mods)
             },
             "dim_gen" = {
               data_model <- data
               mods <- ~g_indicator
               mod <- build_model(data_model, mods = mods)
             },
             "dim_vert" = {
               data_model <- filter(data, tax == "vertebrate")
               mods <- ~g_indicator
               mod <- build_model(data_model, mods = mods)
             },
             "dim_invert" = {
               data_model <- filter(data, tax == "invertebrate")
               mods <- ~g_indicator
               mod <- build_model(data_model, mods = mods)
             },
             "dim_Araneae" = {
               data_model <- subset(data, orden == "Araneae")
               mods <- ~g_indicator
               mod <- build_model(data_model, mods = mods)
             },
             "dim_Coleoptera" = {
               data_model <- subset(data, orden == "Coleoptera")
               mods <- ~g_indicator
               mod <- build_model(data_model, mods = mods)
             },
             "dim_Hemiptera" = {
               data_model <- subset(data, orden == "Hemiptera")
               mods <- ~g_indicator
               mod <- build_model(data_model, mods = mods)
             },
             "dim_Lepidoptera" = {
               data_model <- subset(data, orden == "Lepidoptera")
               mods <- ~g_indicator
               mod <- build_model(data_model, mods = mods)
             },
             "dim_Hymenoptera" = {
               data_model <- subset(data, orden == "Hymenoptera")
               mods <- ~g_indicator
               mod <- build_model(data_model, mods = mods)
             },
             "dim_Odonata" = {
               data_model <- subset(data, orden == "Odonata")
               mods <- ~g_indicator
               mod <- build_model(data_model, mods = mods)
             },
             "dim_Anura" = {
               data_model <- subset(data, orden == "Anura")
               mods <- NULL
               mod <- build_model(data_model, mods = mods)
             },
             "dim_Bird" = {
               data_model <- subset(data, orden == "Bird")
               mods <- NULL
               mod <- build_model(data_model, mods = mods)
             },
             "land_abuden" = {
               data_model <- filter(data, g_desbio == "abu_den")
               mods <- NULL
               mod <- build_model(data_model, mods = mods)
             },
             "land_richness" = {
               data_model <- filter(data, g_desbio == "richness")
               mods <- NULL
               mod <- build_model(data_model, mods = mods)
             },
             "land_shannon" = {
               data_model <- filter(data, g_desbio == "shannon")
               mods <- NULL
               mod <- build_model(data_model, mods = mods)
             }
      )
    
################################################################################
    
    # 1. Model summary
    cat("Model Summary for Analysis Label:", LB, "\n")
    summary(mod)
    
    # 2. Coefficients table
      # Pearson r values
      coef <- data.frame(
        Term = rownames(mod$b),  
        Estimate = as.numeric(mod$b),
        CI.Lower = as.numeric(mod$ci.lb),
        CI.Upper = as.numeric(mod$ci.ub),
        pval = as.numeric(mod$pval)
      )
      
      # Transformed Pearson r to Fisher's z values
      coef_bt <- data.frame(
        Term = coef$Term,
        Estimate = transf.ztor(coef$Estimate),
        CI.Lower = transf.ztor(coef$CI.Lower),
        CI.Upper = transf.ztor(coef$CI.Upper),
        pval = coef$pval
      )

  
    # 3. Forest plot
    forest(mod,
           addpred=TRUE, 
           header=TRUE, 
           atransf=transf.ztor, 
           annotate=TRUE, 
           order="obs", 
           xlim=c(-1.5, 1.5), 
           addfit=TRUE, 
           slab=data_model$id_r)  # Use the subset data
    
    # 4. I2 for multilevel and multivariate models
     # From: https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
     # Overall I2 -> indicated how much of the total variance can be attributed 
     # to the total amount of heterogeneity (which is the sum of between- and 
     # within-cluster heterogeneity).
      W <- diag(1/mod$vi)
      X <- model.matrix(mod)
      P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
      i2 <- 100 * sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P)))
      cat("Overall I2 for", LB, ":", i2, "\n")
      
     # Partial I2 -> how much of the total variance can be attributed to between- 
     # and within-cluster heterogeneity separately  
      i2p <- 100 * mod$sigma2 / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))
      cat("Between- and within-cluster I2 for", LB, ":", i2p, "\n")
    
    # 5. Graphical display
    plot <- ggplot(coef_bt, aes(x = Estimate, y = Term)) +
      geom_point(size = 3, color = "purple") +
      geom_errorbarh(aes(xmin = CI.Lower, xmax = CI.Upper), height = 0.1, color = "purple", linewidth = 1) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
      theme_minimal() +
      labs(title = paste(LB), x = "Pearson's r", y = NULL)
    
    print(plot)
    
    # 6. Outliers analysis
    # NOTE: outliers analysis cannot be performed for models with a small number 
    # of observations (i.e., Lepidoptera and Odonata)
    
      # calculate Cook's distance for each effect size
      x <- cooks.distance(mod, parallel = "snow", ncpus=40)  # this can take a few minutes
      plot(x, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance", las=1) 
    
      # 6.1. Outlier identification  
        treshold <- 4 / length(x)                          # set a treshold         
        outliers <- data_model$id_r[which(x > treshold)]   # identify outliers
        abline(h = treshold, col = "purple", lwd = 2)      # add an horizontal line
        outlier_rows <- data_model[which(x > treshold), ]  # identify dataset rows with outliers
        
      # 6.2. Outliers removal
        # exclude identified outliers
        data_out <- data_model %>% filter(!id_r %in% outliers) 
        
        # repeat model without outliers
        mod_out <- build_model(data = data_out, mods = mods) 
        summary(mod_out)
        

      # 6.3. Compare analysis with and without outliers
        # extract model results
        coef_out <- data.frame(
          Term = rownames(mod_out$b),           # Extract term names
          Estimate = as.numeric(mod_out$b),     # Coefficients
          CI.Lower = as.numeric(mod_out$ci.lb),  # Lower CI bounds
          CI.Upper = as.numeric(mod_out$ci.ub),  # Upper CI bounds
          pval = as.numeric(mod_out$pval)       # p-values
        )
        
        # back-transform coefficients
        coef_bt_out <- data.frame(
          Term = coef_out$Term,
          Estimate = transf.ztor(coef_out$Estimate), 
          CI.Lower = transf.ztor(coef_out$CI.Lower),
          CI.Upper = transf.ztor(coef_out$CI.Upper),
          pval = coef_out$pval)
        
        # compare back-transformed coefficients with and without outliers
        print(coef_bt)
        print(coef_bt_out)
        
      
      
      
      
      
      
      
      
      
      
      
    
    
