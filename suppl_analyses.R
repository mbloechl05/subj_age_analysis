# ===============================================================
# Supplementary analysis for 
# "Psychological benefits of subjective age bias" 
# ==============================================================


library(car)

# ----------------------------------------------------
# Re-run analysis excluding influential observations
# ----------------------------------------------------

# We choose three strategies to identify influential observations: 
# 1) Cook's distance and 2) df fits. Influential observations are identified as 
# data points that are deemed influential using both strategies. 


# 1.) Identify influential observations
# --------------------------------------

# 1.a) Cooks distance

### fit linear model
model <- lm(ls ~ ca.s + ca.s2 + sa.s + sa.s2 + sa.ca, data = df)

### calculate Cook's distance
cooksd <- cooks.distance(model)

### Which are influential observations?
infl.cd <- as.numeric(names(cooksd)[(cooksd > 100*mean(cooksd, na.rm=T))])
infl.cd # row numbers 287, 2495, 8299

### Plot 
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance") # plot cd
abline(h = 100*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, 
     y=cooksd, 
     labels=ifelse(cooksd > 100*mean(cooksd, na.rm=T), names(cooksd),""), 
     col="red")  # add labels


# 1.b) dffits

### calculate df fits
df_fits <- dffits(model)

### Which are influential observations?
infl.df <- as.numeric(names(df_fits)[(df_fits > 0.5)])
infl.df # 287, 8299

### Plot
plot(df_fits)



# 2.) Exclude influential observations 
# --------------------------------------

### Add variable that contains row numbers to identify influential obs.
df$rn <- as.numeric(rownames(df))

### Which participants are influential obs.?
df[df$rn == 287, ]
df[df$rn == 8299, ]


### Exclude these from data frame
df_rm_io <-df[!(df$rn == 287 | df$rn == 8299),]



# 3.) Re-run analyses
# ----------------------

# Fit main models

nu.fit      <- sem(nu.model,      data = df_rm_io, se = "robust", estimator = "MLR")
pb.fit      <- sem(pb.model,      data = df_rm_io, se = "robust", estimator = "MLR")
omfr.fit    <- sem(omfr.model,    data = df_rm_io, se = "robust", estimator = "MLR")
omrr.fit    <- sem(omrr.model,    data = df_rm_io, se = "robust", estimator = "MLR")
somfr_d.fit <- sem(somfr_d.model, data = df_rm_io, se = "robust", estimator = "MLR")
somfr_u.fit <- sem(somfr_u.model, data = df_rm_io, se = "robust", estimator = "MLR")
somrr_d.fit <- sem(somrr_d.model, data = df_rm_io, se = "robust", estimator = "MLR")
somrr_u.fit <- sem(somrr_u.model, data = df_rm_io, se = "robust", estimator = "MLR")
fu.fit      <- sem(fu.model,      data = df_rm_io, se = "robust", estimator = "MLR")


# Fit supplementary models

s1.fit <- sem(s1.model, data = df_rm_io, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df_rm_io, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df_rm_io, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df_rm_io, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df_rm_io, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df_rm_io, se = "robust", estimator = "MLR")
s7.fit <- sem(s7.model, data = df_rm_io, se = "robust", estimator = "MLR")


# Evaluations full model set

### Create list of all models
all.models <- list(fu.fit, 
                   pb.fit, 
                   omfr.fit, 
                   omrr.fit, 
                   somfr_d.fit, 
                   somrr_d.fit, 
                   s1.fit, 
                   s2.fit, 
                   s3.fit, 
                   s4.fit, 
                   s5.fit,   
                   s6.fit, 
                   s7.fit, 
                   nu.fit)

### Set names of models in list
names(all.models) <- c("full", "pb", "omfr", "omrr", "somfr", "somrr", "s1", 
                       "s2", "s3", "s4", "s5", "s6", "s7", "null")

### Get AICs and other information
aics.3 <- as.data.frame(AICcmodavg::aictab(all.models, 
                                           modnames = names(all.models), 
                                           second.ord = F))

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
frm(aics.3)
### Full model, s1, s3, s5, and s7 are redundant

### Exclude redundant models from list of models
all.models.r <- all.models[-c(1, 11, 9, 13, 7)]

### Get AICs and other fit information for final model set
aics.3 <- as.data.frame(AICcmodavg::aictab(all.models.r, 
                                           modnames = names(all.models.r), 
                                           second.ord = F))
aics.3 # print results


# # -------------------------------------------------------------
# # Re-un analyses excluding outliers and people aged < 40 years
# # -------------------------------------------------------------
# 
# # mark people with age < 40 y as missing
# df$ca[df$ca < 40] <- NA
# 
# # 3. Exclude cases with missing data (i.e. outliers)
# df <- na.omit(df)
# 
# # Get descriptive stats after exclusion
# describe(df$ca)
# describe(df$sa)
# 
# # Standardise predictors again
# grandmean <- mean(c(df$sa, df$ca), na.rm = T)
# pooledsd  <- sqrt(((sd(df$sa)^2)+(sd(df$ca)^2))/2)
# 
# df$sa.s  <- (df$sa-grandmean)/pooledsd
# df$ca.s  <- (df$ca-grandmean)/pooledsd
# 
# # Add squared and interaction terms of the predictors again
# df$sa.s2 <- df$sa.s^2
# df$ca.s2 <- df$ca.s^2
# df$sa.ca <- df$sa.s*df$ca.s
# 
# # Refit models
# nu.fit <- sem(nu.model, data = df, se = "robust", estimator = "MLR")
# fu.fit <- sem(fu.model, data = df, se = "robust", estimator = "MLR")
# pb.fit <- sem(pb.model, data = df, se = "robust", estimator = "MLR")
# om.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# # s1.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# # s2.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# # s3.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# # s4.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# # s5.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# # s6.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# # s7.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# 
# # Get AICs of all models in the full model set
# aic.results <- 
#   AIC(nu.fit, fu.fit, 
#       pb.fit, om.fit, som.fit, 
#       s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, s7.fit)
# 
# # Order results according to AIC (models with low AIC first)
# aic.results <- aic.results[order(aic.results$AIC),]
# 
# # Add rownames as column to data frame
# aic.results <- add_rownames(aic.results, var = "rowname")
# 
# # Make a vector from AICs 
# aics.vector <- pull(aic.results, AIC)
# 
# # Calculate Akaike weights and attach to results of model comparisons
# weights <- akaike.weights(aics.vector)
# aic.results$weights <- weights$weights



# -------------------------------------------------------------
# Re-run analysis using FIML (instead of list-wise deletion)
# -------------------------------------------------------------

# Re-do data preparation
# ------------------------



# Re-fit models using full data and FIML
# -----------------------------------------

# select relevant variables
df.full <- df.raw[,c("idauniq", "DhSex", "scold", "dhager", "sclifea", "sclifeb", 
                     "sclifec", "sclifed", "sclifee")]

# 1.) Rename ageing variables
df.full$sa <- df.full$scold # subjective age
df.full$ca <- df.full$dhager # real age


# 2.) Replace missing values with NAs

### Recode values for missing data (-9, -1) in whole dataset as NA
df.full <- df.full %>% mutate_all(funs(na_if(., -9)))
df.full <- df.full %>% mutate_all(funs(na_if(., -1)))

### Exclude old age values (99)
df.full.2 <- df.full
df.full.2$exclude <- "no"

df.full.2$exclude[df.full.2$ca == 99] <- "yes"
df.full.2 <- df.full.2 %>% filter(exclude == "no")

# 3.) Recode sex variable
df.full.2$sex[df.full.2$DhSex == 2] <- "female"
df.full.2$sex[df.full.2$DhSex == 1] <- "male"


# 4.) Prepare outcome variable life satisfaction

### Recode life satisfaction (ls) items so that higher values indicate higher ls
df.full.2$sclifea_r <- 8-df.full.2$sclifea
df.full.2$sclifeb_r <- 8-df.full.2$sclifeb
df.full.2$sclifec_r <- 8-df.full.2$sclifec
df.full.2$sclifed_r <- 8-df.full.2$sclifed
df.full.2$sclifee_r <- 8-df.full.2$sclifee

### Calculate mean score on life satistfaction scale
df.full.2$ls <- rowMeans(df.full.2[,c("sclifea_r",  "sclifeb_r",  "sclifec_r", 
                                  "sclifed_r", "sclifee_r")], na.rm = T)



# 5.) Remove outliers chronological age

### Calculate +/- 3 SD
out <- 3*sd(df.full.2$ca, na.rm = T)

### How many outliers
sum(df.full.2$ca < mean(df.full.2$ca, na.rm = T) - out, na.rm = T) # -3 SD
sum(df.full.2$ca > mean(df.full.2$ca, na.rm = T) + out, na.rm = T) # +3 SD

### Mark people with chronological age +/- 3 SD as missing
df.full.2$ca[df.full.2$ca > mean(df.full.2$ca, na.rm = T) + out | 
             df.full.2$ca < mean(df.full.2$ca, na.rm = T) - out] <- 999


# 6.) Remove outliers subjective age

### Calculate +/- 3 SD
out <- 3*sd(df.full.2$sa, na.rm = T)

### How many outliers 
sum(df.full.2$sa < mean(df.full.2$sa, na.rm = T) - out, na.rm = T) # - 3 SD
sum(df.full.2$sa > mean(df.full.2$sa, na.rm = T) + out, na.rm = T) # + 3 SD


### Mark people with chronological age +/- 3 SD as missing
df.full.2$sa[df.full.2$sa > mean(df.full.2$sa, na.rm = T) + out | 
             df.full.2$sa < mean(df.full.2$sa, na.rm = T) - out] <- 999


### Exclude these from data frame

df.full.2$exclude[df.full.2$sa == 999] <- "yes"
df.full.2$exclude[df.full.2$ca == 999] <- "yes"

### Exclude these from data frame
df.full.2 <- df.full.2 %>% filter(exclude == "no")


# 1.) Standardise predictors
grandmean <- mean(c(df.full.2$sa, df.full.2$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df.full.2$sa,  na.rm = T)^2)+(sd(df.full.2$ca,  na.rm = T)^2))/2)

df.full.2$sa.s  <- (df.full.2$sa-grandmean)/pooledsd
df.full.2$ca.s  <- (df.full.2$ca-grandmean)/pooledsd


# 2.) Add squared and interaction terms of the predictors
df.full.2$sa.s2 <- df.full.2$sa.s^2
df.full.2$ca.s2 <- df.full.2$ca.s^2
df.full.2$sa.ca <- df.full.2$sa.s*df.full.2$ca.s


# Fit main models

nu.fit      <- sem(nu.model,      data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
pb.fit      <- sem(pb.model,      data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
omfr.fit    <- sem(omfr.model,    data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
omrr.fit    <- sem(omrr.model,    data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
somfr_d.fit <- sem(somfr_d.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
somfr_u.fit <- sem(somfr_u.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
somrr_d.fit <- sem(somrr_d.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
somrr_u.fit <- sem(somrr_u.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
fu.fit      <- sem(fu.model,      data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)


# Fit supplementary models

s1.fit <- sem(s1.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
s2.fit <- sem(s2.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
s3.fit <- sem(s3.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
s4.fit <- sem(s4.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
s5.fit <- sem(s5.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
s6.fit <- sem(s6.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)
s7.fit <- sem(s7.model, data = df.full.2, se = "robust", missing = "FIML", fixed.x = T)


# Evaluations full model set

### Create list of all models
all.models <- list(fu.fit, 
                   pb.fit, 
                   omfr.fit, 
                   omrr.fit, 
                   somfr_d.fit, 
                   somrr_d.fit, 
                   s1.fit, 
                   s2.fit, 
                   s3.fit, 
                   s4.fit, 
                   s5.fit,   
                   s6.fit, 
                   s7.fit, 
                   nu.fit)

### Set names of models in list
names(all.models) <- c("full", "pb", "omfr", "omrr", "somfr", "somrr", "s1", 
                       "s2", "s3", "s4", "s5", "s6", "s7", "null")

### Get AICs and other information
aics.4 <- as.data.frame(AICcmodavg::aictab(all.models, 
                                           modnames = names(all.models), 
                                           second.ord = F))

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
frm(aics.4)
### Full model, s1, s3, s5, and s7 are redundant

### Exclude redundant models from list of models
all.models.r <- all.models[-c(1, 7, 9, 11, 13)]

### Get AICs and other fit information for final model set
aics.4 <- as.data.frame(AICcmodavg::aictab(all.models.r, 
                                           modnames = names(all.models.r), 
                                           second.ord = F))
aics.4 # print results


# Re-fit full model using RSA function
r.fu  <- RSA(ls  ~ sa.s*ca.s, df.full)

# modified seaborne
plot(r.fu, model = "full", axes = c("LOC", "LOIC", "PA1"), 
     project = c("LOC", "LOIC", "PA1"))

# # --------------------------------------------------
# # Re-run analyses with full sample (pre-registered)
# # --------------------------------------------------
# 
# ### 
# ### CHECK AGAIN IF THIS IS RIGHT!! 
# ###
# 
# # Re-name full data frame to df
# df <- df.full
# 
# # mark people with age < 40 y as missing
# df$ca[df$ca < 40] <- NA
# 
# # 3. Exclude cases with missing data (i.e. outliers)
# df <- na.omit(df)
# 
# # Get descriptive stats after exclusion
# describe(df$ca)
# describe(df$sa)
# 
# # Standardise predictors again
# grandmean <- mean(c(df$sa, df$ca), na.rm = T)
# pooledsd  <- sqrt(((sd(df$sa)^2)+(sd(df$ca)^2))/2)
# 
# df$sa.s  <- (df$sa-grandmean)/pooledsd
# df$ca.s  <- (df$ca-grandmean)/pooledsd
# 
# # Add squared and interaction terms of the predictors again
# df$sa.s2 <- df$sa.s^2
# df$ca.s2 <- df$ca.s^2
# df$sa.ca <- df$sa.s*df$ca.s
# 
# # Refit models
# nu.fit <- sem(nu.model, data = df, se = "robust", estimator = "MLR")
# fu.fit <- sem(fu.model, data = df, se = "robust", estimator = "MLR")
# pb.fit <- sem(pb.model, data = df, se = "robust", estimator = "MLR")
# om.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s1.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s2.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s3.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s4.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s5.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s6.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s7.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# 
# # Get AICs of all models in the full model set
# aic.results <- 
#   AIC(nu.fit, fu.fit, 
#       pb.fit, om.fit, som.fit, 
#       s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, s7.fit)
# 
# # Order results according to AIC (models with low AIC first)
# aic.results <- aic.results[order(aic.results$AIC),]
# 
# # Add rownames as column to data frame
# aic.results <- add_rownames(aic.results, var = "rowname")
# 
# # Make a vector from AICs 
# aics.vector <- pull(aic.results, AIC)
# 
# # Calculate Akaike weights and attach to results of model comparisons
# weights <- akaike.weights(aics.vector)
# aic.results$weights <- weights$weights


# ----------------------------------------------
# Rescale 1 - center predictors to their mean
# ---------------------------------------------
r.fu  <- RSA(ls  ~ sa*ca, df)
getPar(r.fu, "coef", model="SRR")

### Center predictors to their respective mean (instead of grandmean)
df$sa.c  <- scale(df$sa, center = T, scale = F)
df$ca.c  <- scale(df$ca, center = T, scale = F)

r.fu  <- RSA(ls  ~ sa.c*ca.c, df)
plot(r.fu, model = "full")
getPar(r.fu, "coef", model="SRR")


# --------------------------------------------------
# Rescale 2 - center predictors to their grandmean
# --------------------------------------------------

# 
grandmean <- mean(c(df$sa, df$ca), na.rm = T)

df$sa.c  <- (df$sa-grandmean)
df$ca.c  <- (df$ca-grandmean)

r.fu  <- RSA(ls  ~ sa.c*ca.c, df)
plot(r.fu, model = "full", project = c("PA1", "LOC"))
plot(r.fu, model = "full", type = "contour", project = c("PA1", "LOC"))
getPar(r.fu, "coef", model="SRR")


#### Values of C change due to granmean centering 
#### The coefficient is about 17y years if predictors are centered to their 
#### respective mean. The value for C changes to 26y if the values are 
#### grandmean centered

