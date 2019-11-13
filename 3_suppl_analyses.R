# ===============================================================
# Supplementary analysis for 
# "Psychological benefits of subjective age bias" 
# ==============================================================


# -----------------------------------------------------
# 1) Re-run analyses with full sample (pre-registered)
# ------------------------------------------------------

# Re-name full data frame to df
df <- df.eligible

# Standardise predictors again
grandmean <- mean(c(df$sa, df$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df$sa)^2)+(sd(df$ca)^2))/2)

df$sa.s  <- (df$sa-grandmean)/pooledsd
df$ca.s  <- (df$ca-grandmean)/pooledsd

# Add squared and interaction terms of the predictors again
df$sa.s2 <- df$sa.s^2
df$ca.s2 <- df$ca.s^2
df$sa.ca <- df$sa.s*df$ca.s

# Refit models
nu.fit      <- sem(nu.model,      data = df, se = "robust", estimator = "MLR")
pb.fit      <- sem(pb.model,      data = df, se = "robust", estimator = "MLR")
omfr.fit    <- sem(omfr.model,    data = df, se = "robust", estimator = "MLR")
omrr.fit    <- sem(omrr.model,    data = df, se = "robust", estimator = "MLR")
somfr_d.fit <- sem(somfr_d.model, data = df, se = "robust", estimator = "MLR")
somfr_u.fit <- sem(somfr_u.model, data = df, se = "robust", estimator = "MLR")
somrr_d.fit <- sem(somrr_d.model, data = df, se = "robust", estimator = "MLR")
somrr_u.fit <- sem(somrr_u.model, data = df, se = "robust", estimator = "MLR")
fu.fit      <- sem(fu.model,      data = df, se = "robust", estimator = "MLR")
s1.fit <- sem(s1.model, data = df, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df, se = "robust", estimator = "MLR")
s7.fit <- sem(s7.model, data = df, se = "robust", estimator = "MLR")

### Create list of all models
models <- list(fu.fit, pb.fit, omfr.fit,  omrr.fit, somfr_d.fit, somrr_d.fit, 
               s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, s7.fit, nu.fit)

### Set names of models in list
names(models) <- c("full", "pb","omfr", "omrr", "somfr", "somrr", "s1", "s2", 
                   "s3", "s4", "s5", "s6", "s7", "null")


### Get AICs and other information
aics <- as.data.frame(
  AICcmodavg::aictab(models, modnames = names(models), second.ord = F)
)

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
aics
frm(aics)


### Full model, omrr, s1, s3, s5, and s7 are redundant --> exclude from list of models
models_r <- models[-c( 1, 4, 11, 9, 13, 7)]

### Get AICs and other fit information for final model set
aics <- as.data.frame(
  AICcmodavg::aictab(models_r, modnames = names(models_r), second.ord = F)
)

aics 


# ----------------------------------------------------
# 2) Re-run analysis excluding influential observations
# ----------------------------------------------------

# Data preparation as in main analyses

df <- df.eligible

### Trim variable chronological age
out <- 3*sd(df$ca, na.rm = T)

df$ca[df$ca > mean(df$ca, na.rm = T) + out | 
        df$ca < mean(df$ca, na.rm = T) - out] <- NA

### Trim variable subjective age
out <- 3*sd(df$sa, na.rm = T)

df$sa[df$sa > mean(df$sa, na.rm = T) + out | 
        df$sa < mean(df$sa, na.rm = T) - out] <- NA

### Exclude cases 
df <- na.omit(df)

### Standardise predictors again
grandmean <- mean(c(df$sa, df$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df$sa)^2)+(sd(df$ca)^2))/2)

df$sa.s  <- (df$sa-grandmean)/pooledsd
df$ca.s  <- (df$ca-grandmean)/pooledsd

### Add squared and interaction terms of the predictors again
df$sa.s2 <- df$sa.s^2
df$ca.s2 <- df$ca.s^2
df$sa.ca <- df$sa.s*df$ca.s

# We choose three strategies to identify influential observations: 
# 1) Cook's distance and 2) df fits. Influential observations are identified as 
# data points that are deemed influential using both strategies. 


# 2.1) Identify influential observations
# --------------------------------------

# 2.1.1) Cooks distance

### fit linear model
model <- lm(ls ~ ca.s + ca.s2 + sa.s + sa.s2 + sa.ca, data = df)

### calculate Cook's distance
cooksd <- cooks.distance(model)

### Which are influential observations?
infl.cd <- as.numeric(names(cooksd)[(cooksd > 100*mean(cooksd, na.rm=T))])
infl.cd # row numbers 2495, 8299

### Plot 
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance") # plot cd
abline(h = 100*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, 
     y=cooksd, 
     labels=ifelse(cooksd > 100*mean(cooksd, na.rm=T), names(cooksd),""), 
     col="red")  # add labels


# 2.1.2) dffits

### calculate df fits
df_fits <- dffits(model)

### Which are influential observations?
infl.df <- as.numeric(names(df_fits)[(df_fits > 0.05)])
infl.df # 2495, 8299

### Plot
plot(df_fits)



# 2.2) Exclude influential observations 
# --------------------------------------

### Add variable that contains row numbers to identify influential obs.
df$rn <- as.numeric(rownames(df))

### Which participants are influential obs.?
df[df$rn == 2495, ]
df[df$rn == 8299, ]

### Exclude these from data frame
df_rm_io <-df[!(df$rn == 2495 | df$rn == 8299),]


# 2.3) Re-run analyses
# ----------------------

# Refit models
nu.fit      <- sem(nu.model,      data = df_rm_io, se = "robust", estimator = "MLR")
pb.fit      <- sem(pb.model,      data = df_rm_io, se = "robust", estimator = "MLR")
omfr.fit    <- sem(omfr.model,    data = df_rm_io, se = "robust", estimator = "MLR")
omrr.fit    <- sem(omrr.model,    data = df_rm_io, se = "robust", estimator = "MLR")
somfr_d.fit <- sem(somfr_d.model, data = df_rm_io, se = "robust", estimator = "MLR")
somfr_u.fit <- sem(somfr_u.model, data = df_rm_io, se = "robust", estimator = "MLR")
somrr_d.fit <- sem(somrr_d.model, data = df_rm_io, se = "robust", estimator = "MLR")
somrr_u.fit <- sem(somrr_u.model, data = df_rm_io, se = "robust", estimator = "MLR")
fu.fit      <- sem(fu.model,      data = df_rm_io, se = "robust", estimator = "MLR")
s1.fit <- sem(s1.model, data = df_rm_io, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df_rm_io, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df_rm_io, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df_rm_io, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df_rm_io, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df_rm_io, se = "robust", estimator = "MLR")
s7.fit <- sem(s7.model, data = df_rm_io, se = "robust", estimator = "MLR")

### Create list of all models
models <- list(fu.fit, pb.fit, omfr.fit,  omrr.fit, somfr_d.fit, somrr_d.fit, 
               s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, s7.fit, nu.fit)

### Set names of models in list
names(models) <- c("full", "pb","omfr", "omrr", "somfr", "somrr", "s1", "s2", 
                   "s3", "s4", "s5", "s6", "s7", "null")

### Get AICs and other information
aics <- as.data.frame(
  AICcmodavg::aictab(models, modnames = names(models), second.ord = F)
)

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
aics
frm(aics)

### Full model, s1, s3, s5, and s7 are redundant --> exclude from list of models
models_r <- models[-c( 1, 11, 9, 13, 7)]

### Get AICs and other fit information for final model set
aics <- as.data.frame(
  AICcmodavg::aictab(models_r, modnames = names(models_r), second.ord = F)
)

aics 


# ------------------------------------------------------------------------------
# 3) Re-un analyses with trimmed variables and excluding people aged < 40 years
# ------------------------------------------------------------------------------

# 3.1) Exclude people with +/- 3 SD and < 40 years

df <- df.eligible

### Trim variable chronological age
out <- 3*sd(df$ca, na.rm = T)

df$ca[df$ca > mean(df$ca, na.rm = T) + out | 
        df$ca < mean(df$ca, na.rm = T) - out] <- NA


### Trim variable subjective age
out <- 3*sd(df$sa, na.rm = T)

df$sa[df$sa > mean(df$sa, na.rm = T) + out | 
        df$sa < mean(df$sa, na.rm = T) - out] <- NA

### Mark people with age < 40 y. as missing
df$ca[df$ca < 40] <- NA

### Exclude cases 
df <- na.omit(df)


# 3.2) Re-run analyses 

### Standardise predictors again
grandmean <- mean(c(df$sa, df$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df$sa)^2)+(sd(df$ca)^2))/2)

df$sa.s  <- (df$sa-grandmean)/pooledsd
df$ca.s  <- (df$ca-grandmean)/pooledsd

### Add squared and interaction terms of the predictors again
df$sa.s2 <- df$sa.s^2
df$ca.s2 <- df$ca.s^2
df$sa.ca <- df$sa.s*df$ca.s

### Re-fit models
nu.fit      <- sem(nu.model,      data = df, se = "robust", estimator = "MLR")
pb.fit      <- sem(pb.model,      data = df, se = "robust", estimator = "MLR")
omfr.fit    <- sem(omfr.model,    data = df, se = "robust", estimator = "MLR")
omrr.fit    <- sem(omrr.model,    data = df, se = "robust", estimator = "MLR")
somfr_d.fit <- sem(somfr_d.model, data = df, se = "robust", estimator = "MLR")
somfr_u.fit <- sem(somfr_u.model, data = df, se = "robust", estimator = "MLR")
somrr_d.fit <- sem(somrr_d.model, data = df, se = "robust", estimator = "MLR")
somrr_u.fit <- sem(somrr_u.model, data = df, se = "robust", estimator = "MLR")
fu.fit      <- sem(fu.model,      data = df, se = "robust", estimator = "MLR")
s1.fit <- sem(s1.model, data = df, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df, se = "robust", estimator = "MLR")
s7.fit <- sem(s7.model, data = df, se = "robust", estimator = "MLR")

### Create list of all models
models <- list(fu.fit, pb.fit, omfr.fit,  omrr.fit, somfr_d.fit, somrr_d.fit, 
               s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, s7.fit, nu.fit)

### Set names of models in list
names(models) <- c("full", "pb","omfr", "omrr", "somfr", "somrr", "s1", "s2", 
                   "s3", "s4", "s5", "s6", "s7", "null")

### Get AICs and other information
aics <- as.data.frame(
  AICcmodavg::aictab(models, modnames = names(models), second.ord = F)
)

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
aics
frm(aics)

### Full model, s1, s3, s5, and s7 are redundant --> exclude from list of models
models_r <- models[-c( 1, 11, 9, 13, 7)]

### Get AICs and other fit information for final model set
aics <- as.data.frame(
  AICcmodavg::aictab(models_r, modnames = names(models_r), second.ord = F)
)

aics 


# -------------------------------------------------------------
# 4) Re-run analysis using FIML (instead of list-wise deletion)
# -------------------------------------------------------------

# 4.1) Re-do data preparation

df <- df.full

### Trim variable chronological age
out <- 3*sd(df$ca, na.rm = T)

df$ca[df$ca > mean(df$ca, na.rm = T) + out | 
        df$ca < mean(df$ca, na.rm = T) - out] <- -99


#### Trim variable subjective age
out <- 3*sd(df$sa, na.rm = T)

df$sa[df$sa > mean(df$sa, na.rm = T) + out | 
        df$sa < mean(df$sa, na.rm = T) - out] <- -99

### Exclude cases 
df <- subset(df, sa != -99 | is.na(sa))
df <- subset(df, ca != -99)

### Standardise predictors
grandmean <- mean(c(df$sa, df$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df$sa, na.rm = T)^2)+(sd(df$ca, na.rm = T)^2))/2)

df$sa.s  <- (df$sa-grandmean)/pooledsd
df$ca.s  <- (df$ca-grandmean)/pooledsd

### Add squared and interaction terms of the predictors
df$sa.s2 <- df$sa.s^2
df$ca.s2 <- df$ca.s^2
df$sa.ca <- df$sa.s*df$ca.s

# 4.2) Re-run analyses 

### Refit models
nu.fit      <- sem(nu.model,      data = df, se = "robust", missing = "fiml", fixed.x = T)
pb.fit      <- sem(pb.model,      data = df, se = "robust", missing = "fiml", fixed.x = T)
omfr.fit    <- sem(omfr.model,    data = df, se = "robust", missing = "fiml", fixed.x = T)
omrr.fit    <- sem(omrr.model,    data = df, se = "robust", missing = "fiml", fixed.x = T)
somfr_d.fit <- sem(somfr_d.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
somfr_u.fit <- sem(somfr_u.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
somrr_d.fit <- sem(somrr_d.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
somrr_u.fit <- sem(somrr_u.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
fu.fit      <- sem(fu.model,      data = df, se = "robust", missing = "fiml", fixed.x = T)
s1.fit <- sem(s1.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s2.fit <- sem(s2.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s3.fit <- sem(s3.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s4.fit <- sem(s4.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s5.fit <- sem(s5.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s6.fit <- sem(s6.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s7.fit <- sem(s7.model, data = df, se = "robust", missing = "fiml", fixed.x = T)

### Create list of all models
models <- list(fu.fit, pb.fit, omfr.fit,  omrr.fit, somfr_d.fit, somrr_d.fit, 
               s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, s7.fit, nu.fit)

### Set names of models in list
names(models) <- c("full", "pb","omfr", "omrr", "somfr", "somrr", "s1", "s2", 
                   "s3", "s4", "s5", "s6", "s7", "null")

### Get AICs and other information
aics <- as.data.frame(
  AICcmodavg::aictab(models, modnames = names(models), second.ord = F)
)

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
aics
frm(aics)

### Full model, s1, s3, s5, and s7 are redundant --> exclude from list of models
models_r <- models[-c( 1, 11, 9, 13, 7)]

### Get AICs and other fit information for final model set
aics <- as.data.frame(
  AICcmodavg::aictab(models_r, modnames = names(models_r), second.ord = F)
)

aics 






