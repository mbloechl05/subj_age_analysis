# ===============================================================
# Subjective age bias and life satisfaction
# Script 3: Supplementary analyses
# ==============================================================

# This script contains code to run all additional or supplementary 
# analayses. Script 1 and script 2 have to be run before running 
# this script since they contain model definitions and data, which 
# which are are called upon here. 

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
nu.fit <- sem(nu.model, data = df, se = "robust", estimator = "MLR")
fu.fit <- sem(fu.model, data = df, se = "robust", estimator = "MLR")
a1.fit <- sem(a1.model, data = df, se = "robust", estimator = "MLR")
b1.fit <- sem(b1.model, data = df, se = "robust", estimator = "MLR")
a2.fit <- sem(a2.model, data = df, se = "robust", estimator = "MLR")
b2.fit <- sem(b2.model, data = df, se = "robust", estimator = "MLR")
a3.fit <- sem(a3.model, data = df, se = "robust", estimator = "MLR")
b3.fit <- sem(b3.model, data = df, se = "robust", estimator = "MLR")
s1.fit <- sem(s1.model, data = df, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df, se = "robust", estimator = "MLR")

### Create list of all models
models <- list(fu.fit, a1.fit, b1.fit, a2.fit, b2.fit, a3.fit, b3.fit, 
               s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, nu.fit)

### Set names of models in list
names(models) <- c("full", "model-1a", "model-1b", "model-2a", "model-2b", 
                   "model-3a", "model-3b", "model-s1", "model-s2", "model-s3", 
                   "model-s4", "model-s5", "model-s6", "null")

### Get AICs and other information
aics <- as.data.frame(
  AICcmodavg::aictab(models, modnames = names(models), second.ord = F)
)

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
aics
frm(aics)

### full, model-2b, model-3a, model-s5, model-s3, model-s6, model-s1 are redundant 
### --> exclude from list of models
models_r <- models[-c( 1, 5, 6, 12, 10, 13, 8)]

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
nu.fit <- sem(nu.model, data = df_rm_io, se = "robust", estimator = "MLR")
fu.fit <- sem(fu.model, data = df_rm_io, se = "robust", estimator = "MLR")
a1.fit <- sem(a1.model, data = df_rm_io, se = "robust", estimator = "MLR")
b1.fit <- sem(b1.model, data = df_rm_io, se = "robust", estimator = "MLR")
a2.fit <- sem(a2.model, data = df_rm_io, se = "robust", estimator = "MLR")
b2.fit <- sem(b2.model, data = df_rm_io, se = "robust", estimator = "MLR")
a3.fit <- sem(a3.model, data = df_rm_io, se = "robust", estimator = "MLR")
b3.fit <- sem(b3.model, data = df_rm_io, se = "robust", estimator = "MLR")
s1.fit <- sem(s1.model, data = df_rm_io, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df_rm_io, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df_rm_io, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df_rm_io, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df_rm_io, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df_rm_io, se = "robust", estimator = "MLR")

### Create list of all models
models <- list(fu.fit, a1.fit, b1.fit, a2.fit, b2.fit, a3.fit, b3.fit, 
               s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, nu.fit)

### Set names of models in list
names(models) <- c("full", "model-1a", "model-1b", "model-2a", "model-2b", 
                   "model-3a", "model-3b", "model-s1", "model-s2", "model-s3", 
                   "model-s4", "model-s5", "model-s6", "null")

### Get AICs and other information
aics <- as.data.frame(
  AICcmodavg::aictab(models, modnames = names(models), second.ord = F)
)

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
aics
frm(aics)

### full, model-3a, model-s5, model-s3, model-s6, model-s1 are redundant --> 
### exclude from list of models
models_r <- models[-c( 1, 6, 12, 10, 13, 8)]

### Get AICs and other fit information for final model set
aics <- as.data.frame(
  AICcmodavg::aictab(models_r, modnames = names(models_r), second.ord = F)
)

aics 



# -------------------------------------------------------------
# 3) Re-run analysis using FIML (instead of list-wise deletion)
# -------------------------------------------------------------

# 3.1) Re-do data preparation

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

# 3.2) Re-run analyses 

# Refit models
nu.fit <- sem(nu.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
fu.fit <- sem(fu.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
a1.fit <- sem(a1.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
b1.fit <- sem(b1.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
a2.fit <- sem(a2.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
b2.fit <- sem(b2.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
a3.fit <- sem(a3.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
b3.fit <- sem(b3.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s1.fit <- sem(s1.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s2.fit <- sem(s2.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s3.fit <- sem(s3.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s4.fit <- sem(s4.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s5.fit <- sem(s5.model, data = df, se = "robust", missing = "fiml", fixed.x = T)
s6.fit <- sem(s6.model, data = df, se = "robust", missing = "fiml", fixed.x = T)

### Create list of all models
models <- list(fu.fit, a1.fit, b1.fit, a2.fit, b2.fit, a3.fit, b3.fit, 
               s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, nu.fit)

### Set names of models in list
names(models) <- c("full", "model-1a", "model-1b", "model-2a", "model-2b", 
                   "model-3a", "model-3b", "model-s1", "model-s2", "model-s3", 
                   "model-s4", "model-s5", "model-s6", "null")

### Get AICs and other information
aics <- as.data.frame(
  AICcmodavg::aictab(models, modnames = names(models), second.ord = F)
)

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
aics
frm(aics)

### full, model-3a, model-s5, model-s3, model-s6, model-s1 are redundant --> 
### exclude from list of models
models_r <- models[-c( 1, 6, 12, 10, 13, 8)]

### Get AICs and other fit information for final model set
aics <- as.data.frame(
  AICcmodavg::aictab(models_r, modnames = names(models_r), second.ord = F)
)

aics 


# -------------------------------
# 4) Explore sex differences
# -------------------------------

# 4.1.) Fit models with and without equivalence constrains and compare
# ----------------------------------------------------------------------

# In the following, we take the best-fitting model (3.b) and fit this model model 
# without constraints (i.e. for both sexes, parameters are estimated freely implicating 
# sex differences), then  with constraints (i.e. for both sexes, parameters are constrained
# to be equal, implicating no sex differences) and then compare both models. 

# Model without constraints

b3.model.uc <- "
    ls ~ 1 + c(b1m,b1f)*sa.s + c(b2m,b2f)*ca.s + c(b3m,b3f)*sa.s2 + c(b4m,b4f)*sa.ca + c(b5m,b5f)*ca.s2
         
    # parameter constraints male
    b3m   <  -0.000001
    b5m   <  -0.000001
    b4m^2 == 4*b3m*b5m
    
    # parameter constraints female:
    b3f   <  -0.000001
    b5f   <  -0.000001
    b4f^2 == 4*b3f*b5f

    # parameters male
    X0m  := (b2m*b4m - 2*b1m*b5m) / (4*b3m*b5m - b4m^2)
    Y0m  := (b1m*b4m - 2*b2m*b3m) / (4*b3m*b5m - b4m^2)
    p11m := (b5m - b3m + sqrt(((b3m - b5m)^2) + (b4m^2))) / b4m
    p10m := Y0m - p11m*X0m

    # parameters female
    X0f  := (b2f*b4f - 2*b1f*b5f) / (4*b3f*b5f - b4f^2)
    Y0f  := (b1f*b4f - 2*b2f*b3f) / (4*b3f*b5f - b4f^2)
    p11f := (b5f - b3f + sqrt(((b3f - b5f)^2) + (b4f^2))) / b4f
    p10f := Y0f - p11f*X0f"

fit.uc <- sem(b3.model.uc, data = df, group = "sex", se = "robust", estimator = "MLR")
summary(fit.uc)


# Model with constraints

b3.model.c <- "
    ls ~ 1 + c(b1,b1)*sa.s + c(b2,b2)*ca.s + c(b3,b3)*sa.s2 + c(b4,b4)*sa.ca + c(b5,b5)*ca.s2
    
    # parameter constraints (equal for males and females)
    b3   <  -0.000001
    b5   <  -0.000001
    b4^2 == 4*b3*b5
    "
fit.c <- sem(b3.model.c, data = df, group = "sex", se = "robust", estimator = "MLR")
summary(fit.c)


# Model comparison

anova(fit.uc, fit.c)


# 4.2.) Optimal margins for men and women
# -----------------------------------------

# Prepare data
ca <- c(40,50,60,70,80,90) # Select chronological ages and put them in a data frame
ca.male <- data.frame(ca)

ca.male$ca.sd <- (ca.male$ca - grandmean)/pooledsd # Standardised chron. age

# Extract PA1 parameters from model fit
p10 <- fit.uc@ParTable$est[64] # p10 
p11 <- fit.uc@ParTable$est[63] # p11

# Run calculations 
ca.male$sa.sd <- ((ca.male$ca.sd - p10)/p11) # Using PA1, calculate corresponding standardised subj. ages 
ca.male$sa    <- grandmean + (ca.male$sa.sd*pooledsd) # Unstandard. subjective ages
ca.male$bi    <- ca.male$ca-ca.male$sa # Calculate difference (i.e. optimal subj. age bias)
ca.male$bi.2  <- ca.male$sa-ca.male$ca # negative values indicate feeling younger
ca.male$bi.p  <- (ca.male$sa-ca.male$ca)/ca.male$ca # proportinal values (divided by age)

# Show results
ca.male


# Prepare data
ca <- c(40,50,60,70,80,90) # Select chronological ages and put them in a data frame
ca.female <- data.frame(ca)

ca.female$ca.sd <- (ca.female$ca - grandmean)/pooledsd # Standardised chron. age

# Extract PA1 parameters from model fit
p10 <- fit.uc@ParTable$est[68] # p10
p11 <- fit.uc@ParTable$est[67] # p11

# Run calculations 
ca.female$sa.sd <- ((ca.female$ca.sd - p10)/p11) # Using PA1, calculate corresponding standardised subj. ages 
ca.female$sa    <- grandmean + (ca.female$sa.sd*pooledsd) # Unstandard. subjective ages
ca.female$bi    <- ca.female$ca-ca.female$sa # Calculate difference (i.e. optimal subj. age bias)
ca.female$bi.2  <- ca.female$sa-ca.female$ca # negative values indicate feeling younger
ca.female$bi.p  <- (ca.female$sa-ca.female$ca)/ca.female$ca # proportinal values (divided by age)

# Show results
ca.female


# ------------------------------------------------------------
# 5) "Classical" difference score model and model comparison
# ------------------------------------------------------------

# Model S1B: Low subj. age and high chronological age model

trad.diff.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == -b2
b3 == 0
b4 == 0
b5 == 0"

fit.trad.diff <- sem(trad.diff.model, data = df, se = "robust", estimator = "MLR")
summary(fit.trad.diff)

fitmeasures(fit.trad.diff, fit.measures = c("AIC", "CFI", "SRMR", "RMSEA"))
inspect(fit.trad.diff, 'r2')

