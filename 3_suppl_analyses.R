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
# 4) Explore gender differences
# -------------------------------

# 4.1.) Model for males 
# -----------------------

# 4.1.1) Prepare data

# Create male data set
df.male <- subset(df.eligible, DhSex == 1)

# Trim variable chronological age
out <- 3*sd(df.male$ca, na.rm = T) # Calculate +/- 3 SD

sum(df.male$ca < mean(df.male$ca, na.rm = T) - out, na.rm = T) # -3 SD: number of excluded people
sum(df.male$ca > mean(df.male$ca, na.rm = T) + out, na.rm = T) # +3 SD: number of excluded people

df.male$ca[df.male$ca > mean(df.male$ca, na.rm = T) + out | 
             df.male$ca < mean(df.male$ca, na.rm = T) - out] <- NA # Mark people as NA

# Trim variable subjective age
out <- 3*sd(df.male$sa, na.rm = T) # Calculate +/- 3 SD

sum(df.male$sa < mean(df.male$sa, na.rm = T) - out, na.rm = T) # - 3 SD: number of excluded people
sum(df.male$sa > mean(df.male$sa, na.rm = T) + out, na.rm = T) # + 3 SD: number of excluded people

df.male$sa[df.male$sa > mean(df.male$sa, na.rm = T) + out | 
        df.male$sa < mean(df.male$sa, na.rm = T) - out] <- NA # Mark people as NA

# Exclude cases from trimmed variables
df.male <- na.omit(df.male)

# Standardise predictors
grandmean <- mean(c(df.male$sa, df.male$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df.male$sa)^2)+(sd(df.male$ca)^2))/2)

df.male$sa.s  <- (df.male$sa-grandmean)/pooledsd
df.male$ca.s  <- (df.male$ca-grandmean)/pooledsd

# Add squared and interaction terms of the predictors
df.male$sa.s2 <- df.male$sa.s^2
df.male$ca.s2 <- df.male$ca.s^2
df.male$sa.ca <- df.male$sa.s*df.male$ca.s

# Descriptive statistics 
df.male$sa.bias <- df.male$sa - df.male$ca # substract chron. from subj. age
describe(df.male[,c("ca", "sa", "ls", "sa.bias")])
hist(df.male$ca, xlim = c(30, 90),  main = "Men: chron. age")
hist(df.male$sa, xlim = c(15, 100), main = "Men: subj. age")


# 4.1.2.) Model fitting and comparison

# Refit models
nu.fit <- sem(nu.model, data = df.male, se = "robust", estimator = "MLR")
fu.fit <- sem(fu.model, data = df.male, se = "robust", estimator = "MLR")
a1.fit <- sem(a1.model, data = df.male, se = "robust", estimator = "MLR")
b1.fit <- sem(b1.model, data = df.male, se = "robust", estimator = "MLR")
a2.fit <- sem(a2.model, data = df.male, se = "robust", estimator = "MLR")
b2.fit <- sem(b2.model, data = df.male, se = "robust", estimator = "MLR")
a3.fit <- sem(a3.model, data = df.male, se = "robust", estimator = "MLR")
b3.fit <- sem(b3.model, data = df.male, se = "robust", estimator = "MLR")
s1.fit <- sem(s1.model, data = df.male, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df.male, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df.male, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df.male, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df.male, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df.male, se = "robust", estimator = "MLR")

# Create list of all models
models <- list(fu.fit, a1.fit, b1.fit, a2.fit, b2.fit, a3.fit, b3.fit, 
               s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, nu.fit)

# Set names of models in list
names(models) <- c("full", "model-1a", "model-1b", "model-2a", "model-2b", 
                   "model-3a", "model-3b", "model-s1", "model-s2", "model-s3", 
                   "model-s4", "model-s5", "model-s6", "null")

# Get AICs and other information
aics.male.1 <- as.data.frame(AICcmodavg::aictab(models, modnames = names(models), second.ord = F))
aics.male.1

# Identify models with essentially the same log-likelihood as a simpler model
frm(aics.male.1) # full, model-s3, model-s5, model-s6, model-s1 (model-2b, model-3a not excluded)
models.r.male <- models[-c( 1, 10, 12, 13, 8)] # Exlude redundant models: 

# Get AICs and other fit information for final model set
aics.male.2 <- as.data.frame(AICcmodavg::aictab(models.r.male, 
                                                modnames = names(models.r.male), second.ord = F))
aics.male.2 


# 4.1.3) Inspect parameters and plot model

# Inspect best fitting model
b3.fit.male <- sem(b3.model, data = df.male, se = "robust", estimator = "MLR")
summary(b3.fit.male, ci = T, fit.measures = T)

### 3D RSA plot
r.fu.male  <- RSA(ls  ~ sa.s*ca.s, df.male)

plot(r.fu.male, model = "full",
     axes = c("LOC", "LOIC", "PA1"),
     axesStyles = list(LOC  = list(lty = "solid", lwd = 2, col = "black"),
                       LOIC = list(lty = "solid", lwd = 2, col = "black"),
                       PA1  = list(lty = "dotted", lwd = 2, col = "grey40")),
     xlab = "Subjective age", 
     ylab = "Chronlogical age", 
     zlab ="Life satisfaction",
     cex.tickLabel = 2, cex.axesLabel = 2,
     rotation = list(x = -48, y = 26, z = 20),
     label.rotation = list(x = 22, y = -51, z = 92),
     points = list(show = F), hull = T, legend = T,
     project = c("LOC", "LOIC", "PA1", "points"), param = F,
     pal = colorRampPalette(c("#24353f", "#385061", "#50748D", "#6b8a9f", 
                              "#7d98ab", "#9bafbe", "#bbc8d3", "#cad4dc", 
                              "#e8ecef", "#f0f2f4", "#ffffff"))(14))

# Examine how many points lay below and above first principal axis

# First plot PA1, LOC and dots with quadrants
plot(df.male$ca.s ~ jitter(df.male$sa.s, 3), pch = 16, xlab = "Subjective age",
     ylab = "Chronological age", col = rgb(0,0,0, alpha = 0.3))
abline(a = 3.276, b = 1.390, col = "mediumvioletred", lwd = 3)
abline(a = 0, b = 1, col = "green4", lwd = 3)
abline(v = 0, lty = "dotted", lwd = 3)
abline(h = 0, lty = "dotted", lwd = 3)

# Now calculate how many values are below and above the PA1
fpa.ca.fit <- 3.276 + 1.390 * df.male$sa.s
resi <- df.male$ca.s - fpa.ca.fit

sum(resi < 0) # below
sum(resi > 0) # above

## Check also how many values are below and above the LOC
loc.ca.fit <- 0 + 1 * df.male$sa.s
resi <- df.male$ca.s - loc.ca.fit

sum(resi < 0) # below
sum(resi > 0) # above


# 4.1.4.) Quantify optimal margin

# Prepare data
ca.male <- c(40,50,60,70,80,90) # Select chronological ages and put them in a data frame
ca.male <- data.frame(ca.male)

ca.male$ca.sd <- (ca.male$ca - grandmean)/pooledsd # Standardised chron. age

# Extract PA1 parameters from model fit
p10 <- b3.fit.male@ParTable$est[38] # p10
p11 <- b3.fit.male@ParTable$est[37] # p11

# Run calculations 
ca.male$sa.sd <- ((ca.male$ca.sd - p10)/p11) # Using PA1, calculate corresponding standardised subj. ages 
ca.male$sa    <- grandmean + (ca.male$sa.sd*pooledsd) # Unstandard. subjective ages
ca.male$bi    <- ca.male$ca-ca.male$sa # Calculate difference (i.e. optimal subj. age bias)
ca.male$bi.2  <- ca.male$sa-ca.male$ca # negative values indicate feeling younger
ca.male$bi.p  <- (ca.male$sa-ca.male$ca)/ca.male$ca # proportinal values (divided by age)

# Show results
ca.male


# 4.2.) Model for females
# -------------------------

# 4.2.1) Prepare data

# Create male data set
df.female <- subset(df.eligible, DhSex == 2)

# Trim variable chronological age
out <- 3*sd(df.female$ca, na.rm = T) # Calculate +/- 3 SD

sum(df.female$ca < mean(df.female$ca, na.rm = T) - out, na.rm = T) # -3 SD: number of excluded people
sum(df.female$ca > mean(df.female$ca, na.rm = T) + out, na.rm = T) # +3 SD: number of excluded people

df.female$ca[df.female$ca > mean(df.female$ca, na.rm = T) + out | 
             df.female$ca < mean(df.female$ca, na.rm = T) - out] <- NA # Mark people as NA

# Trim variable subjective age
out <- 3*sd(df.female$sa, na.rm = T) # Calculate +/- 3 SD

sum(df.female$sa < mean(df.female$sa, na.rm = T) - out, na.rm = T) # - 3 SD: number of excluded people
sum(df.female$sa > mean(df.female$sa, na.rm = T) + out, na.rm = T) # + 3 SD: number of excluded people

df.female$sa[df.female$sa > mean(df.female$sa, na.rm = T) + out | 
             df.female$sa < mean(df.female$sa, na.rm = T) - out] <- NA # Mark people as NA

# Exclude cases from trimmed variables
df.female <- na.omit(df.female)

# Standardise predictors
grandmean <- mean(c(df.female$sa, df.female$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df.female$sa)^2)+(sd(df.female$ca)^2))/2)

df.female$sa.s  <- (df.female$sa-grandmean)/pooledsd
df.female$ca.s  <- (df.female$ca-grandmean)/pooledsd

# Add squared and interaction terms of the predictors
df.female$sa.s2 <- df.female$sa.s^2
df.female$ca.s2 <- df.female$ca.s^2
df.female$sa.ca <- df.female$sa.s*df.female$ca.s

# Descriptive statistics 
df.female$sa.bias <- df.female$sa - df.female$ca # substract chron. from subj. age
describe(df.female[,c("ca", "sa", "ls", "sa.bias")])
hist(df.female$ca, xlim = c(30, 90), main = "Women: chron. age")
hist(df.female$sa, xlim = c(15, 100), main = "Women: subj. age")


# 4.2.2) Model fitting and comparison

# Refit models
nu.fit <- sem(nu.model, data = df.female, se = "robust", estimator = "MLR")
fu.fit <- sem(fu.model, data = df.female, se = "robust", estimator = "MLR")
a1.fit <- sem(a1.model, data = df.female, se = "robust", estimator = "MLR")
b1.fit <- sem(b1.model, data = df.female, se = "robust", estimator = "MLR")
a2.fit <- sem(a2.model, data = df.female, se = "robust", estimator = "MLR")
b2.fit <- sem(b2.model, data = df.female, se = "robust", estimator = "MLR")
a3.fit <- sem(a3.model, data = df.female, se = "robust", estimator = "MLR")
b3.fit <- sem(b3.model, data = df.female, se = "robust", estimator = "MLR")
s1.fit <- sem(s1.model, data = df.female, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df.female, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df.female, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df.female, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df.female, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df.female, se = "robust", estimator = "MLR")

# Create list of all models
models <- list(fu.fit, a1.fit, b1.fit, a2.fit, b2.fit, a3.fit, b3.fit, 
               s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, nu.fit)

# Set names of models in list
names(models) <- c("full", "model-1a", "model-1b", "model-2a", "model-2b", 
                   "model-3a", "model-3b", "model-s1", "model-s2", "model-s3", 
                   "model-s4", "model-s5", "model-s6", "null")

# Get AICs and other information
aics.female.1 <- as.data.frame(AICcmodavg::aictab(models, modnames = names(models), second.ord = F))
aics.female.1

# Identify models with essentially the same log-likelihood as a simpler model
frm(aics.female.1) # model-3b(!), full, model-s5, model-s1, model-s2, model-s6, model-s3
models.r.female <- models[-c( 7, 1, 12, 8, 9, 13, 10)] # Exlude redundant models

# Get AICs and other fit information for final model set
aics.female.2 <- as.data.frame(AICcmodavg::aictab(models.r.female, 
                                                  modnames = names(models.r.female), second.ord = F))
aics.female.2 


# 4.2.3) Inspect parameters and plot model

# Inspect model
b3.fit.female <- sem(b3.model, data = df.female, se = "robust", estimator = "MLR") 
summary(b3.fit.female, ci = T, fit.measures = T) 
# STILL FITTING B3 HERE BECAUSE I HAVEN'T DEFINED P10 AND P11 IN MODEL 3A
# BUT SHOULD HOPEFULLY BE SIMILAR VALUES

### 3D RSA plot
r.fu.female  <- RSA(ls  ~ sa.s*ca.s, df.female)

plot(r.fu.female, model = "full",
     axes = c("LOC", "LOIC", "PA1"),
     axesStyles = list(LOC  = list(lty = "solid", lwd = 2, col = "black"),
                       LOIC = list(lty = "solid", lwd = 2, col = "black"),
                       PA1  = list(lty = "dotted", lwd = 2, col = "grey40")),
     xlab = "Subjective age", 
     ylab = "Chronlogical age", 
     zlab ="Life satisfaction",
     cex.tickLabel = 2, cex.axesLabel = 2,
     rotation = list(x = -48, y = 26, z = 20),
     label.rotation=list(x = 22, y = -51, z = 92),
     points = list(show = F), hull = T, legend = T,
     project = c("LOC", "LOIC", "PA1", "points"), param = F,
     pal = colorRampPalette(c("#24353f", "#385061", "#50748D", "#6b8a9f", 
                              "#7d98ab", "#9bafbe", "#bbc8d3", "#cad4dc", 
                              "#e8ecef", "#f0f2f4", "#ffffff"))(14))

# First plot PA1, LOC and dots with quadrants
plot(df.female$ca.s ~ jitter(df.female$sa.s, 3), pch = 16, xlab = "Subjective age",
     ylab = "Chronological age", col = rgb(0,0,0, alpha = 0.3))
abline(a = 2.344, b = 1.144, col = "mediumvioletred", lwd = 3)
abline(a = 0, b = 1, col = "green4", lwd = 3)
abline(v = 0, lty = "dotted", lwd = 3)
abline(h = 0, lty = "dotted", lwd = 3)

# Now calculate how many values are below and above the PA1
fpa.ca.fit <- 2.344 + 1.144 * df.female$sa.s
resi <- df.female$ca.s - fpa.ca.fit

sum(resi < 0) # below
sum(resi > 0) # above

## Check also how many values are below and above the LOC
loc.ca.fit <- 0 + 1 * df.female$sa.s
resi <- df.female$ca.s - loc.ca.fit

sum(resi < 0) # below
sum(resi > 0) # above


# 4.2.4) Quantify optimal margin

# Prepare data
ca.female <- c(40,50,60,70,80,90) # Select chronological ages and put them in a data frame
ca.female <- data.frame(ca.female)

ca.female$ca.sd <- (ca.female$ca - grandmean)/pooledsd # Standardised chron. age

# Extract PA1 parameters from model fit
p10 <- b3.fit.female@ParTable$est[38] # p10
p11 <- b3.fit.female@ParTable$est[37] # p11

# Run calculations 
ca.female$sa.sd <- ((ca.female$ca.sd - p10)/p11) # Using PA1, calculate corresponding standardised subj. ages 
ca.female$sa    <- grandmean + (ca.female$sa.sd*pooledsd) # Unstandard. subjective ages
ca.female$bi    <- ca.female$ca-ca.female$sa # Calculate difference (i.e. optimal subj. age bias)
ca.female$bi.2  <- ca.female$sa-ca.female$ca # negative values indicate feeling younger
ca.female$bi.p  <- (ca.female$sa-ca.female$ca)/ca.female$ca # proportinal values (divided by age)

# Show results
ca.female
