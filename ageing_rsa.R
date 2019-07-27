# ===============================================================
# Analysis for 
# "Psychological benefits of subjective age bias" 
# ==============================================================

# load packages 
library(dplyr)
library(lavaan)
library(qpcR)
library(RSA)
library(psych)
library(ggplot2)

source("frm.R")

# load data  
df.raw  <- read.table("data/elsa/raw/tab/wave_2_core_data_v4.tab", sep = "\t", 
                      header = T)

# select relevant variables
df.full <- df.raw[,c("idauniq", "DhSex", "scold", "dhager", "sclifea", "sclifeb", 
                     "sclifec", "sclifed", "sclifee")]


#-------------------
# Data preparation
#-------------------

# 1.) Rename ageing variables
df.full$sa <- df.full$scold # subjective age
df.full$ca <- df.full$dhager # real age


# 2.) Replace missing values with NAs

### Recode values for missing data (-9, -1) in whole dataset as NA
df.full <- df.full %>% mutate_all(funs(na_if(., -9)))
df.full <- df.full %>% mutate_all(funs(na_if(., -1)))

### Recode old age values (99) as missing
df.full$ca[df.full$ca == 99] <- NA


# 3.) Recode sex variable
df.full$sex[df.full$DhSex == 2] <- "female"
df.full$sex[df.full$DhSex == 1] <- "male"


# 4.) Prepare outcome variable life satisfaction

### Recode life satisfaction (ls) items so that higher values indicate higher ls
df.full$sclifea_r <- 8-df.full$sclifea
df.full$sclifeb_r <- 8-df.full$sclifeb
df.full$sclifec_r <- 8-df.full$sclifec
df.full$sclifed_r <- 8-df.full$sclifed
df.full$sclifee_r <- 8-df.full$sclifee

### Calculate mean score on life satistfaction scale
df.full$ls <- rowMeans(df.full[,c("sclifea_r",  "sclifeb_r",  "sclifec_r", 
                                  "sclifed_r", "sclifee_r")], na.rm = T)


# 5.) Remove outliers chronological age

### Calculate +/- 3 SD
out <- 3*sd(df.full$ca, na.rm = T)

### How many outliers
sum(df.full$ca < mean(df.full$ca, na.rm = T) - out, na.rm = T) # -3 SD
sum(df.full$ca > mean(df.full$ca, na.rm = T) + out, na.rm = T) # +3 SD

### Mark people with chronological age +/- 3 SD as missing
df.full$ca[df.full$ca > mean(df.full$ca, na.rm = T) + out | 
             df.full$ca < mean(df.full$ca, na.rm = T) - out] <- NA


# 6.) Remove outliers subjective age

### Calculate +/- 3 SD
out <- 3*sd(df.full$sa, na.rm = T)

### How many outliers 
sum(df.full$sa < mean(df.full$sa, na.rm = T) - out, na.rm = T) # - 3 SD
sum(df.full$sa > mean(df.full$sa, na.rm = T) + out, na.rm = T) # + 3 SD

### Mark people with chronological age +/- 3 SD as missing
df.full$sa[df.full$sa > mean(df.full$sa, na.rm = T) + out | 
             df.full$sa < mean(df.full$sa, na.rm = T) - out] <- NA


# 7.) Exclude cases with missing data
df <- na.omit(df.full)


# 8.) Standardise predictors
grandmean <- mean(c(df$sa, df$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df$sa)^2)+(sd(df$ca)^2))/2)

df$sa.s  <- (df$sa-grandmean)/pooledsd
df$ca.s  <- (df$ca-grandmean)/pooledsd


# 9.) Add squared and interaction terms of the predictors
df$sa.s2 <- df$sa.s^2
df$ca.s2 <- df$ca.s^2
df$sa.ca <- df$sa.s*df$ca.s



# -------------------------
# Descriptive statistics
# -------------------------

# 1.) Descriptive stats before excluding missing data

### Continuous variables
describe(df.full[,c("ls", "sa", "ca")])

### Categorical variables
table(df.full$sex)
prop.table(table(df.full$sex))


# 2.) Descriptive stats after excluding missing data

### Continuous variables
describe(df[,c("sa", "ca", "ls")])

### Categorical variables
table(df$sex)
prop.table(table(df$sex))


# 3.) Zero-order correlation matrix
cor(df[,c("DhSex", "sa", "ca", "ls")])


# 4.) Cronbachs alpha life satisfaction scale
psych::alpha(df[,c("sclifea_r", "sclifeb_r", "sclifec_r", "sclifed_r", 
                   "sclifee_r")])



# ------------------------------------------
# Specify all polynomial regression models
# ------------------------------------------


# 0.) Null model 
# ------------------

nu.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 == 0
b2 == 0
b3 == 0
b4 == 0
b5 == 0"


# 1.) Positivity bias model
# ---------------------------

pb.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 <  0
b2 == 0
b3 == 0
b4 == 0
b5 == 0"


# 2.) Optimal margin models
# ---------------------------

# 2.1.) Optimal margin flat ridge

omfr.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 == -b2
b3 <  0
b3 == b5
b3 + b4 + b5 == 0"


# 2.2.) Optimal margin rising ridge

omrr.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2
# parameter constraints
b3 == b5
b3 + b4 + b5 == 0"


# 3.) Shifting optimal margin models
# ------------------------------------


# 3.1.) Shifting optimal margin flat ridge

### Flat ridge down

somfr_down.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2
# parameter constraints
b1 == (b2*b4)/2*b5
b3 < -0.000001
b5 < -0.000001
b4^2 == 4*b5*b3"


### Flat ridge up

somfr_up.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2
# parameter constraints
b1 == (b2*b4)/2*b5
b3 > 0.000001
b5 > 0.000001
b4^2 == 4*b5*b3"


# 3.2.) Shifting optimal margin rising ridge

### Rising ridge down

somrr_down.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b3 < -0.000001
b5 < -0.000001
b4^2 == 4*b3*b5"


### Rising ridge up

somrr_up.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b3 > 0.000001
b5 > 0.000001
b4^2 == 4*b3*b5"


# 4.) Full model
# -------------------

fu.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2"

fu.fit <- sem(fu.model, data = df, se = "robust", estimator = "MLR")
summary(fu.fit, fit.measures = T)



# 2.2.4.) Supplementary models
# -------------------------------

# Model S1: Chronological age model I
# negative effect of chronological age

s1.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 == 0
b2 <  0
b3 == 0
b4 == 0 
b5 == 0"


# Model S2: Chronological age model II
# positive effect of chronological age

s2.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 == 0
b2 >  0
b3 == 0
b4 == 0 
b5 == 0"


# Model S3: Curvelinear chronological age model
# first positive effect of age that decreases in older adults

s3.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 == 0
b3 == 0
b4 == 0
b5 <  0"


# Model S4: Curvelinear subjective age model
# positive effect of mean subjective age 

s4.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b2 == 0
b3 <  0
b4 == 0
b5 == 0"


# Model S5: Main effects model I
# positive effect younger subj age, and younger chronol age

s5.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 <  0
b2 <  0
b3 == 0
b4 == 0
b5 == 0"


# Model S6: Main effects model II
# positive effect younger subj age, and older chronol age

s6.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 <  0
b2 >  0
b3 == 0
b4 == 0
b5 == 0"


# Model S7: Congruency model
# positive effect of correct self-evaluation

s7.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 == 0
b2 == b1
b3 <  0
b4 == -2*b3
b5 == b3"



############################
############################
############################


# 2.3.) Fit all models


nu.fit <- sem(nu.model, data = df, se = "robust", estimator = "MLR")
summary(nu.fit, fit.measures = T)

pb.fit <- sem(pb.model, data = df, se = "robust", estimator = "MLR")
summary(pb.fit)

omfr.fit <- sem(omfr.model, data = df, se = "robust", estimator = "MLR")
summary(omfr.fit)


omrr.fit <- sem(omrr.model, data = df, se = "robust", estimator = "MLR")
summary(omrr.fit)

somfr_down.fit <- sem(somfr_down.model, se = "robust", data = df, 
                      estimator = "MLR")
summary(somfr_down.fit)

somfr_up.fit <- sem(somfr_up.model, se = "robust", data = df, estimator = "MLR")
summary(somfr_up.fit)

somrr_down.fit <- sem(somrr_down.model, data = df, se = "robust", 
                      estimator = "MLR")
summary(somrr_down.fit)


somrr_up.fit <- sem(somrr_up.model, data = df, se = "robust", estimator = "MLR")
summary(somrr_up.fit)

s1.fit <- sem(s1.model, data = df, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df, se = "robust", estimator = "MLR")
s7.fit <- sem(s7.model, data = df, se = "robust", estimator = "MLR")

############################
############################
############################




# 2.3.) Model comparisons
# -------------------------


# 2.3.1.) Exclude redundant models from set

### Create list of main models
main.models <- list(fu.fit, pb.fit, omfr.fit, omrr.fit, somfr_down.fit, 
                  somrr_down.fit, nu.fit)
names(main.models) <- c("full", "pb", "omfr", "omrr", "somfr", "somrr", "null")

### Get AICs and other information
aics.main <- as.data.frame(AICcmodavg::aictab(main.models, 
                                              modnames = names(main.models), 
                                              second.ord = F))

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
frm(aics.main)


# 2.3.2.) Model comparisons

### Update Modellist 
main.models.reduced <- main.models[-c(1)] # without full

aics.main <- as.data.frame(AICcmodavg::aictab(main.models.reduced, 
                                              modnames = names(main.models.reduced), 
                                              second.ord = F))
aics.main

# 2.3.1.) Exclude redundant models from set

### Create list of all models
modellist <- list(fu.fit, pb.fit, omfr.fit, omrr.fit, somfr_down.fit, 
                  somrr_down.fit, s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, 
                  s6.fit, s7.fit, nu.fit)
names(modellist) <- c("full", "pb", "omfr", "omrr", "somfr", "somrr",  "s1", 
                      "s2", "s3", "s4", "s5", "s6", "s7", "null")

### Get AICs and other information
aic.results <- as.data.frame(AICcmodavg::aictab(modellist, 
                                                modnames = names(modellist), 
                                                second.ord = F))

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
frm(aic.results)


# 2.3.2.) Model comparisons

### Update Modellist 
modellist <- modellist[-c(1,  11, 9, 13, 7)] # without full, s5, s3, s7, s1

aic.results <- as.data.frame(AICcmodavg::aictab(modellist, 
                                                modnames = names(modellist), 
                                                second.ord = F))
aic.results

### Get CFIs for all models
fitmeasures(srrr_down.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(srr.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(om.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))

# 2.3.3.) Statistical models
anova(fu.fit, srrr_down.fit)
anova(srrr_down.fit, srr.fit)
anova(srr.fit, om.fit)


# --------
# Plots 
# --------

# 1.) Plot AIC results all models
ggplot(data = aic.results, aes(x = rowname, y = AIC)) +
  geom_bar(stat = "identity", fill = "plum") +
  scale_x_discrete(limits=c("som.fit", "fu.fit", "om.fit", "s6.fit", "s4.fit", 
                            "pb.fit", "s5.fit", "s2.fit", "s3.fit", "nu.fit",
                            "s7.fit", "s1.fit")) +
  coord_cartesian(ylim = c(23500, 24200))+
  theme(axis.text  = element_text(colour = "black", size = 20), 
        axis.title = element_text(colour = "black", size = 23), 
        legend.position = "none") 


# 2.) Plot AIC results for 3 models
aic.results.hyp <- aic.results %>% filter(
  rowname1 == "som.fit" | rowname == "om.fit" |rowname1 == "pb.fit"
)

ggplot(data = aic.results.hyp, aes(x = rowname, y = AIC)) +
  geom_bar(stat = "identity", fill = c("#183442", "#91A5AE","#A6A3AA"), color = "black") +
  scale_x_discrete(limits = c("som.fit", "om.fit", "pb.fit"), 
                   labels = c("Shifting optimal \n margin model", 
                              "Optimal margin \n model", 
                              "Positivity bias \n model")) +
  coord_cartesian(ylim = c(23800, 24100))+
  theme_classic(base_size = 20) +
  labs(x = "", y = "Akaike's Information Criterion") +
  theme(axis.text  = element_text(colour = "black", size = 20), 
        axis.title = element_text(colour = "black", size = 23), 
        legend.position = "none") 


# 3.) Plot Akaike weight result
ggplot(data = aic.results, aes(x = rowname, y = weights)) +
  geom_bar(stat = "identity", fill = "plum") +
  scale_x_discrete(limits=c("som.fit", "fu.fit", "om.fit", "s6.fit", "s4.fit", 
                            "pb.fit", "s5.fit", "s2.fit", "s3.fit", "nu.fit",
                            "s7.fit", "s1.fit")) +
  theme_classic()


# 4.) Plot RSA result full model

# Re-fit full model using RSA function
r.fu  <- RSA(ls  ~ sa.s*ca.s, df, model = "full")
summary(r.fu)

# Plot full model

#png("C:/Users/Maria/Desktop/learn/0_PhD/Projects/ageing_rsa/analysis/results/model.png", 
#    width = 6, height = 6, units = 'in', res = 600)
plot(r.fu, type = "contour",
     axes = c("LOC", "LOIC", "PA1"),
     axesStyles = list(LOC  = list(lty = "solid", lwd = 2, col = "black"),
                       LOIC = list(lty = "solid", lwd = 2, col = "black"),
                       PA1  = list(lty = "dashed", lwd = 2, col = "grey40")),
     xlab = "Subjective age", ylab = "Actual age", zlab ="Life satisfaction",
     cex.main = 1.3, cex.tickLabel = 1.25, cex.axesLabel = 1.25,
     points = list(show=TRUE, jitter = 0.01), hull = T, legend = F,
     project = c("LOC", "LOIC", "PA1", "contour"), param = F,
     pal = colorRampPalette(c("#4c3633", "#795a70","#8f86a6",
                              "#a3b5c1", "#c8dcd2", "#EEF4F1"))(15))
#dev.off()

#png("C:/Users/Maria/Desktop/learn/0_PhD/Projects/ageing_rsa/analysis/results/model.png", 
#    width = 6, height = 6, units = 'in', res = 600)
plot(r.fu,
     axes = c("LOC", "LOIC", "PA1"),
     axesStyles = list(LOC  = list(lty = "solid", lwd = 2, col = "black"),
                       LOIC = list(lty = "solid", lwd = 2, col = "black"),
                       PA1  = list(lty = "dashed", lwd = 2, col = "grey40")),
     xlab = "Subjective age", ylab = "Actual age", zlab ="Life satisfaction",
     cex.main = 1.3, cex.tickLabel = 1.25, cex.axesLabel = 1.25,
     points = list(show=F), hull = T, legend = F,
     project = c("LOC", "LOIC", "PA1", "contour"), param = F,
     pal = colorRampPalette(c("#4c3633", "#795a70","#8f86a6",
                              "#a3b5c1", "#c8dcd2", "#EEF4F1"))(15))
#dev.off()

# ------------------------
#### Additional double check stuff: 
# ------------------------------

# Examine how many points lay below and above first principal axis

## First plot PA1, LOC and dots with quadrants
plot(df$ca.s ~ jitter(df$sa.s, 3), pch = 16, xlab = "Subjective age", 
     ylab = "Chronological age", col = rgb(0,0,0, alpha = 0.3))
abline(a = 2.649, b = 1.170, col = "mediumvioletred", lwd = 3)
abline(a = 0, b = 1, col = "green4", lwd = 3)
abline(v = 0, lty = "dotted", lwd = 3)
abline(h = 0, lty = "dotted", lwd = 3)

## Now calculate how many values are below and above the PA1
fpa.ca.fit <- 2.649 + 1.170 * df$sa.s
resi <- df$ca.s - fpa.ca.fit

sum(resi < 0)
sum(resi > 0)

## Chack also how many values are below and above the LOC
loc.ca.fit <- 0 + 1 * df$sa.s
resi <- df$ca.s - loc.ca.fit

sum(resi < 0)
sum(resi > 0)



# -------------------------------------------------------------
# Re-un analyses excluding outliers and people aged < 40 years
# -------------------------------------------------------------

# mark people with age < 40 y as missing
df$ca[df$ca < 40] <- NA

# 3. Exclude cases with missing data (i.e. outliers)
df <- na.omit(df)

# Get descriptive stats after exclusion
describe(df$ca)
describe(df$sa)

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
pb.fit <- sem(pb.model, data = df, se = "robust", estimator = "MLR")
om.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s1.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s2.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s3.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s4.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s5.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s6.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
# s7.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")

# Get AICs of all models in the full model set
aic.results <- 
  AIC(nu.fit, fu.fit, 
      pb.fit, om.fit, som.fit, 
      s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, s7.fit)

# Order results according to AIC (models with low AIC first)
aic.results <- aic.results[order(aic.results$AIC),]

# Add rownames as column to data frame
aic.results <- add_rownames(aic.results, var = "rowname")

# Make a vector from AICs 
aics.vector <- pull(aic.results, AIC)

# Calculate Akaike weights and attach to results of model comparisons
weights <- akaike.weights(aics.vector)
aic.results$weights <- weights$weights



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

