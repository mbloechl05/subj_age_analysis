# ===============================================================
# Subjective age bias and life satisfaction
# Script 2: Preproc and main analysis 
# ==============================================================

# This script contains code to run the main analyses.  
# Script 1 has to be run before running this script since script 1 
# contains all model definitions, which are are called upon here. 

# load packages 
library(dplyr)
library(lavaan)
library(qpcR)
library(RSA)
library(psych)
library(ggplot2)
library(car) # for vif()
library(reshape)
library(aplpack) # for compute bagplot

# source helper scripts
source("analysis/subj_age_analysis/frm.R")
source("analysis/subj_age_analysis/rsa_plot.R")

# load data  
df.raw  <- read.table("data/wave_2_core_data_v4.tab", sep = "\t", header = T)

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

# 5.) Exclude cases with missing data on variables
df.eligible <- na.omit(df.full)


# 6.) Trim variables

df <- df.eligible

### 6.a) Trim variable chronological age

### Calculate +/- 3 SD
out <- 3*sd(df$ca, na.rm = T)

### How many will be excluded?
sum(df$ca < mean(df$ca, na.rm = T) - out, na.rm = T) # -3 SD
sum(df$ca > mean(df$ca, na.rm = T) + out, na.rm = T) # +3 SD

### Mark people with chronological age +/- 3 SD as missing
df$ca[df$ca > mean(df$ca, na.rm = T) + out | 
             df$ca < mean(df$ca, na.rm = T) - out] <- NA


#### 6.b) Trim variable subjective age

### Calculate +/- 3 SD
out <- 3*sd(df$sa, na.rm = T)

### How many will be excluded?
sum(df$sa < mean(df$sa, na.rm = T) - out, na.rm = T) # - 3 SD
sum(df$sa > mean(df$sa, na.rm = T) + out, na.rm = T) # + 3 SD

### Mark people with chronological age +/- 3 SD as missing
df$sa[df$sa > mean(df$sa, na.rm = T) + out | 
             df$sa < mean(df$sa, na.rm = T) - out] <- NA


### 6.c) Exclude cases from trimmed variables
df <- na.omit(df)


# 7.) Standardise predictors
grandmean <- mean(c(df$sa, df$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df$sa)^2)+(sd(df$ca)^2))/2)

df$sa.s  <- (df$sa-grandmean)/pooledsd
df$ca.s  <- (df$ca-grandmean)/pooledsd


# 8.) Add squared and interaction terms of the predictors
df$sa.s2 <- df$sa.s^2
df$ca.s2 <- df$ca.s^2
df$sa.ca <- df$sa.s*df$ca.s


# -------------------------
# Descriptive statistics
# -------------------------

# 1.) Descriptive stats full sample

### Continuous variables
describe(df.full[,c("ls", "sa", "ca")])

### Categorical variables
table(df.full$sex)
prop.table(table(df.full$sex))


# 2.) Descriptive stats eligible sample

### Continuous variables
describe(df.eligible[,c("ls", "sa", "ca")])

### Categorical variables
table(df.eligible$sex)
prop.table(table(df.eligible$sex))


# 3.) Descriptive stats after trimming data

### Continuous variables
df$sa.bias <- df$sa - df$ca # substract chron. from subj. age
describe(df[,c("sa", "ca", "ls", "DhSex", "sa.bias")])

### Categorical variables
table(df$sex)
prop.table(table(df$sex))


# 4.) Zero-order correlation matrix
cor(df[,c("DhSex", "sa", "ca", "ls", "sa.bias")])

# 5.) Cronbachs alpha life satisfaction scale
psych::alpha(df[,c("sclifea_r", "sclifeb_r", "sclifec_r", 
                   "sclifed_r", "sclifee_r")])

# 6.) Check for multicollinearity
lm <- lm(ls ~ sa.s + ca.s + sa.s2 + sa.ca + ca.s2, data = df)
vif(lm)


# ----------------------------------
# Fit polynomial regression models
# ----------------------------------

# Fit full and null model
nu.fit <- sem(nu.model, data = df, se = "robust", estimator = "MLR")
fu.fit <- sem(fu.model, data = df, se = "robust", estimator = "MLR")

# Fit main models
a1.fit <- sem(a1.model, data = df, se = "robust", estimator = "MLR")
b1.fit <- sem(b1.model, data = df, se = "robust", estimator = "MLR")
a2.fit <- sem(a2.model, data = df, se = "robust", estimator = "MLR")
b2.fit <- sem(b2.model, data = df, se = "robust", estimator = "MLR")
a3.fit <- sem(a3.model, data = df, se = "robust", estimator = "MLR")
b3.fit <- sem(b3.model, data = df, se = "robust", estimator = "MLR")

# Fit supplementary models
s1.fit <- sem(s1.model, data = df, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df, se = "robust", estimator = "MLR")


# ----------------------------------
# Model evaluation and comparisons
# ----------------------------------

### Note that all "down" models provided a better fit than the respective "up"
### models; the "up" models (somfr_u, somrr_u) were therefore excluded from 
### further analyses

### Create list of all models
models <- list(fu.fit,
               a1.fit, 
               b1.fit, 
               a2.fit, 
               b2.fit, 
               a3.fit, 
               b3.fit, 
               s1.fit, 
               s2.fit, 
               s3.fit, 
               s4.fit, 
               s5.fit,   
               s6.fit, 
               nu.fit)

### Set names of models in list
names(models) <- c("full",     #1
                   "model-1a", #2
                   "model-1b", #3
                   "model-2a", #4
                   "model-2b", #5
                   "model-3a", #6
                   "model-3b", #7
                   "model-s1", #8
                   "model-s2", #9
                   "model-s3", #10
                   "model-s4", #11
                   "model-s5", #12
                   "model-s6", #13
                   "null")     #14

### Get AICs and other information
aics <- as.data.frame(
  AICcmodavg::aictab(models, modnames = names(models), second.ord = T)
  )

aics

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
frm(aics)


### full, model-3a, model-s5, model-s3, model-s6, model-s1 are redundant --> 
### exclude from list of models
models_r <- models[-c( 1, 6, 12, 10, 13, 8)]


### Get AICs and other fit information for final model set
aics <- as.data.frame(
  AICcmodavg::aictab(models_r, modnames = names(models_r), second.ord = F)
  )

aics 


### Get CFI, SRMR, and RMSEA 
fitmeasures(b3.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(b2.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(fu.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(a2.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(a3.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(b1.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(a1.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(nu.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))

### Get R squared
inspect(b3.fit, 'r2')
inspect(b2.fit, 'r2')
inspect(fu.fit, 'r2')
inspect(a2.fit, 'r2')
inspect(a3.fit, 'r2')
inspect(b1.fit, 'r2')
inspect(a1.fit, 'r2')
inspect(nu.fit, 'r2')

### Likelihood ratio tests of nested models
anova(fu.fit, b3.fit)
anova(fu.fit, b2.fit)
anova(fu.fit, a2.fit)


### Get parameters SOMRR model
summary(b3.fit, ci = T, fit.measures = T)
summary(b2.fit, ci = T, fit.measures = T)
summary(a2.fit, ci = T, fit.measures = T)
summary(fu.fit, ci = T, fit.measures = T)

# -------------------------
# Quantify optimal margin 
# -------------------------

# In the following, we will calculate the "optimal" subjective age for 
# certain chronological ages. 

# 1.) Select chronological ages and put them in a data frame
ca <- c(40,50,60,70,80,90) 
ca <- data.frame(ca)

# 2.) Standardised chronological ages to grandmean and pooled sd
ca$ca.sd <- (ca$ca - grandmean)/pooledsd 

# 3.) Extract PA1 parameters from model fit
p10 <- b3.fit@ParTable$est[38] # p10
p11 <- b3.fit@ParTable$est[37] # p11

# 4.) Using PA1, calculate corresponding standardised subj. ages 
# i.e. these values lay on the PA1 and at them, life satisfaction is highest
ca$sa.sd <- ((ca$ca.sd - p10)/p11) 

# 5.) Unstandard. subjective ages
ca$sa <- grandmean + (ca$sa.sd*pooledsd) 

# 6.) Calculate difference (i.e. optimal subj. age bias)
ca$bi   <- ca$ca-ca$sa 
ca$bi.2 <- ca$sa-ca$ca # negative values indicate feeling younger
ca$bi.p <- (ca$sa-ca$ca)/ca$ca # proportinal values (divided by age)

# 7.) Show results in data frame
ca 


# --------
# Plots 
# --------


# Figure 2A): Response surface of best fitting model

### Re-fit full model using RSA function

r.fu  <- RSA(ls  ~ sa.s*ca.s, df)

### 3D RSA plot

png("analysis/results/3d_new.png", 
    width = 12, height = 12, units = 'in', res = 600)
plot(r.fu, model = "SRRR",
     axes = c("LOC", "PA1"),
     axesStyles = list(LOC  = list(lty = "solid", lwd = 2, col = "black"),
                       LOIC = list(lty = "solid", lwd = 2, col = "black"),
                       PA1  = list(lty = "dotted", lwd = 2, col = "black")),
     xlab = "Subjective age", 
     ylab = "Chronlogical age", 
     zlab ="Life satisfaction",
     cex.tickLabel = 2, cex.axesLabel = 2,
     rotation = list(x = -48, y = 26, z = 20),
     label.rotation = list(x = 22, y = -51, z = 92),
     points = list(show = F), hull = T, legend = T,
     project = c("LOC", "PA1", "points"), param = F,
     pal = colorRampPalette(c("#24353f", "#385061", "#50748D", "#6b8a9f", 
                              "#7d98ab", "#9bafbe", "#bbc8d3", "#cad4dc", 
                              "#e8ecef", "#f0f2f4", "#ffffff"))(14))
dev.off()


# Figure 2B): Countour plot
# note: please source rsa_plot before (contains modified contour plot function)

SP <- RSA.ST(r.fu, model = "SRRR")

p1 <- plot(r.fu, type = "contour", model = "SRRR",
           axes = NA,
           xlab = "Subjective age", ylab = "Chronological age", 
           zlab ="Life satisfaction",
           cex.main = 1.3, cex.tickLabel = 2, cex.axesLabel = 2,
           points = list(show=F, jitter = 0.00, color = "black"), hull = F, 
           legend = F,
           pal = colorRampPalette(c("#24353f", "#385061", "#50748D", "#6b8a9f", 
                                    "#7d98ab", "#9bafbe", "#bbc8d3", "#cad4dc", 
                                    "#e8ecef", "#f0f2f4", "#ffffff"))(14))

png("analysis/results/contour.png", 
    width = 9, height = 9, units = 'in', res = 600)
p1 + 
  geom_abline(aes(intercept = 0, slope = 1), size = 1.5, color="black") +
  geom_abline(data=data.frame(SP[c("p10", "p11")]), 
              aes_string(intercept="p10", slope = "p11"), 
              linetype = "dashed", color="black", size = 1.5)+
  geom_polygon(data = hull.loop, aes(x = x, y = y), inherit.aes = FALSE, 
               colour = "black", fill = NA, 
               alpha = 0.1, linetype = "dashed", size = 0.9) +
  geom_polygon(data = hull.bag, aes(x = x, y = y), inherit.aes = FALSE, 
               colour = "black", fill = NA, alpha = 0.1, size = 0.9) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(color = "black", size = 23), 
        axis.title.y = element_text(color = "black", size = 23, 
                                    margin = margin(t = 0, r = 30, b = 0, l = 0)), 
        axis.title.x = element_text(color = "black", size = 23,
                                    margin = margin(t = 30, r = 0, b = 0, l = 0)))
dev.off()

