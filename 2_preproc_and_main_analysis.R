# ===============================================================
# Preproc and main analysis for 
# "Psychological benefits of subjective age bias" 
# ==============================================================

# load packages 
library(dplyr)
library(lavaan)
library(qpcR)
library(RSA)
library(psych)
library(ggplot2)
library(car) # for vif()

source("C:/Users/Maria/Desktop/learn/0_PhD/Projects/ageing_rsa/analysis/subj_age_analysis/frm.R")

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

# 5.) Average subjective age bias

describe(df$sa.bias)


# 6.) Cronbachs alpha life satisfaction scale
psych::alpha(df
             [,c("sclifea_r", "sclifeb_r", "sclifec_r", "sclifed_r", "sclifee_r")]
             )

# 7.) Check for multicollinearity
lm <- lm(ls ~ sa.s + ca.s + sa.s2 + sa.ca + ca.s2, data = df)
vif(lm)


# ----------------------------------
# Fit polynomial regression models
# ----------------------------------

# Fit main models

nu.fit      <- sem(nu.model,      data = df, se = "robust", estimator = "MLR")
pb.fit      <- sem(pb.model,      data = df, se = "robust", estimator = "MLR")
omfr.fit    <- sem(omfr.model,    data = df, se = "robust", estimator = "MLR")
omrr.fit    <- sem(omrr.model,    data = df, se = "robust", estimator = "MLR")
somfr_d.fit <- sem(somfr_d.model, data = df, se = "robust", estimator = "MLR")
somfr_u.fit <- sem(somfr_u.model, data = df, se = "robust", estimator = "MLR")
somrr_d.fit <- sem(somrr_d.model, data = df, se = "robust", estimator = "MLR")
somrr_u.fit <- sem(somrr_u.model, data = df, se = "robust", estimator = "MLR")
fu.fit      <- sem(fu.model,      data = df, se = "robust", estimator = "MLR")


# Fit supplementary models

s1.fit <- sem(s1.model, data = df, se = "robust", estimator = "MLR")
s2.fit <- sem(s2.model, data = df, se = "robust", estimator = "MLR")
s3.fit <- sem(s3.model, data = df, se = "robust", estimator = "MLR")
s4.fit <- sem(s4.model, data = df, se = "robust", estimator = "MLR")
s5.fit <- sem(s5.model, data = df, se = "robust", estimator = "MLR")
s6.fit <- sem(s6.model, data = df, se = "robust", estimator = "MLR")
s7.fit <- sem(s7.model, data = df, se = "robust", estimator = "MLR")


# ----------------------------------
# Model evaluation and comparisons
# ----------------------------------

### Note that all "down" models provided a better fit than the respective "up"
### models; the "up" models (somfr_u, somrr_u) were therefore excluded from 
### further analyses

### Create list of all models
models <- list(fu.fit, 
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
names(models) <- c("full", 
                   "pb", 
                   "omfr", 
                   "omrr", 
                   "somfr", 
                   "somrr", 
                   "s1", 
                   "s2", 
                   "s3", 
                   "s4", 
                   "s5", 
                   "s6", 
                   "s7", 
                   "null")

### Get AICs and other information
aics <- as.data.frame(
  AICcmodavg::aictab(models, modnames = names(models), second.ord = F)
  )

### Now, identify models with essentially the same log-likelihood as a simpler 
### model and print warning if there are any
frm(aics)


### Full model, s1, s3, s5, and s7 are redundant --> exclude from list of models
models_r <- models[-c( 1, 11, 9, 13, 7)]

### Get AICs and other fit information for final model set
aics <- as.data.frame(
  AICcmodavg::aictab(models_r, modnames = names(models_r), second.ord = F)
  )

aics 


### Get CFI, SRMR, and RMSEA 
fitmeasures(somrr_d.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(omrr.fit,    fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(fu.fit,      fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(omfr.fit,    fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(pb.fit,      fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(somfr_d.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(nu.fit,      fit.measures = c("CFI", "SRMR", "RMSEA"))

### Get R squared
inspect(somrr_d.fit, 'r2')
inspect(omrr.fit,    'r2')
inspect(fu.fit,      'r2')
inspect(omfr.fit,    'r2')
inspect(pb.fit,      'r2')
inspect(somfr_d.fit, 'r2')
inspect(nu.fit,      'r2')

### Likelihood ratio tests of nested models
anova(fu.fit, somrr_d.fit)
anova(somrr_d.fit, omrr.fit)
anova(omrr.fit, omfr.fit)


### Get parameters SOMRR model
summary(somrr_d.fit, ci = T, fit.measures = T)
summary(omrr.fit, ci = T, fit.measures = T)
summary(omfr.fit, ci = T, fit.measures = T)

### Get a4 for SOMRR model and significance test using RSA function
r.fu  <- RSA(ls  ~ sa.s*ca.s, df) 
getPar(r.fu, model = "SRRR")


# --------
# Plots 
# --------

# We will plot the SOMRR model, one as a 3d and once as a contour plot.

# 1.) Re-fit full model using RSA function
r.fu  <- RSA(ls  ~ sa.s*ca.s, df)

# 2.) 3D plot
png("C:/Users/Maria/Desktop/learn/0_PhD/Projects/ageing_rsa/analysis/results/3d.png", 
    width = 12, height = 12, units = 'in', res = 600)
# modified seaborne
plot(r.fu, model = "SRRR",
     axes = c("LOC", "LOIC", "PA1"),
     axesStyles = list(LOC  = list(lty = "solid", lwd = 2, col = "black"),
                       LOIC = list(lty = "solid", lwd = 2, col = "black"),
                       PA1  = list(lty = "dotted", lwd = 2, col = "grey40")),
     xlab = "Subjective age", ylab = "Chronlogical age", zlab ="Life satisfaction",
     cex.tickLabel = 2, cex.axesLabel = 2,
     rotation = list(x=-48, y=26, z=20),
     label.rotation=list(x=22, y=-51, z=92),
     points = list(show=F), hull = T, legend = T,
     project = c("LOC", "LOIC", "PA1", "points"), param = F,
     pal = colorRampPalette(c("#22303a","#3c5466", "#50748D", "#648096", "#748D9D", 
                              "#96A8B6", "#B9C7D0", "#ccd6de", "#eef1f4", "#f4f6f8",
                              "#ffffff"))(14))
dev.off()


# 3.) Contour plot
p1 <- plot(r.fu, type = "contour", model = "SRRR",
           axes = c("LOC", "PA1"),
           xlab = "Subjective age", ylab = "Chronological age", 
           zlab ="Life satisfaction",
           cex.main = 1.3, cex.tickLabel = 2, cex.axesLabel = 2,
           points = list(show=T, jitter = 0.00, color = "black"), hull = F, 
           legend = F,
           pal = colorRampPalette(c("#22303a","#3c5466", "#50748D", "#648096", "#748D9D", 
                                    "#96A8B6", "#B9C7D0", "#ccd6de", "#eef1f4", "#f4f6f8",
                                    "#ffffff"))(14))

SP <- RSA.ST(r.fu, model = "SRRR")

png("C:/Users/Maria/Desktop/learn/0_PhD/Projects/ageing_rsa/analysis/results/contour.png", 
    width = 9, height = 9, units = 'in', res = 600)
p1 + stat_contour(bins=40, alpha=0.7, color = "grey20") +
  geom_abline(aes(intercept=0, slope=1), size=1, color="black") +
#  geom_abline(aes(intercept=0, slope=-1), size=1, color="black") +
  geom_abline(data=data.frame(SP[c("p10", "p11")]), 
              aes_string(intercept="p10", slope="p11"), 
              linetype="solid", color="#eef1f4", size=1)+
  geom_abline(data=data.frame(SP[c("p10", "p11")]), 
              aes_string(intercept="p10", slope="p11"), 
              linetype="dashed", color="black", size=1)+
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text = element_text(color = "black", size = 23), 
        axis.title.y = element_text(color = "black", size = 23, 
                                    margin = margin(t = 0, r = 30, b = 0, l = 0)), 
        axis.title.x = element_text(color = "black", size = 23,
                                    margin = margin(t = 30, r = 0, b = 0, l = 0)))
dev.off()



