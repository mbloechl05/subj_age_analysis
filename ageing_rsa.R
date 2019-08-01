# ===============================================================
# Main analysis for 
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


# 5.) Check for multicollinearity
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
models_r <- models[-c(1,  11, 9, 13, 7)]

### Get AICs and other fit information for final model set
aics <- as.data.frame(
  AICcmodavg::aictab(models_r, modnames = names(models_r), second.ord = F)
  )

aics 


### Get parameters full model and SOMRR model
summary(somrr_d.fit, ci = T)
summary(fu.fit, ci = T)

### Get CFI, SRMR, and RMSEA 
fitmeasures(somrr_d.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(omrr.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))
fitmeasures(omfr.fit, fit.measures = c("CFI", "SRMR", "RMSEA"))


### Likelihood ratio tests of nested models
anova(fu.fit, somrr_d.fit)
anova(somrr_d.fit, omrr.fit)
anova(omrr.fit, omfr.fit)



# --------
# Plots 
# --------

# # 1.) Plot AIC results all models
# ggplot(data = aic.results, aes(x = rowname, y = AIC)) +
#   geom_bar(stat = "identity", fill = "plum") +
#   scale_x_discrete(limits=c("som.fit", "fu.fit", "om.fit", "s6.fit", "s4.fit", 
#                             "pb.fit", "s5.fit", "s2.fit", "s3.fit", "nu.fit",
#                             "s7.fit", "s1.fit")) +
#   coord_cartesian(ylim = c(23500, 24200))+
#   theme(axis.text  = element_text(colour = "black", size = 20), 
#         axis.title = element_text(colour = "black", size = 23), 
#         legend.position = "none") 
# 
# 
# # 2.) Plot AIC results for 3 models
# aic.results.hyp <- aic.results %>% filter(
#   rowname1 == "som.fit" | rowname == "om.fit" |rowname1 == "pb.fit"
# )
# 
# ggplot(data = aic.results.hyp, aes(x = rowname, y = AIC)) +
#   geom_bar(stat = "identity", fill = c("#183442", "#91A5AE","#A6A3AA"), color = "black") +
#   scale_x_discrete(limits = c("som.fit", "om.fit", "pb.fit"), 
#                    labels = c("Shifting optimal \n margin model", 
#                               "Optimal margin \n model", 
#                               "Positivity bias \n model")) +
#   coord_cartesian(ylim = c(23800, 24100))+
#   theme_classic(base_size = 20) +
#   labs(x = "", y = "Akaike's Information Criterion") +
#   theme(axis.text  = element_text(colour = "black", size = 20), 
#         axis.title = element_text(colour = "black", size = 23), 
#         legend.position = "none") 
# 
# 
# # 3.) Plot Akaike weight result
# ggplot(data = aic.results, aes(x = rowname, y = weights)) +
#   geom_bar(stat = "identity", fill = "plum") +
#   scale_x_discrete(limits=c("som.fit", "fu.fit", "om.fit", "s6.fit", "s4.fit", 
#                             "pb.fit", "s5.fit", "s2.fit", "s3.fit", "nu.fit",
#                             "s7.fit", "s1.fit")) +
#   theme_classic()


# 4.) Plot RSA result full model

# Re-fit full model using RSA function
r.fu  <- RSA(ls  ~ sa.s*ca.s, df)

# Plot full model
png("C:/Users/Maria/Desktop/learn/0_PhD/Projects/ageing_rsa/analysis/results/3d.png", 
    width = 12, height = 12, units = 'in', res = 600)
# modified seaborne
plot(r.fu, model = "SRRR",
     axes = c("LOC", "LOIC", "PA1"),
     axesStyles = list(LOC  = list(lty = "solid", lwd = 2, col = "black"),
                       LOIC = list(lty = "solid", lwd = 2, col = "black"),
                       PA1  = list(lty = "dotted", lwd = 2, col = "grey40")),
     xlab = "Subjective age", ylab = "Actual age", zlab ="Life satisfaction",
     cex.tickLabel = 2, cex.axesLabel = 2,
     points = list(show=F), hull = T, legend = T,
     project = c("LOC", "LOIC", "PA1"), param = F,
     pal = colorRampPalette(c("#3e2b31", "#684852","#916572","#D091A4","#DFB1B4",
                              "#EBD2CB","#f7edea","#f8f0ee", "#FBFFFE"))(15))
dev.off()

png("C:/Users/Maria/Desktop/learn/0_PhD/Projects/ageing_rsa/analysis/results/contour.png", 
    width = 9, height = 9, units = 'in', res = 600)
plot(r.fu, type = "contour", model = "SRRR",
     axes = c("LOC", "PA1"),
     axesStyles = list(LOC  = list(lty = "solid", lwd = 6, col = "black"),
                       PA1 = list(lty = "solid", lwd = 6, col = "black")),
                      # LOIC  = list(lty = "solid", lwd = 6, col = "grey40")),
     xlab = "Subjective age", ylab = "Actual age", zlab ="Life satisfaction",
     cex.main = 1.3, cex.tickLabel = 2, cex.axesLabel = 2,
     points = list(show=T, jitter = 0.005, color = "black"), hull = T, legend = F,
     #project = c("LOC", "PA1"), param = F, 
     pal = colorRampPalette(c("#3e2b31", "#684852","#916572","#D091A4","#DFB1B4",
                              "#EBD2CB","#f7edea","#f8f0ee", "#FBFFFE"))(15))
dev.off()



# # ------------------------
# #### Additional double check stuff: 
# # ------------------------------
# 
# 
# # Examine how many points lay below and above first principal axis
# 
# ## First plot PA1, LOC and dots with quadrants
# plot(df$ca.s ~ jitter(df$sa.s, 3), pch = 16, xlab = "Subjective age", 
#      ylab = "Chronological age", col = rgb(0,0,0, alpha = 0.3))
# abline(a = 2.649, b = 1.170, col = "mediumvioletred", lwd = 3)
# abline(a = 0, b = 1, col = "green4", lwd = 3)
# abline(v = 0, lty = "dotted", lwd = 3)
# abline(h = 0, lty = "dotted", lwd = 3)
# 
# ## Now calculate how many values are below and above the PA1
# fpa.ca.fit <- 2.649 + 1.170 * df$sa.s
# resi <- df$ca.s - fpa.ca.fit
# 
# sum(resi < 0)
# sum(resi > 0)
# 
# ## Chack also how many values are below and above the LOC
# loc.ca.fit <- 0 + 1 * df$sa.s
# resi <- df$ca.s - loc.ca.fit
# 
# sum(resi < 0)
# sum(resi > 0)




