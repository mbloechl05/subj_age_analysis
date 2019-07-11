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

# load data from wave 2 
df.raw  <- read.table("data/elsa/raw/tab/wave_2_core_data_v4.tab", sep = "\t", 
                      header = T)

# select relevant variables
df.full <- df.raw[,c("idauniq", "DhSex", "scold", "dhager", "sclifea", "sclifeb", 
                     "sclifec", "sclifed", "sclifee")]


#-------------------
# Data preparation
#-------------------

# 1. Rename ageing variables
df.full$sa <- df.full$scold # subjective age
df.full$ca <- df.full$dhager # real age


# 2. Replace missing values with NAs

# recode values for missing data (-9, -1) in whole dataset as NA
df.full <- df.full %>% mutate_all(funs(na_if(., -9)))
df.full <- df.full %>% mutate_all(funs(na_if(., -1)))

# recode old age values (99) as missing
df.full$ca[df.full$ca == 99] <- NA

# 3. Recode sex variable
df.full$sex[df.full$DhSex == 2] <- "female"
df.full$sex[df.full$DhSex == 1] <- "male"

# 3. Prepare outcome variable life satisfaction

# recode life satisfaction (ls) items so that higher values indicate higher ls
df.full$sclifea_r <- 8-df.full$sclifea
df.full$sclifeb_r <- 8-df.full$sclifeb
df.full$sclifec_r <- 8-df.full$sclifec
df.full$sclifed_r <- 8-df.full$sclifed
df.full$sclifee_r <- 8-df.full$sclifee

# calculate mean score on life satistfaction scale
df.full$ls <- rowMeans(df.full[,c("sclifea_r",  "sclifeb_r",  "sclifec_r", 
                                  "sclifed_r", "sclifee_r")], na.rm = T)

# 4. Exclude cases with missing data
df <- na.omit(df.full)


# standardise predictors
grandmean <- mean(c(df$sa, df$ca), na.rm = T)
pooledsd  <- sqrt(((sd(df$sa)^2)+(sd(df$ca)^2))/2)

df$sa.s  <- (df$sa-grandmean)/pooledsd
df$ca.s  <- (df$ca-grandmean)/pooledsd


# add squared and interaction terms of the predictors (required to inspect multicollinearity)
df$sa.s2 <- df$sa.s^2
df$ca.s2 <- df$ca.s^2
df$sa.ca <- df$sa.s*df$ca.s


# ------------------------
# Descriptive statistics
# ------------------------

# 1. Descriptive stats before excluding missing data

# Continuous variables
describe(df.full[,c("ls", "sa", "ca")])

# Categorical variables
table(df.full$sex)
prop.table(table(df.full$sex))


# 2. Descriptive stats after excluding missing data

# Continuous variables
describe(df[,c("sa", "ca", "ls")])

# Categorical variables
table(df$sex)
prop.table(table(df$sex))


# 3. Zero-order correlation matrix
cor(df[,c("DhSex", "sa", "ca", "ls")])


# 3. Cronbachs alpha life satisfaction scale
alpha(df[,c("sclifea_r", "sclifeb_r", "sclifec_r", "sclifed_r", "sclifee_r")])


# -------------------------------------
# Fit polynomial regression models
# -------------------------------------

# 2.0) Null model 
# ------------------

nu.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == 0
b2 == 0
b3 == 0
b4 == 0
b5 == 0"

nu.fit <- sem(nu.model, data = df, se = "robust", estimator = "MLR")
summary(nu.fit)


# 2.1) Full model
# -------------------

fu.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2"

fu.fit <- sem(fu.model, data=df)
summary(fu.fit, fit.measures = T)


# 2.2) Hypothesised models
# ----------------------------

# Model 1: Positivity bias

pb.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 <  0
b2 == 0
b3 == 0
b4 == 0
b5 == 0"

pb.fit <- sem(pb.model, data = df, se = "robust", estimator = "MLR")
summary(pb.fit)


# Model 2: Optimal margin

om.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == -b2
b3 < 0
b3 == b5
b3 + b4 + b5 == 0"

om.fit <- sem(om.model, data = df, se = "robust", estimator = "MLR")
summary(om.fit)


# Model 3: Shifting optimal margin

som.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == -b2
b3 < 0
b3 < b5
b3 + b4 + b5 == 0"

som.fit <- sem(som.model, data = df, se = "robust", estimator = "MLR")
summary(som.fit)


# 2.2.4.) Supplementary models
# -------------------------------

# Model S1: Chronological age model I
# negative effect of chronological age

s1.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == 0
b2 <  0
b3 == 0
b4 == 0 
b5 == 0"

s1.fit <- sem(s1.model, data = df, se = "robust", estimator = "MLR")


# Model S2: Chronological age model II
# positive effect of chronological age

s2.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == 0
b2 >  0
b3 == 0
b4 == 0 
b5 == 0"

s2.fit <- sem(s2.model, data = df, se = "robust", estimator = "MLR")


# Model S3: Curvelinear chronological age model
# first positive effect of age that decreases in older adults

s3.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == 0
b3 == 0
b4 == 0
b5 <  0"

s3.fit <- sem(s3.model, data = df, se = "robust", estimator = "MLR")


# Model S4: Curvelinear subjective age model
# positive effect of mean subjective age 

s4.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b2 == 0
b3 <  0
b4 == 0
b5 == 0"

s4.fit <- sem(s4.model, data = df, se = "robust", estimator = "MLR")


# Model S5: Main effects model I
# positive effect younger subj age, and younger chronol age

s5.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 <  0
b2 <  0
b3 == 0
b4 == 0
b5 == 0"

s5.fit <- sem(s5.model, data = df, se = "robust", estimator = "MLR")


# Model S6: Main effects model II
# positive effect younger subj age, and older chronol age

s6.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 <  0
b2 >  0
b3 == 0
b4 == 0
b5 == 0"

s6.fit <- sem(s6.model, data = df, se = "robust", estimator = "MLR")


# Model S7: Congruency model
# positive effect of correct self-evaluation

s7.model <- "
ls ~ 1 + 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == 0
b2 == b1
b3 <  0
b4 == -2*b3
b5 == b3"

s7.fit <- sem(s7.model, data = df, se = "robust", estimator = "MLR")


# 2.3.) Model comparisons
# -------------------------

# 2.3.1.) Model comparisons based on AIC

# Get AICs of all models in the full model set
aic.results <- 
  AIC(nu.fit, fu.fit, 
      pb.fit, om.fit, som.fit, 
      s1.fit, s2.fit, s3.fit, s4.fit, s5.fit, s6.fit, s7.fit)

# Order results according to AIC (models with low AIC first)
aic.results <- aic.results[order(aic.results$AIC),]

# Add rownames as column to data frame
aic.results <- add_rownames(aic.results, var = "rowname")

# Plot result
ggplot(data = aic.results, aes(x = rowname, y = AIC)) +
  geom_bar(stat = "identity", fill = "plum") +
  scale_x_discrete(limits=c("som.fit", "fu.fit", "om.fit", "s6.fit", "s4.fit", 
                            "pb.fit", "s5.fit", "s2.fit", "s3.fit", "nu.fit",
                            "s7.fit", "s1.fit")) +
#  scale_y_reverse() +
#  scale_y_continuous(limits = c(23500, 24200)) +
  coord_cartesian(ylim = c(23500, 24200))+
  theme(axis.text  = element_text(colour = "black", size = 20), 
        axis.title = element_text(colour = "black", size = 23), 
        legend.position = "none") 


# Plot results for 3 models
aic.results.hyp <- aic.results %>% filter(
  rowname1 == "som.fit" | rowname1 == "om.fit" |rowname1 == "pb.fit"
)

ggplot(data = aic.results.hyp, aes(x = rowname1, y = AIC)) +
  geom_bar(stat = "identity", fill = c("#183442", "#91A5AE","#A6A3AA"), color = "black") +
  scale_x_discrete(limits = c("som.fit", "om.fit", "pb.fit"), 
                   labels = c("Shifting optimal \n margin model", 
                              "Optimal margin \n model", 
                              "Positivity bias \n model")) +
  #  scale_y_reverse() +
  #  scale_y_continuous(limits = c(23500, 24200)) +
  coord_cartesian(ylim = c(23800, 24100))+
  theme_classic(base_size = 20) +
  labs(x = "", y = "Akaike's Information Criterion") +
  theme(axis.text  = element_text(colour = "black", size = 20), 
        axis.title = element_text(colour = "black", size = 23), 
        legend.position = "none") 


# 2.3.2.) Calculate Akaike weights

# Make a vector from AICs 
aics.vector <- pull(aic.results, AIC)

# Calculate Akaike weights and attach to results of model comparisons
weights <- akaike.weights(aics.vector)
aic.results$weights <- weights$weights


# Plot result
ggplot(data = aic.results, aes(x = rowname1, y = weights)) +
  geom_bar(stat = "identity", fill = "plum") +
  scale_x_discrete(limits=c("som.fit", "fu.fit", "om.fit", "s6.fit", "s4.fit", 
                            "pb.fit", "s5.fit", "s2.fit", "s3.fit", "nu.fit",
                            "s7.fit", "s1.fit")) +
  #  scale_y_reverse() +
  #  scale_y_continuous(limits = c(23500, 24200)) +
  #coord_cartesian(ylim = c(23500, 24200))+
  theme_classic()


# -----------------------------
# Plot RSA result full model
# -----------------------------

# refit model using RSA function
r.fu  <- RSA(ls  ~ sa.s*ca.s, df, model = "full")
summary(r.fu)

# Plot models and save in directory

png("model.png", width = 6, height = 6, units = 'in', res = 600)
plot(r.fu, 
     axes = c("LOC", "LOIC", "PA1"),
     axesStyles = list(LOC  = list(lty = "solid", lwd = 2, col = "black"),
                       LOIC = list(lty = "solid", lwd = 2, col = "black"),
                       PA1  = list(lty = "dashed", lwd = 2, col = "grey40")),
     xlab = "Subjective age", ylab = "Actual age", zlab ="Life satisfaction",
     cex.main = 1.3, cex.tickLabel = 1.25, cex.axesLabel = 1.25,
     points = list(show=FALSE), hull = F, legend = F,
     project = c("LOC", "LOIC", "PA1", "contour"), param = F,
     pal = colorRampPalette(c("#4c3633", "#795a70","#8f86a6",
                              "#a3b5c1", "#c8dcd2", "#EEF4F1"))(15))
dev.off()


