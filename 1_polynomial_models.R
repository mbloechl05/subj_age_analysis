# ===============================================================
# Models for 
# "Psychological benefits of subjective age bias" 
# ==============================================================

# This script contains code to define all polynomial regression models 
# for the analysis of the subjective age bias. 
# It needs to be run before running the main or supplementary analyses
# since the models are called in these scripts. 

# ------------------------------------------
# Specify all polynomial regression models
# ------------------------------------------


# 0.) Null model 
# ---------------------

nu.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 == 0
b2 == 0
b3 == 0
b4 == 0
b5 == 0"


# 1.) Positivity bias models
# --------------------------------

# Model 1A: Low subjective age-only model

a1.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 <  0
b2 == 0
b3 == 0
b4 == 0
b5 == 0"


# Model S1B: Low subj. age and high chronological age model

b1.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 <  0
b2 >  0
b3 == 0
b4 == 0
b5 == 0"


# 2.) Optimal margin models
# --------------------------------

# Model 2A: Optimal margin flat ridge

a2.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == -b2
b3 <   0
b3 ==  b5
b3 + b4 + b5 == 0 

# surface parameters
C := b1/(2*b3)"


# Model 2B: Optimal margin rising ridge

b2.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2

# parameter constraints
b3 == b5
b3 + b4 + b5 == 0

# surface parameters
C  := (b1-b2)/(4*b3)
bM := b1+b2"


# 3.) Increasing optimal margin models
# ------------------------------------

# We only provide code for the models with the down-ward tilting ridge here. 
# Although we also fitted the model with the upward tilting ridge, these models
# had a poorer fit and were excluded in further analyses.


# Model 3A: Increasing optimal margin flat ridge

# Note that we had to give this model starting values for free parameters
# to find an aedequate solution.

a3.model <- "
ls ~ 1 + b1*sa.s + start(-0.2)*sa.s + b2*ca.s + start(0.2)*ca.s + b3*sa.s2 + 
start(-0.10)*sa.s2 + b4*sa.ca + start(0.2)*sa.ca + b5*ca.s2 + start(-0.1)*ca.s2

# parameter constraints
b1   == (b2*b4)/(2*b5)
b3    < -0.000001
b5    < -0.000001
b4^2 == 4*b5*b3

# surface parameters
C := -0.5*(b2/b5)"


# Model 3B: Increasing optimal margin rising ridge

b3.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b3   <  -0.000001
b5   <  -0.000001
b4^2 == 4*b3*b5

# surface parameters
C   := -(2*b1*b5 + b2*b4)/(4*b4*b5)
S   := -b4/(2*b5)
Sht := S-1 # for testing S against 1
bM  := b1/S + b2
X0  := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)
Y0  := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)
p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4
p10 := Y0 - p11*X0"



# 4.) Full model
# -------------------

fu.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2"



# 5.) Supplementary models
# ----------------------------

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


# Model S6: Congruency model
# positive effect of correct self-evaluation

s6.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + + b3*sa.s2 + b4*sa.ca + b5*ca.s2 

# parameter constraints
b1 == 0
b2 == b1
b3 <  0
b4 == -2*b3
b5 == b3"

