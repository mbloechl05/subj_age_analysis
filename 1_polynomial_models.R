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


# 1.) Positivity bias model
# --------------------------------

pb.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 <  0
b2 == 0
b3 == 0
b4 == 0
b5 == 0"


# 2.) Optimal margin models
# --------------------------------

# 2.1.) Optimal margin flat ridge

omfr.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b1 == -b2
b3 <  0
b3 == b5
b3 + b4 + b5 == 0 
C: = b1/(2*b3)"


# 2.2.) Optimal margin rising ridge

omrr.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2
# parameter constraints
b3 == b5
b3 + b4 + b5 == 0
C := (b2-b1)/(4*b3)"


# 3.) Shifting optimal margin models
# ------------------------------------

# 3.1.) Shifting optimal margin flat ridge

### Flat ridge down

somfr_d.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2
# parameter constraints
b1 == (b2*b4)/2*b5
b3 < -0.000001
b5 < -0.000001
b4^2 == 4*b5*b3"


### Flat ridge up

somfr_u.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2
# parameter constraints
b1 == (b2*b4)/2*b5
b3 > 0.000001
b5 > 0.000001
b4^2 == 4*b5*b3"


# 3.2.) Shifting optimal margin rising ridge

### Rising ridge down

somrr_d.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b3 < -0.000001
b5 < -0.000001
b4^2 == 4*b3*b5
C:= -(2*b1*b5 + b2*b4)/(4*b4*b5)
S:= -b4/(2*b5)
bM := b1/S + b2"


### Rising ridge up

somrr_u.model <- "
ls ~ 1 + b1*sa.s + b2*ca.s + b3*sa.s2 + b4*sa.ca + b5*ca.s2 
# parameter constraints
b3 > 0.000001
b5 > 0.000001
b4^2 == 4*b3*b5
C:= -(2*b1*b5 + b2*b4)/(4*b4*b5)
S:= -b4/(2*b5)
bM := b1/S + b2"


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

