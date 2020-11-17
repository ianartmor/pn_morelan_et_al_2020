# This file contains notes for running a repeated measures analysis on
# JFW1 paper's diversity, climate, and must data. 

# This script was written by Amisha Poret-Peterson and Kerri Steenwerth

library(lme4)
library(lsmeans) # notice that this is now emmeans package 
library(ARTool)
library(emmeans)
library(rcompanion)
library(dplyr)
library(broom)
library(knitr)
library(gtools)
library(emmeans)
library(boxcoxmix)

load("Files/JFW_mixed_model.RData")

##############################################################################

citation("lme4")
#To cite lme4 in publications use:
#  
#  Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using
#  lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.



##############################################################################
# Checking residuals 

bact_rich_reschk <- lmer(Bact_Rich ~ AVA*Year + (1 | Vineyard), data)
plot(bact_rich_reschk) 
qqnorm(residuals(bact_rich_reschk))
qqnorm(ranef(bact_rich_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(bact_rich_reschk))  

fung_rich_reschk <- lmer(Fung_Rich ~ AVA*Year + (1 | Vineyard), data)
plot(fung_rich_reschk) 
qqnorm(residuals(fung_rich_reschk))
qqnorm(ranef(fung_rich_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(fung_rich_reschk)) 

# log transform and retest 
data$log_Fung_Rich <- log(data$Fung_Rich)
log_Fung_Rich_reschk <- lmer(log_Fung_Rich ~ AVA*Year + (1 | Vineyard), data)
plot(log_Fung_Rich_reschk) 
qqnorm(residuals(log_Fung_Rich_reschk))
qqnorm(ranef(log_Fung_Rich_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(log_Fung_Rich_reschk)) 

bact_Shannon_reschk <- lmer(Bact_Shannon ~ AVA*Year + (1 | Vineyard), data)
plot(bact_Shannon_reschk) 
qqnorm(residuals(bact_Shannon_reschk))
qqnorm(ranef(bact_Shannon_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(bact_Shannon_reschk))

# log transform and retest 
data$log_Bact_Shannon <- log(data$Bact_Shannon)
log_bact_Shannon_reschk <- lmer(log_Bact_Shannon ~ AVA*Year + (1 | Vineyard), data)
plot(log_bact_Shannon_reschk) 
qqnorm(residuals(log_bact_Shannon_reschk))
qqnorm(ranef(log_bact_Shannon_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(log_bact_Shannon_reschk)) 

fung_Shannon_reschk <- lmer(Fung_Shannon ~ AVA*Year + (1 | Vineyard), data)
plot(fung_Shannon_reschk) 
qqnorm(residuals(fung_Shannon_reschk))
qqnorm(ranef(fung_Shannon_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(fung_Shannon_reschk))

Bact_InvSimpson_reschk <- lmer(Bact_InvSimp ~ AVA*Year + (1 | Vineyard), data)
plot(Bact_InvSimpson_reschk) 
qqnorm(residuals(Bact_InvSimpson_reschk))
qqnorm(ranef(Bact_InvSimpson_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(Bact_InvSimpson_reschk))  

# log transform and retest 
data$log_Bact_InvSimp <- log(data$Bact_InvSimp)
log_Bact_InvSimp_reschk <- lmer(log_Bact_InvSimp ~ AVA*Year + (1 | Vineyard), data)
plot(log_Bact_InvSimp_reschk) 
qqnorm(residuals(log_Bact_InvSimp_reschk))
qqnorm(ranef(log_Bact_InvSimp_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(log_Bact_InvSimp_reschk)) 

fung_InvSimpson_reschk <- lmer(Fung_InvSimp ~ AVA*Year + (1 | Vineyard), data)
plot(fung_InvSimpson_reschk) 
qqnorm(residuals(fung_InvSimpson_reschk))
qqnorm(ranef(fung_InvSimpson_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(fung_InvSimpson_reschk))  

must_TSS_reschk <- lmer(Must_TSS ~ AVA*Year + (1 | Vineyard), data)
plot(must_TSS_reschk)
qqnorm(residuals(must_TSS_reschk))
qqnorm(ranef(must_TSS_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(must_TSS_reschk)) 

must_pH_reschk <- lmer(Must_pH ~ AVA*Year + (1 | Vineyard), data)
plot(must_pH_reschk) 
qqnorm(residuals(must_pH_reschk))
qqnorm(ranef(must_pH_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(must_pH_reschk)) 

must_acidity_reschk <- lmer(Must_acidity ~ AVA*Year + (1 | Vineyard), data)
plot(must_acidity_reschk) 
qqnorm(residuals(must_acidity_reschk))
qqnorm(ranef(must_acidity_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(must_acidity_reschk))  

GDD_reschk <- lmer(GDD ~ AVA*Year + (1 | Vineyard), data)
plot(GDD_reschk) 
qqnorm(residuals(GDD_reschk))
qqnorm(ranef(GDD_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(GDD_reschk)) 

Precip_reschk <- lmer(Precip ~ AVA*Year + (1 | Vineyard), data)
plot(Precip_reschk) 
qqnorm(residuals(Precip_reschk))
qqnorm(ranef(Precip_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(Precip_reschk)) 

# log transform and retest 
data$log_Precip <- log(data$Precip)
log_Precip_reschk <- lmer(log_Precip ~ AVA*Year + (1 | Vineyard), data)
plot(log_Precip_reschk) 
qqnorm(residuals(log_Precip_reschk))
qqnorm(ranef(log_Precip_reschk)$Vineyard[,1])
plotNormalHistogram(residuals(log_Precip_reschk))




##############################################################################
# Repeated measures tests 

######################################## Bact_Richness




bact_rich_full_model <- lmer(Bact_Rich ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
bact_rich_reduced_interaction <- lmer(Bact_Rich ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)
# received "boundary (singular) fit: see ?isSingular" error 
bact_rich_AVA <- lmer(Bact_Rich ~ Year + (1 | Vineyard), data, REML=FALSE) # null for AVA
bact_rich_Year <- lmer(Bact_Rich ~ AVA + (1 | Vineyard), data, REML=FALSE) # null for Year 
# received "boundary (singular) fit: see ?isSingular" error 


bact_rich_AVA_Year_test <- anova(bact_rich_reduced_interaction, bact_rich_full_model)
bact_rich_AVA_Year_test 

# Data: data
# Models:  bact_rich_reduced_interaction: Bact_Rich ~ AVA + Year + (1 | Vineyard)
#          bact_rich_full_model: Bact_Rich ~ AVA * Year + (1 | Vineyard)
#                                Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# bact_rich_reduced_interaction  8 382.24 392.90 -183.12   366.24                            
# bact_rich_full_model          12 373.77 389.76 -174.89   349.77 16.469      4    0.00245 **
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#### THM: significant interaction betweeen AVA and year (chisq = 16.47, chisq df = 4, p=0.002)

bact_rich_AVA_test <- anova(bact_rich_AVA, bact_rich_reduced_interaction)
bact_rich_AVA_test



# Data: data
# Models:   bact_rich_AVA: Bact_Rich ~ Year + (1 | Vineyard)
#           bact_rich_reduced_interaction: Bact_Rich ~ AVA + Year + (1 | Vineyard)
#                               Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# bact_rich_AVA                  4 390.05 395.38 -191.03   382.05                            
# bact_rich_reduced_interaction  8 382.24 392.90 -183.12   366.24 15.811      4   0.003284 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#### THM: bacterial richness differs significantly by AVA (chisq 15.81, chisq df = 4, p = 0.003)

bact_rich_Year_test <- anova(bact_rich_Year, bact_rich_reduced_interaction)
bact_rich_Year_test
# Data: data
# Models:  bact_rich_Year: Bact_Rich ~ AVA + (1 | Vineyard)
#          bact_rich_reduced_interaction: Bact_Rich ~ AVA + Year + (1 | Vineyard)
#                               Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# bact_rich_Year                 7 382.40 391.73 -184.20   368.40                         
# bact_rich_reduced_interaction  8 382.24 392.90 -183.12   366.24 2.1609      1     0.1416
#### THM: bacterial richness does not differ significanty by year (chisq 2.16, chisq df = 1, p = 0.142)

# bacterial richness posthoc 
bact_rich_AVA_posthoc <- lsmeans(bact_rich_full_model, pairwise ~ AVA | Year) 
bact_rich_AVA_posthoc$contrasts
# Year = 2016:
#   contrast                          estimate  SE   df t.ratio p.value
# Mendocino - Monterey                -747.4 198 43.5 -3.766  0.0043 
# Mendocino - Santa Barbara           -391.0 188 43.5 -2.080  0.2471 
# Mendocino - Sonoma                  -254.9 177 43.5 -1.440  0.6055 
# Mendocino - Willamette Valley       -405.4 233 43.6 -1.736  0.4232 
# Monterey - Santa Barbara             356.3 142 43.5  2.506  0.1079 
# Monterey - Sonoma                    492.5 127 43.5  3.873  0.0031 
# Monterey - Willamette Valley         342.0 198 43.5  1.723  0.4310 
# Santa Barbara - Sonoma               136.2 110 43.5  1.236  0.7304 
# Santa Barbara - Willamette Valley    -14.4 188 43.5 -0.076  1.0000 
# Sonoma - Willamette Valley          -150.5 177 43.5 -0.851  0.9130 

# Year = 2017:
#   contrast                          estimate  SE   df t.ratio p.value
# Mendocino - Monterey                -358.5 156 43.5 -2.302  0.1639 
# Mendocino - Santa Barbara            154.0 142 43.5  1.083  0.8142 
# Mendocino - Sonoma                  -245.2 127 43.5 -1.928  0.3185 
# Mendocino - Willamette Valley       -194.5 156 43.5 -1.249  0.7232 
# Monterey - Santa Barbara             512.5 142 43.5  3.604  0.0068 
# Monterey - Sonoma                    113.3 127 43.5  0.891  0.8986 
# Monterey - Willamette Valley         164.0 156 43.5  1.053  0.8291 
# Santa Barbara - Sonoma              -399.2 110 43.5 -3.624  0.0064 
# Santa Barbara - Willamette Valley   -348.5 142 43.5 -2.451  0.1212 
# Sonoma - Willamette Valley            50.7 127 43.5  0.398  0.9945 

# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 5 estimates 
# lsmeans(bact_rich_full_model, pairwise ~ Year | AVA) 

bact_rich_Year_posthoc <- lsmeans(bact_rich_full_model, pairwise ~ Year | AVA) 
bact_rich_Year_posthoc$contrasts
# AVA = Mendocino:
#  contrast    estimate  SE   df t.ratio p.value
# 2016 - 2017   -105.4 196 28.8 -0.538  0.5948 

# AVA = Monterey:
#   contrast    estimate  SE   df t.ratio p.value
# 2016 - 2017    283.5 152 20.5  1.860  0.0774 

# AVA = Santa Barbara:
#  contrast    estimate  SE   df t.ratio p.value
# 2016 - 2017    439.7 124 20.5  3.532  0.0020 

# AVA = Sonoma:
#  contrast    estimate  SE   df t.ratio p.value
# 2016 - 2017    -95.7  88 20.5 -1.087  0.2897 

# AVA = Willamette Valley:
#  contrast    estimate  SE   df t.ratio p.value
# 2016 - 2017    105.5 196 28.8  0.539  0.5942 

# Degrees-of-freedom method: kenward-roger



#############################################log_Bact_Shannon 

log_bact_Shannon_full_model <- lmer(log_Bact_Shannon ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
log_bact_Shannon_reduced_interaction <- lmer(log_Bact_Shannon ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)

log_bact_Shannon_AVA <- lmer(log_Bact_Shannon ~ Year + (1 | Vineyard), data, REML=FALSE) # null for AVA
log_bact_Shannon_Year <- lmer(log_Bact_Shannon ~ AVA + (1 | Vineyard), data, REML=FALSE) # null for Year 


log_bact_Shannon_AVA_Year_test <- anova(log_bact_Shannon_reduced_interaction, log_bact_Shannon_full_model)
log_bact_Shannon_AVA_Year_test 

log_bact_Shannon_AVA_test <- anova(log_bact_Shannon_AVA, log_bact_Shannon_reduced_interaction)
log_bact_Shannon_AVA_test

log_bact_Shannon_Year_test <- anova(log_bact_Shannon_Year, log_bact_Shannon_reduced_interaction)
log_bact_Shannon_Year_test

log_bact_Shannon_full_model <- lmer(log_Bact_Shannon ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
log_bact_Shannon_reduced_interaction <- lmer(log_Bact_Shannon ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)

log_bact_Shannon_AVA <- lmer(log_Bact_Shannon ~ Year + (1 | Vineyard), data, REML=FALSE) # null for AVA
log_bact_Shannon_Year <- lmer(log_Bact_Shannon ~ AVA + (1 | Vineyard), data, REML=FALSE) # null for Year 



log_bact_Shannon_AVA_Year_test <- anova(log_bact_Shannon_reduced_interaction, log_bact_Shannon_full_model)
log_bact_Shannon_AVA_Year_test 

log_bact_Shannon_AVA_test <- anova(log_bact_Shannon_AVA, log_bact_Shannon_reduced_interaction)
log_bact_Shannon_AVA_test


# bacterial Shannonness posthoc 
log_bact_Shannon_AVA_posthoc <- lsmeans(log_bact_Shannon_full_model, pairwise ~ AVA | Year) 
log_bact_Shannon_AVA_posthoc$contrasts

log_bact_Shannon_Year_posthoc <- lsmeans(log_bact_Shannon_full_model, pairwise ~ Year | AVA) 
log_bact_Shannon_Year_posthoc$contrasts
log_bact_Shannon_AVA_posthoc <- lsmeans(log_bact_Shannon_full_model, pairwise ~ AVA | Year) 
log_bact_Shannon_AVA_posthoc$contrasts


log_bact_Shannon_Year_posthoc <- lsmeans(log_bact_Shannon_full_model, pairwise ~ Year | AVA) 
log_bact_Shannon_Year_posthoc$contrasts






############################################# Fungal Shannon

fung_Shannon_full_model <- lmer(Fung_Shannon ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
# received "boundary (singular) fit: see ?isSingular" error
fung_Shannon_reduced_interaction <- lmer(Fung_Shannon ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)
# received "boundary (singular) fit: see ?isSingular" error 
fung_Shannon_AVA <- lmer(Fung_Shannon ~ Year + (1 | Vineyard), data, REML=FALSE) # null for AVA
# boundary (singular) fit: see ?isSingular
fung_Shannon_Year <- lmer(Fung_Shannon ~ AVA + (1 | Vineyard), data, REML=FALSE) # null for Year 
# received "boundary (singular) fit: see ?isSingular" error 

fung_Shannon_AVA_Year_test <- anova(fung_Shannon_reduced_interaction, fung_Shannon_full_model)
fung_Shannon_AVA_Year_test 
# Data: data
# Models:   fung_Shannon_reduced_interaction: Fung_Shannon ~ AVA + Year + (1 | Vineyard)
#           fung_Shannon_full_model: Fung_Shannon ~ AVA * Year + (1 | Vineyard)
#                                  Df    AIC     BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# fung_Shannon_reduced_interaction  8 90.439 101.097 -37.219   74.439                             
# fung_Shannon_full_model          12 67.236  83.222 -21.618   43.236 31.203      4  2.783e-06 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#### THM: significant interaction b/w AVA and year for fungal Shannon (chisq = 31.20, chisq df=4, p=2.783e-06 or P < 0.001)

fung_Shannon_AVA_test <- anova(fung_Shannon_AVA, fung_Shannon_reduced_interaction)
fung_Shannon_AVA_test
# Data: data
# Models:  fung_Shannon_AVA: Fung_Shannon ~ Year + (1 | Vineyard)
#          fung_Shannon_reduced_interaction: Fung_Shannon ~ AVA + Year + (1 | Vineyard)
#                                  Df    AIC     BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# fung_Shannon_AVA                  4 87.888  93.216 -39.944   79.888                         
# fung_Shannon_reduced_interaction  8 90.439 101.097 -37.219   74.439 5.4486      4     0.2443
#### THM: Fungal Shannon does not differ significantly by AVA (chisq 5.45, chisq df = 4, p = 0.244)

fung_Shannon_Year_test <- anova(fung_Shannon_Year, fung_Shannon_reduced_interaction)
fung_Shannon_Year_test
# Data: data
# Models:  fung_Shannon_Year: Fung_Shannon ~ AVA + (1 | Vineyard)
#          fung_Shannon_reduced_interaction: Fung_Shannon ~ AVA + Year + (1 | Vineyard)
#                                  Df    AIC     BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# fung_Shannon_Year                 7 88.613  97.938 -37.306   74.613                         
# fung_Shannon_reduced_interaction  8 90.439 101.097 -37.219   74.439 0.1736      1     0.6769
#### THM: fungal shannon does not differ signicantly by year (chisq 0.17, chisq df = 1, p = 0.677)

# fungal Shannonness posthoc 
fung_Shannon_AVA_posthoc <- lsmeans(fung_Shannon_full_model, pairwise ~ AVA | Year) 
fung_Shannon_AVA_posthoc$contrasts

fung_Shannon_Year_posthoc <- lsmeans(fung_Shannon_full_model, pairwise ~ Year | AVA) 
fung_Shannon_Year_posthoc$contrasts



########################################  Fungal Inverse Simpson

fung_InvSimp_full_model <- lmer(Fung_InvSimp ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
# boundary (singular) fit: see ?isSingular
fung_InvSimp_reduced_interaction <- lmer(Fung_InvSimp ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)
# received "boundary (singular) fit: see ?isSingular" error 
fung_InvSimp_AVA <- lmer(Fung_InvSimp ~ Year + (1 | Vineyard), data, REML=FALSE) # null for AVA
# boundary (singular) fit: see ?isSingular
fung_InvSimp_Year <- lmer(Fung_InvSimp ~ AVA + (1 | Vineyard), data, REML=FALSE) # null for Year 
# received "boundary (singular) fit: see ?isSingular" error 

fung_InvSimp_AVA_Year_test <- anova(fung_InvSimp_reduced_interaction, fung_InvSimp_full_model)
fung_InvSimp_AVA_Year_test 
# Data: data
# Models:    fung_InvSimp_reduced_interaction: Fung_InvSimp ~ AVA + Year + (1 | Vineyard)
#            fung_InvSimp_full_model: Fung_InvSimp ~ AVA * Year + (1 | Vineyard)
#                                   Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
# fung_InvSimp_reduced_interaction  8 82.549 93.207 -33.275   66.549                            
# fung_InvSimp_full_model          12 64.029 80.016 -20.015   40.029 26.52      4  2.485e-05 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#### THM: significant interaction b/w AVA and year (chisq 26.52, chisq df = 4, p =2.485e-05 or P < 0.001)

fung_InvSimp_AVA_test <- anova(fung_InvSimp_AVA, fung_InvSimp_reduced_interaction)
fung_InvSimp_AVA_test
# Data: data
# Models:     fung_InvSimp_AVA: Fung_InvSimp ~ Year + (1 | Vineyard)
#             fung_InvSimp_reduced_interaction: Fung_InvSimp ~ AVA + Year + (1 | Vineyard)
#                                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# fung_InvSimp_AVA                  4 79.761 85.090 -35.881   71.761                         
# fung_InvSimp_reduced_interaction  8 82.549 93.207 -33.275   66.549 5.2118      4     0.2662
#### THM: Fungal Inverse Simpson - no sig diffs in AVA (chisq 5.21, chisq df = 4, p = 0.266)

fung_InvSimp_Year_test <- anova(fung_InvSimp_Year, fung_InvSimp_reduced_interaction)
fung_InvSimp_Year_test
# Data: data
# Models:       fung_InvSimp_Year: Fung_InvSimp ~ AVA + (1 | Vineyard)
#               fung_InvSimp_reduced_interaction: Fung_InvSimp ~ AVA + Year + (1 | Vineyard)
#                                  Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# fung_InvSimp_Year                 7 80.601 89.926 -33.300   66.601                         
# fung_InvSimp_reduced_interaction  8 82.549 93.207 -33.275   66.549 0.0511      1     0.8211
#### THM: no sig diffs b/w years (chisq = 0.05, chisq df = 1, p = 0.821)

fung_InvSimp_AVA_posthoc <- lsmeans(fung_InvSimp_full_model, pairwise ~ AVA | Year) 
fung_InvSimp_AVA_posthoc$contrasts

fung_InvSimp_Year_posthoc <- lsmeans(fung_InvSimp_full_model, pairwise ~ Year | AVA) 
fung_InvSimp_Year_posthoc$contrasts



########################################  GDD 

GDD_full_model <- lmer(GDD ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
GDD_reduced_interaction <- lmer(GDD ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)
GDD_AVA <- lmer(GDD ~ Year + (1 | Vineyard), data, REML=FALSE) # null for AVA
GDD_Year <- lmer(GDD ~ AVA + (1 | Vineyard), data, REML=FALSE) # null for Year 

GDD_AVA_Year_test <- anova(GDD_reduced_interaction, GDD_full_model)
GDD_AVA_Year_test 
# Data: data
# Models:   GDD_reduced_interaction: GDD ~ AVA + Year + (1 | Vineyard)
#           GDD_full_model: GDD ~ AVA * Year + (1 | Vineyard)
#                           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# GDD_reduced_interaction  8 340.02 350.68 -162.01   324.02                             
# GDD_full_model          12 316.71 332.69 -146.35   292.71 31.316      4   2.64e-06 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##### THM: GDD has a significant interaction b/w AVA and year (chisq 31.31, chisq df = 4, p =2.64e-06, P < 0.001)

GDD_AVA_test <- anova(GDD_AVA, GDD_reduced_interaction)
GDD_AVA_test
# Data: data
# Models:       GDD_AVA: GDD ~ Year + (1 | Vineyard)
#               GDD_reduced_interaction: GDD ~ AVA + Year + (1 | Vineyard)
#                          Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# GDD_AVA                  4 346.13 351.46 -169.07   338.13                            
# GDD_reduced_interaction  8 340.02 350.68 -162.01   324.02 14.113      4   0.006943 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##### THM: GDD sig diff by AVA (chisq 14.11, chisq df = 4, p = 0.007)

GDD_Year_test <- anova(GDD_Year, GDD_reduced_interaction)
GDD_Year_test
# Data: data
# Models:   GDD_Year: GDD ~ AVA + (1 | Vineyard)
#           GDD_reduced_interaction: GDD ~ AVA + Year + (1 | Vineyard)
#                          Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
# GDD_Year                 7 355.57 364.90 -170.79   341.57                            
# GDD_reduced_interaction  8 340.02 350.68 -162.01   324.02 17.55      1  2.799e-05 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##### THM: GDD significantly differs by year (chisq 17.55, chisq df=1, P = 2.799e-05 or P < 0.001)

GDD_AVA_posthoc <- lsmeans(GDD_full_model, pairwise ~ AVA | Year) 
GDD_AVA_posthoc$contrasts

GDD_Year_posthoc <- lsmeans(GDD_full_model, pairwise ~ Year | AVA) 
GDD_Year_posthoc$contrasts


########################################  Must_TSS 

Must_TSS_full_model <- lmer(Must_TSS ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
Must_TSS_reduced_interaction <- lmer(Must_TSS ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)
# boundary (singular) fit: see ?isSingular
Must_TSS_AVA <- lmer(Must_TSS ~ Year + (1 | Vineyard), data, REML=FALSE) # null for AVA
Must_TSS_Year <- lmer(Must_TSS ~ AVA + (1 | Vineyard), data, REML=FALSE) # null for Year 
# boundary (singular) fit: see ?isSingular

Must_TSS_AVA_Year_test <- anova(Must_TSS_reduced_interaction, Must_TSS_full_model)
Must_TSS_AVA_Year_test 
# Data: data
# Models:  Must_TSS_reduced_interaction: Must_TSS ~ AVA + Year + (1 | Vineyard)
#          Must_TSS_full_model: Must_TSS ~ AVA * Year + (1 | Vineyard)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Must_TSS_reduced_interaction  8 92.400 103.06  -38.2   76.400                         
# Must_TSS_full_model          12 94.601 110.59  -35.3   70.601 5.7991      4     0.2147
##### THM: Must_TSS does not have a significant interction b/w AVA and year (chisq = 5.80, chisq df = 4, p=0.215)

Must_TSS_AVA_test <- anova(Must_TSS_AVA, Must_TSS_reduced_interaction)
Must_TSS_AVA_test
# Data: data
# Models:       Must_TSS_AVA: Must_TSS ~ Year + (1 | Vineyard)
#               Must_TSS_reduced_interaction: Must_TSS ~ AVA + Year + (1 | Vineyard)
#                              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# Must_TSS_AVA                  4 94.146  99.475 -43.073   86.146                           
# Must_TSS_reduced_interaction  8 92.400 103.057 -38.200   76.400 9.7459      4    0.04493 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##### THM: Must_TSS sig diff by AVA (chisq 9.75, chisq df = 4, p = 0.045)

Must_TSS_Year_test <- anova(Must_TSS_Year, Must_TSS_reduced_interaction)
Must_TSS_Year_test
# Data: data
# Models:   Must_TSS_Year: Must_TSS ~ AVA + (1 | Vineyard)
#           Must_TSS_reduced_interaction: Must_TSS ~ AVA + Year + (1 | Vineyard)
#                              Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
# Must_TSS_Year                 7 90.406  99.731 -38.203   76.406                        
# Must_TSS_reduced_interaction  8 92.400 103.057 -38.200   76.400 0.006      1     0.9383
##### THM: Must_TSS does not differ by year (chisq 0.006, chisq df=1, P = 0.9383)

Must_TSS_AVA_posthoc <- lsmeans(Must_TSS_full_model, pairwise ~ AVA | Year) 
Must_TSS_AVA_posthoc$contrasts

Must_TSS_Year_posthoc <- lsmeans(Must_TSS_full_model, pairwise ~ Year | AVA) 
Must_TSS_Year_posthoc$contrasts

Must_TSS_AVA_posthoc <- lsmeans(Must_TSS_reduced_interaction, pairwise ~ AVA ) 
Must_TSS_AVA_posthoc$contrasts



########################################  Must_pH 

Must_pH_full_model <- lmer(Must_pH ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
Must_pH_reduced_interaction <- lmer(Must_pH ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)
Must_pH_AVA <- lmer(Must_pH ~ Year + (1 | Vineyard), data, REML=FALSE) # null for AVA
Must_pH_Year <- lmer(Must_pH ~ AVA + (1 | Vineyard), data, REML=FALSE) # null for Year 

Must_pH_AVA_Year_test <- anova(Must_pH_reduced_interaction, Must_pH_full_model)
Must_pH_AVA_Year_test 
# Data: data
# Models:  Must_pH_reduced_interaction: Must_pH ~ AVA + Year + (1 | Vineyard)
#          Must_pH_full_model: Must_pH ~ AVA * Year + (1 | Vineyard)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Must_pH_reduced_interaction  8 -53.980 -43.322 34.990  -69.980                         
# Must_pH_full_model          12 -53.142 -37.156 38.571  -77.142 7.1623      4     0.1276
##### THM: no sig interaction b/w AVA and year for Must_pH (chisq 7.16, chisq df = 4, p = 0.128)

Must_pH_AVA_test <- anova(Must_pH_AVA, Must_pH_reduced_interaction)
Must_pH_AVA_test
# Data: data
# Models:       Must_pH_AVA: Must_pH ~ Year + (1 | Vineyard)
#               Must_pH_reduced_interaction: Must_pH ~ AVA + Year + (1 | Vineyard)
#                             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# Must_pH_AVA                  4 -52.002 -46.673 30.001  -60.002                           
# Must_pH_reduced_interaction  8 -53.980 -43.322 34.990  -69.980 9.9774      4    0.04081 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##### THM: Must_pH sig diff by AVA (chisq 9.98, chisq df = 4, p = 0.04)

Must_pH_Year_test <- anova(Must_pH_Year, Must_pH_reduced_interaction)
Must_pH_Year_test
# Data: data
# Models:   Must_pH_Year: Must_pH ~ AVA + (1 | Vineyard)
#           Must_pH_reduced_interaction: Must_pH ~ AVA + Year + (1 | Vineyard)
#                             Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
# Must_pH_Year                 7 -52.504 -43.179 33.252  -66.504                           
# Must_pH_reduced_interaction  8 -53.980 -43.322 34.990  -69.980 3.4757      1    0.06228 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##### THM: Must_pH does not differ by year or is marginally sig (chisq 3.48, chisq df=1, P = 0.06228)

Must_pH_AVA_posthoc <- lsmeans(Must_pH_full_model, pairwise ~ AVA | Year) 
Must_pH_AVA_posthoc$contrasts

Must_pH_Year_posthoc <- lsmeans(Must_pH_full_model, pairwise ~ Year | AVA) 
Must_pH_Year_posthoc$contrasts

Must_pH_AVA_posthoc <- lsmeans(Must_pH_reduced_interaction, pairwise ~ AVA ) 
Must_pH_AVA_posthoc$contrasts


########################################  Must_acidity - parametric 

Must_acidity_full_model <- lmer(Must_acidity ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
Must_acidity_reduced_interaction <- lmer(Must_acidity ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)
Must_acidity_AVA <- lmer(Must_acidity ~ Year + (1 | Vineyard), data, REML=FALSE) # null for AVA
Must_acidity_Year <- lmer(Must_acidity ~ AVA + (1 | Vineyard), data, REML=FALSE) # null for Year 

Must_acidity_AVA_Year_test <- anova(Must_acidity_reduced_interaction, Must_acidity_full_model)
Must_acidity_AVA_Year_test 
# Data: data
# Models:  Must_acidity_reduced_interaction: Must_acidity ~ AVA + Year + (1 | Vineyard)
#          Must_acidity_full_model: Must_acidity ~ AVA * Year + (1 | Vineyard)
#                                  Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# Must_acidity_reduced_interaction  8 47.628 58.285 -15.814   31.628                         
# Must_acidity_full_model          12 54.041 70.028 -15.021   30.041 1.5865      4     0.8112
##### THM: No sig interaction b/w AVA and year (chisq 1.59, chisq df = 4, p = 0.811)

Must_acidity_AVA_test <- anova(Must_acidity_AVA, Must_acidity_reduced_interaction)
Must_acidity_AVA_test
# Data: data
# Models:       Must_acidity_AVA: Must_acidity ~ Year + (1 | Vineyard)
#               Must_acidity_reduced_interaction: Must_acidity ~ AVA + Year + (1 | Vineyard)
#                                   Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# Must_acidity_AVA                  4 51.455 56.784 -21.727   43.455                           
# Must_acidity_reduced_interaction  8 47.628 58.285 -15.814   31.628 11.828      4    0.01868 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##### THM: Must_acidity sig diff by AVA (chisq 11.82, chisq df = 4, p = 0.01868) 

Must_acidity_Year_test <- anova(Must_acidity_Year, Must_acidity_reduced_interaction)
Must_acidity_Year_test
# Data: data
# Models:   Must_acidity_Year: Must_acidity ~ AVA + (1 | Vineyard)
#           Must_acidity_reduced_interaction: Must_acidity ~ AVA + Year + (1 | Vineyard)
#                                 Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
# Must_acidity_Year                 7 45.639 54.964 -15.819   31.639                        
# Must_acidity_reduced_interaction  8 47.628 58.285 -15.814   31.628 0.011      1     0.9164
##### THM: Must_acidity does not differ by year (chisq 0.011, chisq df=1, P = 0.9164)

Must_acidity_AVA_posthoc <- lsmeans(Must_acidity_full_model, pairwise ~ AVA | Year) 
Must_acidity_AVA_posthoc$contrasts

Must_acidity_Year_posthoc <- lsmeans(Must_acidity_full_model, pairwise ~ Year | AVA) 
Must_acidity_Year_posthoc$contrasts

Must_acidity_AVA_posthoc <- lsmeans(Must_acidity_reduced_interaction, pairwise ~ AVA ) 
Must_acidity_AVA_posthoc$contrasts







############################

########################################  Precip Logged for normality

log_Precip_full_model <- lmer(log_Precip ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
log_Precip_reduced_interaction <- lmer(log_Precip ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)

log_Precip_AVA <- lmer(log_Precip ~ Year + (1 | Vineyard), data, REML=FALSE) 
log_Precip_Year <- lmer(log_Precip ~ AVA + (1 | Vineyard), data, REML=FALSE) 

log_Precip_AVA_Year_test <- anova(log_Precip_reduced_interaction, log_Precip_full_model)
log_Precip_AVA_Year_test 

#Data: data
#Models:
#  log_Precip_reduced_interaction: log_Precip ~ AVA + Year + (1 | Vineyard)
#log_Precip_full_model: log_Precip ~ AVA * Year + (1 | Vineyard)
#npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
#log_Precip_reduced_interaction    8 4.3013 14.959  5.8493  -11.699                       
#log_Precip_full_model            12 0.4096 16.396 11.7952  -23.590 11.892  4    0.01817 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

log_Precip_AVA_test <- anova(log_Precip_AVA, log_Precip_reduced_interaction)
log_Precip_AVA_test

#Data: data
#Models:
#  log_Precip_AVA: log_Precip ~ Year + (1 | Vineyard)
#log_Precip_reduced_interaction: log_Precip ~ AVA + Year + (1 | Vineyard)
#npar    AIC    BIC   logLik deviance  Chisq Df Pr(>Chisq)    
#log_Precip_AVA                    4 35.766 41.095 -13.8829   27.766                         
#log_Precip_reduced_interaction    8  4.301 14.959   5.8493  -11.699 39.465  4  5.585e-08 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

log_Precip_Year_test <- anova(log_Precip_Year, log_Precip_reduced_interaction)
log_Precip_Year_test

#Data: data
#Models:
#  log_Precip_Year: log_Precip ~ AVA + (1 | Vineyard)
#log_Precip_reduced_interaction: log_Precip ~ AVA + Year + (1 | Vineyard)
#npar    AIC    BIC   logLik deviance  Chisq Df Pr(>Chisq)    
#log_Precip_Year                   7 51.292 60.617 -18.6460   37.292                         
#log_Precip_reduced_interaction    8  4.301 14.959   5.8493  -11.699 48.991  1  2.572e-12 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

log_Precip_AVA_posthoc <- lsmeans(log_Precip_full_model, pairwise ~ AVA | Year) 
log_Precip_AVA_posthoc$contrasts

#Year = 2016:
#  contrast                          estimate    SE   df t.ratio p.value
#Mendocino - Monterey                 1.604 0.279 32.2  5.742  <.0001 
#Mendocino - Santa Barbara            1.129 0.259 33.3  4.364  0.0010 
#Mendocino - Sonoma                   0.232 0.236 34.9  0.982  0.8615 
#Mendocino - Willamette Valley       -1.047 0.299 37.1 -3.500  0.0102 
#Monterey - Santa Barbara            -0.475 0.235 26.1 -2.017  0.2863 
#Monterey - Sonoma                   -1.372 0.211 26.1 -6.512  <.0001 
#Monterey - Willamette Valley        -2.650 0.279 32.2 -9.490  <.0001 
#Santa Barbara - Sonoma              -0.897 0.182 26.1 -4.916  0.0004 
#Santa Barbara - Willamette Valley   -2.176 0.259 33.3 -8.411  <.0001 
#Sonoma - Willamette Valley          -1.279 0.236 34.9 -5.413  <.0001 

#Year = 2017:
#  contrast                          estimate    SE   df t.ratio p.value
#Mendocino - Monterey                 1.781 0.258 26.2  6.906  <.0001 
#Mendocino - Santa Barbara            1.184 0.235 26.2  5.029  0.0003 
#Mendocino - Sonoma                   0.153 0.211 26.2  0.726  0.9485 
#Mendocino - Willamette Valley       -0.641 0.258 26.3 -2.484  0.1251 
#Monterey - Santa Barbara            -0.597 0.235 26.1 -2.536  0.1131 
#Monterey - Sonoma                   -1.629 0.211 26.1 -7.732  <.0001 
#Monterey - Willamette Valley        -2.422 0.258 26.2 -9.390  <.0001 
#Santa Barbara - Sonoma              -1.031 0.182 26.1 -5.654  0.0001 
#Santa Barbara - Willamette Valley   -1.825 0.235 26.2 -7.750  <.0001 
#Sonoma - Willamette Valley          -0.794 0.211 26.2 -3.768  0.0069 

#Degrees-of-freedom method: kenward-roger 
#P value adjustment: tukey method for comparing a family of 5 estimates 

log_Precip_Year_posthoc <- lsmeans(log_Precip_full_model, pairwise ~ Year | AVA) 
log_Precip_Year_posthoc$contrasts

#AVA = Mendocino:
 
# contrast    estimate     SE   df t.ratio p.value
#2016 - 2017   -0.881 0.1522 22.5  -5.791 <.0001 

#AVA = Monterey:
#  contrast    estimate     SE   df t.ratio p.value
#2016 - 2017   -0.703 0.1081 20.8  -6.505 <.0001 

#AVA = Santa Barbara:
#  contrast    estimate     SE   df t.ratio p.value
#2016 - 2017   -0.826 0.0883 20.8  -9.352 <.0001 

#AVA = Sonoma:
#  contrast    estimate     SE   df t.ratio p.value
#2016 - 2017   -0.960 0.0624 20.8 -15.380 <.0001 

#AVA = Willamette Valley:
#  contrast    estimate     SE   df t.ratio p.value
#2016 - 2017   -0.475 0.1522 22.5  -3.121 0.0049 

#Degrees-of-freedom method: kenward-roger 


########################################  Bact_Shannon Logged for normality

log_Bact_Shannon_full_model <- lmer(log_Bact_Shannon ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
log_Bact_Shannon_reduced_interaction <- lmer(log_Bact_Shannon ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)

log_Bact_Shannon_AVA <- lmer(log_Bact_Shannon ~ Year + (1 | Vineyard), data, REML=FALSE) 
log_Bact_Shannon_Year <- lmer(log_Bact_Shannon ~ AVA + (1 | Vineyard), data, REML=FALSE) 

log_Bact_Shannon_AVA_Year_test <- anova(log_Bact_Shannon_reduced_interaction, log_Bact_Shannon_full_model)
log_Bact_Shannon_AVA_Year_test 

#Data: data
#Models:
#  log_Bact_Shannon_reduced_interaction: log_Bact_Shannon ~ AVA + Year + (1 | Vineyard)
#log_Bact_Shannon_full_model: log_Bact_Shannon ~ AVA * Year + (1 | Vineyard)
#npar    AIC     BIC  logLik deviance  Chisq Df Pr(>Chisq)   
#log_Bact_Shannon_reduced_interaction    8 91.894 102.552 -37.947   75.894                        
#log_Bact_Shannon_full_model            12 83.856  99.842 -29.928   59.856 16.039  4   0.002968 **
  ---
  #  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

log_Bact_Shannon_AVA_test <- anova(log_Bact_Shannon_AVA, log_Bact_Shannon_reduced_interaction)
log_Bact_Shannon_AVA_test

#Data: data
#Models:
#  log_Bact_Shannon_AVA: log_Bact_Shannon ~ Year + (1 | Vineyard)
#log_Bact_Shannon_reduced_interaction: log_Bact_Shannon ~ AVA + Year + (1 | Vineyard)
#npar     AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#log_Bact_Shannon_AVA                    4 102.965 108.29 -47.483   94.965                         
#log_Bact_Shannon_reduced_interaction    8  91.894 102.55 -37.947   75.894 19.071  4  0.0007612 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

log_Bact_Shannon_Year_test <- anova(log_Bact_Shannon_Year, log_Bact_Shannon_reduced_interaction)
log_Bact_Shannon_Year_test

#Data: data
#Models:
#  log_Bact_Shannon_Year: log_Bact_Shannon ~ AVA + (1 | Vineyard)
#log_Bact_Shannon_reduced_interaction: log_Bact_Shannon ~ AVA + Year + (1 | Vineyard)
#npar    AIC     BIC  logLik deviance  Chisq Df Pr(>Chisq)
#log_Bact_Shannon_Year                   7 90.176  99.501 -38.088   76.176                     
#log_Bact_Shannon_reduced_interaction    8 91.894 102.552 -37.947   75.894 0.2813  1     0.5959

log_Bact_Shannon_AVA_posthoc <- lsmeans(log_Bact_Shannon_full_model, pairwise ~ AVA | Year) 
log_Bact_Shannon_AVA_posthoc$contrasts

#Year = 2016:
#  contrast                          estimate    SE   df t.ratio p.value
#Mendocino - Monterey              -5.23282 1.120 43.5 -4.671  0.0003 
#Mendocino - Santa Barbara         -2.52077 1.061 43.5 -2.377  0.1412 
#Mendocino - Sonoma                -1.80815 0.998 43.5 -1.813  0.3797 
#Mendocino - Willamette Valley     -2.52561 1.315 43.6 -1.920  0.3224 
#Monterey - Santa Barbara           2.71205 0.806 42.5  3.365  0.0134 
#Monterey - Sonoma                  3.42467 0.721 42.5  4.750  0.0002 
#Monterey - Willamette Valley       2.70721 1.120 43.5  2.417  0.1302 
#Santa Barbara - Sonoma             0.71262 0.624 42.5  1.141  0.7838 
#Santa Barbara - Willamette Valley -0.00483 1.061 43.5 -0.005  1.0000 
#Sonoma - Willamette Valley        -0.71745 0.998 43.5 -0.719  0.9509 

#Year = 2017:
#  contrast                          estimate    SE   df t.ratio p.value
#Mendocino - Monterey              -2.25989 0.883 42.6 -2.559  0.0966 
#Mendocino - Santa Barbara          0.36329 0.806 42.6  0.451  0.9912 
#Mendocino - Sonoma                -1.52107 0.721 42.6 -2.110  0.2349 
#Mendocino - Willamette Valley     -1.52492 0.883 42.6 -1.727  0.4289 
#Monterey - Santa Barbara           2.62319 0.806 42.5  3.254  0.0180 
#Monterey - Sonoma                  0.73882 0.721 42.5  1.025  0.8425 
#Monterey - Willamette Valley       0.73497 0.883 42.6  0.832  0.9190 
#Santa Barbara - Sonoma            -1.88436 0.624 42.5 -3.018  0.0330 
#Santa Barbara - Willamette Valley -1.88822 0.806 42.6 -2.342  0.1516 
#Sonoma - Willamette Valley        -0.00385 0.721 42.6 -0.005  1.0000 

#Degrees-of-freedom method: kenward-roger 
#P value adjustment: tukey method for comparing a family of 5 estimates 

log_Bact_Shannon_Year_posthoc <- lsmeans(log_Bact_Shannon_full_model, pairwise ~ Year | AVA) 
log_Bact_Shannon_Year_posthoc$contrasts

#AVA = Mendocino:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017   -1.493 1.070 28.3 -1.396  0.1737 

#AVA = Monterey:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017    1.480 0.818 20.5  1.810  0.0850 

#AVA = Santa Barbara:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017    1.391 0.668 20.5  2.083  0.0500 

#AVA = Sonoma:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017   -1.206 0.472 20.5 -2.553  0.0187 

#AVA = Willamette Valley:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017   -0.492 1.070 28.3 -0.460  0.6490 

#Degrees-of-freedom method: kenward-roger 


########################################  Bact_InvSimp Logged for normality



log_Bact_InvSimp_full_model <- lmer(log_Bact_InvSimp ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
log_Bact_InvSimp_reduced_interaction <- lmer(log_Bact_InvSimp ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)

log_Bact_InvSimp_AVA <- lmer(log_Bact_InvSimp ~ Year + (1 | Vineyard), data, REML=FALSE) 
log_Bact_InvSimp_Year <- lmer(log_Bact_InvSimp ~ AVA + (1 | Vineyard), data, REML=FALSE) 

log_Bact_InvSimp_AVA_Year_test <- anova(log_Bact_InvSimp_reduced_interaction, log_Bact_InvSimp_full_model)
log_Bact_InvSimp_AVA_Year_test 

#Data: data
#Models:
#  log_Bact_InvSimp_reduced_interaction: log_Bact_InvSimp ~ AVA + Year + (1 | Vineyard)
#log_Bact_InvSimp_full_model: log_Bact_InvSimp ~ AVA * Year + (1 | Vineyard)
#npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#log_Bact_InvSimp_reduced_interaction    8 81.885 92.543 -32.943   65.885                         
#log_Bact_InvSimp_full_model            12 64.715 80.701 -20.357   40.715 25.171  4  4.648e-05 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

log_Bact_InvSimp_AVA_test <- anova(log_Bact_InvSimp_AVA, log_Bact_InvSimp_reduced_interaction)
log_Bact_InvSimp_AVA_test

#Data: data
#Models:
#  log_Bact_InvSimp_AVA: log_Bact_InvSimp ~ Year + (1 | Vineyard)
#log_Bact_InvSimp_reduced_interaction: log_Bact_InvSimp ~ AVA + Year + (1 | Vineyard)
#npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
#log_Bact_InvSimp_AVA                    4 90.531 95.860 -41.265   82.531                        
#log_Bact_InvSimp_reduced_interaction    8 81.885 92.543 -32.943   65.885 16.645  4   0.002265 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

log_Bact_InvSimp_Year_test <- anova(log_Bact_InvSimp_Year, log_Bact_InvSimp_reduced_interaction)
log_Bact_InvSimp_Year_test

#Data: data
#Models:
#  log_Bact_InvSimp_Year: log_Bact_InvSimp ~ AVA + (1 | Vineyard)
#log_Bact_InvSimp_reduced_interaction: log_Bact_InvSimp ~ AVA + Year + (1 | Vineyard)
#npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#log_Bact_InvSimp_Year                   7 80.272 89.597 -33.136   66.272                     
#log_Bact_InvSimp_reduced_interaction    8 81.885 92.543 -32.943   65.885 0.3864  1     0.5342


log_Bact_InvSimp_AVA_posthoc <- lsmeans(log_Bact_InvSimp_full_model, pairwise ~ AVA | Year) 
log_Bact_InvSimp_AVA_posthoc$contrasts

#Year = 2016:
#  contrast                           estimate    SE   df t.ratio p.value
#Mendocino - Monterey              -4.454901 0.797 42.9 -5.592  <.0001 
#Mendocino - Santa Barbara         -1.624658 0.753 43.1 -2.159  0.2149 
#Mendocino - Sonoma                -1.128379 0.706 43.4 -1.599  0.5063 
#Mendocino - Willamette Valley     -1.416118 0.927 43.6 -1.527  0.5510 
#Monterey - Santa Barbara           2.830243 0.584 39.0  4.845  0.0002 
#Monterey - Sonoma                  3.326522 0.522 39.0  6.367  <.0001 
#Monterey - Willamette Valley       3.038783 0.797 42.9  3.814  0.0038 
#Santa Barbara - Sonoma             0.496278 0.452 39.0  1.097  0.8071 
#Santa Barbara - Willamette Valley  0.208540 0.753 43.1  0.277  0.9987 
#Sonoma - Willamette Valley        -0.287738 0.706 43.4 -0.408  0.9940 
#
#Year = 2017:
  #  contrast                           estimate    SE   df t.ratio p.value
#Mendocino - Monterey              -1.117651 0.640 39.2 -1.747  0.4183 
#Mendocino - Santa Barbara          0.137553 0.584 39.2  0.235  0.9993 
#Mendocino - Sonoma                -1.105664 0.522 39.3 -2.116  0.2337 
#Mendocino - Willamette Valley     -1.118259 0.640 39.4 -1.748  0.4177 
#Monterey - Santa Barbara           1.255204 0.584 39.0  2.149  0.2207 
#Monterey - Sonoma                  0.011987 0.522 39.0  0.023  1.0000 
#Monterey - Willamette Valley      -0.000607 0.640 39.2 -0.001  1.0000 
#Santa Barbara - Sonoma            -1.243217 0.452 39.0 -2.748  0.0648 
#Santa Barbara - Willamette Valley -1.255812 0.584 39.2 -2.150  0.2201 
#Sonoma - Willamette Valley        -0.012595 0.522 39.3 -0.024  1.0000 

#Degrees-of-freedom method: kenward-roger 
#P value adjustment: tukey method for comparing a family of 5 estimates 

log_Bact_InvSimp_Year_posthoc <- lsmeans(log_Bact_InvSimp_full_model, pairwise ~ Year | AVA) 
log_Bact_InvSimp_Year_posthoc$contrasts

#AVA = Mendocino:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017   -1.124 0.711 27.1 -1.581  0.1254 

#AVA = Monterey:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017    2.213 0.529 20.5  4.182  0.0004 

#AVA = Santa Barbara:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017    0.638 0.432 20.5  1.477  0.1550 

#AVA = Sonoma:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017   -1.101 0.306 20.5 -3.605  0.0017 

#AVA = Willamette Valley:
#  contrast    estimate    SE   df t.ratio p.value
#2016 - 2017   -0.826 0.711 27.1 -1.162  0.2553 

#Degrees-of-freedom method: kenward-roger 



########################################################################Fung_Rich log transformed



log_Fung_Rich_full_model <- lmer(log_Fung_Rich ~ AVA*Year + (1 | Vineyard), data, REML=FALSE)
log_Fung_Rich_reduced_interaction <- lmer(log_Fung_Rich ~ AVA + Year + (1 | Vineyard), data, REML=FALSE)

log_Fung_Rich_AVA <- lmer(log_Fung_Rich ~ Year + (1 | Vineyard), data, REML=FALSE) 
log_Fung_Rich_Year <- lmer(log_Fung_Rich ~ AVA + (1 | Vineyard), data, REML=FALSE) 

log_Fung_Rich_AVA_Year_test <- anova(log_Fung_Rich_reduced_interaction, log_Fung_Rich_full_model)
log_Fung_Rich_AVA_Year_test 

#Data: data
#Models:
#  log_Fung_Rich_reduced_interaction: log_Fung_Rich ~ AVA + Year + (1 | Vineyard)
#log_Fung_Rich_full_model: log_Fung_Rich ~ AVA * Year + (1 | Vineyard)
#npar     AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
#log_Fung_Rich_reduced_interaction    8 -3.1247 7.5330  9.5623  -19.125                        
#log_Fung_Rich_full_model            12 -9.8815 6.1049 16.9408  -33.882 14.757  4   0.005233 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

log_Fung_Rich_AVA_test <- anova(log_Fung_Rich_AVA, log_Fung_Rich_reduced_interaction)
log_Fung_Rich_AVA_test

#Data: data
#Models:
#  log_Fung_Rich_AVA: log_Fung_Rich ~ Year + (1 | Vineyard)
#log_Fung_Rich_reduced_interaction: log_Fung_Rich ~ AVA + Year + (1 | Vineyard)
#npar     AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#log_Fung_Rich_AVA                    4 10.3666 15.695 -1.1833   2.3666                         
#log_Fung_Rich_reduced_interaction    8 -3.1247  7.533  9.5623 -19.1247 21.491  4   0.000253 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 

log_Fung_Rich_Year_test <- anova(log_Fung_Rich_Year, log_Fung_Rich_reduced_interaction)
log_Fung_Rich_Year_test

#Data: data
#Models:
#  log_Fung_Rich_Year: log_Fung_Rich ~ AVA + (1 | Vineyard)
#log_Fung_Rich_reduced_interaction: log_Fung_Rich ~ AVA + Year + (1 | Vineyard)
#npar     AIC    BIC logLik deviance Chisq Df Pr(>Chisq)
#log_Fung_Rich_Year                   7 -4.1817 5.1437 9.0909  -18.182                    
#log_Fung_Rich_reduced_interaction    8 -3.1247 7.5330 9.5623  -19.125 0.943  1     0.3315


log_Fung_Rich_AVA_posthoc <- lsmeans(log_Fung_Rich_full_model, pairwise ~ AVA | Year) 
log_Fung_Rich_AVA_posthoc$contrasts

#Year = 2016:
#  contrast                          estimate    SE   df t.ratio p.value
#Mendocino - Monterey                0.0464 0.213 40.3  0.218  0.9995 
#Mendocino - Santa Barbara          -0.0991 0.200 41.2 -0.496  0.9873 
#Mendocino - Sonoma                 -0.4766 0.186 42.2 -2.565  0.0956 
#Mendocino - Willamette Valley      -1.0621 0.241 43.1 -4.402  0.0006 
#Monterey - Santa Barbara           -0.1455 0.165 32.4 -0.884  0.9009 
#Monterey - Sonoma                  -0.5229 0.147 32.4 -3.551  0.0098 
#Monterey - Willamette Valley       -1.1085 0.213 40.3 -5.204  0.0001 
#Santa Barbara - Sonoma             -0.3774 0.128 32.4 -2.960  0.0423 
#Santa Barbara - Willamette Valley  -0.9630 0.200 41.2 -4.818  0.0002 
#Sonoma - Willamette Valley         -0.5855 0.186 42.2 -3.151  0.0236 

#Year = 2017:
#  contrast                          estimate    SE   df t.ratio p.value
#Mendocino - Monterey               -0.1444 0.180 32.6 -0.801  0.9285 
#Mendocino - Santa Barbara           0.0745 0.165 32.7  0.452  0.9909 
#Mendocino - Sonoma                 -0.2080 0.147 32.8 -1.413  0.6241 
#Mendocino - Willamette Valley      -0.6356 0.180 32.9 -3.524  0.0104 
#Monterey - Santa Barbara            0.2189 0.165 32.4  1.330  0.6752 
#Monterey - Sonoma                  -0.0636 0.147 32.4 -0.432  0.9924 
#Monterey - Willamette Valley       -0.4912 0.180 32.6 -2.724  0.0719 
#Santa Barbara - Sonoma             -0.2825 0.128 32.4 -2.215  0.1997 
#Santa Barbara - Willamette Valley  -0.7101 0.165 32.7 -4.313  0.0012 
#Sonoma - Willamette Valley         -0.4276 0.147 32.8 -2.904  0.0480 


log_Fung_Rich_Year_posthoc <- lsmeans(log_Fung_Rich_full_model, pairwise ~ Year | AVA) 
log_Fung_Rich_Year_posthoc$contrasts

#AVA = Mendocino:
#  contrast    estimate     SE   df t.ratio p.value
#2016 - 2017  -0.1075 0.1641 24.8 -0.655  0.5184 

#AVA = Monterey:
#  contrast    estimate     SE   df t.ratio p.value
#2016 - 2017  -0.2983 0.1187 20.5 -2.513  0.0204 

#AVA = Santa Barbara:
#  contrast    estimate     SE   df t.ratio p.value
#2016 - 2017   0.0661 0.0969 20.5  0.682  0.5027 

#AVA = Sonoma:
#  contrast    estimate     SE   df t.ratio p.value
#2016 - 2017   0.1610 0.0685 20.5  2.350  0.0289 

#AVA = Willamette Valley:
#  contrast    estimate     SE   df t.ratio p.value
#2016 - 2017   0.3190 0.1641 24.8  1.944  0.0633 

#Degrees-of-freedom method: kenward-roger 


####################################################Comparing different transformations of InvSimpson for Fungi vs. untransformed




---
#  title: "Mixed model results tables"
#author: "Ian Morelan"
#date: "11/5/2020"
output: word_document
---
  
  ```{r setup, include=FALSE}
#load(file = "JFW_mixed_model.RData")
library(dplyr)
library(broom)
library(knitr)
library(gtools)
```

### Must characteristics and Climate


Must_acidity_AVA_Year_test %>% kable(caption = "**Must TA AVA/Vintage interaction test**")
Must_acidity_AVA_test %>% kable(caption = "**Must TA AVA test**")
Must_acidity_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**Must TA AVA contrasts**")
Must_acidity_Year_test %>% kable(caption = "**Must TA Vintage test**")

Must_pH_AVA_Year_test %>% kable(caption = "**Must pH AVA/Vintage interaction test**")
Must_pH_AVA_test %>% kable(caption = "**Must pH AVA tes**t")
Must_pH_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**Must pH AVA contrasts**")
Must_pH_Year_test %>% kable(caption = "**Must pH Vintage test**")

Must_TSS_AVA_Year_test %>% kable(caption = "**Must TSS AVA/Vintage interaction test**")
Must_TSS_AVA_test %>% kable(caption = "**Must TSS AVA test**")
Must_TSS_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**Must TSS AVA contrasts**")
Must_TSS_Year_test %>% kable(caption = "**Must TSS Vintage test**")

log_Precip_AVA_Year_test %>% kable(caption = "**Log Precipitation AVA/Vintage interaction test**")
log_Precip_AVA_test %>% kable(caption = "**Log Precipitation AVA test**")
log_Precip_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**Log Precipitation AVA contrasts**")
log_Precip_Year_test %>% kable(caption = "**Log Precipitation Vintage test**")
log_Precip_Year_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(p.value)) %>%
  kable(caption = "**Log Precipitation Vintage contrasts**")

GDD_AVA_Year_test %>% kable(caption = "**GDD AVA/Vintage interaction test**")
GDD_AVA_test %>% kable(caption = "**GDD AVA test**")
GDD_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**GDD AVA contrasts**")
GDD_Year_test %>% kable(caption = "**GDD Vintage test**")
GDD_Year_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(p.value)) %>%
  kable(caption = "**GDD Vintage contrasts**")



### Bacterial diversity


bact_rich_AVA_Year_test %>% kable(caption = "**Bacterial richness AVA/Vintage interaction test**")
bact_rich_AVA_test %>% kable(caption = "**Bacterial richness AVA test**")
bact_rich_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**Bacterial richness AVA contrasts**")
bact_rich_Year_test %>% kable(caption = "**Bacterial richness Vintage test**")


log_Bact_Shannon_AVA_Year_test %>% kable(caption = "**log Bacterial Exponential Shannon AVA/Vintage interaction test**")
log_Bact_Shannon_AVA_test %>% kable(caption = "**log Bacterial Exponential Shannon AVA test**")
log_Bact_Shannon_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**log Bacterial Exponential Shannon AVA contrasts**")
log_Bact_Shannon_Year_test %>% kable(caption = "**log Bacterial Exponential Shannon Vintage test**")


log_Bact_InvSimp_AVA_Year_test %>% kable(caption = "**Log Bacterial Inverse Simpson AVA/Vintage interaction test**")
log_Bact_InvSimp_AVA_test %>% kable(caption = "**Log Inverse Simpson AVA test**")
log_Bact_InvSimp_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**Log Bacterial Inverse Simpson AVA contrasts**")
log_Bact_InvSimp_Year_test %>% kable(caption = "**Log Bacterial Inverse Simpson Vintage test**")
log_Bact_InvSimp_Year_posthoc$contrasts %>% 
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(p.value)) %>%
  kable(caption = "**Log Bacterial Inverse Simpson Vintage contrasts**")


### Fungal diversity


log_Fung_Rich_AVA_Year_test %>% kable(caption = "**log Fungal richness AVA/Vintage interaction test**")
log_Fung_Rich_AVA_test %>% kable(caption = "**log Fungal richness AVA test**")
fung_rich_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**log Fungal richness AVA contrasts**")
log_Fung_Rich_Year_test %>% kable(caption = "**log Fungal richness Vintage test**")


fung_Shannon_AVA_Year_test %>% kable(caption = "**Fungal Exponential Shannon AVA/Vintage interaction test**")
fung_Shannon_AVA_test %>% kable(caption = "**Fungal Exponential Shannon AVA test**")
fung_Shannon_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**Fungal Exponential Shannon AVA contrasts**")
fung_Shannon_Year_test %>% kable(caption = "**Fungal Exponential Shannon Vintage test**")
fung_Shannon_Year_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(p.value)) %>%
  kable(caption = "**Fungal Exponential Shannon Vintage contrasts**")

fung_InvSimp_AVA_Year_test %>% kable(caption = "**Fungal Inverse Simpson AVA/Vintage interaction test**")
fung_InvSimp_AVA_test %>% kable(caption = "**Fungal Inverse Simpson AVA test**")
fung_InvSimp_AVA_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(adj.p.value)) %>%
  kable(caption = "**Fungal Inverse Simpson AVA contrasts**")
fung_InvSimp_Year_test %>% kable(caption = "**Fungal Inverse Simpson Vintage test**")
fung_InvSimp_Year_posthoc$contrasts %>%
  tidy() %>%
  select(!null.value) %>%
  mutate(signif = stars.pval(p.value)) %>%
  kable(caption = "**Fungal Inverse Simpson Vintage contrasts**")


