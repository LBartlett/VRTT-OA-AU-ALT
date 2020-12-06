# Oxalic Acid FFAR study 
# Shop towel data (Auburn): Duck Pond study.
# Analysis by Dr. Lewis J Bartlett & Christian Baker
# lewis.bartlett@uga.edu
# Data: Williams Lab, Auburn University

# GitDir: https://github.com/LBartlett/VRTT-OA-AU-ALT.git

#load packages
library('lme4')
library('afex')
library('car')

# set random for reproducability
set.seed(137880)

# Read in main data
STFull <- read.csv('D:/Google Drive/PostDoc2/Varroa Control/OxalicAcid/DuckPond/ShopTowelOA_AU_Clean.csv', header=T)

# check
head(STFull)

##############################################
#### Let's start with the Varroa mite analysis

# Data cleaning. Strip out everything where we don't have mite numbers.
OAData.M <- STFull[which(!is.na(STFull$Mites)),]

# Data cleaning. Include only colonies which are recorded in all 3 time points.
OAData.M <- OAData.M[which(!is.na(match(as.character(OAData.M$Colony),as.character(OAData.M$Colony[which(OAData.M$Timepoint == 1)])))),]
OAData.M <- OAData.M[which(!is.na(match(as.character(OAData.M$Colony),as.character(OAData.M$Colony[which(OAData.M$Timepoint == 2)])))),]
OAData.M <- OAData.M[which(!is.na(match(as.character(OAData.M$Colony),as.character(OAData.M$Colony[which(OAData.M$Timepoint == 3)])))),]

# Check data frame looks alright
head(OAData.M)

# Inspect distribution of response variable
hist(OAData.M$Mites)

#Poisson modelling requires whole numbers. Round up to nearest whole mite where X = 1 for 0 < X =< 1.
# Distribution otherwise well matches what we expect for this count data & a Poisson-distributed linear modelling approach.
OAData.M$Mites.R <- ceiling(OAData.M$Mites)


## Look for significant effect of treatment (control, shop towel, cardboard)
## We're looking for if treatment changes rate of change of mite loads with time (interaction)
## Error structure reflects that each colony is sampled multiple times (Timepoint).
## Note no 'treatment alone' effect as at time Timepoint 1, no had been yet applied, so we don't expect differences in intercepts beyond controlling for each colony in the random component

mixed(Mites.R ~ Timepoint + Timepoint:Treatment + (1|Timepoint/Colony),
      family = 'poisson', method = 'LRT',
      data = OAData.M)

# No significant effect of treatments compared control, significant effect of time.

## Look at direction of effect sizes to check against significance.
coef(summary(glmer(Mites.R ~ Timepoint + Timepoint:Treatment + (1|Timepoint/Colony),
                   family = 'poisson',
                   data = OAData.M)))[,1:2]

# Mites increase with time (positive effect of Timepoint); positive intercept (sanity check)
# Effect size for both treatments negative compared to control but (from above), not significant.
# We see mite loads increase with time throughout the experiment. 
# Treated colonies showed mite load increases which were not significantly different to control colonies.

# We can make a graph of this.

# Quick transparency function for plotting
Transpa <- function(color, percent) {
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  
  return(t.col)
  
}

# Apply x-axis jitter for better plotting.
OAData.M$TimepointJ <- sapply((OAData.M$Timepoint+5), jitter, factor = 0.5) - 5

# Colour scheme for whole analysis series
ColRef <- data.frame(Treatment = levels(as.factor(OAData.M$Treatment)), Col =  c('red3','purple2','orange3'))

# Make legends easier
ColRef$LongName <- c('Control','Cardboard','Shoptowel')

# Make blank plot
plot(OAData.M$Mites.R ~ OAData.M$TimepointJ, 
     type = 'n', 
     xaxp = c(1,3,2), 
     ylab = 'Percent Mite Intensity', 
     xlab ='Timepoint',
     cex.lab = 1.3)

# Add points showing each colony's mite intensity across time, colour coded by treatment

for(C in unique(OAData.M$Colony)){
  
  CCol <- as.character(ColRef$Col[which(ColRef$Treatment == unique(OAData.M$Treatment[which(OAData.M$Colony == C)]))])
  
  points(x = OAData.M$TimepointJ[which(OAData.M$Colony == C)],
         y = OAData.M$Mites[which(OAData.M$Colony == C)],
         col = Transpa(CCol, 65), type = 'b', lwd = 0.5, pch = 20, lty = 3, cex = 2)
  
}

# Add naive regression lines to roughly represent the mixed-modelling done

for(Trt in 1:NROW(ColRef)){
  
  LinMod <- glm(OAData.M$Mites.R[which(OAData.M$Treatment == ColRef$Treatment[Trt])] ~ OAData.M$Timepoint[which(OAData.M$Treatment == ColRef$Treatment[Trt])],
                family = 'poisson')
  
  V1 <- seq(from = par('usr')[1]*1.15, to = par('usr')[2]*0.955, length.out = 100)
  V2 <- exp(coef(LinMod)[[1]] + (V1*coef(LinMod)[[2]]))
  
  points(V2 ~ V1, type = 'l', lty = 1, lwd = 6, col = ColRef$Col[Trt])

}

legend(x = 1, y = 20,
       legend = ColRef$LongName,
       col = as.character(ColRef$Col),
       pch = 19,
       bty = 'n',
       cex = 1.2)

# The above framework only works for the mite data as it's the only data with three Timepoints. It makes for a good supplementary analysis.
# Ideally, for this sort opf framework, we would have wanted much more dense Timepoint sampling but that is difficult to impossible in this system (destructive sampling).

# We can also simply look at overall change in mites from start to end, essentially just looking at TP1 vs TP3 according to treatment.
# This is in the same fashion as we will do for the colony health and survival data. 
# Calculate and analyse change in mite loads from P1 to P3.

# Assemble necessary data frame
MD <- as.data.frame(unique(OAData.M$Colony[which(OAData.M$Timepoint == 3)]))

colnames(MD) <- 'Colony'
MD$Treatment <- NA
MD$P1M <- NA
MD$P3M <- NA
MD$DeltaMites <- NA

# Populate and calculate DeltaMites - difference in mite loads between Timepoint 1 and Timepoint 3 for each colony
for(C in 1:NROW(MD)){
  
  Col <- MD$Colony[C]
  
  MD$Treatment[C] <- as.character(unique(OAData.M$Treatment[which(OAData.M$Colony==Col)]))
  
  MD$P1M[C] <- OAData.M$Mites[which(OAData.M$Colony==Col & OAData.M$Timepoint == 1)]
  
  MD$P3M[C] <- OAData.M$Mites[which(OAData.M$Colony==Col & OAData.M$Timepoint == 3)]
  
  MD$DeltaMites[C] <- (MD$P3M[C] - MD$P1M[C])
  
}

# Look at distribution. Arguably normal enough bar one outlier, pay some attention to model graphical inspections.
hist(MD$DeltaMites)

# Plot.
par(mar=c(5,5,2,2))

boxplot(MD$DeltaMites ~ MD$Treatment, 
        main = NA, ylab = expression(paste(Delta, ' Percent Mite Intensity',sep='')), xlab = 'Treatment', 
        border = 'transparent', cex.axis = 1.3, cex.lab = 1.5, outline = FALSE, lwd = 1.2,
        boxlty = 0, whisklty = 0, staplelty = 1, boxwex = 0.3, col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 70))
stripchart(MD$DeltaMites ~ MD$Treatment,
           col = as.character(ColRef$Col),
           vertical = T, add = T, pch = c(16,4,8), cex = 1.5, 
           method = 'jitter', lwd = 2)
PlotAg <- aggregate(MD$DeltaMites ~ MD$Treatment, FUN = mean)
stripchart(PlotAg[,2] ~ PlotAg[,1],
           col = as.character(ColRef$Col),
           vertical = T, add = T, pch = 95, cex = 6, lwd = 2)
abline(a = 0, b = 0, lty = 3)

# Test for significance of difference between groups.

# Create model
AVMod.M <- aov(MD$DeltaMites~MD$Treatment)

# Inspect model fit & residuals
qqnorm(resid(AVMod.M))
qqline(resid(AVMod.M), col = "blue1", lwd = 2)
par(mfrow = c(2,2))
plot(AVMod.M)
par(mfrow = c(1,1))
hist(resid(AVMod.M))
shapiro.test(resid(AVMod.M))

# Proceed cautiously with analysis - fit could be better but isn't so bad as to be prohibitive. Field data for you!
Anova(AVMod.M, type = 2)

# p = 0.071
# We can see even in treated colonies, delta mites is mostly positive (mite loads increased despite 'treatment').
# Any hint of an effect is at best an inhibition of mite population growth - little evidence for reduction of mites.

# We can crudely marry this to a modified Henderson-Tilton's efficacy measure.
# Take means across colonies to get a 'plot' measure.
# You could do this MUCH fancier with some complicated bootstrapping but I don't think it's worth it.

CB.Effic <- 1 - (
  (mean(MD$P1M[which(MD$Treatment == 'C')]) * mean(MD$P3M[which(MD$Treatment == 'CB')]))
  /
    (mean(MD$P3M[which(MD$Treatment == 'C')]) * mean(MD$P1M[which(MD$Treatment == 'CB')]))  
)

ST.Effic <- 1 - (
  (mean(MD$P1M[which(MD$Treatment == 'C')]) * mean(MD$P3M[which(MD$Treatment == 'ST')]))
  /
    (mean(MD$P3M[which(MD$Treatment == 'C')]) * mean(MD$P1M[which(MD$Treatment == 'ST')]))  
)

# CB: 61%, ST: 49%

# For completeness, let's do a t-test confirming mite numbers did increase with time like the mixed-modelling says.
# For this, we want to see if the distribution of the delta-mites values is significantly different to averaging zero
# It's an easy two-sided t-test with 1 sample set
t.test(x = MD$DeltaMites, mu = 0)
# We can see mites significantly increased over time, in agreement with our generalized linear mixed modelling.

# We could do this for each specific treatment, if we wanted, so hammer home that treatment didn't reduce mite load
t.test(x = MD$DeltaMites[which(MD$Treatment == 'C')], mu = 0)
t.test(x = MD$DeltaMites[which(MD$Treatment == 'CB')], mu = 0)
t.test(x = MD$DeltaMites[which(MD$Treatment == 'ST')], mu = 0)

# No treatment reduced mite load.

########################################################################################
###### Mites analysis done. Let's look at other colony health metrics.
## Note we only have start / end measures for these, which large simplifies analysis.
## So, let's follow the same framework as the 'Delta-PMI' above for 'colony health' data

## First: adults
# Data cleaning. Strip out everything where we don't have numbers.
OAData.A <- STFull[which(!is.na(STFull$Adults)),]

# Data cleaning. Include only colonies which are recorded in P1 and P3.
OAData.A <- OAData.A[which(!is.na(match(as.character(OAData.A$Colony),as.character(OAData.A$Colony[which(OAData.A$Timepoint == 1)])))),]
OAData.A <- OAData.A[which(!is.na(match(as.character(OAData.A$Colony),as.character(OAData.A$Colony[which(OAData.A$Timepoint == 3)])))),]

# Check data frame looks alright
head(OAData.A)

# Assemble necessary data frame
AD <- as.data.frame(unique(OAData.A$Colony[which(OAData.A$Timepoint == 3)]))
colnames(AD) <- 'Colony'
AD$Treatment <- NA
AD$P1M <- NA
AD$P3M <- NA
AD$DeltaAdults <- NA

# Populate and calculate DeltaAdults - difference in mite loads between Timepoint 1 and Timepoint 3 for each colony
for(C in 1:NROW(AD)){
  
  Col <- AD$Colony[C]
  
  AD$Treatment[C] <- as.character(unique(OAData.A$Treatment[which(OAData.A$Colony==Col)]))
  
  AD$P1M[C] <- OAData.A$Adults[which(OAData.A$Colony==Col & OAData.A$Timepoint == 1)]
  
  AD$P3M[C] <- OAData.A$Adults[which(OAData.A$Colony==Col & OAData.A$Timepoint == 3)]
  
  AD$DeltaAdults[C] <- (AD$P3M[C] - AD$P1M[C])
  
}

# Look at distribution. Not terrible!
hist(AD$DeltaAdults)

# Plot.
par(mar=c(5,5,2,2))

boxplot(AD$DeltaAdults ~ AD$Treatment, 
        main = NA, ylab = expression(paste(Delta, ' Adult Population',sep='')), xlab = 'Treatment', 
        border = 'transparent', cex.axis = 1.3, cex.lab = 1.5, outline = FALSE, lwd = 1.2,
        boxlty = 0, whisklty = 0, staplelty = 1, boxwex = 0.3, col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 70))
stripchart(AD$DeltaAdults ~ AD$Treatment,
           col = as.character(ColRef$Col),
           vertical = T, add = T, pch = c(16,4,8), cex = 1.5, 
           method = 'jitter', lwd = 2)
PlotAg <- aggregate(AD$DeltaAdults ~ AD$Treatment, FUN = mean)
stripchart(PlotAg[,2] ~ PlotAg[,1],
           col = as.character(ColRef$Col),
           vertical = T, add = T, pch = 95, cex = 6, lwd = 2)
abline(a = 0, b = 0, lty = 3)

# Test for significance of difference between groups.

# Create model
AVMod.A <- aov(AD$DeltaAdults~AD$Treatment)

# Inspect model fit & residuals
qqnorm(resid(AVMod.A))
qqline(resid(AVMod.A), col = "blue1", lwd = 2)
par(mfrow = c(2,2))
plot(AVMod.A)
par(mfrow = c(1,1))
hist(resid(AVMod.A))
shapiro.test(resid(AVMod.A))

# Residuals distributed beautifully but qq-plot is a touch concerning at the lower tail...
Anova(AVMod.A, type = 2)
# Less concerning however when we are not claiming a positive result, and graphical interpretation is still most usful for practitioners.

###################
## Exactly the same again just with brood.
# Data cleaning. Strip out everything where we don't have numbers.
OBData.B <- STFull[which(!is.na(STFull$Brood)),]

# Data cleaning. Include only colonies which are recorded in P1 and P3.
OBData.B <- OBData.B[which(!is.na(match(as.character(OBData.B$Colony),as.character(OBData.B$Colony[which(OBData.B$Timepoint == 1)])))),]
OBData.B <- OBData.B[which(!is.na(match(as.character(OBData.B$Colony),as.character(OBData.B$Colony[which(OBData.B$Timepoint == 3)])))),]

# Check data frame looks alright
head(OBData.B)

# Assemble necessary data frame
BD <- as.data.frame(unique(OBData.B$Colony[which(OBData.B$Timepoint == 3)]))
colnames(BD) <- 'Colony'
BD$Treatment <- NA
BD$P1M <- NA
BD$P3M <- NA
BD$DeltaBrood <- NA

# Populate and calculate DeltaBrood - difference in mite loBDs between Timepoint 1 and Timepoint 3 for each colony
for(C in 1:NROW(BD)){
  
  Col <- BD$Colony[C]
  
  BD$Treatment[C] <- as.character(unique(OBData.B$Treatment[which(OBData.B$Colony==Col)]))
  
  BD$P1M[C] <- OBData.B$Brood[which(OBData.B$Colony==Col & OBData.B$Timepoint == 1)]
  
  BD$P3M[C] <- OBData.B$Brood[which(OBData.B$Colony==Col & OBData.B$Timepoint == 3)]
  
  BD$DeltaBrood[C] <- (BD$P3M[C] - BD$P1M[C])
  
}

# Look at distribution. Bit skewed, we'll see what the model looks like.
hist(BD$DeltaBrood)

# Plot.
par(mar=c(5,5,2,2))

boxplot(BD$DeltaBrood ~ BD$Treatment, 
        main = NA, ylab = expression(paste(Delta, ' Brood Population',sep='')), xlab = 'Treatment', 
        border = 'transparent', cex.Bxis = 1.3, cex.lab = 1.5, outline = FALSE, lwd = 1.2,
        boxlty = 0, whisklty = 0, staplelty = 1, boxwex = 0.3, col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 70))
stripchart(BD$DeltaBrood ~ BD$Treatment,
           col = as.character(ColRef$Col),
           vertical = T, add = T, pch = c(16,4,8), cex = 1.5, 
           method = 'jitter', lwd = 2)
PlotAg <- aggregate(BD$DeltaBrood ~ BD$Treatment, FUN = mean)
stripchart(PlotAg[,2] ~ PlotAg[,1],
           col = as.character(ColRef$Col),
           vertical = T, add = T, pch = 95, cex = 6, lwd = 2)
abline(a = 0, b = 0, lty = 3)

# Test for significance of difference between groups.

# Create model
AVMod.B <- aov(BD$DeltaBrood~BD$Treatment)

# Inspect model fit & residuals
qqnorm(resid(AVMod.B))
qqline(resid(AVMod.B), col = "blue1", lwd = 2)
par(mfrow = c(2,2))
plot(AVMod.B)
par(mfrow = c(1,1))
hist(resid(AVMod.B))
shapiro.test(resid(AVMod.B))

# Contrary to before, the qq-plot looks pretty decent but the residuals are modestly skewed.
Anova(AVMod.B, type = 2)
# Again, allover less concerning however when we are not claiming a positive result, and graphical interpretation is still most useful for practitioners.

##############
# Final 'colony health' metric we want to take a look at is surival. Few colonies died, so not expecting much.

OBData.D <- STFull[which(STFull$Timepoint != 2),]

# Exclude colonies which were dead already at TP 1
OBData.D <- OBData.D[which(!(OBData.D$Colony %in% OBData.D$Colony[which(OBData.D$Timepoint == 1 & OBData.D$Dead)])),]

# Simply test if it's dead or alive still at Timepoint 3
# in my opinion, easiest way to do this from an R perspective is a quick binomial glm
# so that's what we're doing. it's effectively just a chi-squared test.

Mod.D <- glm(OBData.D$Dead[which(OBData.D$Timepoint == 3)] ~ OBData.D$Treatment[which(OBData.D$Timepoint == 3)],
             family = 'binomial')
Anova(Mod.D, type = 2)

# No difference in colony mortality according to treatment group.
