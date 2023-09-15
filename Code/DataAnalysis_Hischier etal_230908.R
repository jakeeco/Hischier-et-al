##################################################################################################
######                R Code for Data Analysis Included in the Article:               	  	######
######  					Positive and negative plant-plant interactions 						                ######
######					influence seedling establishment at both high and low elevations		        ###### 
##################################################################################################

##################################################################################################
### Section 1: Information about the script -----
##################################################################################################


############ A. Script information      -----------------------

# Script Content:         This R script explores data of seed germination of 10 plant species in an field experiment with 3 different treatments
#                         at two elevations (1400 m.a.s.l. and 2000 m.a.s.l.) on Calanda Mountain (Switzerland) and related vegetation, 
#                         soil temperature, soil moisture and light data. 
# Needed data:            The data files  HischierEtal_DataForAnalysis_230629.txt
#                                         HOBO_allData_C.Hischier.csv,
#                                         Soil temperature,moisture_C.Hischier.csv,
#                                         Light_C.Hischier.csv                                  
#                         are needed for this script.
# Authors of this script:  Chantal Hischier, Jake Alexander (jake.alexander@usys.ethz.ch) & Evelin Iseli (evelin.iseli@usys.ethz.ch)
# Date created:           06.12.2021
# Last changes:           07.07.2023



############ B. Focal species abbreviations --------------------------------

#### Lowland species:
# Pla.med         Plantago media
# Sca.col         Scabiosa columbaria
# Med.lup         Medicago lupulina
# Bra.pin         Brachypodium pinnatum
# Bro.ere         Bromus erectus

#### Highland species:
# Pla.atr         Plantago atrata
# Sca.luc         Scabiosa lucida
# Lot.alp         Lotus alpinus
# Poa.alp         Poa alpinus
# Ses.cae         Sesleria caerulea



############ C. Sites  ----------

# Low elevation site (at 1400 m.a.s.l.)         Nesselboden
# High elevation site (at 2000 ma.s.l.)         Calanda



##################################################################################################
### Section 2: Package installation, functions and data import -----
##################################################################################################


############ Working directory ############

setwd("/Users/alexander/polybox2/Shared/Paper Hischier et al./Data and code")
setwd("/Users/eviseli/Desktop/Submission 2_code and data") 


############ Packages ############

library("tidyverse")      
library("car")   				#for Anova on glmer objects (establishment models)
library("MuMIn")				#for AICc
library("lme4")					#for glmer
library("lmerTest")   	#for tests based on lmer models (env data)
library("raster"); library("sp") # for climate data handling (CHELSA)


############ Functions ############

#function to produce model-checking plots for the fixed effects of an lmer model
fix.check <- function(mod){
  par(mfrow = c(1,3))
  plot(fitted(mod),resid(mod),main="Scale-location plot") 	#should have no pattern
  print(anova(lm(fitted(mod)~resid(mod)))) 					#should be non-significant
  qqnorm(resid(mod), ylab="Residuals") 						#should be approximately straight line
  qqline(resid(mod))
  plot(density(resid(mod))) 								#should be roughly normally distributed
  rug(resid(mod))}

# see https://github.com/lme4/lme4/issues/220
overdisp_fun <- function(model) {
  ## number of variance parameters in
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

se <- function (x) sqrt(var(x,na.rm=T)/ length(x[is.na(x)==FALSE]) )



seed.germination.max_data <- read.table("HischierEtal_DataForAnalysis_230629.txt", header=T)
head(seed.germination.max_data)
  

##################################################################################################
### Section 3: Data analysis - establishment (main analysis) -----
##################################################################################################

############ Summary Statistics ############


# mean of germination rate per site
germinationrate.mean.site <- tapply(seed.germination.max_data$germinationrate, seed.germination.max_data$site, mean)
germinationrate.mean.site  
germinationrate.mean.site [1] / germinationrate.mean.site  [2]         ##  1.356296     #35.6% higher germination rate in Cal (high site)                             
tapply(seed.germination.max_data$germinationrate, seed.germination.max_data$site, se)


# mean of germination rate per functional group and site
round(tapply(seed.germination.max_data$germinationrate, 
                                              list(seed.germination.max_data$focalgroup,
                                                   seed.germination.max_data$site), mean) ,3)
round(tapply(seed.germination.max_data$germinationrate, 
                                              list(seed.germination.max_data$focalgroup,
                                                   seed.germination.max_data$site), se) ,3)                                      

# mean of germination rate per functional group and site and treatment
round(tapply(seed.germination.max_data$germinationrate, 
                                              list(seed.germination.max_data$focalgroup, seed.germination.max_data$treatment, seed.germination.max_data$site), mean) ,3)
round(tapply(seed.germination.max_data$germinationrate, 
                                              list(seed.germination.max_data$focalgroup, seed.germination.max_data$treatment, seed.germination.max_data$site), se) ,3)

# mean of germination rate per functional group and site, excluding Brachipodium and Bromus
datsub <- seed.germination.max_data[! seed.germination.max_data$species %in% c("Bra.pin", "Bro.ere"),]
round(tapply(datsub$germinationrate, list(datsub$focalgroup, datsub$treatment, datsub$site), mean) ,3)
round(tapply(datsub$germinationrate, list(datsub$focalgroup, datsub$treatment, datsub$site), se) ,3)





# ii) species origin (lowland, alpine)
# mean germination rate per site
germinationrate.mean.site <- tapply(seed.germination.max_data$germinationrate, 
                                    list(seed.germination.max_data$site,
                                         seed.germination.max_data$focalorigin), mean)
germinationrate.mean.site 

tapply(seed.germination.max_data$germinationrate, 
                                    list(seed.germination.max_data$site,
                                         seed.germination.max_data$focalorigin), se)

# mean germination rate per site and treatment
germinationrate.mean.site.treatment <- tapply(seed.germination.max_data$germinationrate, 
                                              list(seed.germination.max_data$site,
                                                   seed.germination.max_data$treatment,
                                                   seed.germination.max_data$focalorigin), mean)
germinationrate.mean.site.treatment 
  
tapply(seed.germination.max_data$germinationrate, 
                                              list(seed.germination.max_data$site,
                                                   seed.germination.max_data$treatment,
                                                   seed.germination.max_data$focalorigin), se)



############ Statistical Analysis ############     

# i) elevation origin * site

germination.speciestog_glmm <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site 
                                     + (1|species) + (1|plot) + (1|block),family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))  
summary(germination.speciestog_glmm)
isSingular(germination.speciestog_glmm, tol = 1e-4)                                 # to test model for singularity - is singular, drop block (0 variance)

germination.speciestog_glmm <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site 
                                     + (1|species) + (1|plot),family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))  # GLMM
summary(germination.speciestog_glmm)	
overdisp_fun(germination.speciestog_glmm)	#model is overdispersed
isSingular(germination.speciestog_glmm, tol = 1e-4)                                 

## add an observation level random effect (toothpick) https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
germination.speciestog_glmm <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site 
                                     + (1|species) + (1|plot) + (1|toothpick),family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))  # GLMM
summary(germination.speciestog_glmm) #model converges
Anova(germination.speciestog_glmm)

# Analysis of Deviance Table (Type II Wald chisquare tests)

# Response: cbind(seedlings, meanseeds - seedlings)
                # Chisq Df Pr(>Chisq)    
# treatment      46.308  2  8.797e-11 ***
# site           15.767  1  7.163e-05 ***
# treatment:site 26.892  2  1.447e-06 ***

AICc(germination.speciestog_glmm)	# 5724.927
aovtab <- as.data.frame(Anova(germination.speciestog_glmm))
aovtab$Chisq <- round(aovtab$Chisq, 2)
aovtab[,3] <- round(aovtab[,3], 3)
aovtab[,3][aovtab[,3] < 0.001] <- "<0.001"

aovtab.base <- aovtab


#same model but with grasses exluded - we see general pattern of facilitation, with slightly larger effects of artificial veg in both sites
germination.speciestog_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site 
                                     + (1|species) + (1|plot) + (1|toothpick),family = "binomial", subset=focalgroup!="grass", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))  # GLMM
summary(germination.speciestog_glmm2) #model converges
Anova(germination.speciestog_glmm2)

# Response: cbind(seedlings, meanseeds - seedlings)
                 # Chisq Df Pr(>Chisq)    
# treatment      44.3803  2  2.306e-10 ***
# site            0.0956  1     0.7572    
# treatment:site  8.0122  2     0.0182 * 




# ii) species origin (lowland, alpine)
germination.focalorigin_glmm <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                      + (1|species) + (1|plot), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.focalorigin_glmm) 
overdisp_fun(germination.focalorigin_glmm)
isSingular(germination.focalorigin_glmm, tol = 1e-4)                             # to test model for singularity
#model is overdispersed, add observation level random effect

germination.focalorigin_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                      + (1|species) + (1|plot) + (1|toothpick), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.focalorigin_glmm2)  # converges
isSingular(germination.focalorigin_glmm2, tol = 1e-4)                             # to test model for singularity
overdisp_fun(germination.focalorigin_glmm2)
Anova(germination.focalorigin_glmm2)

# Analysis of Deviance Table (Type II Wald chisquare tests)

# Response: cbind(seedlings, meanseeds - seedlings)
                             # Chisq Df Pr(>Chisq)    
# treatment                  45.5562  2  1.281e-10 ***
# site                       15.6150  1  7.764e-05 ***
# focalorigin                 2.3832  1    0.12265    
# treatment:site             26.4347  2  1.819e-06 ***
# treatment:focalorigin      28.8728  2  5.375e-07 ***
# site:focalorigin            1.0275  1    0.31074    
# treatment:site:focalorigin  7.7974  2    0.02027 *  

AICc(germination.focalorigin_glmm2)	#5699.066
aovtab <- as.data.frame(Anova(germination.focalorigin_glmm2))
aovtab$Chisq <- round(aovtab$Chisq, 2)
aovtab[,3] <- round(aovtab[,3], 3)
aovtab[,3][aovtab[,3] < 0.001] <- "<0.001"

aovtab.origin <- aovtab



## ASIDE: fit of this model using default optimiser returned a convergence warning. I use "allFit" to compare different optimisers: (?convergence ; https://joshua-nugent.github.io/allFit/)
## this led me to use the "bobyqa" optimiser for all models
germination.focalorigin_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                      + (1|species) + (1|plot) + (1|toothpick), family = "binomial", data = seed.germination.max_data)
library(dfoptim)
diff_optims <- allFit(germination.focalorigin_glmm2, maxfun=1e5)
diff_optims_OK <- diff_optims[sapply(diff_optims, is, "merMod")]
lapply(diff_optims_OK, function(x) x@optinfo$conv$lme4$messages)
germination.focalorigin_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                      + (1|species) + (1|plot) + (1|toothpick), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
Anova(germination.focalorigin_glmm2)


# iii) mean seed mass

plot(seed.germination.max_data$meanseedweight  , seed.germination.max_data$meanseedweightstand)

germination.seedm_glmm <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * meanseedweight 
                                      + (1|species) + (1|plot), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.seedm_glmm)  #model is overdispersed
overdisp_fun(germination.seedm_glmm)
isSingular(germination.seedm_glmm, tol = 1e-4)                             # to test model for singularity

germination.seedm_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * meanseedweight 
                                      + (1|species) + (1|plot) + (1|toothpick), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.seedm_glmm2)
overdisp_fun(germination.seedm_glmm2)
isSingular(germination.seedm_glmm2, tol = 1e-4)                             # to test model for singularity
Anova(germination.seedm_glmm2)

# Analysis of Deviance Table (Type II Wald chisquare tests)

# Response: cbind(seedlings, meanseeds - seedlings)
                                # Chisq Df Pr(>Chisq)    
# treatment                     45.2632  2  1.483e-10 ***
# site                          15.6561  1  7.597e-05 ***
# meanseedweight                 2.0418  1    0.15303    
# treatment:site                26.0833  2  2.168e-06 ***
# treatment:meanseedweight      34.2708  2  3.616e-08 ***
# site:meanseedweight            5.5224  1    0.01877 *  
# treatment:site:meanseedweight 40.3471  2  1.733e-09 ***

AICc(germination.seedm_glmm2)	#5657.158
aovtab <- as.data.frame(Anova(germination.seedm_glmm2))
aovtab$Chisq <- round(aovtab$Chisq, 2)
aovtab[,3] <- round(aovtab[,3], 3)
aovtab[,3][aovtab[,3] < 0.001] <- "<0.001"

aovtab.seeds <- aovtab


# iv) functional type
germination.funt_glmm <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalgroup 
                                      + (1|species) + (1|plot), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.funt_glmm)  #model is overdispersed
overdisp_fun(germination.funt_glmm)
isSingular(germination.funt_glmm, tol = 1e-4)                             # to test model for singularity

germination.funt_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalgroup 
                                      + (1|species) + (1|plot) + (1|toothpick), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.funt_glmm2)  
overdisp_fun(germination.funt_glmm2)
isSingular(germination.funt_glmm2, tol = 1e-4)                             # to test model for singularity
Anova(germination.funt_glmm2)

# Analysis of Deviance Table (Type II Wald chisquare tests)

# Response: cbind(seedlings, meanseeds - seedlings)
                            # Chisq Df Pr(>Chisq)    
# treatment                 41.7437  2  8.619e-10 ***
# site                      14.2956  1  0.0001562 ***
# focalgroup                 2.9606  2  0.2275677    
# treatment:site            22.0256  2  1.649e-05 ***
# treatment:focalgroup      46.6352  4  1.816e-09 ***
# site:focalgroup           46.4004  2  8.400e-11 ***
# treatment:site:focalgroup 64.0906  4  3.999e-13 ***

AICc(germination.funt_glmm2)	#5598.497
aovtab <- as.data.frame(Anova(germination.funt_glmm2))
aovtab$Chisq <- round(aovtab$Chisq, 2)
aovtab[,3] <- round(aovtab[,3], 3)
aovtab[,3][aovtab[,3] < 0.001] <- "<0.001"

aovtab.funtype <- aovtab


### Brachypodium and Bromus seem to strongly influence the results - refit model after dropping them...
germination.funt_glmm3 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalgroup 
                                      + (1|species) + (1|plot), family = "binomial", subset= ! species %in% c("Bra.pin", "Bro.ere"), data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))

# germination.funt_glmm3 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalgroup 
                                      # + (1|species) + (1|plot), family = "binomial", subset= focalgroup != "grass", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.funt_glmm3)  
overdisp_fun(germination.funt_glmm3)	#model is overdispersed
isSingular(germination.funt_glmm3, tol = 1e-4)                             # to test model for singularity

germination.funt_glmm4 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalgroup 
                                      + (1|species) + (1|plot) + (1|toothpick), family = "binomial", subset= ! species %in% c("Bra.pin", "Bro.ere"), data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.funt_glmm4)  
overdisp_fun(germination.funt_glmm4)
isSingular(germination.funt_glmm4, tol = 1e-4)                             # to test model for singularity
Anova(germination.funt_glmm4)

# Analysis of Deviance Table (Type II Wald chisquare tests)

# Response: cbind(seedlings, meanseeds - seedlings)
                            # Chisq Df Pr(>Chisq)    
# treatment                 62.6463  2  2.492e-14 ***
# site                       3.4776  1    0.06220 .  
# focalgroup                 0.1354  2    0.93455    
# treatment:site             6.2034  2    0.04497 *  
# treatment:focalgroup      24.8186  4  5.472e-05 ***
# site:focalgroup           36.3364  2  1.287e-08 ***
# treatment:site:focalgroup 24.6490  4  5.918e-05 ***

AICc(germination.funt_glmm4)	#3870.426


############ Anova Tables ############ 


aovtab.funtype
aovtab.seeds
aovtab.origin
aovtab.base

aovtab.all <- rbind(aovtab.base, aovtab.origin, aovtab.seeds , aovtab.funtype)
aovtab.all$model <- c( rep("base",3), rep("origin",7), rep("seeds",7), rep("fun.type",7))
aovtab.all$AIC <- round(c( rep(AIC(germination.speciestog_glmm)	,3), rep(AIC(germination.focalorigin_glmm2)	,7), rep(AIC(germination.seedm_glmm2),7), rep(AIC(germination.funt_glmm2),7)) ,2)
aovtab.all$fixed.effect <- rownames(aovtab.all)
rownames(aovtab.all) <- NULL
aovtab.all <- aovtab.all[,c(4,6,5,1,2,3)]

#write.table(aovtab.all, "anova_table_all.txt")



##################################################################################################
### Section 4: Figures -----
##################################################################################################


############ Labels & Functions ############

# Label graphs with "Low elevation" instead of "Nes" and with "High elevation" instead of "Cal"
site <- list(
  'Nes'="Low elevation",
  'Cal'="High elevation"
)
site_labeller <- function(variable,value){
  return(site[value])
}

makeTransparent85<-function(someColor, alpha=85) # transparent equivalent of specific colour
{newColor<-col2rgb(someColor)
apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                            blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}


bars <- function(x,y,z,c,l){for (k in 1:length(y)) for (i in c(-1, 1)) arrows(x[k], y[k], x[k], y[k]+i*z[k], angle=90, length=0, col=c, lwd=l)}
bars2 <- function(x,y,z1,z2,c, l){for (k in 1:length(y)){ arrows(x[k], y[k], x[k], z1[k], angle=90, length=0, col=c, lwd=l)
                                                       arrows(x[k], y[k], x[k], z2[k], angle=90, length=0, col=c, lwd=l)}}



### Plot including elevation, focal species origin and species seperately (i.e. including i) elevation and ii) origin)
#### Germination rate per species
## calculate mean, min., max. and standard error of max. seed germination data 

seed.germination.max.plot_data <- seed.germination.max_data %>% 
  group_by(species,treatment,site, focalorigin, focalgroup) %>%
  summarize(seedsize = mean(meanseedweight),
  			germinationratemean = mean(germinationrate),
            germinationratemin = min(germinationrate),
            germinationratemax = max(germinationrate),
            germinationratesd = sd(germinationrate),
            germinationratese = se(germinationrate))
seed.germination.max.plot_data <- data.frame(seed.germination.max.plot_data)
 # View(seed.germination.max.plot_data)

seed.germination.max.plot_data$treatment <- as.factor(seed.germination.max.plot_data $treatment)
seed.germination.max.plot_data$treat.num <- as.numeric(seed.germination.max.plot_data$treatment)

############ Figure 2 ############

pdf(file="Fig_establishment by origin_221212.pdf",width= 8,height=5.2, useDingbats=FALSE)


par(mfrow=c(1,2),xaxs="i",yaxs="i",tck=0.02, mar=c(2,2,2,0.5), oma=c(2,2,0,0),bty="l", cex=0.9)


dat <- subset(seed.germination.max.plot_data, seed.germination.max.plot_data$site=="Nes")
	x <- dat$treat.num 
	y <- dat$germinationratemean
	origins <- dat$focalorigin
	xlimit <- c(min(x)-(max(x)-min(x))/20,max(x)+(max(x)-min(x))/20)
  	ylimit <- c(-0.05, 0.8)

 
plot(x ,y, type="n", xlim=xlimit, ylim=ylimit, frame.plot=T, axes=F, main="Low elevation site")
axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = c(1,2,3), cex.axis=0.8, labels = FALSE)
text(x = 1:3 -0.2,  y = ylimit[1] - 0.07, labels = c("Bare soil", "Artificial veg.", "Natural veg."), xpd = NA, srt = 35, cex = 1.1) # Rotate the labels by 35 degrees



lowspp <- unique(dat$species[dat$focalorigin=="lowland"])
for(i in 1:5){
points(c(1:3), y[dat$species== lowspp[i]], type="l" , col=makeTransparent85("brown3"), lty=i) #lowland species
points(c(1:3), y[dat$species== lowspp[i]], pch=16 , col=makeTransparent85("brown3"), cex=0.8)
points(c(1:3), y[dat$species== lowspp[i]], pch=16 , col=makeTransparent85("brown3"), cex=0.8)
bars(c(1:3), y[dat$species== lowspp[i]], dat$germinationratese[dat$species== lowspp[i]], makeTransparent85("brown3"), 1)
}

alpspp <- unique(dat$species[dat$focalorigin=="alpine"])
for(i in c(1:5)){
points(c(1.05:3.05), y[dat$species== alpspp[i]], type="l" , col=makeTransparent85("dodgerblue"), lty=i) #lowland species
points(c(1.05:3.05), y[dat$species== alpspp[i]], pch=16 , col=makeTransparent85("dodgerblue"), cex=0.8)
points(c(1.05:3.05), y[dat$species== alpspp[i]], pch=16 , col=makeTransparent85("dodgerblue"), cex=0.8)
bars(c(1.05:3.05), y[dat$species== alpspp[i]], dat$germinationratese[dat$species== alpspp[i]], makeTransparent85("dodgerblue"), 1)
}

bars(c(1.01:3.01), tapply(y,list(x, origins),mean)[,1], tapply(y,list(x, origins),se)[,1], "brown3", 2)
bars(c(1.06:3.06), tapply(y,list(x, origins),mean)[,2], tapply(y,list(x, origins),se)[,2], "dodgerblue",2)
points(c(1.01:3.01), tapply(y,list(x, origins),mean)[,1], type="l", col="brown3", lwd=4) #lowland species
points(c(1.06:3.06), tapply(y,list(x, origins),mean)[,2], type="l", col="dodgerblue", lwd=4)	#alpine species
points(c(1.01:3.01), tapply(y,list(x, origins),mean)[,1], pch=21, bg="brown3") #lowland species
points(c(1.06:3.06), tapply(y,list(x, origins),mean)[,2], pch=21, bg="dodgerblue")	#alpine species


legend(1,0.8, lty=c(rep(1:5,2),1,1), col=c(rep("brown3",5),rep("dodgerblue",5),"brown3","dodgerblue"), lwd=c(rep(1,10),3,3),
 c(as.expression(bquote(italic("Plantago media"))), 
 as.expression(bquote(italic("Scabiosa columbaria"))), 
 as.expression(bquote(italic("Bromus erectus"))),  
 as.expression(bquote(italic("Brachypodium pinnatum"))), 
 as.expression(bquote(italic("Medicago lupulina"))),
   as.expression(bquote(italic("Plantago atrata"))), 
   as.expression(bquote(italic("Scabiosa lucida"))),
   as.expression(bquote(italic("Poa alpina"))), 
   as.expression(bquote(italic("Sesleria caerulea"))),  
   as.expression(bquote(italic("Lotus alpinus"))),
   "Lowland species", "Alpine species"), cex=0.8, bty="n")



dat <- subset(seed.germination.max.plot_data, seed.germination.max.plot_data$site=="Cal")
	x <- dat$treat.num 
	y <- dat$germinationratemean
	origins <- dat$focalorigin
	xlimit <- c(min(x)-(max(x)-min(x))/20,max(x)+(max(x)-min(x))/20)
  	ylimit <- c(-0.05, 0.8)

 
plot(x ,y, type="n", xlim=xlimit, ylim=ylimit, frame.plot=T, axes=F, main="High elevation site")
axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = c(1,2,3), cex.axis=0.8, labels = FALSE)
text(x = 1:3 -0.2,  y = ylimit[1] - 0.07, labels = c("Bare soil", "Artificial veg.", "Natural veg."), xpd = NA, srt = 35, cex = 1.1) # Rotate the labels by 35 degrees

lowspp <- unique(dat$species[dat$focalorigin=="lowland"])
for(i in 1:5){
points(c(1:3), y[dat$species== lowspp[i]], type="l" , col=makeTransparent85("brown3"), lty=i) #lowland species
points(c(1:3), y[dat$species== lowspp[i]], pch=16 , col=makeTransparent85("brown3"), cex=0.8)
points(c(1:3), y[dat$species== lowspp[i]], pch=16 , col=makeTransparent85("brown3"), cex=0.8)
bars(c(1:3), y[dat$species== lowspp[i]], dat$germinationratese[dat$species== lowspp[i]], makeTransparent85("brown3"), 1)
}

alpspp <- unique(dat$species[dat$focalorigin=="alpine"])
for(i in c(1:5)){
points(c(1.05:3.05), y[dat$species== alpspp[i]], type="l" , col=makeTransparent85("dodgerblue"), lty=i) #lowland species
points(c(1.05:3.05), y[dat$species== alpspp[i]], pch=16 , col=makeTransparent85("dodgerblue"), cex=0.8)
points(c(1.05:3.05), y[dat$species== alpspp[i]], pch=16 , col=makeTransparent85("dodgerblue"), cex=0.8)
bars(c(1.05:3.05), y[dat$species== alpspp[i]], dat$germinationratese[dat$species== alpspp[i]], makeTransparent85("dodgerblue"), 1)
}

bars(c(1.01:3.01), tapply(y,list(x, origins),mean)[,1], tapply(y,list(x, origins),se)[,1], "brown3", 2)
bars(c(1.06:3.06), tapply(y,list(x, origins),mean)[,2], tapply(y,list(x, origins),se)[,2], "dodgerblue",2)
points(c(1.01:3.01), tapply(y,list(x, origins),mean)[,1], type="l", col="brown3", lwd=4) #lowland species
points(c(1.06:3.06), tapply(y,list(x, origins),mean)[,2], type="l", col="dodgerblue", lwd=4)	#alpine species
points(c(1.01:3.01), tapply(y,list(x, origins),mean)[,1], pch=21, bg="brown3") #lowland species
points(c(1.06:3.06), tapply(y,list(x, origins),mean)[,2], pch=21, bg="dodgerblue")	#alpine species


mtext("Maximum establishment rate", side=2, line=0, cex=1.2, outer=TRUE)

dev.off()



############ Figure 4 ############

pdf(file="Fig_establishment by functional type_221212.pdf",width= 8,height=5.2, useDingbats=FALSE)


par(mfrow=c(1,2),xaxs="i",yaxs="i",tck=0.02, mar=c(2,2,2,0.5), oma=c(2,2,0,0),bty="l", cex=0.9)


dat <- subset(seed.germination.max.plot_data, seed.germination.max.plot_data$site=="Nes")
	x <- dat$treat.num 
	y <- dat$germinationratemean
	funtypes <- dat$focalgroup
	xlimit <- c(min(x)-(max(x)-min(x))/20,max(x)+(max(x)-min(x))/20)
  	ylimit <- c(-0.05, 0.8)

 
plot(x ,y, type="n", xlim=xlimit, ylim=ylimit, frame.plot=T, axes=F, main="Low elevation site")
axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = c(1,2,3), cex.axis=0.8, labels = FALSE)
#text(1, line=0, at = c(1,2,3), cex.axis=0.8, labels = c("Soil", "Artificial veg.", "Natural veg."))
text(x = 1:3 -0.2,  y = ylimit[1] - 0.07, labels = c("Bare soil", "Artificial veg.", "Natural veg."), xpd = NA, srt = 35, cex = 1.1) # Rotate the labels by 35 degrees


forbs <- unique(dat$species[dat$focalgroup =="forb"])
for(i in 1:4){
points(c(1:3), y[dat$species== forbs[i]], type="l" , col=makeTransparent85("green"), lty=i) #lowland species
points(c(1:3), y[dat$species== forbs[i]], pch=16 , col=makeTransparent85("green"), cex=0.8)
points(c(1:3), y[dat$species== forbs[i]], pch=16 , col=makeTransparent85("green"), cex=0.8)
bars(c(1:3), y[dat$species== forbs[i]], dat$germinationratese[dat$species== forbs[i]], makeTransparent85("lime green"), 1)
}
grasses <- unique(dat$species[dat$focalgroup =="grass"])
for(i in 1:4){
points(c(1:3), y[dat$species== grasses[i]], type="l" , col=makeTransparent85("orange"), lty=i) 
points(c(1:3), y[dat$species== grasses[i]], pch=16 , col=makeTransparent85("orange"), cex=0.8)
points(c(1:3), y[dat$species== grasses[i]], pch=16 , col=makeTransparent85("orange"), cex=0.8)
bars(c(1:3), y[dat$species== grasses[i]], dat$germinationratese[dat$species== grasses[i]], makeTransparent85("orange"), 1)
}
legumes <- unique(dat$species[dat$focalgroup =="legume"])
for(i in 1:2){
points(c(1:3), y[dat$species== legumes[i]], type="l" , col=makeTransparent85("blue"), lty=i) #lowland species
points(c(1:3), y[dat$species== legumes[i]], pch=16 , col=makeTransparent85("blue"), cex=0.8)
points(c(1:3), y[dat$species== legumes[i]], pch=16 , col=makeTransparent85("blue"), cex=0.8)
bars(c(1:3), y[dat$species== legumes[i]], dat$germinationratese[dat$species== legumes[i]], makeTransparent85("blue"), 1)
}



bars(c(1.01:3.01), tapply(y,list(x, funtypes),mean)[,1], tapply(y,list(x, funtypes),se)[,1], "green", 2)
bars(c(1.06:3.06), tapply(y,list(x, funtypes),mean)[,2], tapply(y,list(x, funtypes),se)[,2], "orange",2)
bars(c(0.96:2.96), tapply(y,list(x, funtypes),mean)[,3], tapply(y,list(x, funtypes),se)[,3], "blue",2)

points(c(1.01:3.01), tapply(y,list(x, funtypes),mean)[,1], type="l", col="green", lwd=4) 
points(c(1.06:3.06), tapply(y,list(x, funtypes),mean)[,2], type="l", col="orange", lwd=4)	
points(c(0.96:2.96), tapply(y,list(x, funtypes),mean)[,3], type="l", col="blue",lwd=4)

points(c(1.01:3.01), tapply(y,list(x, funtypes),mean)[,1], pch=21, bg="green")
points(c(1.06:3.06), tapply(y,list(x, funtypes),mean)[,2], pch=21, bg="orange")
points(c(0.96:2.96), tapply(y,list(x, funtypes),mean)[,3], pch=21, bg="blue")


legend(1,0.8, lty=c(rep(1:4,2),1,2,1,1,1), col=c(rep("green",4),rep("orange",4),rep("blue",2),"green","orange","blue" ), lwd=c(rep(1,10),rep(3,3)),
 c(as.expression(bquote(italic("Plantago media"))), 
    as.expression(bquote(italic("Plantago atrata"))), 
 	as.expression(bquote(italic("Scabiosa columbaria"))), 
   	as.expression(bquote(italic("Scabiosa lucida"))),
	as.expression(bquote(italic("Bromus erectus"))),  
    as.expression(bquote(italic("Poa alpina"))), 
	as.expression(bquote(italic("Brachypodium pinnatum"))), 
    as.expression(bquote(italic("Sesleria caerulea"))),  
 	as.expression(bquote(italic("Medicago lupulina"))),
   	as.expression(bquote(italic("Lotus alpinus"))),
   "Forbs", "Grasses", "Legumes"), cex=0.8, bty="n")


dat <- subset(seed.germination.max.plot_data, seed.germination.max.plot_data$site=="Cal")
	x <- dat$treat.num 
	y <- dat$germinationratemean
	funtypes <- dat$focalgroup
	xlimit <- c(min(x)-(max(x)-min(x))/20,max(x)+(max(x)-min(x))/20)
  	ylimit <- c(-0.05, 0.8)

 
plot(x ,y, type="n", xlim=xlimit, ylim=ylimit, frame.plot=T, axes=F, main="High elevation site")
axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = c(1,2,3), cex.axis=0.8, labels = FALSE)
text(x = 1:3 -0.2,  y = ylimit[1] - 0.07, labels = c("Bare soil", "Artificial veg.", "Natural veg."), xpd = NA, srt = 35, cex = 1.1) # Rotate the labels by 35 degrees

forbs <- unique(dat$species[dat$focalgroup =="forb"])
for(i in 1:4){
points(c(1:3), y[dat$species== forbs[i]], type="l" , col=makeTransparent85("green"), lty=i) #lowland species
points(c(1:3), y[dat$species== forbs[i]], pch=16 , col=makeTransparent85("green"), cex=0.8)
points(c(1:3), y[dat$species== forbs[i]], pch=16 , col=makeTransparent85("green"), cex=0.8)
bars(c(1:3), y[dat$species== forbs[i]], dat$germinationratese[dat$species== forbs[i]], makeTransparent85("lime green"), 1)
}
grasses <- unique(dat$species[dat$focalgroup =="grass"])
for(i in 1:4){
points(c(1:3), y[dat$species== grasses[i]], type="l" , col=makeTransparent85("orange"), lty=i) 
points(c(1:3), y[dat$species== grasses[i]], pch=16 , col=makeTransparent85("orange"), cex=0.8)
points(c(1:3), y[dat$species== grasses[i]], pch=16 , col=makeTransparent85("orange"), cex=0.8)
bars(c(1:3), y[dat$species== grasses[i]], dat$germinationratese[dat$species== grasses[i]], makeTransparent85("orange"), 1)
}
legumes <- unique(dat$species[dat$focalgroup =="legume"])
for(i in 1:2){
points(c(1:3), y[dat$species== legumes[i]], type="l" , col=makeTransparent85("blue"), lty=i) #lowland species
points(c(1:3), y[dat$species== legumes[i]], pch=16 , col=makeTransparent85("blue"), cex=0.8)
points(c(1:3), y[dat$species== legumes[i]], pch=16 , col=makeTransparent85("blue"), cex=0.8)
bars(c(1:3), y[dat$species== legumes[i]], dat$germinationratese[dat$species== legumes[i]], makeTransparent85("blue"), 1)
}



bars(c(1.01:3.01), tapply(y,list(x, funtypes),mean)[,1], tapply(y,list(x, funtypes),se)[,1], "green", 2)
bars(c(1.06:3.06), tapply(y,list(x, funtypes),mean)[,2], tapply(y,list(x, funtypes),se)[,2], "orange",2)
bars(c(0.96:2.96), tapply(y,list(x, funtypes),mean)[,3], tapply(y,list(x, funtypes),se)[,3], "blue",2)

points(c(1.01:3.01), tapply(y,list(x, funtypes),mean)[,1], type="l", col="green", lwd=4) 
points(c(1.06:3.06), tapply(y,list(x, funtypes),mean)[,2], type="l", col="orange", lwd=4)	
points(c(0.96:2.96), tapply(y,list(x, funtypes),mean)[,3], type="l", col="blue",lwd=4)

points(c(1.01:3.01), tapply(y,list(x, funtypes),mean)[,1], pch=21, bg="green")
points(c(1.06:3.06), tapply(y,list(x, funtypes),mean)[,2], pch=21, bg="orange")
points(c(0.96:2.96), tapply(y,list(x, funtypes),mean)[,3], pch=21, bg="blue")

mtext("Maximum establishment rate", side=2, line=0, cex=1.2, outer=TRUE)

dev.off()

############ Figure 3 ############

pdf(file="Fig_establishment by seed size_230629.pdf",width= 8,height=5, useDingbats=FALSE)


par(mfrow= c(1,2),xaxs="i",yaxs="i",tck=0.02, mar=c(2,2,2,0.5), oma=c(2,2,0,0),bty="l", cex=0.9)


dat <- subset(seed.germination.max.plot_data, seed.germination.max.plot_data$site=="Nes")
	x <- dat$seedsize 
	y <- dat$germinationratemean
	treats <- dat$treatment
	xlimit <- c(min(x)-(max(x)-min(x))/20,max(x)+(max(x)-min(x))/20)
  	ylimit <- c(-0.05, 0.8)


###labels = c("bare soil", "artificial vegetation", "natural vegetation"), values=c("burlywood4", "orange", "yellowgreen")) 
 
plot(x ,y, type="n", xlim=xlimit, ylim=ylimit, frame.plot=T, axes=F, main="Low elevation site")
axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = c(0:5), cex.axis=0.8)

mod <- lm(y ~ x*treats, dat)
xs <- seq(xlimit[1], xlimit[2], 0.1)
datnew <- data.frame( x = rep(xs,3), treats = sort(rep(levels(treats), length(xs))))
datnew$y <- predict(mod, newdata=datnew, type="response")
conf_interval <- predict(mod, newdata= datnew, interval="confidence", level = 0.95)

polygon(x = c(datnew$x[datnew$treats=="baresoil"], rev(datnew$x[datnew$treats=="baresoil"])),
        y = c(conf_interval[datnew$treats=="baresoil",2], 
              rev(conf_interval[datnew$treats=="baresoil",3])),
        col =  adjustcolor("burlywood4", alpha.f = 0.10), border = NA)

polygon(x = c(datnew$x[datnew$treats=="fakevege"], rev(datnew$x[datnew$treats=="fakevege"])),
        y = c(conf_interval[datnew$treats=="fakevege",2], 
              rev(conf_interval[datnew$treats=="fakevege",3])),
        col =  adjustcolor("orange", alpha.f = 0.10), border = NA)

polygon(x = c(datnew$x[datnew$treats=="natuvege"], rev(datnew$x[datnew$treats=="natuvege"])),
        y = c(conf_interval[datnew$treats=="natuvege",2], 
              rev(conf_interval[datnew$treats=="natuvege",3])),
        col =  adjustcolor("yellowgreen", alpha.f = 0.10), border = NA)

points(x[dat$treatment =="baresoil"], y[dat$treatment =="baresoil"], pch=21 , bg="burlywood4", lty=i)
points(x[dat$treatment =="fakevege"], y[dat$treatment =="fakevege"], pch=21 , bg="orange", lty=i)
points(x[dat$treatment =="natuvege"], y[dat$treatment =="natuvege"], pch=21 , bg="yellowgreen", lty=i)

lines(datnew$x[datnew$treats=="baresoil"], datnew$y[datnew$treats=="baresoil"], col="burlywood4")
lines(datnew$x[datnew$treats=="fakevege"], datnew$y[datnew$treats=="fakevege"], col="orange")
lines(datnew$x[datnew$treats=="natuvege"], datnew$y[datnew$treats=="natuvege"], col="yellowgreen")

legend(0.5,0.8, lty=rep(1,3), col=c("burlywood4", "orange", "yellowgreen"), lwd=rep(2,3), c("Bare soil", "Artificial vegetation", "Natural vegetation"), cex=0.8, bty="n")



dat <- subset(seed.germination.max.plot_data, seed.germination.max.plot_data$site=="Cal")
	x <- dat$seedsize 
	y <- dat$germinationratemean
	treats <- dat$treatment
	xlimit <- c(min(x)-(max(x)-min(x))/20,max(x)+(max(x)-min(x))/20)
  	ylimit <- c(-0.05, 0.8)

 
plot(x ,y, type="n", xlim=xlimit, ylim=ylimit, frame.plot=T, axes=F, main="High elevation site")
axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = c(0:5), cex.axis=0.8)



mod <- lm(y ~ x*treats, dat)
xs <- seq(xlimit[1], xlimit[2], 0.1)
datnew <- data.frame( x = rep(xs,3), treats = sort(rep(levels(treats), length(xs))))
datnew$y <- predict(mod, newdata=datnew, type="response")
conf_interval <- predict(mod, newdata= datnew, interval="confidence", level = 0.95)

polygon(x = c(datnew$x[datnew$treats=="baresoil"], rev(datnew$x[datnew$treats=="baresoil"])),
        y = c(conf_interval[datnew$treats=="baresoil",2], 
              rev(conf_interval[datnew$treats=="baresoil",3])),
        col =  adjustcolor("burlywood4", alpha.f = 0.10), border = NA)

polygon(x = c(datnew$x[datnew$treats=="fakevege"], rev(datnew$x[datnew$treats=="fakevege"])),
        y = c(conf_interval[datnew$treats=="fakevege",2], 
              rev(conf_interval[datnew$treats=="fakevege",3])),
        col =  adjustcolor("orange", alpha.f = 0.10), border = NA)

polygon(x = c(datnew$x[datnew$treats=="natuvege"], rev(datnew$x[datnew$treats=="natuvege"])),
        y = c(conf_interval[datnew$treats=="natuvege",2], 
              rev(conf_interval[datnew$treats=="natuvege",3])),
        col =  adjustcolor("yellowgreen", alpha.f = 0.10), border = NA)


points(x[dat$treatment =="baresoil"], y[dat$treatment =="baresoil"], pch=21 , bg="burlywood4", lty=i)
points(x[dat$treatment =="fakevege"], y[dat$treatment =="fakevege"], pch=21 , bg="orange", lty=i)
points(x[dat$treatment =="natuvege"], y[dat$treatment =="natuvege"], pch=21 , bg="yellowgreen", lty=i)

lines(datnew$x[datnew$treats=="baresoil"], datnew$y[datnew$treats=="baresoil"], col="burlywood4")
lines(datnew$x[datnew$treats=="fakevege"], datnew$y[datnew$treats=="fakevege"], col="orange")
lines(datnew$x[datnew$treats=="natuvege"], datnew$y[datnew$treats=="natuvege"], col="yellowgreen")


mtext("Maximum establishment rate", side=2, line=0, cex=1.2, outer=TRUE)
mtext("Seed size (mg)", side=1, line=0.5, cex=1.2, outer=TRUE)

dev.off()


##################################################################################################
############ Section 5: Environmental data -----
##################################################################################################


############ Load Data ############

###### HOBO soil temperature data
HOBO_data <- read.csv("HOBO_allData_C.Hischier.csv", sep=";")
  # View(HOBO_data)

###### Moisture data of one day (20.08.2021)
soil.temperature.moisture_data <- read.csv("Soil temperature,moisture_C.Hischier.csv", sep=";")
  # View(soil.temperature.moisture_data)


###### Light data of one day (01.09.2021)
light_data <- read.csv("Light_C.Hischier.csv", sep=";")
  # View(light_data)


############ Soil Temperature ############
    
HOBO_data$site <- factor(HOBO_data$site,levels=c("Nes", "Cal"))                  # to define the order of "site" factors in plots


# define date as type "Date" (so that you can change the order of day, month, year also later in plots)
HOBO_data <- HOBO_data %>%
  mutate(date = as.Date(date, format = "%d.%m.%Y")) %>%
  mutate(block = paste(site, block, sep="."))
   # View(HOBO_data)

# get rid of one NA value in Cal fakevege
HOBO_data <- HOBO_data %>%
  drop_na(soiltemp)

# calculate mean, min., max. and standard error data (date, treatment, site: per day --> without date: for the entire vegetation period)
HOBO.mean.min.max.se_data <- HOBO_data %>% 
  group_by(treatment,site) %>%
  summarize(soiltempmean = mean(soiltemp),
            soiltempmin = min(soiltemp),
            soiltempmax = max(soiltemp),
            soiltempsd = sd(soiltemp),
            soiltempse = se(soiltemp),
            date_start = min(date),
            date_end = max(date))
HOBO.mean.min.max.se_data <- data.frame(HOBO.mean.min.max.se_data)
  # View(HOBO.mean.min.max.se_data)

with(HOBO.mean.min.max.se_data, tapply(soiltempmean ,site,mean))  #difference is 3.64638C. 15.16290 - 11.51652

#### test for significance of differences in soil temperatures between sites and treatments
hobomod <- lmer(soiltemp ~ treatment * site + (1|hobo), data = HOBO_data) #logger ID accounts for block effects (2 loggers per treatment)
fix.check(hobomod) #uneven variance
hobomod <- lmer(sqrt(soiltemp) ~ treatment * site + (1|hobo), data = HOBO_data) #logger ID accounts for block effects (2 loggers per treatment)
fix.check(hobomod) #variance more even
                                      
summary(hobomod) 
Anova(hobomod)

# Analysis of Deviance Table (Type II Wald chisquare tests)

# Response: sqrt(soiltemp)
#                   Chisq     Df   Pr(>Chisq)   
#
# treatment         18.0264   2    0.0001218 ***
#  site             379.9418  1    < 2.2e-16 ***
#  treatment:site   4.4314    2    0.1090775     


aovtab <- as.data.frame(Anova(hobomod))
aovtab$Chisq <- round(aovtab$Chisq, 2)
aovtab[,3] <- round(aovtab[,3], 3)
aovtab[,3][aovtab[,3] < 0.001] <- "<0.001"

aovtab.hobo <- aovtab



############ Soil Moisture ############

soil.temp.moist.02.09.2021_data <- subset(soil.temperature.moisture_data, date == "02.09.2021")    # subset data of 02.09.2021
# View(soil.temp.moist.02.09.2021_data)



## calculate mean soil moisture
# calculate mean, min., max. and standard error data per day 
soil.temp.moist.mean.02.09.2021_data <- soil.temp.moist.02.09.2021_data %>% 
  group_by(treatment,site) %>%
  summarize(soilmoistmean = mean(soilmoist),
            soilmoistmin = min(soilmoist),
            soilmoistmax = max(soilmoist),
            soilmoistsd = sd(soilmoist),
            soilmoistse = se(soilmoist))
soil.temp.moist.mean.02.09.2021_data <- data.frame(soil.temp.moist.mean.02.09.2021_data)
# View(soil.temp.moist.mean.02.09.2021_data)

soil.temp.moist.02.09.2021_data$block <- as.factor(paste(soil.temp.moist.02.09.2021_data$site, soil.temp.moist.02.09.2021_data$plot, sep="."))
soil.temp.moist.02.09.2021_data$plot <- as.factor(paste(soil.temp.moist.02.09.2021_data$site, soil.temp.moist.02.09.2021_data$plot, substr(soil.temp.moist.02.09.2021_data$treatment,1,3), sep="."))


#### test for significance of differences in soil moisture between sites and treatments
moistmod <- lmer(soilmoist ~ treatment * site + (1|plot) , data = soil.temp.moist.02.09.2021_data) #can't fit block due to singularity
fix.check(moistmod) #looks OK
summary(moistmod) 
Anova(moistmod)

# Analysis of Deviance Table (Type II Wald chisquare tests)

# Response: soilmoist
# Chisq Df Pr(>Chisq)    
# treatment      146.113  2  < 2.2e-16 ***
# site            40.008  1  2.529e-10 ***
# treatment:site  11.711  2   0.002863 ** 


aovtab <- as.data.frame(Anova(moistmod))
aovtab$Chisq <- round(aovtab$Chisq, 2)
aovtab[,3] <- round(aovtab[,3], 3)
aovtab[,3][aovtab[,3] < 0.001] <- "<0.001"

aovtab.moist <- aovtab


############ Light ############

light.01.09.2021_data <- subset(light_data, date == "01.09.2021")    # subset data of 01.09.2021

light.01.09.2021_data <- light.01.09.2021_data %>%
  mutate(block = paste(site, plot, sep=".")) %>%
  mutate(plot = paste(site, plot, substr(treatment,1,3),sep=".")) %>%
  mutate(lightreduction_pc = (lightground/lightaboveground)*100)
# View(HOBO_data)

# calculate mean, min., max. and standard error data per day 
light.mean.01.09.2021_data <- light.01.09.2021_data %>% 
  group_by(treatment,site) %>%
  summarize(lightground_pc_mean = mean(lightreduction_pc),
            lightground_pc_se = se(lightreduction_pc),
            lightgroundmean = mean(lightground),
            lightgroundmin = min(lightground),
            lightgroundmax = max(lightground),
            lightgroundsd = sd(lightground),
            lightgroundse = se(lightground))
light.mean.01.09.2021_data <- data.frame(light.mean.01.09.2021_data)

#### test for significance of differences in light interception between sites and treatments
lightmod <- lmer(lightreduction_pc ~ treatment * site + (1|plot) + (1|block) , data = light.01.09.2021_data, subset=treatment != "baresoil") 
fix.check(lightmod) #uneven variance
lightmod <- lmer(sqrt(lightreduction_pc) ~ treatment * site + (1|plot) + (1|block) , data = light.01.09.2021_data, subset=treatment != "baresoil") 
fix.check(lightmod) #fixes model

summary(lightmod) 
Anova(lightmod)

# Analysis of Deviance Table (Type II Wald chisquare tests)

# Response: sqrt(lightreduction_pc)
# Chisq Df Pr(>Chisq)    
# treatment      57.3723  1  3.607e-14 ***
# site            1.4653  1     0.2261    
# treatment:site  0.4144  1     0.5197    

aovtab <- as.data.frame(Anova(lightmod))
aovtab$Chisq <- round(aovtab$Chisq, 2)
aovtab[,3] <- round(aovtab[,3], 3)
aovtab[,3][aovtab[,3] < 0.001] <- "<0.001"
aovtab.light <- aovtab


############ Figure S1 ############

# calculate mean, min., max. and standard error data for each date during the entire vegetation period for temperature
HOBO.dailymean.min.max.se_data <- HOBO_data %>% 
  group_by(date, treatment,site) %>%
  summarize(soiltempmean = mean(soiltemp),
            soiltempmin = min(soiltemp),
            soiltempmax = max(soiltemp),
            soiltempsd = sd(soiltemp),
            soiltempse = se(soiltemp))
HOBO.dailymean.min.max.se_data <- data.frame(HOBO.dailymean.min.max.se_data)
# View(HOBO.mean.min.max.se_data)

# reorder factor levels for all data frames
HOBO.dailymean.min.max.se_data$site <- factor(HOBO.dailymean.min.max.se_data$site , levels=c("Nes", "Cal"))
soil.temp.moist.02.09.2021_data$site <- factor(soil.temp.moist.02.09.2021_data$site , levels=c("Nes", "Cal"))
light.01.09.2021_data$site <- factor(light.01.09.2021_data$site , levels=c("Nes", "Cal"))

# define beginnings of month for axis ticks & labels
xax <- data.frame(start_date = as.Date(c("2021-07-01", "2021-08-01", "2021-09-01", "2021-10-01")),
                  label = c("July", "August", "September", "October"))


# 4-PANEL VERSION
pdf(file="FigS1_temp_moist_230707.pdf",width= 10,height=5, useDingbats=FALSE)

layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)

layout(mat = layout.matrix,
       heights = c(2, 2), # Heights of the two rows
       widths = c(2, 1.2)) # Widths of the two columns


par(xaxs="i",yaxs="i",tck=0.02, mar=c(2,2,2,0.5), oma=c(2,2,0,0),bty="l", cex=0.9) # mfrow= c(2,1),

# plot temperature data
dat <- subset(HOBO.dailymean.min.max.se_data,site == "Nes")
x <- dat$date 
y <- dat$soiltempmean
treats <- dat$treatment
xlimit <- as.Date(c("2021-06-28", "2021-10-28"))
ylimit <- c(2, 27)


plot(x ,y, type="n", ylim=ylimit, xlim=xlimit, frame.plot=T, axes=F, main="Low elevation site")
axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = xax$start_date, cex.axis=0.8, labels=xax$label)
#text(tck, par("usr")[3], labels=xax$label, srt=60,
#     xpd=TRUE, cex=0.8, adj = 1.2)

lines(x[dat$treatment=="baresoil"], y[dat$treatment=="baresoil"], col="burlywood4", lwd = 2)
lines(x[dat$treatment=="fakevege"], y[dat$treatment=="fakevege"], col="orange", lwd = 2)
lines(x[dat$treatment=="natuvege"], y[dat$treatment=="natuvege"], col="yellowgreen", lwd = 2)

legend("topright", col=c("burlywood4", "orange", "yellowgreen"), c("Bare soil", "Artificial vegetation", "Natural vegetation"), cex=0.9, pch = 19, bty="n", xpd = "NA") # lty=rep(1,3), lwd=rep(2,3)
mtext("(a)", 2, adj=2, las=1, padj=-7)

dat <- subset(HOBO.dailymean.min.max.se_data, site == "Cal") 
x <- dat$date 
y <- dat$soiltempmean
treats <- dat$treatment
xlimit <- as.Date(c("2021-06-28", "2021-10-28"))
ylimit <- c(2, 27)


plot(x ,y, type="n", ylim=ylimit, xlim = xlimit, frame.plot=T, axes=F, main="High elevation site")
axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = xax$start_date, cex.axis=0.8, labels=xax$label)

lines(x[dat$treatment=="baresoil"], y[dat$treatment=="baresoil"], col="burlywood4", lwd = 2)
lines(x[dat$treatment=="fakevege"], y[dat$treatment=="fakevege"], col="orange", lwd = 2)
lines(x[dat$treatment=="natuvege"], y[dat$treatment=="natuvege"], col="yellowgreen", lwd = 2)

mtext("Soil Temperature [°C]", side=2, line=-1, cex=0.9, outer=TRUE)
mtext("Month of Year 2021", side=1, line=2, cex=0.9, )

# plot average temperature data

boxplot(soiltempmean ~ treatment*site, data = HOBO.dailymean.min.max.se_data, col = c("burlywood4", "orange", "yellowgreen"), yaxt = "n", xaxt = "n", ylimit = c(2, 36))

axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = c(1:3, 4:6), cex.axis=0.8, labels=c("", "Low elevation\n site", "", "", "High elevation\n site", ""), mgp = c(1.1,1,1))

mtext(side=2, line=0.8, "Soil Temperature [°C]", cex = 0.9)
mtext("(b)", side = 3, line = 0, adj = 1)


# plot soil moisture data

boxplot(soilmoist ~ treatment*site, data = soil.temp.moist.02.09.2021_data , col = c("burlywood4", "orange", "yellowgreen"), yaxt = "n", xaxt = "n", ylimit = c(10, 55))

axis(2, line=0, mgp = c(1.1,0,0), cex.axis=0.8)
axis(1, line=0, at = c(1:3, 4:6), cex.axis=0.8, labels=c("", "Low elevation\n site", "", "", "High elevation\n site", ""), mgp = c(1.1,1,0))

mtext(side=2, line=0.8, "Volumetric Water Content [%]", cex = 0.9)
mtext("(c)", side = 3, line = 0, adj = 1)


dev.off()

############ Summary ############

aovtab.all <- rbind(aovtab.hobo, aovtab.moist, aovtab.light)
aovtab.all$model <- c( rep("Temperature (?C)",3), rep("Volumetric water content (%)",3), rep("Light interception (%)",3))
aovtab.all$fixed.effect <- rep(c("Treatmemnt (T)", "Site (S)", "T x S"),3)
rownames(aovtab.all) <- NULL
aovtab.all <- aovtab.all[,c(4,5,1,2,3)]
names(aovtab.all) <- c("Response variable", "Fixed effect", "Chisq", "df", "P")

write.table(aovtab.all, "anova_table_env.txt")




##################################################################################################
############ Section 6: Climate data -----
##################################################################################################


############ Preparation ############ 

# save coordinates of the high and low elevation site
coord <- data.frame(lon = c(9.490098, 9.489510), lat = c(46.869266, 46.887824))


# turning this into a spatial points dataframe to extract data CHELSA at each coordinate
Coordinates_spatial<-SpatialPoints(coord, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


############ 2000 - 2016 Temperature ############ 

# read timeseries data for tmax
setwd("~/Desktop/Downloaded Data/CHELSA Cruts/tmax") # CHELSA Cruts data can be downloaded on https://chelsa-climate.org/chelsacruts/

# read tifs for timeseries data
rasValue_late_tmax<-c()
rasValue_annual_late_tmax<-c()
for (j in 1:17) #for a period of 17 Years, starting in 1999+1 = 2000
{
  print(j)
  for (i in 1:12) {
    s <- raster(paste(paste("CHELSAcruts_tmax",i, 1999+j,"V.1.0",sep="_"),"tif",sep="."))
    rasValue_late_tmax=cbind(rasValue_late_tmax,raster::extract(s, Coordinates_spatial)) # extract monthly means
  } 
  rasValue_annual_late_tmax=cbind(rasValue_annual_late_tmax,rowMeans(rasValue_late_tmax[,(ncol(rasValue_late_tmax)-11):ncol(rasValue_late_tmax)])) #calculating annual means
}

# plot the temperature data of the last year & the sampling coordinates to check whether they are correct
plot(s) 
points(Coordinates_spatial)

# convert to degrees, as data is saved a bit unusual (°C/10)
rasValue_late_tmax<-rasValue_late_tmax/10
rasValue_annual_late_tmax<-rasValue_annual_late_tmax/10


# add coordinates and site to temperature data (rasValue)
combinePointValue_late_tmax=cbind("site" = c("Nes", "Cal"),coord,rasValue_late_tmax) 

# create sensible column names
column_names<-c()
for (j in 1:17) {
  for (i in 1:12) {
    column_names<-c(column_names,paste("CHELSA",1999+j,i,sep="_"))
  }
}
column_names<-c("site","lon","lat",column_names)
colnames(combinePointValue_late_tmax)<-column_names

# add coordinates and site to yearly temperature data (rasValue_annual)
combinePointValue_annual_late_tmax=cbind("site" = c("Nes", "Cal"),as.data.frame(cbind(coord,rasValue_annual_late_tmax)))

# create sensible column names
column_names<-c()
for (j in 1:17) {
  column_names<-c(column_names,paste("CHELSA",1999+j,sep="_"))
}

column_names<-c("site","lon","lat",column_names)
colnames(combinePointValue_annual_late_tmax)<-column_names

# transform into data frame
Calanda_annual_late_tmax<-as.data.frame(combinePointValue_annual_late_tmax)

# same for tmin
setwd("~/Desktop/Downloaded Data/CHELSA Cruts/tmin") # CHELSA Cruts data can be downloaded on https://chelsa-climate.org/chelsacruts/

# read tifs for timeseries data
rasValue_late_tmin<-c()
rasValue_annual_late_tmin<-c()
for (j in 1:17) #for a period of 17 Years, starting in 1999+1 = 2000
{
  print(j)
  for (i in 1:12) {
    s <- raster(paste(paste("CHELSAcruts_tmin",i, 1999+j,"V.1.0",sep="_"),"tif",sep="."))
    rasValue_late_tmin=cbind(rasValue_late_tmin,raster::extract(s, Coordinates_spatial))
  } 
  rasValue_annual_late_tmin=cbind(rasValue_annual_late_tmin,rowMeans(rasValue_late_tmin[,(ncol(rasValue_late_tmin)-11):ncol(rasValue_late_tmin)])) #calculating annual means
}

# plot the temperature data of the last year & the sampling coordinates to check whether they are correct
plot(s) 
points(Coordinates_spatial)

# convert to degrees, as data is saved a bit unusual (°C/10)
rasValue_late_tmin<-rasValue_late_tmin/10
rasValue_annual_late_tmin<-rasValue_annual_late_tmin/10


# add coordinates and pidn to temperature data (rasValue)
combinePointValue_late_tmin=cbind("site" = c("Nes", "Cal"),coord,rasValue_late_tmin) 

# create sensible column names
column_names<-c()
for (j in 1:17) {
  for (i in 1:12) {
    column_names<-c(column_names,paste("CHELSA",1999+j,i,"tmin", sep="_"))
  }
}
column_names<-c("site","lon","lat",column_names)
colnames(combinePointValue_late_tmin)<-column_names

# add coordinates and pidn to temperature data (rasValue_annual)
combinePointValue_annual_late_tmin=cbind("site" = c("Nes", "Cal"),as.data.frame(cbind(coord,rasValue_annual_late_tmin)))

# create sensible column names
column_names<-c()
for (j in 1:17) {
  column_names<-c(column_names,paste("CHELSA","tmin",1999+j, sep="_"))
}
column_names<-c("site","lon","lat",column_names)
colnames(combinePointValue_annual_late_tmin)<-column_names

# transform into data frame
Calanda_annual_late_tmin<-as.data.frame(combinePointValue_annual_late_tmin)

# merge tmin and tmax by coordinates and site names to calculate mean
Calanda_annual_late <- left_join(Calanda_annual_late_tmax, Calanda_annual_late_tmin, by = c("site", "lon", "lat"))

# take the mean of tmin and tmax
aa <- split.default(Calanda_annual_late[,-c(1:3)], str_remove(names(Calanda_annual_late[,-c(1:3)]), "\\D+")) %>% map_dfc(rowMeans)
Calanda_annual_late_temp <- bind_cols(Calanda_annual_late, aa) %>%
  dplyr::select(!c(4:37)) 


############ 2000 - 2016 Precipitation ############ 

# read timeseries data for precipitation
setwd("~/Desktop/Downloaded Data/CHELSA Cruts/prec") # CHELSA Cruts data can be downloaded on https://chelsa-climate.org/chelsacruts/

# read tifs for timeseries data
rasValue_late_prec<-c()
rasValue_annual_late_prec<-c()
for (j in 1:17) #for a period of 17 Years, starting in 1999+1 = 2000
{
  print(j)
  for (i in 1:12) {
    s <- raster(paste(paste("CHELSAcruts_prec",i, 1999+j,"V.1.0",sep="_"),"tif",sep="."))
    rasValue_late_prec=cbind(rasValue_late_prec,raster::extract(s, Coordinates_spatial)) # extract monthly sums
  } 
  rasValue_annual_late_prec=cbind(rasValue_annual_late_prec,rowSums(rasValue_late_prec[,(ncol(rasValue_late_prec)-11):ncol(rasValue_late_prec)])) # calculating annual sums
}

# plot the temperature data of the last year & the sampling coordinates to check whether they are correct
plot(s)
points(Coordinates_spatial)

# add coordinates and site to temperature data (rasValue)
combinePointValue_late_prec=cbind("site" = c("Nes", "Cal"),coord,rasValue_late_prec) 

# create sensible column names
column_names<-c()
for (j in 1:17) {
  for (i in 1:12) {
    column_names<-c(column_names,paste("CHELSA",1999+j,i,sep="_"))
  }
}
column_names<-c("site","lon","lat",column_names)
colnames(combinePointValue_late_prec)<-column_names

# add coordinates and site to yearly temperature data (rasValue_annual)
combinePointValue_annual_late_prec=cbind("site" = c("Nes", "Cal"),as.data.frame(cbind(coord,rasValue_annual_late_prec)))

# create sensible column names
column_names<-c()
for (j in 1:17) {
  column_names<-c(column_names,paste("CHELSA",1999+j,sep="_"))
}

column_names<-c("site","lon","lat",column_names)
colnames(combinePointValue_annual_late_prec)<-column_names

# transform into data frame
Calanda_annual_late_prec<-as.data.frame(combinePointValue_annual_late_prec)



############ Long-term Sum/ Mean ############ 

# save annual means/ sums as vectors
temp_Nes <- unlist(Calanda_annual_late_temp[Calanda_annual_late_temp$site == "Nes", c(4:20)], use.names = FALSE)
prec_Nes <- unlist(Calanda_annual_late_prec[Calanda_annual_late_prec$site == "Nes", c(4:20)], use.names = FALSE)


temp_Cal <- unlist(Calanda_annual_late_temp[Calanda_annual_late_temp$site == "Cal", c(4:20)], use.names = FALSE)
prec_Cal <- unlist(Calanda_annual_late_prec[Calanda_annual_late_prec$site == "Cal", c(4:20)], use.names = FALSE)

# get mean +- SD for mean annual temp. and annual precipitation for both sites
mean(temp_Nes); sd(temp_Nes)
mean(prec_Nes); sd(prec_Nes)
mean(temp_Cal); sd(temp_Cal)
mean(prec_Cal); sd(prec_Cal)
