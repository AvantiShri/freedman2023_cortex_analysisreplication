library(lme4)
library(pillar)
library(lmerTest)
library(effectsize)


### TABLE 3 ###

dat1 <- read.csv(file="C://...//TMS_TRIAL_DATA_2023JUL20.CSV",header=TRUE)

### Define variables to be used in model:

dat1$Stim <- ifelse(dat1$tmsgroup==1,"Right",ifelse(dat1$tmsgroup==2,"Left",ifelse(dat1$tmsgroup==3,"Sham","Whale")))
dat1$Stim <- factor(dat1$Stim, levels=c("Sham","Right","Left"))
contrasts(dat1$Stim) <- contr.treatment(3, base=3)

dat1$Direct <- ifelse(dat1$condition==1,"Right",ifelse(dat1$condition==2,"Right",ifelse(dat1$condition==3,"Left",ifelse(dat1$condition==4,"Left","Whale"))))
dat1$Direct <- factor(dat1$Direct, levels=c("Right","Left"))
contrasts(dat1$Direct) <- contr.treatment(2, base=1)

dat1$Intend <- ifelse(dat1$condition==1,"Experimental",ifelse(dat1$condition==2,"Control",ifelse(dat1$condition==3,"Experimental",ifelse(dat1$condition==4,"Control","Whale"))))
dat1$Intend <- factor(dat1$Intend, levels=c("Experimental","Control"))
contrasts(dat1$Intend) <- contr.sum(2)


### Accommodate time between various events:  (i) offset of rTMS, (ii) onset of "experimental" trials, (iii) break between blocks of 100 experimental trials, (iv) break between experimental and control trials.

onset <- 150 # time between TMS administration and start of the first "experimental" session approx 2-3 minutes
iibi <- 12 # Between each block of 100 intention trials is a gap of approx 10 - 15 seconds while we ask if they are ok to continue
icgap <- 60 # Approximately 1 minute between end of intention trials and start of control trials
icbi <- 22 # Between the two blocks of control trials there is a gap of approx 15 - 30 seconds

# Intention trials
dat1$sec[dat1$randomization==1&dat1$condition==3] <- dat1$trial[dat1$randomization==1&dat1$condition==3]-1+onset+I(trunc((dat1$trial[dat1$randomization==1&dat1$condition==3]-1)/100,0))*iibi
dat1$sec[dat1$randomization==2&dat1$condition==1] <- dat1$trial[dat1$randomization==2&dat1$condition==1]-1+onset+I(trunc((dat1$trial[dat1$randomization==2&dat1$condition==1]-1)/100,0))*iibi
dat1$sec[dat1$randomization==2&dat1$condition==3] <- dat1$trial[dat1$randomization==2&dat1$condition==3]-1+onset+500+5*iibi+I(trunc((dat1$trial[dat1$randomization==2&dat1$condition==3]-1)/100,0))*iibi
dat1$sec[dat1$randomization==1&dat1$condition==1] <- dat1$trial[dat1$randomization==1&dat1$condition==1]-1+onset+500+5*iibi+I(trunc((dat1$trial[dat1$randomization==1&dat1$condition==1]-1)/100,0))*iibi

dat1$sec[dat1$randomization==3&dat1$condition==3] <- dat1$trial[dat1$randomization==3&dat1$condition==3]-1+onset+I(trunc((dat1$trial[dat1$randomization==3&dat1$condition==3]-1)/100,0))*iibi
dat1$sec[dat1$randomization==4&dat1$condition==1] <- dat1$trial[dat1$randomization==4&dat1$condition==1]-1+onset+I(trunc((dat1$trial[dat1$randomization==4&dat1$condition==1]-1)/100,0))*iibi
dat1$sec[dat1$randomization==4&dat1$condition==3] <- dat1$trial[dat1$randomization==4&dat1$condition==3]-1+onset+500+5*iibi+I(trunc((dat1$trial[dat1$randomization==4&dat1$condition==3]-1)/100,0))*iibi
dat1$sec[dat1$randomization==3&dat1$condition==1] <- dat1$trial[dat1$randomization==3&dat1$condition==1]-1+onset+500+5*iibi+I(trunc((dat1$trial[dat1$randomization==3&dat1$condition==1]-1)/100,0))*iibi

dat1$sec[dat1$randomization==5&dat1$condition==3] <- dat1$trial[dat1$randomization==5&dat1$condition==3]-1+onset+I(trunc((dat1$trial[dat1$randomization==5&dat1$condition==3]-1)/100,0))*iibi
dat1$sec[dat1$randomization==6&dat1$condition==1] <- dat1$trial[dat1$randomization==6&dat1$condition==1]-1+onset+I(trunc((dat1$trial[dat1$randomization==6&dat1$condition==1]-1)/100,0))*iibi
dat1$sec[dat1$randomization==6&dat1$condition==3] <- dat1$trial[dat1$randomization==6&dat1$condition==3]-1+onset+500+5*iibi+I(trunc((dat1$trial[dat1$randomization==6&dat1$condition==3]-1)/100,0))*iibi
dat1$sec[dat1$randomization==5&dat1$condition==1] <- dat1$trial[dat1$randomization==5&dat1$condition==1]-1+onset+500+5*iibi+I(trunc((dat1$trial[dat1$randomization==5&dat1$condition==1]-1)/100,0))*iibi

# Control trials 
dat1$sec[dat1$randomization==1&dat1$condition==4] <- dat1$trial[dat1$randomization==1&dat1$condition==3]-1+onset+1000+9*iibi+icgap
dat1$sec[dat1$randomization==2&dat1$condition==2] <- dat1$trial[dat1$randomization==2&dat1$condition==1]-1+onset+1000+9*iibi+icgap
dat1$sec[dat1$randomization==2&dat1$condition==4] <- dat1$trial[dat1$randomization==2&dat1$condition==3]-1+onset+1000+9*iibi+icgap+500+icbi
dat1$sec[dat1$randomization==1&dat1$condition==2] <- dat1$trial[dat1$randomization==1&dat1$condition==1]-1+onset+1000+9*iibi+icgap+500+icbi

dat1$sec[dat1$randomization==3&dat1$condition==4] <- dat1$trial[dat1$randomization==3&dat1$condition==3]-1+onset+1000+9*iibi+icgap
dat1$sec[dat1$randomization==4&dat1$condition==2] <- dat1$trial[dat1$randomization==4&dat1$condition==1]-1+onset+1000+9*iibi+icgap
dat1$sec[dat1$randomization==4&dat1$condition==4] <- dat1$trial[dat1$randomization==4&dat1$condition==3]-1+onset+1000+9*iibi+icgap+500+icbi
dat1$sec[dat1$randomization==3&dat1$condition==2] <- dat1$trial[dat1$randomization==3&dat1$condition==1]-1+onset+1000+9*iibi+icgap+500+icbi

dat1$sec[dat1$randomization==5&dat1$condition==4] <- dat1$trial[dat1$randomization==5&dat1$condition==3]-1+onset+1000+9*iibi+icgap
dat1$sec[dat1$randomization==6&dat1$condition==2] <- dat1$trial[dat1$randomization==6&dat1$condition==1]-1+onset+1000+9*iibi+icgap
dat1$sec[dat1$randomization==6&dat1$condition==4] <- dat1$trial[dat1$randomization==6&dat1$condition==3]-1+onset+1000+9*iibi+icgap+500+icbi
dat1$sec[dat1$randomization==5&dat1$condition==2] <- dat1$trial[dat1$randomization==5&dat1$condition==1]-1+onset+1000+9*iibi+icgap+500+icbi


### Linear models

nlmerControl(optimizer = "Nelder_Mead", tolPwrss = 1e-10, optCtrl = list())
lmerControl(check.conv.grad     = .makeCC("warning", tol = 3e-3, relTol = NULL))

contrasts(dat1$Stim)
contrasts(dat1$Intend)
contrasts(dat1$Direct)


### Weighted linear model with random intercept and slope:

# mid-taper about 300th trial or about 480 seconds
wt <- pmax(1*(1+exp((dat1$sec-480)/60))^(-1),2*atan((dat1$sec-1270)/20)/pi)
z2 <- lmer(I(reg-100) ~ Stim*Intend*Direct + (1 + Intend |id), weight = wt, data=dat1, REML = FALSE, control=lmerControl(optimizer = "Nelder_Mead",optCtrl=list(maxfun=40000)))
summary(z2)

# mid-taper about 500th trial or about 704 seconds
wt <- pmax(1*(1+exp((dat1$sec-704)/60))^(-1),2*atan((dat1$sec-1270)/20)/pi)
z2 <- lmer(I(reg-100) ~ Stim*Intend*Direct + (1 + Intend |id), weight = wt, data=dat1, REML = FALSE, control=lmerControl(optimizer = "Nelder_Mead",optCtrl=list(maxfun=40000)))
summary(z2)

# mid-taper about 700th trial or about 928 seconds
wt <- pmax(1*(1+exp((dat1$sec-928)/60))^(-1),2*atan((dat1$sec-1270)/20)/pi)
z2 <- lmer(I(reg-100) ~ Stim*Intend*Direct + (1 + Intend |id), weight = wt, data=dat1, REML = FALSE, control=lmerControl(optimizer = "Nelder_Mead",optCtrl=list(maxfun=50000)))
summary(z2)


### Unweighted linear model with random intercept and slope:

z2 <- lmer(I(reg-100) ~ Stim*Intend*Direct + (1 + Intend |id), data=dat1, REML = FALSE, control=lmerControl(optimizer = "Nelder_Mead",optCtrl=list(maxfun=50000)))
summary(z2)




### TABLE 2 ###

part1 <- read.csv(file="C://...//TMS_PARTICIPANT_DATA_2023JUL20.csv",header=TRUE)

fun1 <- function(ex){
print("MEAN,   STDEV,  MIN, MAX, MISSING")
print(cbind(by(ex,part1$randomization,function(x) mean(x,na.rm=TRUE)),
      by(ex,part1$randomization,function(x) sd(x,na.rm=TRUE)),
      by(ex,part1$randomization,function(x) min(x,na.rm=TRUE)),
      by(ex,part1$randomization,function(x) max(x,na.rm=TRUE)),
      by(ex,part1$randomization,function(x) sum(ifelse(is.na(x)=="TRUE",1,0)))))
print("EFFECT SIZE")
z2 <- aov(ex ~ as.factor(part1$randomization))
print(effectsize(z2, type="eta"))
}

fun1(part1$amt)
fun1(part1$intensity)

