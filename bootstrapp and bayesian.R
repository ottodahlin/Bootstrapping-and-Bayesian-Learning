library(Lock5Data)
library(boot)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(latex2exp)


###############################################################
# Bootstrapping 
###############################################################

# Non-parametric 95% CI for correlation between Bodyfat and Weight.
data("BodyFat")
BodyFat

# Extract variables: Bodyfat and Weight.
Bodyfat <- BodyFat[,1]
Weight <- BodyFat[,3]

# variables:
Bodyfat
Weight

data = data.frame(Bodyfat, Weight)
data

bootcorr <- boot(data, statistic = function(data, i){
  cor(data[i, "Bodyfat"], data[i, "Weight"], method="pearson")  
}, R=1000)

# corr
bootcorr
summary(bootcorr)


# bootstrapped CI 95% different types:
boot.ci(bootcorr, type=c("basic"))

# plot
plot(density(bootcorr$t)) # corr 0.59
#abline(v=0, lty="dashed", col="red")
hist(bootcorr$t)

# corr 0.59
cor.test(data$Bodyfat, data$Weight)


## Version 2: får samma resultat som ovan.
## På annat sätt få fram 95% KI av korrelationen mellan Weight och Bodyfat
results <- boot(data=data, statistic = cor, R=1000)
results

# i = sample, data=data
corboot <- function(data, i){
  cor(data[i, ])[1, 2]
}

results2 <- boot(data=data, statistic = corboot, R=1000)
results2

# 95% CI
boot.ci(results2, type=c("basic"))
# får samma resultat som ovan 



###############################################################
###############################################################

# Predicera med linjär regression BodyFat givet de andra variablerna

data("BodyFat")
BodyFat
ggplot(data, aes(x=Bodyfat, y=Weight)) + geom_point()


# Fitting a linear reg modell
df.BodyFat <- data.frame(BodyFat)
df.BodyFat
colnames(df.BodyFat) <- c("Bodyfat", "Age", "Weight", "Height", "Neck", "Chest", "Abdomen", "Ankle", "Biceps", "Wrist")
reg = summary(lm(Bodyfat ~ ., data=df.BodyFat))
reg

###############################################################
###############################################################

# Beräknar med bootstrap ett 95% KI för Weights effekten

reg
data <- df.BodyFat

boot.body <- function(data, indices) {
  data <- data[indices,] # select obs in bootstrap sample
  mod <- lm(Bodyfat ~ Age + Weight+ Height+ Neck+ Chest+Abdomen+Ankle+Biceps+Wrist, data=data)
  coefficients(mod) # return coefficient vector
}
boot.body

system.time(boot.bodyfat2 <- boot(data, boot.body, R=1000))
boot.bodyfat2 # bootstrapped statistikor

bodyfat.array <- boot.array(boot.bodyfat2)
bodyfat.array [1:2,]

# Plotting:  Weight Coefficient
# OBS streckade linjen som skär igenom koefficient värdet.
plot(boot.bodyfat2, index=3) 


# Weight koefficienten 95% KI
boot.ci(boot.bodyfat2, index=3, type=c("basic"))


###########################################################################
###########################################################################

# Jämför KI ovan med det man får under vanliga normalfördelningsantagandet.
boot.ci(boot.bodyfat2, index=3, type=c("norm"))
boot.ci(boot.bodyfat2, index=3, type=c("basic"))

# BOOT CI:
# Lower: -0.2315  
# Upper: 0.0591


# FRÅN REG OUPUT:
# Coef:  -0.083322
# std error: 0.084706
# z-krit = 1.96 (alpha = 0.05)

# Upper CI: -0.083322 + 1.96*(0.084706)
# Lower CI: -0.083322 - 1.96*(0.084706)

Upper_CI = -0.083322 + 1.96*(0.084706)
Lower_CI = -0.083322 - 1.96*(0.084706)

Upper_CI
Lower_CI

CI = cbind(Lower_CI, Upper_CI)
CI
# Lower: -0.2493458 
# Upper: 0.08270176 


###########################################################################
###########################################################################


# Genes data

gener <- read.csv("gener.csv", header=TRUE)
gener
n <- gener$X
genotype <- gener$genotype
n
genotype

# Final data:
data.df <- data.frame(genotype,n)
data.df

# med ID dubblett kolumn ID
data <- data.frame(genotype)
data

n.AA <- sum(data.df$genotype == "AA")
n.Aa <- sum(data.df$genotype == "Aa")
n.aa <- sum(data.df$genotype == "aa")

# Frequencies:
n.AA # 27 obs of AA
n.Aa # 42 obs of Aa
n.aa # 31 of obs aa

sum(total = n.AA + n.Aa + n.aa)


###########################################################################
###########################################################################

# log likelihooden för p och f givet observationer.

# Frekvens variable (f)
n.AA <- sum(data.df$genotype == "AA")
n.Aa <- sum(data.df$genotype == "Aa")
n.aa <- sum(data.df$genotype == "aa")

# function counts eller data?
# Log Likelihood:
logLFcn <- function(param){
  p <- param[1]
  f <- param[2]
  if ((p>=0 & p<=1) & (f>=0 & f<=1)) {
    n.AA*log( f*p + (1-f)* p^2) + n.Aa*log((1-f)*2*p*(1-p)) + n.aa*log(f*(1-p)+(1-f)*(1-p)^2)
  }
  else{-Inf}
}

count.df <- count(data.df, genotype)
count.df

logLn <- logLFcn(count.df)
logLn


################################################################
################################################################

# Väljer apriori-fördelningen för p och f. Följer Beta fördelningen:

# Alter values
a = 1
b = 1

# Priors för p och f (apriori tätheten)
logPrior <- function(param){
  p <-param[1]
  f <- param[2]
  if ((p>=0 & p<=1) & (f>=0 & f<=1)){
    (a-1)*log(p)+(b-1)*log(1-p) + (a-1)*log(f)+(b-1)*log(1-f)
  }
  else{-Inf}
}

# Posterior (aposteriori tätheten)
logPosterior <- function(param){logLFcn(param) + logPrior(param)}

# tätheten för en Beta fördelning:
# Beta pdf: x^(a-1)*(1-x)^(beta-1)
# (a-1)*log(x)+(beta-1)*log(1-x)


# Våra priors för Beta dist:
# (a-1)*log(p)+(beta-1)*log(1-p)
# (a-1)*log(f)+(beta-1)*log(1-f)


################################################################
################################################################

# Beräknar aposteriorifördelningen

# Markov Chain Monte Carlo (MCMC) function:
mcmc.iter <- function(x, logDensity, sigma, n.iter){
  #Random walk Metropolis Hastings MCMC
  
  res <- matrix(NA, n.iter+1, length(x)) #Create empty matrix
  res[1,] <- x
  logD <- logDensity(x)
  
  accProb <- 0 #keep track of the proportion of proposals that are accepted
  
  for (i in seq_len(n.iter)){
    #New proposal
    xProp <- x + rnorm(length(x), 0, sigma)
    #Log density of proposal
    logDProp <- logDensity(xProp)
    #Acceptance probability
    r <- min( c(1, exp(logDProp - logD) ) )
    
    if(r>runif(1)){ #Accept with probability r, else keep old x
      x <- xProp
      logD <- logDProp
      accProb <- accProb + 1
    }
    res[i+1,] <- x
  }
  list(sample = res, accProb = accProb/n.iter)
}

nIter <- 10000
start_param <- c(0.5,0.5) # start värde på parametern, kan vara vilket värde som helst
set.seed(35)

# Error med "n.AA
# 0.08 = variansen på dragninen.
param.mcmc <- mcmc.iter(start_param,
                        logPosterior,
                        0.05, # variansen i förslagsfördelningen
                        1000)


param.mcmc <- mcmc.iter(param.mcmc$sample[nrow(param.mcmc$sample),],
                        logPosterior,
                        0.05, # variansen i förslagsfördelningen.
                        nIter)


################################################################
################################################################

# Histogram över aposteriori-fördelningen
param.mcmc

df.1 <- data.frame(param.mcmc)
p <- df.1$sample.1
hist(p)

f <- df.1$sample.2
hist(f)

par(mfrow=c(1,2))
hist(df.1$sample.1, xlab="p", main="p")
hist(df.1$sample.2, xlab="f", main="f")

mean(p)
mean(f)

# Figurer, 1000 första dragningarna.

# Plotten för parametern "p"
ggplot(data.frame(p = param.mcmc$sample[1:1000]),
       aes(x = seq(1,1000), y = p )) +
  geom_line() +
  theme_minimal()


ggplot(data.frame(f = param.mcmc$sample[1:1000]),
       aes(x = seq(1,1000), y = f )) +
  geom_line() +
  theme_minimal()



# stabilitets grafen av MCMC
# y-axeln = kumulativa medelvärdet
# x-axeln = iterative sekvenslängden.

# Plottar Random Walk Metropolis Hasting MCMC algoritm nedan:

# PLOT 1. för Parametern "p"
ggplot(data.frame(p = param.mcmc$sample),
       aes(x = seq(0,nIter),
           y = cummean(p))) +
  geom_line() +
  theme_minimal() +
  ylab(TeX("$\\p$"))



# PLOT 2.  for parametern "f"
ggplot(data.frame(f = param.mcmc$sample),
       aes(x = seq(0,nIter),
           y = cummean(f))) +
  geom_line() +
  theme_minimal() + 
  ylab(TeX("$\\f$"))



# alla dragingar syns i nedan kommando:
param.mcmc$sample

# punktskattning nedan:
apply(param.mcmc$sample, 2 , mean)
mean(param.mcmc$sample, 2, mean)


# sd:
sd(param.mcmc$sample)


# 95% Kredibilitetsintervall
apply(param.mcmc$sample, 2, quantile, probs=c(0.05, 0.95))


################################################################
################################################################

# Diagnostik, andel förslag som accepteras.
param.mcmc$accProb

# point estimate
apply(param.mcmc$sample, 2, mean)


# 95% kredibilitetsintervall
apply(param.mcmc$sample, 2, quantile, probs=c(0.05, 0.95))
# column 1 = p
# column 2 = f


################################################################
# END #
################################################################
