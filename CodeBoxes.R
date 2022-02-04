
## BOX ONE
# Data generation
set.seed(7777)
n <- 1000
y <- runif(n, 0, 1)
theoretical_mu <- 0.5 
empirical_mu <- mean(y)
# Functional delta-method: influence curve for the sample mean (first derivative=1(constant))
IF <- 1 * (y - empirical_mu)
mean(IF) #zero by definition
# Plug-in estimation of the sample mean
Yhat <- y + IF # Plug-in estimator
mean(Yhat)
# Geometry of the IF
plot(y, IF)
# Standard Error: Influence Fuction
varYhat.IF <- var(IF) / n # hatvar(hatIF) =1/n sum(yi-bary)^2
seIF <- sqrt(varYhat.IF);seIF
# 0.009161893
# Asymptotic linear inference 95% Confidence Intervals
Yhat_95CI <- c(mean(Yhat) - qnorm(0.975) * sqrt(varYhat.IF), mean(Yhat) + qnorm(0.975) * sqrt(varYhat.IF)); 
mean(Yhat)
## [1] 0.508518
Yhat_95CI
## [1] 0.490561 0.526475

# Check with implemented delta-method in library msm 
library(msm)
se <- deltamethod(g = ~ x1, mean = empirical_mu, cov = varYhat.IF)
se
## [1] 0.009161893

# Check 95%CI delta-method computed by hand with delta-method implemented in RcmdrMisc library
library(RcmdrMisc)
DeltaMethod(lm(y ~ 1), "b0")
#       Estimate SE           2.5%   97.5%
#    b0 0.508518 0.009161893 0.490561 0.526475

## BOX TWO: ratio two mean
# Data generation
library(mvtnorm)
set.seed(123)
sample  <- as.data.frame(rmvnorm(1000, c(3,4), matrix(c(1,0.3,0.3,2), ncol = 2)))
colnames(sample) <- c("X","Y")
sample <- as.data.frame(sample)
# SE estimation for the ratio 
attach(sample)
ratio <- X/Y;mean(ratio)
n <- 1000
a <-    (1 / (mean(Y))^2) * var(X) 
b <-   ((mean(X))^2 / (mean(Y))^4) * var(Y) 
c <-    2 * ((mean(X)) / (mean(Y))^3) * cov(X,Y)
var.IF <- 1/n *(a+b-c); var.IF
SE <- sqrt(var.IF); SE
CI = c(mean(ratio)-qnorm(0.975)*SE,mean(ratio)+qnorm(0.975)*SE); mean(ratio); CI

## BOX THREE
install.packages("epitools")
library(epitools)
RRtable <- matrix(c(60,40,40,60),nrow = 2, ncol = 2)
RRtable
# The next line asks R to compute the RR and 95% confidence interval
riskratio.wald(RRtable)
p1 <- 0.6
p2 <- 0.4
N1 <- 100
N2 <- 100
ratio <- 0.6 / 0.4; ratio
var.IF <- (1 / (p1)^2 * (p1 * (1 - p1)/ N1)) + (1 / (p2)^2 * (p2 * (1 - p2)/ N2));var.IF
SE <- sqrt(var.IF); SE
CI = c(log(ratio)-qnorm(.975)*SE,log(ratio)+qnorm(.975)*SE); ratio; exp(CI)

## BOX FOUR
#This should point to **your** Python path as explained in last section
Sys.setenv(RETICULATE_PYTHON = #"/usr/local/Caskroom/miniconda/base/envs/DeltaMethod/bin/python")
               "/Users/MALF/opt/miniconda3/envs/DeltaMethod/bin/python")
library(caracas)
library(reticulate)
library(MASS)

#Parser for Sympy (you need sympy version 1.10 or 1.9 development)
#check sympy_version() to see you have the right one
sympy_version()

#Create parsers to find the functions check sympy's documentation
#and substitute dots for $ in https://docs.sympy.org/latest/index.html
sympy            <- get_sympy()
Symbol           <- sympy$Symbol
Derivative       <- sympy$derive_by_array
Taylor           <- sympy$series
LaTeX            <- sympy$latex
#Get the variables
x      <- Symbol('x', positive=T)
#Tenth order Taylor
Taylor("(exp(1/x) - 1)/exp(1/x)", x0 = 0, n = 10)$removeO()


## BOX FIVE
# Delta-method for the SE of the correlation between two vectors U and V based on the IF.
#install.packages("psychometric")
#install.packages("MASS")
library(psychometric)
library('MASS')

# Generate the data
samples = 1000
R = 0.83
library('MASS')
data = mvrnorm(n=samples, mu=c(0, 0), Sigma=matrix(c(1, R, R, 1), nrow=2), empirical=TRUE)
X = data[, 1]  # standard normal (mu=0, sd=1)
Y = data[, 2]  # standard normal (mu=0, sd=1)
# Assess that it works
cor(U, V)  # r = 0.83

mu1 = mean(X*Y)
mu2 = mean(X)
mu3 = mean(Y)
mu4 = mean(X^2)
mu5 = mean(Y^2) 

IF1 = X*Y-mu1 
IF2 = X-mu2
IF3 = Y-mu3 
IF4 = X^2-mu4 
IF5 = Y^2-mu5

IF = 
    (sqrt(mu4-mu2^2)*sqrt(mu5-mu3^2))^(-1)*IF1+ 
    (-mu3*sqrt(mu4-mu2^2)*sqrt(mu5-mu3^2)+(mu1-mu2*mu3)*mu2*sqrt(mu5-mu3^2)/sqrt(mu4-mu2^2))/((mu4-mu2^2)*(mu5-mu3^2))*IF2+ 
    (-mu2*sqrt(mu4-mu2^2)*sqrt(mu5-mu3^2)+(mu1-mu2*mu3)* mu3*sqrt(mu4-mu2^2)/sqrt(mu5-mu3^2))/((mu4-mu2^2)*(mu5-mu3^2))*IF3+ 
    (-mu1+mu2*mu3)/(2*(mu4-mu2^2)^1.5*(mu5-mu3^2)^.5)*IF4+ 
    (-mu1+mu2*mu3)/(2*(mu4-mu2^2)^.5*(mu5-mu3^2)^1.5)*IF5

SE = sd(IC)/sqrt(1000); SE

rho_hat = (mu1-mu2*mu3)/(sqrt(mu4-mu2^2)*sqrt(mu5-mu3^2)); rho_hat
CI = c(rho_hat-qnorm(0.975)*SE,rho_hat+qnorm(0.975)*SE); CI
## CI [1] 0.810 0.849

## Using the function CIr from the psychometric package in R: 
CIr(rho_hat,1000,.95)
## [1] 0.8096 0.8483

## BOX SIX: Data generation (simulated example) to apply the Delta-method in a multiple regression 
# Data generation
set.seed(05061972)
N <- 1000
# Age (1: > 65; 0: <= 65)
age <- rbinom(N,1,0.6)                                    
# SES (1: middle class and higher)
SES <- rbinom(N,1,plogis(0.97 + 0.03*age))  
# comorbidities (1: any, 0: none)
morbidity <- rbinom(N,1,plogis(0.7 + 0.010*age - 0.15*SES))     
# stage (1: advanced, TNS=3 or 4; 0: TNS= 1 or 2)
stage <- rbinom(N,1,plogis(-1.0 + 0.1*age - 0.05*SES + 0.2*morbidity)) 
# treatment (1: dual; 0=mono)
treat <- rbinom(N,1,plogis(0.35 + 0.015*stage - 0.15*age + 0.015*SES - 0.3*morbidity)) 
# Counterfactual outcome under A=1 and A=0 respectively
death.1 <-  rbinom(N,1,plogis(-2 - 1*1 + 0.65*age - 0.025*SES + 0.25*morbidity + 0.75*stage + 0.025*1*SES - 0.35*1*morbidity))
death.0 <-  rbinom(N,1,plogis(-2 - 1*0 + 0.65*age - 0.025*SES + 0.25*morbidity + 0.75*stage + 0.025*0*SES - 0.35*0*morbidity))
# Observed outcome: mortality at 1 year after treatment initiation (1: death)
death <- death.1*treat + death.0*(1 - treat)
# OR -> if treated, lower prob of death
exp(coefficients(glm(death~treat,family=binomial)))
# RR
exp(coefficients(glm(death ~ treat + age + SES + morbidity + stage,family=poisson)))
# risk differences
mean(death.1-death.0)
## -0.124

## BOX SEVEN
# Delta method to derive the SE for the conditional RR
data  <- as.data.frame(cbind(death , treat , age , SES , morbidity , stage))
m1 <- glm(death ~ age + treat, family = binomial, data = data)
pMono <- predict(m1, newdata = data.frame(age = 1, treat = 0), type = "response")
pDual <- predict(m1, newdata = data.frame(age = 0, treat = 1), type = "response")
rr <- pMono / pDual
cat("Conditional risk ratio: ", rr)
# Conditional Risk Ratio:  4.406118

# The partial derivative are computed in R as follows:
x1 <- 1
x2 <- 0
x3 <- 1
x4 <- 0
b0 <- coef(m1)[1]
b1 <- coef(m1)[2]
b2 <- coef(m1)[3]
e1 <- exp(- b0 - 1 * b1 - 0 * b2)
e2 <- exp(- b0 - 0 * b1 - 1 * b2)
p1 <- 1 / (1 + e1)
p2 <- 1 / (1 + e2)
dfdb0 <- - e2 * p1 + (1 + e2) * p1 * (1 - p1)
dfdb1 <- - x2 * e2 * p1 + (1 + e2) * x1 * p1 * (1 - p1)
dfdb2 <- - x4 * e2 * p1 + (1 + e2) * x3 * p1 * (1 - p1)
grad <- c(dfdb0, dfdb1, dfdb2)
vG <- t(grad) %*% vcov(m1) %*% (grad)
se_rr <- c(sqrt(vG))
se_rr
## [1] 1.038303

# Check with implemented delta-method in library msm 
library(msm)
se_rr_delta <- deltamethod( ~ (1 + exp(- x1 - 0 * x2 - 1 * x3)) /
                                (1 + exp(- x1 - 1 * x2 - 0 * x3)), 
                            c(b0, b1, b2), 
                            vcov(m1)
)
se_rr_delta
## [1] 1.038303
# We obtain the same results for the SE of the RR computed before

# Finally, we compute the type Wald 95% CI
lb <- rr - 1.96 * sqrt(vG)
ub <- rr + 1.96 * sqrt(vG)
# Conditional Risk Ratio (95%CI)
c(lb, ub)
## [1] 2.386165 6.426072
