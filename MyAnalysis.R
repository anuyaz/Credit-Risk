#Install packages
library(timeDate)
library(timeSeries)
library(fBasics)
library(fMultivar)
library(Matrix)
library(mvtnorm)
library(numDeriv)
library(gsl)
library(copula)
library(fCopulae)
library(QRM)
library(MASS)

library(abind)
library(cvar)
library(FactorCopula)
library(fBasics)
library(fMultivar)
library(ggplot2)

library(gumbel)
library(Kendall)
library(MatrixModels)
library(plot3D)
library(polycor)
library(graphics)
library(Matrix)
library(stats)
library(e1071)
library(VineCopula)
dim = 3 # Dimension
# Import the data set
#CreditRisk <- read.table("c:\mydata.csv", header=TRUE,sep=",", row.names="id")
val.ln <- cbind(CreditRisk[,4],CreditRisk[,5],CreditRisk[,6])
val.ln<- as.matrix(val.ln)
n<-nrow(val.ln)
head(val.ln)
summary(val.ln) # descriptive statistics
sd(val.ln[,1]) # Standard deviation
sd(val.ln[,2])
sd(val.ln[,3])
skewness(val.ln[,1]) # Asymmetry
skewness(val.ln[,2])
skewness(val.ln[,3])
par(mfrow=c(1,2)) #Histogram
lo<-layout(matrix(c(1,1,2,2,3,3),byrow=T,ncol=2))
hist(val.ln[,1],50, main="Variation of LGD")
hist(val.ln[,2],50, main="Variation of PD")
hist(val.ln[,3],50, main="Variation of EAD")
Udata <- pobs(val.ln) #Calculation of pseudo observations
#Transformation of the original variables to uniform variables
ro<-cor(val.ln,method = "spearman") # Definition of initial parameters
ro2<-cor(val.ln,method = "kendall") 
rotau<-Spearman(val.ln)
#Define the copula object
#Fit - Gaussian Copula
norm.cop <- normalCopula(c(ro[2], ro[3],ro[6]),dim=3,dispstr = "un")
NormCopEst<- fitCopula(norm.cop, Udata, method="mpl")
logLik(NormCopEst) # Value of the likelihood function at the optimum
AIC(NormCopEst) # AIC value
BIC(NormCopEst) # BIC value
# Tail dependence - Dependency on the queue (symmetric)
norm.cop_est<-normalCopula(NormCopEst@estimate,dim,dispstr ="un")
TailDep<-lambda(norm.cop)
# Fit - t-Student Copula
t.cop <- tCopula(c(ro2[2], ro2[3],ro2[6]),dim=3,dispstr ="un") 
TCopEst<- fitCopula(t.cop, Udata, method = "mpl", estimate.variance=TRUE)
logLik(TCopEst) # Value of the likelihood function at the optimum
AIC(TCopEst) # AIC value
BIC(TCopEst) # BIC value
# Tail dependence - Dependency on the queue (symmetric)
t.cop_est<-tCopula(c(TCopEst@estimate[1],TCopEst@estimate[2],TCopEst@estimate[3]),dim,dispstr = "un",df=TCopEst@estimate[4])
TailDep<-lambda(t.cop_est)
#Estimation of the marginals
# Normal Marginals - fitting Gauss distribution
mod.GAUSS1 <- fitdistr(val.ln[,1], "normal")
AIC(mod.GAUSS1)
BIC(mod.GAUSS1)
mod.GAUSS2 <- fitdistr(val.ln[,2], "normal")
AIC(mod.GAUSS2)
BIC(mod.GAUSS2)
mod.GAUSS3 <- fitdistr(val.ln[,3], "normal")
AIC(mod.GAUSS3)
BIC(mod.GAUSS3)
# Marginal t-Student - fitting Student's t distribution
mod.t1 <- fitdistr(val.ln[,1], "t")
AIC(mod.t1)
BIC(mod.t1)
mod.t2 <- fitdistr(val.ln[,2],"t")
AIC(mod.t2)
BIC(mod.t2)
mod.t3 <- fitdistr(val.ln[,3],"t")
AIC(mod.t3)
BIC(mod.t3)
#Estimate the risk using Monte Carlo Simulation with each copula and each marginal
set.seed (123) # Fixed a seed
r <-100000 # Total number of values to simulate
# Simulated values
Sim_val.ln_normC<-rCopula(r,norm.cop_est)# Fit - Gaussian Copula
Sim_val.ln_tC<-rCopula(r,t.cop_est)# Fit - t-Student Copula
# Simulation of risk factors and calculation of risk
alfa<-c(0.99, 0.995,0.999)
w <- cbind(1,1,1) # Weights of each risk
# GAUSSIAN COPULA
# Normal Marginals
sim.val.ln1 <- qnorm(Sim_val.ln_normC[,1], mean=mod.GAUSS1$estimate[[1]],sd=mod.GAUSS1$estimate[[2]])
sim.val.ln2 <- qnorm(Sim_val.ln_normC[,2], mean=mod.GAUSS2$estimate[[1]],sd=mod.GAUSS2$estimate[[2]])
sim.val.ln3 <- qnorm(Sim_val.ln_normC[,3], mean=mod.GAUSS3$estimate[[1]],sd=mod.GAUSS3$estimate[[2]])
MC.data<-cbind(sim.val.ln1, sim.val.ln2, sim.val.ln3)
MC.Lsim1 <- -(MC.data%*%t(w))
quantile(MC.Lsim1, alfa)#VaR
# Marginals t-Student
sim.val.ln1t <- qst(Sim_val.ln_normC[,1], mu=mod.t1$estimate[[1]],sd=mod.t1$estimate[[2]],df=mod.t1$estimate[[3]])
sim.val.ln2t <- qst(Sim_val.ln_normC[,2], mu=mod.t2$estimate[[1]],sd=mod.t2$estimate[[2]],df=mod.t2$estimate[[3]])
sim.val.ln3t <- qst(Sim_val.ln_normC[,3], mu=mod.t3$estimate[[1]],sd=mod.t3$estimate[[2]],df=mod.t3$estimate[[3]])
MC.data<-cbind(sim.val.ln1t, sim.val.ln2t, sim.val.ln3t)
MC.Lsim2 <- -(MC.data%*%t(w))
quantile(MC.Lsim2, alfa)#VaR
# COPULA t-STUDENT
# Normal Marginals
t.sim.val.ln1 <- qnorm(Sim_val.ln_tC[,1], mean=mod.GAUSS1$estimate[1],sd=mod.GAUSS1$estimate[2])
t.sim.val.ln2 <- qnorm(Sim_val.ln_tC[,2], mean=mod.GAUSS2$estimate[1],sd=mod.GAUSS2$estimate[2])
t.sim.val.ln3 <- qnorm(Sim_val.ln_tC[,3], mean=mod.GAUSS3$estimate[1],sd=mod.GAUSS3$estimate[2])
MC.data<- cbind(t.sim.val.ln1, t.sim.val.ln2,t.sim.val.ln3)
MC.Lsim3 <- -(MC.data%*%t(w))
quantile(MC.Lsim3, alfa)#VaR
# Marginals t-Student
t.sim.val.ln1t <- qst(Sim_val.ln_tC[,1], mu=mod.t1$estimate[[1]],sd=mod.t1$estimate[[2]],df=mod.t1$estimate[[3]])
t.sim.val.ln2t <- qst(Sim_val.ln_tC[,2], mu=mod.t2$estimate[[1]],sd=mod.t2$estimate[[2]],df=mod.t2$estimate[[3]])
t.sim.val.ln3t <- qst(Sim_val.ln_tC[,3], mu=mod.t3$estimate[[1]],sd=mod.t3$estimate[[2]],df=mod.t3$estimate[[3]])
MC.data<- cbind(t.sim.val.ln1t, t.sim.val.ln2t,t.sim.val.ln3t)
MC.Lsim4 <- -(MC.data%*%t(w))
quantile(MC.Lsim4, alfa)#VaR
