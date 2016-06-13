
# -------------------------------------------------------------------------------
# Published in:     Statistical Programming Languages
# -------------------------------------------------------------------------------
# Description:      The following code estimates copula and margins, simulates
#                   multivariate meta-distributions, calculates VaR and 
#                   Expected Shortfall for different models and creates the
#                   corresponding plots. 
# -------------------------------------------------------------------------------
# Keywords:         copula, Gaussian, Student t, Clayton, VaR, 
#                   Expected Shortfall 
# -------------------------------------------------------------------------------
# Author:           Georg Keilbar, Karolina Stanczak

library(VineCopula)
library(copula)
library(scatterplot3d)
db = read.csv("/Users/YinanWu/Downloads/a1.csv", sep= ',')
cb = read.csv("/Users/YinanWu/Downloads/g1.csv", sep = ',')

data        = data.frame(db$Date, as.numeric(db$Adj.Close), as.numeric(cb$Adj.Close))
names(data) = c("Date", "DB" , "CB")
date        = as.Date(data$Date[2:755],"%Y-%m-%d")

# Calculate log returns
db_lr = diff(log(data$DB))
cb_lr = diff(log(data$CB))
lr    = data.frame(db_lr, cb_lr)


# Plot the returns together
plot(db_lr, cb_lr, pch = ".", xlab = "DB", ylab = "CB")
abline(lm(cb_lr~db_lr),col='red',lwd=1)
cor(db_lr,cb_lr,method='spearman')

# Calculate mean and variance
db_mu = mean(db_lr)
db_sd = sd(db_lr)
cb_mu = mean(cb_lr)
cb_sd = sd(cb_lr)

# Plot data with histograms
xhist = hist(lr[,1], breaks=50, plot=FALSE)
yhist = hist(lr[,2], breaks=50, plot=FALSE) 
top = max(c(xhist$counts, yhist$counts)) 
xrange = c(-0.3,0.3)
yrange = c(-0.3,0.3)
nf = layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE) 
par(mar=c(3,3,1,1)) 
plot(lr[,1], lr[,2], xlim=xrange, ylim=yrange, xlab="", ylab="") 
par(mar=c(0,3,1,1)) 
barplot(xhist$counts, col = "blue", axes=FALSE, ylim=c(0, top), space=0) 
par(mar=c(3,0,1,1)) 
barplot(yhist$counts, col = "yellow", axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)



dev.off()

#Convert log-returns into pseudo observations
m = pobs(lr)

######################### Gaussian Copula #################################

# Estimate and simulate Gaussian copula
normcop = normalCopula(dim=2)
set.seed(123)
normfit = fitCopula(normcop,m,method="ml")

# Find and save the coefficients
normrho = coef(normfit)

# Create a contour plot of a gaussian copula
contour(normalCopula(dim=2,normrho),dCopula, main = "Gaussian")

# Create a perspective plot and a scatter plot for the Gaussian Copula
persp(normalCopula(dim=2,normrho),dCopula, phi = 15, theta = 30)
u1 = rCopula(4135, normalCopula(dim=2,normrho))
plot(u1[,1],u1[,2],pch='.',col='blue', main = "Gaussian",xlab = "u1", ylab = "u2")
cor(u1,method='spearman') # Okay - correlation between the logreturn was about 0.6

######################### Student t Copula #################################

# Estimate and simulate t Copula
tcop = tCopula(dim=2)
set.seed(123)
tfit = fitCopula(tcop,m,method='ml')

# Find and save the coefficients
coef(tfit)
trho = coef(tfit)[1]
tdf  = coef(tfit)[2]

# Create a contour plot of a gaussian copula
contour(tCopula(trho, dim = 2, df = tdf),dCopula, main = "Student t")

# Create a perspective plot and a scatter plot for the t Copula
persp(tCopula(dim=2,trho,df=tdf),dCopula, phi = 15, theta = 30)
u2 = rCopula(4135, tCopula(dim=2,trho,df=tdf))
plot(u2[,1],u2[,2],pch='.',col='blue', main = "Student t",xlab = "u1", ylab = "u2") # Looks good - our correlation was about 0.6
cor(u2,method='spearman') # Okay - correlation between the logret was about 0.6

######################### Clayton Copula #################################

# Estimate and simulate Clayton copula
claytoncop = claytonCopula(dim=2)
clayfit    = fitCopula(claytoncop,m,method="ml")

# Find and save the coefficients
phi = coef(clayfit)

# Create a contour plot of a gaussian copula
contour(claytonCopula(phi, dim=2),dCopula, main = "Clayton")

# Create a perspective plot and a scatter plot for the Clayton Copula
persp(claytonCopula(dim=2,phi),dCopula, phi = 15, theta = 30)
u3 = rCopula(4135, claytonCopula(dim=2,phi))
plot(u3[,1],u3[,2],pch=".",col="blue", main = "Clayton",xlab = "u1", ylab = "u2")
cor(u3,method='spearman')


########################### Model for Margins ####################################

#Normal margins

############################ Gaussian Copula ######################################
normcopula_dist = mvdc(copula=normalCopula(normrho,dim=2), margins=c("norm","norm"),
                    paramMargins=list(list(mean=db_mu, sd=db_sd),
                                      list(mean=cb_mu, sd=cb_sd)))

normsim = rMvdc(10000, normcopula_dist)

# Create a contour plot
contour(normcopula_dist, dMvdc, main='Gaussian', xlim = c(-0.07, 0.07), ylim = c(-0.07, 0.07))

# Plot the observed and copula based simulated returns
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
plot(db_lr,cb_lr,main='Gaussian', pch = ".", xlab = "DB", ylab = "CB")
points(normsim[,1],normsim[,2],col=alpha('red', alpha = 0.4), pch = ".")
legend("topright", inset = c(-1,0), legend=c('Observed','Simulated'), col = c("black", "red"),pch = 20, ncol=1,cex=1, bty="o")

############################ Student t Copula #####################################
tcopula_dist = mvdc(copula=tCopula(trho,dim=2,df=tdf), margins=c("norm","norm"),
                    paramMargins=list(list(mean=db_mu, sd=db_sd),
                                      list(mean=cb_mu, sd=cb_sd)))
tsim         = rMvdc(10000, tcopula_dist)

# Create a contour plot
par(mar=c(5.1, 4.1, 4.1, 2.1), xpd=TRUE)
contour(tcopula_dist, dMvdc, main='Student t', xlim = c(-0.07, 0.07), ylim = c(-0.07, 0.07))

# Plot the observed and copula based simulated returns
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
plot(db_lr,cb_lr,main='Student t', pch = ".", xlab = "DB", ylab = "CB")
points(tsim[,1],tsim[,2],col=alpha('red', alpha = 0.4), pch = ".")
legend("topright", inset = c(-1,0), legend=c('Observed','Simulated'), col = c("black", "red"),pch = 20, ncol=1,cex=1, bty="o")

############################ Clayton Copula ######################################
claycopula_dist = mvdc(copula=claytonCopula(phi,dim=2), margins=c("norm","norm"),
                     paramMargins=list(list(mean=db_mu, sd=db_sd),
                                       list(mean=cb_mu, sd=cb_sd)))
claysim = rMvdc(10000, claycopula_dist)

# Create a contour plot
par(mar=c(5.1, 4.1, 4.1, 2.1), xpd=TRUE)
contour(claycopula_dist, dMvdc, main = "Clayton",xlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1))

# Plot the observed and copula based simulated returns
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
plot(db_lr,cb_lr,main='Clayton', pch = ".", xlab = "DB", ylab = "CB")
points(claysim[,1],claysim[,2],col=alpha('red', alpha = 0.4), pch = ".")
legend("topright", inset = c(-1,0), legend=c('Observed','Simulated'), col = c("black", "red"),pch = 20, ncol=1,cex=1, bty="o")

# t margins

# Standardize margins
strets = scale(lr)
uret = pnorm(strets)

# Estimate df for t-distribution
l = function(k, x){
  -sum(dt(x, df = k, log = TRUE))
}

# Finding DF for t distributed margins
db_df = optimize(f = l, interval = c(0, 30), x = strets[,1])$minimum 
cb_df = optimize(f = l, interval = c(0, 30), x = strets[,2])$minimum
uret  = cbind(pt(strets[,1], df = db_df), pt(strets[,2], df = cb_df))

# Fit copula

# Gaussian Copula 
normfit_tmar    = fitCopula(normalCopula(normrho, dim=2), uret, method = "ml") 

# Student t Copula 
tfit_tmar = fitCopula(tCopula(trho, dim=2, df = tdf, df.fixed = FALSE), uret, method = "ml")

# Clayton Copula 
clayfit_tmar = fitCopula(claytonCopula(phi, dim=2), uret, method = "ml")

# AIC new
fit = c(fit, normfit_tmar, tfit_tmar, clayfit_tmar)
minAIC(fit, c(1,2,1,1,2,1)) # Student t with t margins!

###################################################################################
#                         Profit and Loss Function                                #
###################################################################################


# Create a portfolio of 0.5 DB and 0.5 CB from the Copula portfolios and the 
# historical portfolio
 
normpf = (normsim[,1]+normsim[,2])/2
tpf    = (tsim[,1]+tsim[,2])/2
claypf = (claysim[,1]+claysim[,2])/2
histpf = (lr[,1]+lr[,2])/2

pf = list(normpf, tpf, claypf, histpf)

# Use Silvermn's Rule to determine bandwith for Kernel Density Estimator 
h = vector(length = 4)
J = vector(length = 4)
for (i in 1: length(h)){
  J[i] = length(pf[[i]])  
} 
for (i in 1:length(h)){
  h[i] = 2.62 * 1.06 * sd(pf[[i]], na.rm = T) * J[i]^(-1/5)
}

# Establish the Densities 
pfdens = list()
for (i in 1:length(h)) {
  pfdens[[i]] = density(pf[[i]], bw = h[i], na.rm = T, kernel = "biweight")
}

# Plot all the densities in one graphic
par(mar = c(6,4,4,2))
col = c("red", "#A7A7A7", "dodgerblue", "forestgreen")
plot(pfdens[[1]], xlim = c(-0.07, 0.07), ylim = c(0, 18), col = col[1], 
     lwd = 3, xlab = "Log-Returns", main = "")

# Add plots for the next maturities
for (i in 2:length(h)) {
  lines(pfdens[[i]], col = col[i], lwd = 3)
}

# 3D plot of risk neutral densities for all maturities
logret_3d = c(pfdens[[1]]$x, pfdens[[2]]$x, pfdens[[3]]$x, pfdens[[4]]$x)
pfdens_3d = c(pfdens[[1]]$y, pfdens[[2]]$y, pfdens[[3]]$y, pfdens[[4]]$y)
port_3d   = c(rep(10, length(pfdens[[1]]$x)), rep(20, length(pfdens[[2]]$x)), 
              rep(30, length(pfdens[[3]]$x)), rep(40, length(pfdens[[4]]$x)))
df_3d     = cbind(port_3d, logret_3d, pfdens_3d)
all_3d    = scatterplot3d(df_3d, type = "h", color = alpha("darkgrey", alpha = 0.7), 
                          main = "Density for different portfolios", 
                          tick.marks = TRUE,  font.main = 2, xlab = "Portfolios", 
                          ylab = "Log-Returns", zlab = "Density", 
                          x.ticklabs = c("Gaussian", "", "Student t",  "","Clayton",  "","Historical"), zlim = c(0, 20), 
                          ylim = c(-0.45, 0.45))


VaR = data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c("VaR5", "VaR1", "VaR0.1"))))
for (i in 1:length(h)){
  VaR[i,] = quantile(pf[[i]], probs = c(0.05,0.01,0.001))  
}



# Plot the Portfolio Distributions with VaR 
port_names = c("Gaussian", "Student t", "Clayton", "Historical")
par(mfrow = c(2,2))
for (i in 1:length(h)){
  plot(pfdens[[i]],lwd=2,xlim=c(-0.10,0.02),ylim=c(0,20), main = paste(port_names[i], "Portfolio"), xlab ="")
  abline(v=VaR[i,2],col="blue",lwd=2)
  abline(v=ES[i,2],col="red",lwd=2)
  text(-0.05, 12,"VaR",cex=0.8,col="blue",font=2)
}
