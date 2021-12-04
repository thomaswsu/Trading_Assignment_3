

rm(list=ls())
setwd("/Users/thomassu/Documents/GitHub/Trading_Assignment_3")

library(quantmod)
library(ggplot2)
library(reshape)

# —————————————————————————-------------------------
tickers <- c(read.csv("tickers.csv"))$Ticker
# tickers <- c("IVW", "SGOV")
P.list <- lapply(tickers, function(x) get(getSymbols(x, from = "2016-12-31")))

sapply(P.list,nrow)

n <- length(P.list)
for(i in 1:n)
{
  print(tickers[i])
  print(head(P.list[[i]], 4))
  print("-----------------------------------------")
}

# Get the adjusted prices into a single object
P.adj <- lapply(P.list, function(p) p[,6])

# Merge the elements of the list
P <- Reduce(merge, P.adj)
names(P) <- tickers

#head(P,10)

#tail(P,10)

# Use mean return to compute expected return 

# Calculate daily return 
returns <- P / lag(P) - 1

# Calculate mean return for all non-empty observations (We remove the the NA)
mu <- apply(returns, 2, function(x) mean(x, na.rm = TRUE)) 

# Calculate the covariance matrix. (Thee built in functions make this so easy...)
sigma <- var(returns, use = "pairwise")
n <- length(mu)

# Annualize the data 
mu <- (1 + mu) ** 252 - 1
sigma <- 252 * sigma

#------------------------------------------------
# Report the mean and covariance Matrix for problem 1 
head(mu, 11)
print(sigma)

muPortfolio <- function(w , mu = mu)
  return(t(w) %*% mu)

sigmaPortfolio <- function(w, sigma = sigma)
  return(t(w) %*% sigma %*% w)

randomPortfolio <- function(mu, sigma)
{
  w <- rep(NA, length = length(mu))
  w[1] <- runif(1)
  for(i in  2:length(mu))
    w[i] <- runif(1, min = 0, max = 1 - sum(w, na.rm = TRUE))
  w <- w / sum(w)
  w <- sample(w)
  mu0 <- muPortfolio(w, mu)
  sigma0 <- sigmaPortfolio(w, sigma)
  return(list("mu" = mu0, "sigma" = sigma0, "w" = w))
}

N <- n * 1000
portfolios <- as.data.frame(matrix(NA, ncol = 2 + n, nrow = N))

for(i in 1:N)
{
  randPortfolio <- randomPortfolio(mu, sigma)
  portfolios[i, ] <- matrix(c(randPortfolio$mu, sqrt(randPortfolio$sigma), t(randPortfolio$w)))
}

# For when we optimize the portfolio
x0 <- rep(0, length = n)
x0[which.max(mu)] <- 1
portfolios <- rbind(portfolios, c(max(mu), sqrt(sigmaPortfolio(x0, sigma)), x0))

# Plot 
p <- ggplot(data = portfolios, aes(x = V2, y = V1)) + geom_point() + 
  ylab(label = "mu Porfolio") + xlab(label = "sigma Portfolio") + ggtitle("Random Portfolios") + 
  expand_limits(x = 0, y = 0)


# Plot all the portfolios
for (i in 1:n){
  wtmp <- rep(0,length=n)
  wtmp[i] = 1
  p <- p + geom_point(x = c(sqrt(sigmaPortfolio(wtmp, sigma))),
                      y = c(muPortfolio(wtmp, mu)),color = "green")
} 
#print(p)

# Set up linear inequality constraints for minimum variance portfolio
ui <- rbind(diag(rep(1, length = n)), rep(1, length = n), rep(-1, length = n))
epsilon = 1e-7
ci <- c(rep(0, length = n), 1, -1) - epsilon

# Objective function to minimize variance 
eval_f <- function(x)
  return(t(x) %*% sigma %*% x)

# -----SET PARAMTERS:
mu0 <- 1e-6
outereps0 <- 1e-8
reltol0 <- 1e-16
eps <- 1e-7
maxN <- 5000
maxNouter <- 1000
abstol0 <- 1e-16 #.Machine$double.eps (We need this value to be a double)

#
x0 <- rep(1/n,length=n)
tmp <- constrOptim(x0, eval_f, ui=ui, ci=ci, mu = mu0,
                   outer.iterations = maxNouter, outer.eps = outereps0, method = "Nelder-Mead", 
                   control=list("abstol"=abstol0,"reltol"=reltol0, "maxit"=maxN))

Rmin <- muPortfolio(tmp$par, mu)
Smin <- sqrt(sigmaPortfolio(tmp$par, sigma))
message(sprintf("Minimum variance portfolio:\n\tR = %5.4f\n\tSigma=%5.4f",Rmin,Smin))

# Plot the minimum variance portfolio (orange):
p <- p + geom_point(x=c(Smin),y=c(Rmin),color="orange") + annotate("segment", x = 0.25, xend = Smin,
                                                                   y = 0.1, yend = Rmin, colour = "orange",
                                                                   size=1, alpha=0.6, arrow=arrow()) + 
  annotate("text", x = c(0.25), y = c(0.1),
           label = c("Min Variance") , color="black", size=4 , 
           angle=0)
#  The plot the maximum return portfolio:
wmax <- rep(0,length=n) 
wmax[which.max(mu)] <- 1
Rmax <- muPortfolio(wmax, mu)
Smax <- sqrt(sigmaPortfolio(wmax, sigma))

p <- p + geom_point(x=c(Smax),y=c(Rmax),color="orange") + annotate("segment", x = 0.5, xend = Smax,
                                                                   y = 0.35, yend = Rmax, colour = "orange",
                                                                   size=1, alpha=0.6, arrow=arrow()) + 
  annotate("text", x = c(0.5), y = c(0.35),
           label = c("Max Return") , color="black", 
           size=4 , angle=0)
#print(p)
message(sprintf("Maximum return portfolio:\n\tR = %5.4f\n\tSigma=%5.4f",Rmax,Smax))


# Plot the efficient frontier 
# No short selling (wi >= 0)

ui <- rbind(diag(rep(1, length = n)), rep(1, length = n), rep(-1, length = n))
ci <- c(rep(0, length = n), 1 - epsilon, -1 - epsilon)

ui0 <- rbind(ui, t(mu))

# Set number of portfolios along [muMin, muMax]
numberOfPortfolios <- 100
soln <- as.data.frame(matrix(NA, ncol = 2 + n, nrow = numberOfPortfolios))
muConst <- seq(Rmin, max(mu), length = numberOfPortfolios)

Rmax <- which.max(mu)
x0 <- rep(0, length = n)
x0[which.max(mu)] <- 1

# Same objective function as previous

# Loop through to solve for efficient frontier 
outMatrix <- as.data.frame(matrix(NA,ncol = 9,nrow = numberOfPortfolios))
outWeights <- as.data.frame(matrix(NA,ncol = (n+1),nrow = numberOfPortfolios))

x0 <- wmax
for(i in 1:numberOfPortfolios)
{
  ci0 <- c(ci,muConst[i]) - epsilon
  tmp <- constrOptim(x0, eval_f, ui=ui0, ci=ci0, mu = mu0,
                     outer.iterations = maxNouter, outer.eps = outereps0, method = "Nelder-Mead", 
                     control=list("abstol"=abstol0, "reltol"=reltol0, "maxit"=maxN))
  wtmp <- tmp$par
  soln[i,2] <- sqrt(sigmaPortfolio(wtmp, sigma))
  soln[i,1] <- muPortfolio(wtmp, mu)
  soln[i,c(3:(n+2))] <- wtmp
  message(sprintf("%d. %6.5f\t%6.5f\t%6.5f\t%d\t%6.5f\t%8.7f\t%d\t%d",
                  i,muConst[i],soln[i,1],soln[i,2], tmp$outer.iterations,sum(wtmp), muConst[i] - soln[i,1], tmp$convergence,
                  tmp$counts[1]))
  outMatrix[i,] <- c(i, muConst[i],soln[i,1],soln[i,2],
                     tmp$outer.iterations,sum(wtmp), muConst[i] - soln[i,1],
                     tmp$convergence,
                     tmp$counts[1])
  outWeights[i,] <- c(i,tmp$par)
}

p <- p + geom_point(data=soln, aes(x=soln[,2],y=soln[,1]),color="blue") 
print(p)

#--------------------------------------
# Plot the expected utility 

# Constrain wi >= 0
ui1 <- cbind(diag(rep(1, length = (n-1))), matrix(rep(0, length=(n-1)), ncol=1))

ui <- rbind(ui1, rep(1,length=n), rep(-1, length=n))

ci <- c(rep(0, length = (n-1)), 1 - epsilon, -1 - epsilon)

lambda <- 1.5 
solution <- as.data.frame(matrix(NA, ncol = 2, nrow = length(lambda)))
x0 <- rep(1 / n, length = n)

eval_f <- function(x)
  return(-t(x) %*% mu + lambda / 2 * (t(x) %*% sigma %*% x))

ci0 <- ci - epsilon

tmp <- constrOptim(x0, eval_f, theta = theta, ui=ui1, ci=ci, mu = mu0,
                   outer.iterations = maxNouter, outer.eps = outereps0, method = "Nelder-Mead", 
                   control=list("abstol"=abstol0, "reltol"=reltol0,"maxit"=maxN))
wtmp <- tmp$par
soln[1,2] <- sqrt(sigmaPortfolio(wtmp, sigma))
soln[1,1] <- muPortfolio(wtmp, mu)
message(sprintf("lambda = %3.1f: R = %5.4f, S = %5.4f", lambda, soln[1,1], soln[1,2]))
print(soln)

# Add results to the plot:
p <- p + geom_point(data=soln, aes(x=soln[,2],y=soln[,1]), shape=9,color="black")
print(p)
