context("samMCMC")

# clear the workspace
rm(list=ls())

# Define some functions used in testing
# For calling samMCMC with the direct approach
normLik <- function(x,mu=0,Sig=1,invSig=NA) {
  # To speed computation, invSig (the inverse covariance matrix) can be input
  # instead of Sig (the covariance matrix). Sig is ignored if invSig is input.
  if(all(is.na(invSig))) {
    invSig <- solve(Sig)
  }
  v <- as.matrix(x-mu)
  return(exp(-0.5*t(v) %*% invSig %*% v))
}

# Parameters for testing scalar function
mu_scalar <- -1.5
Sig_scalar <- 0.2
x0_scalar <- 0

# Parameters for testing vector function
mu_vector <- c(-1.5,0.25)
Sig_vector <- diag(c(1,2))
invSig_vector <- solve(Sig_vector)
x0_vector <- c(0,0)


# Check sampling a scalar
test_that("Can sample scalar", {
  expect_is(samMCMC(normLik,x0_scalar,mu=mu_scalar,Sig=Sig_scalar),"sam")
})

test_that("Can sample vector", {
  expect_is(samMCMC(normLik,x0_vector,mu=mu_vector,invSig=invSig_vector),"sam")
})

# Andy: Is it possible to raise particular errors and ensure that we get those particular errors here?
test_that("Expect error if t0 > numSampBurn", {
  expect_error(samMCMC(normLik,x0,mu=mu_scalar,Sig=Sig_scalar,contro=list(t0=100,numSampBurn=10)))
})

test_that("Expect error if temp is NA [default] and direct is FALSE", {
  expect_error(samMCMC(normLik,x0,mu=mu_scalar,Sig=Sig_scalar,contro=list(direct=FALSE)))
})

test_that("Expect error if temp is given and direct is TRUE [default]", {
  expect_error(samMCMC(normLik,x0,mu=mu_scalar,Sig=Sig_scalar,contro=list(temp=10)))
})

# Check that sampling from a scalar normal distribution yields normal samples
# Set the random number seed (from random.org between 1 and 1,000,000)
set.seed(697111) # from random.org
mu_samp <- -.5
Sig_samp <- 1
sd_samp <- sqrt(Sig_samp)
numSamp <- 1000000
numSampBurn <- 1000
thinning <- 500
out <- samMCMC(normLik,x0,mu=mu_samp,Sig=Sig_samp,control=list(numSamp=numSamp,numSampBurn=numSampBurn,thinning=thinning))
xthin <- out$X_mat[seq(numSampBurn+1,numSamp+numSampBurn,by=thinning)]
tol <- 1e-2

test_that("Ensure normalcy of sample ", {
  expect_true(shapiro.test(xthin)$p.value > 0.05)
})

test_that("Check mean of sample ", {
  xmean <- sd(out$X_mat)
  expect_true((mu_samp - abs(mu_samp)*tol < xmean) && (mu_samp + abs(mu_samp)*tol < -.49))
})

test_that("Check standard deviation of sample ", {
  xsd <- sd(out$X_mat)
  expect_true((sd_samp - abs(sd_samp)*tol < xsd) && (xsd < sd_samp + abs(sd_samp)*tol ))
})