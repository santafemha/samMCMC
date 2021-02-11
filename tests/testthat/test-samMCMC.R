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

# Check functioning of sampsToAdd when explicitly given
out <- samMCMC(normLik,x0_scalar,mu=mu_scalar,Sig=Sig_scalar,control=list(sampsToAdd=93))
test_that("sampsToAdd samples are made [new chain; sampsToAdd given]", {
  expect_equal(length(out$X_mat),93)
})

test_that("Input value is used for output control [new chain; sampsToAdd given]", {
  expect_equal(out$control$sampsToAdd,93)
})

out <- samMCMC(normLik,out,mu=mu_scalar,Sig=Sig_scalar)
test_that("sampsToAdd new samples are made [extended chain; sampsToAdd given]", {
  expect_equal(length(out$X_mat),2*93)
})

test_that("Original input value is used for output control [new chain; sampsToAdd given]", {
  expect_equal(out$control$sampsToAdd,93)
})

out <- samMCMC(normLik,out,mu=mu_scalar,Sig=Sig_scalar,control=list(sampsToAdd=44))
test_that("Test that changing sampsToAdd is correctly done for an existing chain", {
  expect_equal(length(out$X_mat),2*93+44)
  expect_equal(out$control$sampsToAdd,44)
})


# Check functioning of sampsToAdd when not explicitly given
out <- samMCMC(normLik,x0_scalar,mu=mu_scalar,Sig=Sig_scalar,control=list(numSamp=205,numSampBurn=108))
test_that("sampsToAdd samples are made [new chain; sampsToAdd not given]", {
  expect_equal(length(out$X_mat),205+108)
})

test_that("numSamp, not numSamp + numSampBurn, is used for output control [new chain; sampsToAdd not given]", {
  expect_equal(out$control$sampsToAdd,205)
})

out <- samMCMC(normLik,x0_vector,mu=mu_vector,invSig=invSig_vector,control=list(numSamp=404,numSampBurn=77))
test_that("Can sample vector", {
  expect_equal(nrow(out$X_mat),2)
  expect_equal(ncol(out$X_mat),404+77)
  expect_is(out,"sam")
})

# test that the *correct* errors are given for the *correct* reasons
test_that("t0 is not a recognized control variable", { # the correct control variable is t_0
  expect_error(samMCMC(normLik,x0_scalar,mu=mu_scalar,Sig=Sig_scalar,control=list(t0=100)))
  badControl <- try(samMCMC(normLik,x0_scalar,mu=mu_scalar,Sig=Sig_scalar,control=list(direct=FALSE)), 
                       silent = TRUE)
  expect_true(grepl('temp is NA but direct is FALSE', badControl))
})


test_that("Expect error if temp is NA [default] and direct is FALSE", {
  expect_error(samMCMC(normLik,x0_scalar,mu=mu_scalar,Sig=Sig_scalar,control=list(direct=FALSE)))
  tempNAdirectF <- try(samMCMC(normLik,x0_scalar,mu=mu_scalar,Sig=Sig_scalar,control=list(direct=FALSE)), 
                       silent = TRUE)
  expect_true(grepl('temp is NA but direct is FALSE', tempNAdirectF))
})

test_that("Expect error if temp is given and direct is TRUE [default]", {
  expect_error(samMCMC(normLik,x0_scalar,mu=mu_scalar,Sig=Sig_scalar,control=list(temp=10)))
  tempXdirectT <- try(samMCMC(normLik,x0_scalar,mu=mu_scalar,Sig=Sig_scalar,control=list(temp=10)), 
                      silent = TRUE)
  expect_true(grepl('temp is given but direct is TRUE', tempXdirectT))
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
out <- samMCMC(normLik,x0_scalar,mu=mu_samp,Sig=Sig_samp,control=list(numSamp=numSamp,numSampBurn=numSampBurn,thinning=thinning))
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
