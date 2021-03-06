#' @title Adaptive Metropolis sampling
#'
#' @description This function samples using the Metropolis algorithm described by
#' Harrio, Saksman, and Tamminen (2001). \code{samMCMC} works in one of two ways:
#'
#' (1) [default] The function to be sampled is directly specified
#' (2)           The negative log density of the function to be sampled is specified,
#'               along with a temperature
#'
#' @details \code{sampFunc} must be a function taking
#'  the variable to be sampled as its first input. Additional inputs to sampmFunc can be
#'  given as named variable inputs.
#'
#' The \code{control} argument is a list that can supply any of the following (otherwise 
#' defaults, in brackets, are used):
#' \itemize{
#'   \item{\code{direct}} Is the sampling function specified directly or as a negative log with associated temperature  ["direct"]
#'   \item{\code{temp}} The temperature for sampling [NA]
#'   \item{\code{numSamp}} Number of samples after initialization [1000]
#'   \item{\code{t_0}} Number of samples using initial proporsal ditribution [100]
#'   \item{\code{C_0}} The covariance to use for sampling for t <= t_0 [diag(c(rep(1e-6,length(X_0))))]
#'   \item{\code{s_d}} The rescaling of the sampling covariance relative to the data covariance [(2.4)^2 / length(X_0)]
#'   \item{\code{epsilon}} The multiplier of the identity matrix that is added to the data covariance to ensure numerical stability [1e-12]
#'   \item{\code{verbose}} Whether to print out information as the run proceeds [F]
#'   \item{\code{numSampBurn}} Number of burn-in samples prior [1000]
#'   \item{\code{thinning}} Thinning of samples to apply [1]
#'   \item{\code{sampsToAdd}} Number of samples to add (for continuing runs) [Either given or set equal to numSamp for a new run]
#' }
#'
#' @param sampFunc The sampling function (a negative log-likelihood in case 2)
#' @param init The starting point for sampling, X_0, or the output from a previous call 
#' to samMCMC
#' @param ... Further arguments to be passed to sampFunc
#' @param control A list of control parameters. See Details.
#'
#' @return An object of class \code{chain} that is a list containing the samples along 
#' with summary information
#'
#' @author Michael Holton Price <MichaelHoltonPrice@@gmail.com>
#' 
#' @references Harrio, Saksman, and Tamminen (2001), An adaptive Metropolis algorithm, 
#' in Bernoulli, Volume 7, Number 2, pages 223 through 242.
#' 
#' @export

samMCMC <- function(sampFunc,init,...,control=list()) {
  # First, handle inputs. Each varible can be specified in one of three ways:
  # the default, the value from a prevous chain, and the value in control.
  # Preference is given to the value in control, followed by the value from a
  # previous chain, followed by the default. That is, the value from the
  # previous chain overwrites the default and the value from control overwrites
  # the value from a previous chain.

  # init is either a vector (X_0) or the result of a previous call to samMCMC
  haveChain <- 'sam' %in% class(init)

  # Check that all variables in the input control are recognized
  # See:
  #     test_that("Expect error if temp is NA [default] and direct is FALSE")
  #     Context: samMCMC
  if(length(control) > 0) {
    allControlVar <- c("direct",
                       "temp",
                       "numSamp",
                       "verbose",
                       "C_0",
                       "t_0",
                       "s_d",
                       "epsilon",
                       "numSampBurn",
                       "thinning",
                       "sampsToAdd"
                      )
    for(controlVar in names(control)) {
      if( !( controlVar %in% allControlVar ) ) {
        stop(paste0(controlVar," is not a recognized control variable"))
      }
    }
  }

  # Save the input control as inputControl
  inputControl <- control
  # Create a new control object
  control <- list()
  # Set prevControl (if applicable; if not, set to NA)
  if(haveChain) {
    prevControl <- init$control
  } else {
    prevControl <- NA
  }

 
  # Set X_0
  if(!haveChain) {
    X_0 <- init
  } else {
    X_0 <- as.vector(init$X_mat[,ncol(init$X_mat)])
  }

  # Set control values by calling the "helper" function chooseControlValue
  # [error handling is done below]
  control$direct <- chooseControlValue('direct',T,inputControl,haveChain,prevControl)
  control$temp <- chooseControlValue('temp',NA,inputControl,haveChain,prevControl)
  control$numSamp <- chooseControlValue('numSamp',1000,inputControl,haveChain,prevControl)
  control$verbose <- chooseControlValue('verbose',F,inputControl,haveChain,prevControl)

  # Must handle the special case where X_0 is a scalar
  if(length(X_0) == 1) {
    control$C_0 <- chooseControlValue('C_0',1e-6,inputControl,haveChain,prevControl)
  } else { # X_0 is a vector
    control$C_0 <- chooseControlValue('C_0',diag(c(rep(1e-6,length(X_0)))),inputControl,haveChain,prevControl)
  }

  control$t_0 <- chooseControlValue('t_0',100,inputControl,haveChain,prevControl)
  control$s_d <- chooseControlValue('s_d',(2.4)^2 / length(X_0),inputControl,haveChain,prevControl)
  control$epsilon <- chooseControlValue('epsilon',1e-12,inputControl,haveChain,prevControl)
  control$numSampBurn <- chooseControlValue('numSampBurn',1000,inputControl,haveChain,prevControl)
  control$thinning <- chooseControlValue('thinning',1,inputControl,haveChain,prevControl)
  
  # Determine the number of samples to make
  # For new chains this is either explicitly given or set to
  # numSamp + numSampBurn. However, if it is not expliclty given it is assumed
  # that future chains should sample numSamp additional times.
  #
  # If this is a continued chain, ues the value stored in the chain's control
  # variable unless it is overridden by user input. If the latter is true,
  # replace the value of sampsToAdd in the chain.
  if(!haveChain) {
    if('sampsToAdd' %in% names(inputControl)) {
      sampsToAdd <- inputControl$sampsToAdd
      control$sampsToAdd <- inputControl$sampsToAdd
    } else { # sampsToAdd not given in control
      sampsToAdd <- control$numSamp + control$numSampBurn
      control$sampsToAdd <- control$numSamp # Future calls with this chain should not do the burn in
    }
  } else { # have chain
    if('sampsToAdd' %in% names(inputControl)) {
      sampsToAdd <- inputControl$sampsToAdd
      control$sampsToAdd <- inputControl$sampsToAdd
    } else { # sampsToAdd not in input control
      sampsToAdd <- prevControl$sampsToAdd
      control$sampsToAdd <- prevControl$sampsToAdd
    }
  }

  # Handle errors
  #
  # Direct must be TRUE if temp is NA
  # See:
  #     test_that("Expect error if temp is NA [default] and direct is FALSE")
  #     Context: samMCMC
  if(is.na(control$temp) && !control$direct) {
    stop('temp is NA but direct is FALSE')
  }

  # temp should be NA if direct is TRUE
  # See:
  #     test_that("Expect error if temp is given and direct is TRUE [default]")
  #     Context: samMCMC
  if(!is.na(control$temp) && control$direct) {
    stop('temp is given but direct is TRUE')
  }

  I_d <- diag(length(X_0))
  sf <- function(X) sampFunc(X,...)

  # Initialize variables
  X_t <- X_0
  sampFunc_t <- sf(X_t)

  if(!haveChain) {
    covObj <- updateCov(X_t)
    ttOffset <- 0
  } else {
    covObj <- init$covObj
    ttOffset <- ncol(init$X_mat)
  }

  X_mat <- matrix(NA,length(X_0),sampsToAdd)
  sampFuncVect <- vector()
  acceptVect <- vector()
  for(ii in 1:sampsToAdd) {
    tt <- ttOffset + ii ## tt because t is transpose
    if(tt <= control$t_0) {    
      if(control$verbose) {
        print('-- C_0 --')
      }

      if(length(X_t)==1) {
        X_tp1 <- X_t +  rnorm(1,sd=sqrt(control$C_0))
      } else {
        X_tp1 <- mvrnormR(1,X_t,control$C_0)
      }
    } else { # tt > t_0
      if(control$verbose) {
        print('-- C_t --')
      }

      C_t <- control$s_d * (covObj$cov + control$epsilon*I_d)
      if(length(X_t)==1) {
        X_tp1 <- X_t +  rnorm(1,sd=sqrt(C_t))
      } else {
        X_tp1 <- mvrnormR(1,X_t,C_t)
      }
    }
    if(control$verbose) {
      print(tt)
      if(!control$direct) {
        print(control$temp)
      }
    }
 
    # Calculate the acceptance ratio, alpha, accounting for direct
    # By construction, the proposal distribution is symmetric.
    sampFunc_tp1 <- sf(X_tp1)
    if(!is.finite(sampFunc_tp1)) {
      accept <- F
    } else {
      if(control$direct) {
        alpha <- min(1,sampFunc_tp1/sampFunc_t)
      } else {
        alpha <- min(1,exp(-(sampFunc_tp1-sampFunc_t)/control$temp))
      }
      accept <- runif(1) < alpha
    }

    acceptVect[ii] <- accept
    if(!accept) {
      X_tp1 <- X_t
      sampFunc_tp1 <- sampFunc_t
    }
    X_mat[,ii] <- X_tp1
    sampFuncVect[ii] <- sampFunc_tp1

    # Get ready for next sample
    X_t <- X_tp1
    sampFunc_t <- sampFunc_tp1
    covObj <- updateCov(X_t,covObj)
    if(control$verbose) {
      print(sampFunc_t)
      print(accept)
    }
  } # end main loop
  if(haveChain) {
    X_mat <- cbind(init$X_mat,X_mat)
    sampFuncVect <- c(init$sampFuncVect,sampFuncVect)
    acceptVect <- c(init$acceptVect,acceptVect)
  }
  returnList <- list(X_mat=X_mat,sampFuncVect=sampFuncVect,acceptVect=acceptVect,covObj=covObj,control=control,temp=control$temp)
  if(!haveChain) {
    returnList$firstX <- X_0
  } else {
    returnList$firstX <- init$firstX
  }
  
  class(returnList) <- c('sam', 'mcmc')
  return(returnList)
}

# See https://gallery.rcpp.org/articles/simulate-multivariate-normal/
# Especially for large dimensions, the following code is much faster than
# MASS::mvrnorm
mvrnormR <- function(n, mu, sigma) {
  ncols <- ncol(sigma)
  mu <- rep(mu, each = n) ## not obliged to use a matrix (recycling)
  return(as.vector(mu + matrix(rnorm(n * ncols), ncol = ncols) %*% chol(sigma)))
}
