#' @title Choose control variable value
#'
#' @description This "helper" function chooses the value for a control variable that can either be (in order of increasing preference) the default, the value fro a previous run, and a directly input value.
#'
#' @param varName The variable name
#' @param defaultVal The default value
#' @param inputControl The user specified control list
#' @param haveChain A boolean indicating whether the values from a previous run are available
#' @param prevControl The control list from a previous run
#'
#' @return The preferred value for the variable
#'
#' @author Michael Holton Price <MichaelHoltonPrice@@gmail.com>
#' 
#' @export
chooseControlValue <- function(varName,defaultVal,inputControl,haveChain,prevControl=NA) {
  # Each input is set in one of three ways with the following preference ordering:
  #
  # default < prevControl < input control

  outputVal <- defaultVal # (a) default

  if(haveChain) {
    outputVal <- prevControl[[varName]] # (b) from the input chain
  }

  if(varName %in% names(inputControl)) {
    outputVal <- inputControl[[varName]] # (c) from the input control
  }

  return(outputVal)
}
