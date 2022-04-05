

## Non-exported function to gracefully fail when external dependencies are missing.
available.external <- function(what){
  chr <- NULL
  try(chr <- system(str_c(what, '--help'), intern = TRUE), silent = TRUE)
  if(is.null(chr)){
    stop(paste("vsearch was not found,please executable command:", what))
    return(FALSE)
  } else {
    return(TRUE)
  }
}
