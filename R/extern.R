

## Non-exported function to gracefully fail when external dependencies are missing.
available.external <- function(what){
  if(what == "vsearch"){
    chr <- NULL
    try(chr <- system('vsearch -help', intern = TRUE), silent = TRUE)
    if(is.null(chr)){
      stop( paste('vsearch was not found by R.',
                  'Please install vsearch from: https://github.com/torognes/vsearch',
                  'After installation, re-start R and make sure vsearch can be run from R by',
                  'the command \'system("vsearch -help")\'.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}
