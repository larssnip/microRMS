

## Non-exported function to gracefully fail when external dependencies are missing.
available.external <- function(what){
  if(what=="vsearch"){
    chr <- NULL
    try(chr <- system('vearch -help', intern=TRUE), silent=TRUE)
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
  # } else if( what=="hmmer" ){
  #   chr <- NULL
  #   try( chr <- system('hmmscan -h', intern=TRUE), silent=TRUE )
  #   if( is.null( chr ) ){
  #     stop( paste('hmmer was not found by R.',
  #                 'Please install hmmer from: http://hmmer.org/download.html',
  #                 'After installation, re-start R and make sure the hmmer softwares can be run from R by',
  #                 'the command \'system("hmmscan -h")\'.', sep = '\n'))
  #     return(FALSE)
  #   } else {
  #     return(TRUE)
  #   }
  # } else if( what=="blast+" ){
  #   chr <- NULL
  #   try( chr <- system('makeblastdb -help', intern=TRUE), silent=TRUE )
  #   if( is.null( chr ) ){
  #     stop( paste('blast+ was not found by R.',
  #                 'Please install blast+ from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/',
  #                 'After installation, re-start R and make sure the blast+ softwares can be run from R by',
  #                 'the command \'system("makeblastdb -help")\'.', sep = '\n'))
  #     return(FALSE)
  #   } else {
  #     return(TRUE)
  #   }
  # }
}
