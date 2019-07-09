#' @name getRMSfragments
#' @title Retrieving RMS fragments from genomes
#'
#' @description Retrieves a set of fragments from a genome, given restriction enzyme cutting motifs.
#'
#' @param genome A \code{\link{Fasta}} object with genome data.
#' @param genome.ID Unique text identifier for each genome, will be added to FASTA-headers.
#' @param min.length Minimum amplicon length (bases).
#' @param max.length Maximum amplicon length (bases).
#' @param verbose Logical to turn on/off output text during processing.
#' @param left Text with first, long, restriction enzyme cut motif. Deafult is the EcoRI.
#' @param right Text with second, short, restriction enzyme cut motif. Deafult is the MseI.
#' @param trim Logical indicating if the restriction motifs above should be trimmed off the ends of the fragments. Default is \code{TRUE}
#'
#' @details This function is used to find and retrieve all RMS fragments from a genome. A \code{\link{Fasta}}-object
#' with the genome sequence(s) is required. A \code{genome.ID} may/should be supplied, and will be added to the header-lines of the output
#' \code{\link{Fasta}} object to identify the origin of the fragments. All retrieved fragments will then get
#' a FASTA-header starting with the token <genome.ID>_RMSx, where x is an integer (1,2,...,). This first token is followed
#' by a blank. This ensures that all first tokens are unique and that the genome of its origin is indicated.
#'
#' The default restriction enzymes are EcoRI and MseI,
#' with cutting motifs \code{"GAATTC"} and \code{"TTAA"}, respectively. Change cutting motifs accordingly if
#' you use other restriction enzymes.
#'
#' @return A \code{\link{Fasta}} object with all fragment sequences (5'-3').
#'
#' @author Lars Snipen.
#'
#' @seealso more here.
#'
#' @importFrom microseq readFasta
#' @importFrom micropan gff2fasta
#' @importFrom stringr str_c
#' @importFrom dplyr mutate filter
#'
#' @examples more here.
#'
#' @export getRMSfragments
#'
getRMSfragments <- function(genome, genome.ID = NULL, min.length = 30, max.length = 500,
                            left = "GAATTC", right = "TTAA", verbose = TRUE, trim = TRUE){
  if(verbose) cat("getRMSfragments: ")
  gff <- getRMS(genome, left, right, trim)
  if(nrow(gff)>0){
    gff %>%
      mutate(Length = abs(Start - End) + 1) %>%
      filter((Length >= min.length) & (Length <= max.length)) -> gff
    if(verbose) cat("found", nrow(gff), "RMS-fragments\n")
    if(nrow(gff) > 0){
      fsa <- gff2fasta(gff, genome)
      if(!is.null(genome.ID)){
        fsa$Header <- str_c(str_c(genome.ID, str_c("RMS", 1:nrow(fsa)), sep = "_"), fsa$Header, sep = " ")
      }
    } else {
      if(verbose) cat("found no RMS-fragments within min and max length!\n")
    }
  } else {
    if(verbose) cat("found no RMS-fragments!\n")
    fsa <- data.frame(Header = NULL, Sequence = NULL, stringsAsFactors = F)
  }
  class(fsa) <- c("Fasta", "data.frame")
  return(fsa)
}



### Local function
getRMS <- function(genome, long = "GAATTC", short = "TTAA", trim = TRUE){
  require(stringr)
  gff <- data.frame(Seqid      = NULL,
                    Source     = NULL,
                    Type       = NULL,
                    Start      = NULL,
                    End        = NULL,
                    Score      = NULL,
                    Strand     = NULL,
                    Phase      = NULL,
                    Attributes = NULL,
                    stringsAsFactors = F)
  ### looping over genome-sequences
  for(i in 1:nrow(genome)){
    long.m <- str_locate_all(genome$Sequence[i], pattern = long)[[1]]
    short.m <- str_locate_all(genome$Sequence[i], pattern = short)[[1]]
    if(length(long.m) > 0 & length(short.m) > 0){
      ### positive strand
      if(trim){
        lng <- long.m[,2] + 1
        shrt <- short.m[,1] - 1
      } else {
        lng <- long.m[,1]
        shrt <- short.m[,2]
      }
      lft.rght <- t(sapply(lng, function(ll){
        ll2 <- suppressWarnings(min(lng[lng > ll]))
        ss <- suppressWarnings(min(shrt[shrt > ll]))
        return(c(ll, ifelse(ss < ll2, ss, Inf)))
      }))
      nf <- nrow(lft.rght)
      if(nf > 0){
        gff <- rbind(gff, data.frame(Seqid      = rep(word(genome$Header[i], 1, 1), nf),
                                     Source     = NA,
                                     Type       = rep("RMS_fragment", nf),
                                     Start      = lft.rght[,1],
                                     End        = lft.rght[,2],
                                     Score      = NA,
                                     Strand     = rep("+", nf),
                                     Phase      = NA,
                                     Attributes = NA,
                                     stringsAsFactors = F))
      }
      ### negative strand
      if(trim){
        lng <- long.m[,1] - 1
        shrt <- short.m[,2] + 1
      } else {
        lng <- long.m[,2]
        shrt <- short.m[,1]
      }
      lft.rght <- t(sapply(lng, function(ll){
        ll2 <- suppressWarnings(max(lng[lng < ll]))
        ss <- suppressWarnings(max(shrt[shrt < ll]))
        return(c(ifelse(ss > ll2, ss, -Inf), ll))
      }))
      nn <- nrow(lft.rght)
      if(nn > 0){
        gff <- rbind(gff, data.frame(Seqid      = rep(word(genome$Header[i], 1, 1), nn),
                                     Source     = NA,
                                     Type       = rep("RMS_fragment", nf),
                                     Start      = lft.rght[,1],
                                     End        = lft.rght[,2],
                                     Score      = NA,
                                     Strand     = rep("-", nf),
                                     Phase      = NA,
                                     Attributes = NA,
                                     stringsAsFactors = F))
      } # end if
    } # end if
  }  # end for
  return(gff)
}

