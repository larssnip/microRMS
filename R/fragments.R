#' @name getRMSfragments
#' @title Retrieving RMS fragments from genomes
#'
#' @description Retrieves a set of fragments from a genome, given restriction enzyme cutting motifs.
#'
#' @param genome A \code{\link{Fasta}} object with genome data.
#' @param genome.ID Unique identifier for each genome, will be added to FASTA-headers (text).
#' @param min.length Minimum fragment length (integer).
#' @param max.length Maximum fragment length (integer).
#' @param verbose Turn on/off output text during processing (logical).
#' @param left Text with first, long, restriction enzyme cut motif (text).
#' @param right Text with second, short, restriction enzyme cut motif (text).
#'
#' @details This function is used to find and retrieve all RMS fragments from a genome.
#' A \code{\link{Fasta}}-object with the genome sequence(s) is required. A \code{genome.ID}, if
#' supplied, will be added to the FASTA headers of the output, which is useful if you want to trace
#' the origin of the fragments. All retrieved fragments will then get a FASTA-header starting with the
#' token <genome.ID>_RMSx, where x is an integer (1,2,...,). This first token is followed by a blank.
#' This ensures that all first tokens are unique and that the genome of its origin is indicated.
#'
#' The default restriction enzymes are EcoRI and MseI, with cutting motifs \code{"G|AATTC"} and
#' \code{"T|TAA"}, respectively. The vertical bar indicates where in the motif the enzyme cuts. The
#' forward primers are ligated to the left
#'
#' @return A \code{\link{Fasta}} object with all fragment sequences (5'-3').
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{RMSdbase}}.
#'
#' @importFrom microseq readFasta gff2fasta
#' @importFrom stringr str_c str_remove_all
#' @importFrom dplyr mutate filter
#'
#' @examples
#' # A small genome in this package
#' xpth <- file.path(path.package("microrms"),"extdata")
#' genome.file <- file.path(xpth,"GCF_000009605.1_ASM960v1_genomic.fna")
#' 
#' # Read genome, find fragments
#' gnm <- readFasta(genome.file)
#' frg <- getRMSfragments(gnm, genome.ID = "my.ID")
#'
#' # Write to file with writeFasta(frg, out.file = <filename>)
#' 
#' @export getRMSfragments
#'
getRMSfragments <- function(genome, genome.ID = NULL, min.length = 30, max.length = 500,
                            left = "G|AATTC", right = "T|TAA", verbose = TRUE){
  if(verbose) cat("getRMSfragments: ")
  lft <- str_remove_all(left, "\\|")
  rght <- str_remove_all(right, "\\|")
  gff <- getRMS(genome, lft, rght)
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
      ct.lft <- str_length(str_remove(left, "\\|.+"))
      ct.rght <- str_length(str_remove(right, "\\|.+"))
      fsa %>% 
        mutate(Sequence = str_sub(Sequence, ct.lft+1, -(ct.rght+1))) -> fsa
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
getRMS <- function(genome, left, right){
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
    left.m <- str_locate_all(genome$Sequence[i], pattern = left)[[1]]
    right.m <- str_locate_all(genome$Sequence[i], pattern = right)[[1]]
    if(length(left.m) > 0 & length(right.m) > 0){
      ### positive strand
      lft <- left.m[,1]
      rght <- right.m[,2]
      lft.rght <- t(sapply(lft, function(ll){
        ll2 <- suppressWarnings(min(lft[lft > ll]))
        ss <- suppressWarnings(min(rght[rght > ll]))
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
      lft <- left.m[,2]
      rght <- right.m[,1]
      lft.rght <- t(sapply(lft, function(ll){
        ll2 <- suppressWarnings(max(lft[lft < ll]))
        ss <- suppressWarnings(max(rght[rght < ll]))
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

