##' Extract and write old ptt and rnt
##'
##' \code{ExtractPtt()}:
##'
##'
##' @title
##' @param featureTable
##' @return
##' @examples
##' @author Yulong Niu \email{niuylscu@@gmail.com}
ExtractPtt <- function(featureTable) {

  geneTable <- featureTable[featureTable[, 1] == 'gene' &
                            featureTable[, 2] == 'protein_coding', ]

  ptt <- data.frame(Location = paste(geneTable[, 8], geneTable[, 9], sep = '..'),
                    Strand = geneTable[, 10],
                    Length = geneTable[, 18] / 3 - 1,
                    PID = '-',
                    Gene = geneTable[, 15],
                    Synonym = geneTable[, 17],
                    Code = '-',
                    COG = '-',
                    Product = '-')
  return(ptt)
}


ExtractRnt <- function(featureTable) {

  rnaTable <- featureTable[featureTable[, 1] != 'gene' &
                           featureTable[, 1] != 'CDS', ]

  rnt <- data.frame(Location = paste(geneTable[, 8], geneTable[, 9], sep = '..'),
                    Strand = geneTable[, 10],
                    Length = geneTable[, 18],
                    PID = '-',
                    Gene = geneTable[, 15],
                    Synonym = geneTable[, 17],
                    Code = '-',
                    COG = '-',
                    Product = '-')
  return(rnt)
}
