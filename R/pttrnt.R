##' Extract and write old ptt and rnt
##'
##' \code{ExtractPtt()}: Extract ptt from the featureTable.
##'
##' \code{ExtractRnt()}: Extract rnt from the featureTable.
##'
##' \code{write.ptt()}: Write the ptt file.
##'
##' \code{write.rnt()}: Write the rnt file.
##' @title ptt and rnt
##' @param featureTable A \code{matrix} represents the feature table retrieved from the NCBI ftp.
##' @return
##' \code{ExtractPtt()}: Extract ptt from the featureTable.
##'
##' \code{ExtractRnt()}: Extract rnt from the featureTable.
##' @examples
##' gzPath <- system.file('extdata', 'eco.feature_table.txt.gz', package = 'ProGenome')
##' ft <- read.gff(gzPath)
##'
##' ptt <- ExtractPtt(ft)
##' rnt <- ExtarctRnt(rnt)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname pttrnt
##' @export
##'
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


##' @inheritParams ExtractPtt
##' @rdname pttrnt
##' @export
##'
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
