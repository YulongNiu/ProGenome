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
##' rnt <- ExtractRnt(ft)
##'
##' \dontrun{
##' write.ptt(ptt, 'eco.ptt')
##' write.rnt(rnt, 'eco.rnt')
##' }
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

  rnt <- data.frame(Location = paste(rnaTable[, 8], rnaTable[, 9], sep = '..'),
                    Strand = rnaTable[, 10],
                    Length = rnaTable[, 18],
                    PID = '-',
                    Gene = rnaTable[, 15],
                    Synonym = rnaTable[, 17],
                    Code = '-',
                    COG = '-',
                    Product = '-')
  return(rnt)
}


##' @param ptt A \code{matrix} indicating the ptt matrix.
##' @param file A \code{character vector} indicating the file path.
##' @importFrom utils write.table
##' @rdname pttrnt
##' @export
##'
write.ptt <- function(ptt, file) {

  f <- file(file, "w")

  writeLines('ptt header for certain species.', f)
  writeLines(paste(nrow(ptt), 'proteins'), f)
  write.table(ptt, sep = '\t', row.names = FALSE, quote = FALSE, f)

  close(f)
}

##' @param rnt A \code{matrix} indicating the rnt matrix
##' @inheritParams write.ptt
##' @importFrom utils write.table
##' @rdname pttrnt
##' @export
##'
write.rnt <- function(rnt, file) {
  f <- file(file, "w")

  writeLines('rnt header for certain species.', f)
  writeLines(paste(nrow(rnt), 'RNAs'), f)
  write.table(rnt, sep = '\t', row.names = FALSE, quote = FALSE, f)

  close(f)
}
