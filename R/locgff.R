##' Get genomic genes location
##'
##' GetLocsfgff(): get genomic gene location information from gff matrix. The whole genomic gene are divided by chromosomes and plasmids (if the organism has).
##'
##' GetLocsfKEGGSpe(): get genes locus of KEGG species from NCBI gff files. It trys the RefSeq database at first; if RefSeq is not found, then changes to the database to GenBank. The prefix of locus name is the abbreviation of KEGG genomes.
##'
##' ExtractLocs(): extract gene location from gff raw file.
##'
##' download.SpeAnno(): download gff/feature_table and md5sum check files. If the md5sum check fails, download the files again until it passes.
##'
##' @title Get genomic gene locations from the gff file
##' @param gffRawMat raw gff matrix.
##' @param genePrefix prefix to locus gene names.
##' @return GetLocsfgff(): a list of genomes containing gene location information. The locus_tags is used for the gene names.
##' @importFrom stringr str_extract_all str_sub
##' @examples
##' ## read in the dra (Deinococcus radiodurans R1) gff gz file in local disk
##' gzPath <- system.file("extdata", "dra.gff.gz", package = "ProGenome")
##' dragff <- read.gff(gzPath, isurl = FALSE, isgz = TRUE)
##'
##' ## only extract whole locus without genome division
##' locsRawMat <- ExtractLocs(dragff)
##'
##' ## whole locus divided by genomes and plasmids
##' locusList <- GetLocsfgff(dragff, genePrefix = 'dra:')
##'
##' ## get dra genomic locus through FTP URL
##' draLocs <- GetLocsfKEGGSpe('dra')
##'
##' ## two locus names for aac (Alicyclobacillus acidocaldarius subsp. acidocaldarius DSM 446)
##' gzPath <- system.file("extdata", "aac.gff.gz", package = "ProGenome")
##' aacgff <- read.gff(gzPath, isurl = FALSE, isgz = TRUE)
##' locusList <- GetLocsfgff(aacgff, genePrefix = 'aac:')
##'
##' \dontrun{
##' hxaLocs <- GetLocsfKEGGSpe('hxa')
##' csuLocs <- GetLocsfKEGGSpe('csu')}
##' @rdname locsFromgff
##' @seealso read.gff read in gff files
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##'
GetLocsfgff <- function(gffRawMat, genePrefix = character(0)) {

  ## identify using feature "Is_circular" to select genomes and plasmids
  gTagIdx <- which(grepl('Is_circular=', gffRawMat[, 9]))
  gTag <- gffRawMat[gTagIdx, 9]
  pLoc <- grepl('plasmid', gTag)
  if (sum(pLoc) > 0) {
    ## must have at least one genome
    gNames <- character(length(pLoc))
    ## first genome names
    gNames[!pLoc] <- paste0('genome', 1:sum(!pLoc))
    ## second plasmid names
    pNames <- unlist(str_extract_all(gTag, 'plasmid-name=.+?;'))
    pNames <- str_sub(pNames, start = 14, end = -2)
    gNames[pLoc] <- pNames
  } else {
    gNames <- paste0('genome', 1:sum(!pLoc))
  }

  ## split genomes
  ## length(pLoc) == nrow(gMat) is TRUE
  regionLoc <- c(gTagIdx, nrow(gffRawMat))
  startVec <- regionLoc[-length(regionLoc)] + 1
  endVec <- regionLoc[-1] - 1
  gMat <- cbind(startVec, endVec)

  gList <- vector('list', length(gNames))
  for (i in 1:nrow(gMat)) {
    eachLocMat <- ExtractLocs(gffRawMat[gMat[i, 1] : gMat[i, 2], ])
    eachLocMat[, 1] <- paste0(genePrefix, eachLocMat[, 1])
    gList[[i]] <- eachLocMat
  }
  names(gList) <- gNames

  return(gList)
}




##' @inheritParams getGenomicGenes
##' @return GetLocsfKEGGSpe():  a list of genomes containing gene location information. The locus_tags is used for the gene names.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname locsFromgff
##' @export
##'
##'
GetLocsfKEGGSpe <- function(KEGGSpe) {

  ##read in url gff gz files
  speUrl <- AutoSpeFtpUrl(KEGGSpe)
  speFiles <- ListFileFtpUrl(speUrl)
  gffUrl <- speFiles[grepl('gff', speFiles)]
  spegff <- read.gff(gffUrl, isurl = TRUE, isgz = TRUE)

  gList <- GetLocsfgff(spegff, genePrefix = paste0(KEGGSpe, ':'))

  return(gList)

}


##' @param annoStr character strings of gff annotation, which is sperated by ';'.
##' @return a matrix. Locus_tags, old_locus_tags will also be return if provided.
##' @rdname locsFromgff
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##'
##'
GetLocsTag <- function(annoStr) {

  ## split the annotation
  locusSplit <- strsplit(annoStr, split = ';', fixed = TRUE)
  locusTag <- lapply(locusSplit, function(x) {
    ## locus name
    eachLocus <- x[grepl('locus_tag=', x)]
    eachLocusSplit <- strsplit(eachLocus, split = '=', fixed = TRUE)
    eachLocus <- str_trim(sapply(eachLocusSplit, '[[', 2))
    names(eachLocus) <- str_trim(sapply(eachLocusSplit, '[[', 1))
    return(eachLocus)
  })
  locusTag <- do.call(rbind, locusTag)

  return(locusTag)
}




##' @inheritParams GetLocsfgff
##' @return ExtractLocs(): 4 or 5-column matrix, 1st (or 1st and 2ed) is the locus_tags, the last three are start position, end postion, and DNA strand.
##' @rdname locsFromgff
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##'
ExtractLocs <- function(gffRawMat) {

  geneLogic <- gffRawMat[, 3] == 'gene'
  geneMat <- gffRawMat[geneLogic, ]

  locMat <- cbind(GetLocsTag(geneMat[, 9]),
                  geneMat[, 4],
                  geneMat[, 5],
                  geneMat[, 7],
                  deparse.level = 0)

  ## names to the last three column
  last3Idx <- (ncol(locMat) - 2) : ncol(locMat)
  colnames(locMat)[last3Idx] <- c('Start', 'End', 'Strand')
  rownames(locMat) <- NULL

  return(locMat)
}


##' @inheritParams getGenomicGenes
##' @param pattern A \code{character string} whether "gff" or "feature_table"
##' @param saveFolder A folder to save gff and md5sum check files. If the folder does not exist, then creat a one at first.
##' @return download.SpeAnno(): download the gff, feature_table, or md5sum files.
##' @rdname locsFromgff
##' @examples
##' \dontrun{
##' download.SpeAnno('eco', 'gff', 'tmpEco')
##' }
##' @importFrom tools md5sum
##' @importFrom utils download.file
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##'
download.SpeAnno <- function(KEGGSpe, pattern, saveFolder){

  ## check folder
  if (!dir.exists(saveFolder)) {
    dir.create(saveFolder)
  } else {}

  ## FTP urls
  speUrl <- AutoSpeFtpUrl(KEGGSpe)
  speFiles <- ListFileFtpUrl(speUrl)

  ## select and download gff and md5sum check files
  sumPat <- paste(pattern, 'md5checksums', sep = '|')
  fileUrls <- speFiles[grepl(sumPat, speFiles)]
  splitNames <- strsplit(fileUrls, split = '/', fixed = TRUE)
  fileNames <- sapply(splitNames, function(x) {
    eachLast <- x[length(x)]
    return(eachLast)
  })

  ## check md5
  while(TRUE) {
    for(i in 1:length(fileUrls)) {
      download.file(fileUrls[i], file.path(saveFolder, fileNames[i]))
    }

    md5File <- file.path(saveFolder, fileNames[grepl('md5checksums', fileNames)])
    gffFile <- file.path(saveFolder, fileNames[grepl(pattern, fileNames)])

    if (md5sum(gffFile) == ExtractMd5(md5File, pattern)) {
      break
    } else {
      print('Md5sum check failed. Try to download files again.')
    }

  }
}


##' Extract md5sum for a given patterned files
##'
##' md5sum for the "gff" and "feature_table" files
##'
##' @title Extract md5sum
##' @param mdfile The md5sum file path.
##' @inheritParams download.SpeAnno
##' @return A \code{Character} md5sum string.
##' @importFrom utils read.table
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
ExtractMd5 <- function(mdfilePath, pattern) {

  mdMat <- read.table(mdfilePath, stringsAsFactors = FALSE)
  mdStr <- mdMat[grepl(pattern, mdMat[, 2]), 1]

  return(mdStr)
}

