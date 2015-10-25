##' Get genomic genes location
##'
##' GetLocsfgff(): get genomic gene location information from gff matrix. The whole genomic gene are divided by chromosomes and plasmids (if the organism has).
##'
##' GetLocsfKEGGSpe(): get genes locus of KEGG species from NCBI gff files. It trys the RefSeq database at first; if RefSeq is not found, then changes to the database to GenBank. The prefix of locus name is the abbreviation of KEGG genomes.
##'
##' ExtractLocs(): extract gene location from gff raw file.
##' 
##' download.Spegff(): download gff and md5sum check files. If the md5sum check fails, download the files again until it passes.
##' 
##' @title Get genomic gene locations from the gff file
##' @param gffRawMat raw gff matrix.
##' @param genePrefix prefix to locus gene names.
##' @return GetLocsfgff(): a list of genomes containing gene location information. The locus_tags is used for the gene names.
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

  ## identify using feature "Is_circular"
  gTagIdx <- which(grepl('Is_circular=', gffRawMat[, 9]))
  gTag <- gffRawMat[gTagIdx, 9]
  pLoc <- grepl('plasmid', gTag)
  if (sum(pLoc) > 0) {
    ## must have one genome
    gNames <- c(paste0('genome', 1:sum(!pLoc)),
                paste0('plasmid', 1:sum(pLoc)))
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
##' @return locus_tags
##' @rdname locsFromgff
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##'
##' 
GetLocsTag <- function(annoStr) {

  locusSplit <- strsplit(annoStr, split = ';')
  locusTag <- sapply(locusSplit, function(x) {
    eachLocus <- x[grepl('^locus_tag=', x)]
    eachLocus <- sapply(strsplit(eachLocus, split = '=', fixed = TRUE), '[[', 2)

    return(eachLocus)
  })
  
  return(locusTag)
}




##' @inheritParams GetLocsfgff
##' @return ExtractLocs(): 4-column matrix, 1st is the locus_tags, 2ed is the from locus, 3rd is the to locus, 4th is the strand.
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
  
  colnames(locMat) <- c('GeneName', 'From', 'To', 'Strand')
  rownames(locMat) <- NULL

  return(locMat)
}


##' @inheritParams getGenomicGenes
##' @param saveFolder A folder to save gff and md5sum check files. If the folder does not exist, then creat a one at first.
##' @return download.Spegff(): download the gff and md5sum files.
##' @rdname locsFromgff
##' @examples
##' \dontrun{
##' download.Spegff('eco', 'tmpEco')
##' }
##' @importFrom tools md5sum
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##' 
download.Spegff <- function(KEGGSpe, saveFolder){

  ExtractGffMd5 <- function(mdfile) {
    ## USE: extract md5sum string for gff file
    ## INPUT: 'mdfile' is the md5sum file path
    ## OUTPUT: md5sum string

    mdMat <- read.table(mdfile, stringsAsFactors = FALSE)
    mdStr <- mdMat[grepl('gff|md5checksums', mdMat[, 2]), 1]

    return(mdStr)
  }

  ## check folder
  if (!dir.exists(saveFolder)) {
    dir.create(saveFolder)
  } else {}
  
  ## FTP urls
  speUrl <- AutoSpeFtpUrl(KEGGSpe)
  speFiles <- ListFileFtpUrl(speUrl)

  ## select and download gff and md5sum check files
  fileUrls <- speFiles[grepl('gff|md5checksums', speFiles)]
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
    gffFile <- file.path(saveFolder, fileNames[grepl('gff', fileNames)])

    if (md5sum(gffFile) == ExtractGffMd5(md5File)) {
      break
    } else {
      print('Md5sum check failed. Try to download files again.')
    }
    
  }
}
  

