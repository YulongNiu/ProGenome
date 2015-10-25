##' Get genomic genes location of KEGG species from NCBI gff files
##'
##' getLocsfgff(): get genomic gene location information from gff file. It trys the RefSeq database at first; if RefSeq is not found, then changes to the database to GenBank.
##'
##' download.Spegff(): download gff and md5sum check files. If the md5sum check fails, download the files again until it passes.
##' @title Get genomic gene locations from the gff file
##' @inheritParams getGenomicGenes
##' @return GetLocsfgff(): a list of genomes containing gene location information. The locus_tags is used for the gene names.
##' @examples
##' draLocs <- GetLocsfgff('dra')
##' \dontrun{
##' hxaLocs <- GetLocsfgff('hxa')
##' csuLocs <- GetLocsfgff('csu')}
##' @rdname locsFromgff
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##' 
GetLocsfgff <- function(KEGGSpe) {


  getLocsTag <- function(annoStr) {
    ## USE: extract locus_tags from a gff annotation information.
    ## INPUT: character strings of gff annotation, which is sperated by ';'.
    ## OUTPU: locus_tags
    ## EXAMPLE:
    ## testStr <- c('ID=gene1;Dbxref=GeneID:10795663;Name=HALXA_RS00010;gbkey=Gene;gene_biotype=protein_coding;locus_tag=HALXA_RS00010;old_locus_tag=Halxa_0687', 'ID=gene0;Dbxref=GeneID:10799347;Name=HALXA_RS00005;gbkey=Gene;gene_biotype=protein_coding;locus_tag=HALXA_RS00005;old_locus_tag=Halxa_0686')
    ## testLocs <- getLocsTag(testStr)

    locusSplit <- strsplit(annoStr, split = ';')
    locusTag <- sapply(locusSplit, function(x) {
      eachLocus <- x[grepl('^locus_tag=', x)]
      eachLocus <- sapply(strsplit(eachLocus, split = '=', fixed = TRUE), '[[', 2)

      return(eachLocus)
    })
    
    return(locusTag)
  }

 ExtractLocs <- function(gffMat) {
    ## USE: extract gene location from gff raw file.
    ## INPUT: 'gffMat' is the gff matrix.
    ## OUTPUT: 4-column matrix, 1st is the locus_tags, 2ed is the from locus, 3rd is the to locus, 4th is the strand.

    geneLogic <- gffMat[, 3] == 'gene'
    geneMat <- gffMat[geneLogic, ]

    locMat <- cbind(getLocsTag(geneMat[, 9]),
                    geneMat[, 4],
                    geneMat[, 5],
                    geneMat[, 7],
                    deparse.level = 0)
    
    colnames(locMat) <- c('GeneName', 'From', 'To', 'Strand')
    rownames(locMat) <- NULL

    return(locMat)
  }

  
  ##~~~~~~~~~~~~~~~~~~~~read in gff files~~~~~~~~~~~~~~~~~~
  speUrl <- AutoSpeFtpUrl(KEGGSpe)
  speFiles <- ListFileFtpUrl(speUrl)
  gffUrl <- speFiles[grepl('gff', speFiles)]
  spegff <- read.gff(gffUrl, isurl = TRUE, isgz = TRUE)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~get gene locs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## identify using feature "Is_circular"
  gTagIdx <- which(grepl('Is_circular=', spegff[, 9]))
  gTag <- spegff[gTagIdx, 9]
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
  regionLoc <- c(gTagIdx, nrow(spegff))
  startVec <- regionLoc[-length(regionLoc)] + 1
  endVec <- regionLoc[-1] - 1
  gMat <- cbind(startVec, endVec)

  gList <- vector('list', length(gNames))
  for (i in 1:nrow(gMat)) {
    eachLocMat <- ExtractLocs(spegff[gMat[i, 1] : gMat[i, 2], ])
    eachLocMat[, 1] <- paste(KEGGSpe, eachLocMat[, 1], sep = ':')
    gList[[i]] <- eachLocMat
  }
  names(gList) <- gNames
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  return(gList)
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
  

