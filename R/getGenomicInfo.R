##' Get genomic information NCBI FTP url
##'
##' getFtpUrl(): get the FTP url storing genome information from GenBank or RefSeq database. This function is mainly used to get the species genome assembly information.
##' listFtpFileUrl(): list the file download FTP urls.
##' read.gffUrl(): read in raw gff file from a FTP url.
##' getLocsfgff(): get genomic gene location information from gff file. It trys the RefSeq database at first; if RefSeq is not found, then changes to the database to GenBank.
##' @title Retrieve genome FTP URL
##' @inheritParams getGenomicGenes
##' @param database "GenBank", "RefSeq", or c("GenBank", "RefSeq").
##' @return FTP url.
##' @examples
##' ## E. coli genome GenBank FTP
##' ecoliUrl <- getSpeFtpUrl('eco', c('GenBank', 'RefSeq'))
##' ## list E. coli FTP files
##' ecoliFiles <- listFileFtpUrl(ecoliUrl[1])
##' ## read in gff file
##' gffUrl <- ecoliFiles[grepl('gff', ecoliFiles)]
##' ecoligff <- read.gffUrl(gffUrl)
##' \dontrun{
##' draliUrl <- getSpeFtpUrl('dra')
##' draliFiles <- listFileFtpUrl(draliUrl)
##' }
##' @importFrom KEGGAPI getKEGGSpeInfo
##' @importFrom stringr str_extract str_detect
##' @importFrom RCurl getURL
##' @importFrom xml2 read_xml xml_attr xml_text
##' @importFrom foreach foreach %do%
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname genomeFTP
##' @export
##'
##' 
getSpeFtpUrl <- function(KEGGSpe, database = 'GenBank') {

  ##~~~~~~~~~~~~~~assembly ID~~~~~~~~~~~~~~~~~~~~~
  speInfo <- getKEGGSpeInfo(KEGGSpe)
  sourceEle <- speInfo$`Data source`

  ## Assembly: GCF_000005845.2
  assNum <- str_extract(sourceEle, 'Assembly: \\w+[\\.\\w]*')
  assNum <- sapply(strsplit(assNum, split = ': ', fixed = TRUE), '[[', 2)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~~~~~~~FTP URL~~~~~~~~~~~~~~~~~~~
  ## collapse database
  database <- paste(database, collapse = '|')
  
  ## NCBI assembly url
  assBaseUrl <- 'http://www.ncbi.nlm.nih.gov/assembly/'
  assUrl <- paste0(assBaseUrl, assNum)

  webStr <- getURL(assUrl)
  webVec <- unlist(strsplit(webStr, split = '\n', fixed = TRUE))
  matchPattern <- paste0('(', database, ')', ' FTP site')
  ftpLine <- webVec[str_detect(webVec, matchPattern)]

  ## find url
  ftpUrl <- foreach(i = 1:length(ftpLine), .combine = c) %do% {
    eachFtpXml <- read_xml(ftpLine[i])
    eachFtpUrl <- xml_attr(eachFtpXml, 'href')
    eachName <- sapply(strsplit(xml_text(eachFtpXml), ' ', fixed = TRUE),
                       '[[',
                       1)
    names(eachFtpUrl) <- eachName
    return(eachFtpUrl)
  }
  
  ## ftp url often point to a folder
  ftpUrlNames <- names(ftpUrl)
  ftpUrl <- paste0(ftpUrl, '/')
  names(ftpUrl) <- ftpUrlNames
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  return(ftpUrl)
}




##' @param ftpUrl FTP url that can be the output of function  getFtpUrl()
##' @return A list of download FTP urls
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @rdname genomeFTP
##' @export
##'
##' 
listFileFtpUrl <- function(ftpUrl) {
  
  fileStr <- getURL(ftpUrl, dirlistonly = TRUE)
  fileVec <- unlist(strsplit(fileStr, split = '\n', fixed = TRUE))

  fileUrl <- paste0(ftpUrl, fileVec)

  return(fileUrl)
}



##' @param gffUrl FTP url of one gff file.
##' @return A table of raw gff file
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname genomeFTP
##' @references How to read in gz files into R from FTP URL \url{https://stackoverflow.com/questions/9548630/read-gzipped-csv-directly-from-a-url-in-r}
##' @export
##'
##' 
read.gffUrl <- function(gffUrl) {

  ## read gz files from stream
  conStream <- gzcon(url(gffUrl))
  eachLine <- readLines(conStream)
  gffMat <- read.table(textConnection(eachLine),
                       sep = '\t',
                       quote = '',
                       stringsAsFactors = FALSE)

  ## remove rows and columns names
  row.names <- NULL
  col.names <- NULL

  return(gffMat)
  
}


##' @inheritParams getGenomicGenes
##' @return A list of genomes containing gene location information. The locus_tags is used for the gene names.
##' @examples
##' draLocs <- getLocsfgff('dra')
##' \dontrun{
##' hxaLocs <- getLocsfgff('hxa')
##' csuLocs <- getLocsfgff('csu')}
##' @rdname genomeFTP
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##' 
getLocsfgff <- function(KEGGSpe) {


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
  speUrl <- getSpeFtpUrl(KEGGSpe, database = c('GenBank', 'RefSeq'))
  if (length(speUrl) == 2) {
    ## choose RefSeq
    speUrl <- speUrl[names(speUrl) == 'RefSeq']
  } else {}
  speFiles <- listFileFtpUrl(speUrl)
  gffUrl <- speFiles[grepl('gff', speFiles)]
  spegff <- read.gffUrl(gffUrl)
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
  
