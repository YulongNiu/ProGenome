##' Get genomic information NCBI FTP url
##'
##' GetSpeFtpUrl(): get the FTP url storing genome information from GenBank or RefSeq database. This function is mainly used to get the species genome assembly information.
##'
##' AutoSpeFtpUrl(): automatically choose the FTP url. This function trys the RefSeq at first and then GenBank.
##' 
##' ListFtpFileUrl(): list the file download FTP urls.
##' 
##' read.gff(): read in a raw gff file (or gz file) from the local disk or web url.
##' 
##' @title Retrieve genome FTP URL
##' @inheritParams getGenomicGenes
##' @param database "GenBank", "RefSeq", or c("GenBank", "RefSeq").
##' @return FTP url.
##' @examples
##' ## E. coli genome GenBank and RefSeq FTP
##' ecoliUrl <- GetSpeFtpUrl('eco', c('GenBank', 'RefSeq'))
##' ## automatically choose FTP
##' ecoliAutoUrl <- AutoSpeFtpUrl('eco')
##' ## list E. coli FTP files
##' ecoliFiles <- ListFileFtpUrl(ecoliAutoUrl)
##' ## read in gff file
##' gffUrl <- ecoliFiles[grepl('gff', ecoliFiles)]
##' ecoligff <- read.gff(gffUrl, isurl = TRUE, isgz = TRUE)
##' \dontrun{
##' draliUrl <- GetSpeFtpUrl('dra')
##' draliFiles <- ListFileFtpUrl(draliUrl)
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
GetSpeFtpUrl <- function(KEGGSpe, database = 'GenBank') {

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
    eachName <- sapply(strsplit(xml_text(eachFtpXml),
                                ' ',
                                fixed = TRUE),
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


##' @inheritParams getGenomicGenes
##' @return AutoSpeFtpUrl(): RefSeq or GenBank FTP url
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname genomeFTP
##' @export
##'
##' 
AutoSpeFtpUrl <- function(KEGGSpe) {
  speUrl <- GetSpeFtpUrl(KEGGSpe, database = c('GenBank', 'RefSeq'))

  if (length(speUrl) == 2) {
    ## choose RefSeq
    speUrl <- speUrl[names(speUrl) == 'RefSeq']
  } else {}

  return(speUrl)
}


##' @param ftpUrl FTP url that can be the output of function  getFtpUrl()
##' @return ListFileFtpUrl(): a list of download FTP urls
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @rdname genomeFTP
##' @export
##'
##' 
ListFileFtpUrl <- function(ftpUrl) {
  
  fileStr <- getURL(ftpUrl, dirlistonly = TRUE)
  fileVec <- unlist(strsplit(fileStr, split = '\n', fixed = TRUE))

  fileUrl <- paste0(ftpUrl, fileVec)

  return(fileUrl)
}



##' @param filePath A local file path or a web url.
##' @param isurl Whether a url (TRUE) or not (FALSE).
##' @param isgz Whether a gzfile (TRUE) or not (FALSE).
##' @return read.gff(): a table of raw gff file
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname genomeFTP
##' @references How to read in gz files into R from FTP URL \url{https://stackoverflow.com/questions/9548630/read-gzipped-csv-directly-from-a-url-in-r}
##' @export
##'
##' 
read.gff <- function(filePath, isurl = FALSE, isgz = FALSE) {

  if (isurl && isgz) {
    ## read gz gff file from url
    conStream <- gzcon(url(filePath))
  }
  else if (isurl && !isgz) {
    ## read gff file from url
    conStream <- url(filePath)
  }
  else if (!isurl && isgz) {
    ## read gz gff file from local disk
    conStream <- gzcon(gzfile(filePath))
  }
  else {
    ## read gff file from local disk
    conStream <- file(filePath)
  }
  
  eachLine <- textConnection(readLines(conStream))
  close(conStream)
  gffMat <- read.table(eachLine,
                       sep = '\t',
                       quote = '',
                       stringsAsFactors = FALSE)
  close(eachLine)

  ## remove rows and columns names
  row.names <- NULL
  col.names <- NULL

  return(gffMat)
}

