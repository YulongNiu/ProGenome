##' Get genomic information NCBI FTP url
##'
##' getFtpUrl(): get the FTP url storing genome information from GenBank or RefSeq database. This function is mainly used to get the species genome assembly information.
##' listFtpFileUrl(): list the file download FTP urls.
##' read.gffUrl(): read in raw gff file from a FTP url.
##' @title Retrieve genome FTP URL
##' @inheritParams getGenomicGenes
##' @param database "GenBank" or "RefSeq".
##' @return FTP url.
##' @examples
##' ## E. coli genome GenBank FTP
##' ecoliUrl <- getSpeFtpUrl('eco')
##' ## list E. coli FTP files
##' ecoliFiles <- listFileFtpUrl(ecoliUrl)
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
##' @importFrom xml2 read_xml xml_attr
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
  ## NCBI assembly url
  assBaseUrl <- 'http://www.ncbi.nlm.nih.gov/assembly/'
  assUrl <- paste0(assBaseUrl, assNum)

  webStr <- getURL(assUrl)
  webVec <- unlist(strsplit(webStr, split = '\n', fixed = TRUE))
  matchPattern <- paste0(database, ' FTP site')
  ftpLine <- webVec[str_detect(webVec, matchPattern)]

  ## find url
  ftpXml <- read_xml(ftpLine)
  ftpUrl <- xml_attr(ftpXml, 'href')

  ## ftp url often point to a folder
  ftpUrl <- paste0(ftpUrl, '/')
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  return(ftpUrl)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  
  conStream <- gzcon(url(gffUrl))
  
  eachLine <- readLines(conStream)
  
  gffMat <- read.table(textConnection(eachLine), sep = '\t', quote = '', stringsAsFactors = FALSE)

  return(gffMat)
  
}

