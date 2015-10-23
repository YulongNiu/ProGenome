##' Get genome FTP url
##'
##' getFtpUrl(): get the FTP url storing genome information from GenBank or RefSeq database. This function is mainly used to get the species genome assembly information.
##' listFtpFileUrl(): list the file download FTP urls.
##' @title Retrieve genome FTP URL
##' @inheritParams getGenomicGenes
##' @param database "GenBank" or "RefSeq".
##' @return FTP url.
##' @examples
##' ecoliUrl <- getSpeFtpUrl('eco')
##' ecoliFiles <- listFileFtpUrl(ecoliUrl)
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
