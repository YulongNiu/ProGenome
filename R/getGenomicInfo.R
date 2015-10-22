##' Get genome FTP url
##'
##' Get the FTP url storing genome information from GenBank or RefSeq database. This function is mainly used to get the species genome assembly information.
##' @title Retrieve genome FTP URL
##' @inheritParams getGenomicGenes
##' @param database "GenBank" or "RefSeq".
##' @return FTP url.
##' @examples
##' ecoliUrl <- getFTPUrl('eco')
##' @importFrom KEGGAPI getKEGGSpeInfo
##' @importFrom stringr str_extract str_detect
##' @importFrom RCurl getURL
##' @importFrom xml2 read_xml xml_attr
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##' 
getFTPUrl <- function(KEGGSpe, database = 'GenBank') {

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
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  return(ftpUrl)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
