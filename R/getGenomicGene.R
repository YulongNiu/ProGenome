##' Download genomic genes.
##'
##' Download the whole genomic gene from a prokaryotic species. The function returns gene location in all genomes, including the cellular genoms and plasmids.
##' @title Retrieve genomic genes
##' @param KEGGSpe A KEGG species ID
##' @return A list. Each element represents a genome, and the gene IDs are in KEGG format.
##' @examples
##' \dontrun{
##' ## from GenBank
##' ## one genome with four plasmids
##' hxaGenes <- getGenomicGenes('hxa')
##' 
##' ## from RefSeq
##' ## two genomes with two plasmids
##' draGenes <- getGenomicGenes('dra')}
##' @importFrom KEGGAPI getKEGGSpeInfo convKEGG
##' @importFrom NCBIAPI singleGenomeAnno getNCBIGenesInfo
##' @importFrom stringr str_extract
##' @importFrom foreach foreach %do%
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname genomicGenes
##' @export
##'
##' 
getGenomicGenes <- function(KEGGSpe) {

  ##~~~~~~~~~~~~~~~~~~~~source of genomes~~~~~~~~~~~~~~~~
  speInfo <- getKEGGSpeInfo(KEGGSpe)
  sourceEle <- speInfo$`Data source`

  genLogic <- grepl('GenBank', sourceEle)
  refLogic <- grepl('RefSeq', sourceEle)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (genLogic && !refLogic) {
    ## source from GenBank
    ## extract GB number
    genomeVec <- c(speInfo$Chromosomes[, 2])
    ## the length of genomeNum must > 0
    genomeNum <- str_extract(genomeVec, 'GB: \\w+')
    genomeName <- paste0('genome', 1:length(genomeNum))


    ## NULL if not found
    plasmidVec <- c(speInfo$Plasmids[, 2])
    ## the length of plamidNum may be 0, then return character(0)
    plasmidNum <- str_extract(plasmidVec, 'GB: \\w+')
    if (length(plasmidNum) == 0) {
      plasmidName <- character(0)
    } else {
      plasmidName <- paste0('plasmid', 1:length(plasmidNum))
    }
    
    ## whole genomes number and names
    wholeNum <- c(genomeNum, plasmidNum)
    wholeNames <- c(genomeName, plasmidName)

    wholeGenes <- foreach (i = 1:length(wholeNum)) %do% {
      eachAnno <- singleGenomeAnno(wholeNum[i], type = 'CDS')
      eachNames <- names(eachAnno)
      eachLocs <- lapply(eachAnno, function(x) {
        geneLoc <- x$GBInterval
        return(geneLoc)
      })
      eachLocs <- do.call(rbind, eachLocs)
      eachMat <- cbind(paste(KEGGSpe, eachNames, sep = ':'), eachLocs)
      colnames(eachMat) <- c('geneName', 'From', 'To')

      return(eachMat)
    }

    names(wholeGenes) <- wholeNames
    
  }
  else if (!genLogic && refLogic) {
    ## source from RefSeq
    wholeGenes <- getGenesfGeneids(KEGGSpe)
  }

  return(wholeGenes)
}



##' @inheritParams getGenomicGenes
##' @param ... Inherited parameters from singleGenomeAnno() of NCBIAPI package
##' @return  A list. Each element represents a genome, and the gene IDs are in KEGG format.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname genomicGenes
##' @importFrom KEGGAPI getKEGGSpeInfo
##' @importFrom NCBIAPI singleGenomeAnno
##' @importFrom stringr str_extract str_trim
##' @importFrom foreach foreach %do%
##' @keywords internal
##' @param KEGGSpe 
##'
##' 
getGenesfGenomes <- function(KEGGSpe, ...) {
  
  ##~~~~~~~~~~~~~~~~~~~~source of genomes~~~~~~~~~~~~~~~~
  speInfo <- getKEGGSpeInfo(KEGGSpe)
  sourceEle <- speInfo$`Data source`

  genLogic <- grepl('GenBank', sourceEle)
  refLogic <- grepl('RefSeq', sourceEle)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## source from GenBank
  ## extract GB number
  genomeVec <- c(speInfo$Chromosomes[, 2])
  ## the length of genomeNum must > 0
  genomeNum <- str_extract(genomeVec, 'GB: \\w+')
  genomeName <- paste0('genome', 1:length(genomeNum))


  ## NULL if not found
  plasmidVec <- c(speInfo$Plasmids[, 2])
  ## the length of plamidNum may be 0, then return character(0)
  plasmidNum <- str_extract(plasmidVec, 'GB: \\w+')
  if (length(plasmidNum) == 0) {
    plasmidName <- character(0)
  } else {
    plasmidName <- paste0('plasmid', 1:length(plasmidNum))
  }
  
  ## whole genomes number and names
  wholeNum <- c(genomeNum, plasmidNum)
  wholeNum <- str_trim(wholeNum)
  wholeNum <- sapply(strsplit(wholeNum, split = ':', fixed = TRUE), '[[', 2)
  wholeNames <- c(genomeName, plasmidName)

  wholeGenes <- foreach (i = 1:length(wholeNum)) %do% {
    
    eachAnno <- singleGenomeAnno(wholeNum[i], type = 'CDS', n = 4)
    
    eachNames <- names(eachAnno)
    
    eachLocs <- lapply(eachAnno, function(x) {
      ## may have multiple gene locations for one gene
      geneLoc <- x$GBInterval
      return(geneLoc)
    })
    eachRepNum <- sapply(eachLocs, nrow)
    eachLocs <- do.call(rbind, eachLocs)
    
    eachMat <- cbind(paste(KEGGSpe,
                           rep(eachNames, each = eachRepNum),
                           sep = ':'),
                     eachLocs)
    colnames(eachMat) <- c('geneName', 'From', 'To')

    return(eachMat)
  }

  names(wholeGenes) <- wholeNames

  return(wholeGenes)
}



##' @inheritParams getGenomicGenes
##' @return  A list. Each element represents a genome, and the gene IDs are in KEGG format.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname genomicGenes
##' @importFrom KEGGAPI convKEGG
##' @importFrom NCBIAPI getNCBIGenesInfo
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##'
##' 
getGenesfGeneids <- function(KEGGSpe) {

  speGenes <- convKEGG(KEGGSpe, 'ncbi-geneid')
  speGeneNames <- sapply(strsplit(speGenes[, 1], split = ':', fixed = TRUE), '[[', 2)
  speInfo <- getNCBIGenesInfo(speGeneNames, n = 4, maxEach = 10000)
  names(speInfo) <- speGenes[, 2]

  ## extract genomic
  genomeType <- sapply(speInfo, function(x) {
    return(x$GeneticSource)
  })
  ##~~~~~~~~~~~~~~transform genome type~~~~~~~~~~~~~~~~
  ## c('genomic', 'plasmid CP1', 'plasmid CP1', 'plasmid MP1') -->
  ## c('genome', 'plasmid1', 'plasmid2')
  ## must have 1 and only 1 'genomic'
  ## may have 0, 1, or more plasmids
  genomeType <- factor(genomeType)
  gLevels <- levels(genomeType)
  gLog <- gLevels == 'genomic'
  levels(genomeType)[gLog] <- 'genome'
  
  if (length(gLevels) > 1) {
    ## has plasmid
    levels(genomeType)[!gLog] <- paste0('plasmid',
                                       1:sum(!gLog))
  } else {}
  
  genomeType <- as.character(genomeType)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  genomeSeq <- sapply(speInfo, function(x) {
    return(x$Chromosome)
  })
  
  genomeLoc <- sapply(speInfo, function(x) {
    ## may have multiple gene location, here return the first
    return(x$GenomicInfo[1, 3:4])
  })
  genomeLoc <- t(genomeLoc)
  
  genomeMat <- cbind(names(speInfo),
                     genomeLoc,
                     paste0(genomeType, genomeSeq),
                     deparse.level = 0)
  colnames(genomeMat)[1:3] <- c('geneName', 'From', 'To')
  rownames(genomeMat) <- NULL

  wholeGenesIdx <- split(1:nrow(genomeMat), genomeMat[, 4])
  wholeGenes <- lapply(wholeGenesIdx, function(x) {
    return(genomeMat[x, 1:3])
  })

  return(wholeGenes)
}
