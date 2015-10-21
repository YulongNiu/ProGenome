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
    speGenes <- convKEGG(KEGGSpe, 'ncbi-geneid')
    speGeneNames <- sapply(strsplit(speGenes[, 1], split = ':', fixed = TRUE), '[[', 2)
    speInfo <- getNCBIGenesInfo(speGeneNames, n = 4)
    names(speInfo) <- speGenes[, 2]

    ## extract genomic
    genomeType <- sapply(speInfo, function(x) {
      return(x$GeneticSource)
    })
    ##~~~~~~~~~~~~~~transform genome type~~~~~~~~~~~~~~~~
    ## c('genomic', 'plasmid CP1', 'plasmid CP1', 'plasmid MP1') -->
    ## c('genome', 'plasmid1', 'plasmid2')
    genomicLogic <- genomeType == 'genomic'
    plasmidName <- unique(genomeType[!genomicLogic])
    genomeType[genomicLogic] <- 'genome'
    if (length(plasmidName) > 0) {
      ## has plasmid
      genomeType[!genomicLogic] <- as.character(
        factor(genomeType[!genomicLogic],
               labels = paste0('plasmid',
                               1:length(plasmidName))))
    } else {}
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
  }

  return(wholeGenes)
}
