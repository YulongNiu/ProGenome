##' Retrieve genomic genes.
##'
##' Retrieve the whole genomic gene from a prokaryotic species. The function returns gene location in all genomes, including the cellular genoms and plasmids. For some gene, multiple location may return. If the genome has no gene featurs, "NULL" will be returned.
##' getGenesfGenomes(): retrieve whole genomic genes from NCBI genomes.
##' getGenesfGeneids(): retrieve whole genomic genes from NCBI gene IDs.
##' getGenomicGenes(): a wrapped of these top two function, which trys getGenesfGeneids() at first and then getGenesfGenomes() if the genome is from GenBank.
##' @title Retrieve genomic genes
##' @param KEGGSpe A KEGG species ID
##' @return A list. Each element represents a genome, and the gene IDs are in KEGG format.
##' @examples
##' ## no genomic genes
##' noGenes <- getGenomicGenes('csu')
##' \dontrun{
##' ## from GenBank
##' ## one genome with three plasmids
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
    wholeGenes <- getGenesfGenomes(KEGGSpe)
  }
  else if (!genLogic && refLogic) {
    ## source from RefSeq
    wholeGenes <- getGenesfGeneids(KEGGSpe)
  }

  return(wholeGenes)
}



##' @inheritParams getGenomicGenes
##' @return  A list. Each element represents a genome, and the gene IDs are in KEGG format.
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname genomicGenes
##' @importFrom KEGGAPI getKEGGSpeInfo
##' @importFrom NCBIAPI singleGenomeAnno
##' @importFrom stringr str_extract str_trim
##' @importFrom foreach foreach %do%
##' @export
##'
##' 
getGenesfGenomes <- function(KEGGSpe) {
  
  speInfo <- getKEGGSpeInfo(KEGGSpe)
  
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
  wholeNum <- sapply(strsplit(wholeNum, split = ':', fixed = TRUE), '[[', 2)
  wholeNum <- str_trim(wholeNum)
  wholeNames <- c(genomeName, plasmidName)

  wholeGenes <- foreach (i = 1:length(wholeNum)) %do% {
    
    eachAnno <- singleGenomeAnno(wholeNum[i], type = 'gene', n = 4)

    if (is.null(eachAnno)) {
      ## no genomic genes were found
      eachMat <- NULL
    } else {
      eachNames <- names(eachAnno)
      eachLocs <- lapply(eachAnno, function(x) {
        ## may have multiple gene locations for one gene
        geneLoc <- x$GBInterval
        return(geneLoc)
      })
      eachRepNum <- sapply(eachLocs, nrow)
      eachLocs <- do.call(rbind, eachLocs)
      
      eachMat <- cbind(paste(KEGGSpe,
                             rep(eachNames, eachRepNum),
                             sep = ':'),
                       eachLocs)
      colnames(eachMat) <- c('geneName', 'From', 'To')
    }

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
##' @export
##'
##' 
getGenesfGeneids <- function(KEGGSpe) {

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
  genomeSeq <- paste0(genomeType, genomeSeq)
  
  genomeLoc <- lapply(speInfo, function(x) {
    ## one gene may have multiple location
    return(x$GenomicInfo[, 3:4, drop = FALSE])
  })
  repNum <- sapply(genomeLoc, nrow)
  genomeLoc <- do.call(rbind, genomeLoc)

  
  genomeMat <- cbind(rep(names(speInfo), repNum),
                     genomeLoc,
                     rep(genomeSeq, repNum),
                     deparse.level = 0)
  colnames(genomeMat)[1:3] <- c('geneName', 'From', 'To')
  rownames(genomeMat) <- NULL

  wholeGenesIdx <- split(1:nrow(genomeMat), genomeMat[, 4])
  wholeGenes <- lapply(wholeGenesIdx, function(x) {
    return(genomeMat[x, 1:3])
  })

  return(wholeGenes)
}
