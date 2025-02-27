# EDIT: Added more details to the man pages EDIT END
#' Find the target genes of the regulatory, condition-specific clusters
#' @description
#' This functions aims to connect the differential enhancers identified in the prior step (getDynamics) to the potential target genes.
#' It exploits normalized gene expression counts for every samples, derived from RNA-seq experiments, and correlates the gene activities 
#' with predicted enhancer probabilities over each replicate of all conditions. If there is high correlation between a gene and an enhancer,
#' the respective gene is considered a target gene of the enhancer. The number of comparisons is narrowed down by only comparing enhancers to genes within the same 
#' topologically associating domain (TAD). Alternatively, one can also just check the nearest gene for every enhancer.
#' 
#' @details
#' This functions depends on the presence of gene expression counts for every every sample. These can 
#' be created by using e.g. DESeq2. For an example code of how to get the counts, check /inst/script/crupR_example_files.txt. 
#' 
#' This functions compares the gene expression counts of the candidate target genes for a differential enhancer and correlates
#' them to the enhancer probabilities as computed in getEnhancers. If the correlation surpasses a threshold (default: 0.9), the candidate gene is considered
#' a target gene. There are two ways to define candidate target genes: Either by using the topologically associating domains (TADs) as boundaries for potential interactions or
#' in case no TAD annotations are available, the nearest gene of the respective differential enhancer is considered. 
#' 
#' @param data condition-specific clusters in a GRanges object (output of getDynamics())
#' @param expr SummarizedExperiment object containing the gene expression counts of RNA-seq experiments for each condition and its replicates
#' @param genome Genome used in the .bam files of the RNA-seq experiments. Possible options are 'mm9', 'mm10', 'hg19' and 'hg38'.
#' @param TAD.file Path to the TAD file to use for finding the target genes. If set to NULL, the default file is used (only if the 'mm10' genome was used)
#' @param cutoff cut-off for correlation between cluster and gene. Default is 0.9.
#' @param nearest [LOGICAL] If set, the nearest gene is taken to build the regulatory regions.
#' @param BPPARAM An object of class SerialParam that is used as input for the BiocParallel functions.
#' @return GRanges object containing the dynamic regulatory units
#' @examples
#' 
#' #first get the output of crupR::getDynamics so skip this
#' files <- c(system.file('extdata', 'Condition1.H3K4me1.bam', package='crupR'),
#'            system.file('extdata', 'Condition1.H3K4me3.bam', package='crupR'),
#'            system.file('extdata', 'Condition1.H3K27ac.bam', package='crupR'),
#'            system.file('extdata', 'Condition2.H3K4me1.bam', package='crupR'),
#'            system.file('extdata', 'Condition2.H3K4me3.bam', package='crupR'),
#'            system.file('extdata', 'Condition2.H3K27ac.bam', package='crupR'))
#' inputs <- rep(system.file('extdata', 'Condition1.Input.bam', package='crupR'), 3)
#' inputs2 <- rep(system.file('extdata', 'Condition2.Input.bam', package='crupR'), 3)  
#' metaData <- data.frame(HM = rep(c('H3K4me1', 'H3K4me3', 'H3K27ac'),2),
#'                        condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
#'                        bamFile = files, inputFile = c(inputs, inputs2))
#' clusters <- readRDS(system.file('extdata', 'differential_enhancers.rds', package='crupR'))
#' S4Vectors::metadata(clusters) <- metaData
#' #load your SummarizedExperiments object containing the gene expressions counts and run the function
#' expr <- readRDS(system.file('extdata', 'expressions.rds', package='crupR'))
#' getTargets(data = clusters, expr = expr, genome = 'mm10')
#'
#' @export
#' @importFrom GenomicRanges mcols GRanges makeGRangesFromDataFrame start end promoters strand seqnames nearest distance
#' @importFrom GenomicFeatures genes
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @importFrom IRanges IRanges subsetByOverlaps %within%
#' @import TxDb.Mmusculus.UCSC.mm9.knownGene
#' @import TxDb.Mmusculus.UCSC.mm10.knownGene
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom utils read.table
#' @importFrom SummarizedExperiment rownames rowRanges assay
#' @importFrom S4Vectors metadata
#' @importFrom BiocParallel SerialParam

getTargets <- function(data, expr = NULL, genome, TAD.file = NULL, cutoff = 0.9,
    nearest = FALSE, BPPARAM = BiocParallel::SerialParam()) {
    start_time <- Sys.time()
    message("\n")

    # EDIT: adjusted these steps to the new input type
    metaData <- S4Vectors::metadata(data)
    conds <- unique(metaData$condition)
    gr <- data
    GenomeInfoDb::seqlevels(gr) <- paste0("chr", gsub("chr|Chr", "", GenomeInfoDb::seqlevels(gr)))
    GenomeInfoDb::genome(gr) <- genome

    IDs <- list()
    for (i in seq_along(conds)) {
        sub <- subset(metaData, condition == conds[i])
        if (is.numeric(conds[i]))
            IDs[[i]] <- paste0("cond", conds[i], "_", unique(sub$replicate)) else IDs[[i]] <- paste0(conds[i], "_", unique(sub$replicate))
    }

    if (!(genome %in% genome_values))
        stop(paste0("Genome ", genome, " currently not provided. Choose one of:",
            paste(genome_values, collapse = ",")))
    if (cutoff < 0.5 | cutoff > 1)
        stop(paste0(cutoff, " is not in range [0.5,1]."))

    if (genome == "mm10") {
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    } else if (genome == "mm9") {
        txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
    } else if (genome == "hg19") {
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if (genome == "hg38") {
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }

    if (nearest == FALSE) {
        if (is.null(TAD.file) & genome == "mm10") {
            TAD.file <- system.file("extdata", "mESC_mapq30_KR_all_TADs.bed", package = "crupR")
        } else if (is.null(TAD.file) & genome != "mm10") {
            stop("You have to provide your own file with TAD domains (fitting to the genome of choice).")
        }
    }

    if (is.null(expr)) {
        warning("No gene expression counts provided... setting nearest to TRUE")
        nearest <- TRUE
    }

    TAD <- GenomicRanges::GRanges()
    if (nearest == FALSE) {
        TAD <- read.table(TAD.file, col.names = GR_header_short)
        TAD <- GenomicRanges::makeGRangesFromDataFrame(TAD[which((TAD$end - TAD$start) >
            0), ])
        GenomeInfoDb::seqlevels(TAD) <- paste0("chr", gsub("chr|Chr", "", GenomeInfoDb::seqlevels(TAD)))
    }
    units <- get_units(gr, expr, TAD, unlist(IDs), BPPARAM = BPPARAM, cutoff, txdb,
        nearest)
    message(paste0("time: ", format(Sys.time() - start_time), "\n"))
    # EDIT: adjusted the output to a single GRanges object with meta data
    S4Vectors::metadata(units) <- metaData
    return(units)
}
