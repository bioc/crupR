# EDIT: Added more details to the man pages EDIT END
#' Normalization of ChIP-seq counts by input
#' 
#' @description
#' This function normalizes the ChIP-seq counts for the three histone modifications (H3K4me3, H3K4me1, H2K27ac) using an input experiment, 
#' as they are needed for the enhancer prediction in the next step. This step is optional, but recommended.
#' @details
#' The function normalizes the ChIP-seq profiles of the three histone modifications using the counts of the provided control/input experiment:
#' First, the binned counts (bin size: 100bp) are computed for every target (H3K4me3, H3K4me1, H2K27ac, Input) with low-quality reads being fitlered.
#' Next, the counts of the histone modifications are normalized by input by forming the log2 fold change of the counts:
#' counts_norm(bin) = log2((counts_raw(bin) + 1)/(counts_input(bin) + 1)).
#' In case input experiments are not available, the input free mode will calculate the binned counts, 
#' but won't further normalize them.
#' @param metaData A data frame containing all important information about the ChIP-seq experiments, i.e, for which histone modification they were conducted, 
#' path to the file, condition or replicate, path to control/input experiments.
#' @param condition The condition of the sample that is supposed to be normalized
#' @param replicate The replicate number of the sample that is supposed to be normalized
#' @param genome Reference genome. Either a character (one of: mm9, mm10, hg19, hg38) or a SeqInfo object
#' @param mapq Minimum mapping quality of the reads that should be included. Default is 10.
#' @param sequencing The type of sequencing that was used. Possible options are paired and single.
#' @param input.free [LOGICAL] Whether or not to run the function in the input free mode, in case there is no input experiment available. Default is FALSE.
#' @param chroms A vector of strings. Specify the relevant chromosomes if needed. Then, only the reads on these chromosomes are considered for normalization and used in further analysis. Default is 'all'.
#' @param BPPARAM An object of class SerialParam that is used as input for the BiocParallel functions.
#' @return GRanges object containing the normalized counts
#' @examples
#'
#'files <- c(system.file('extdata', 'Condition1.H3K4me1.bam', package='crupR'),
#'           system.file('extdata', 'Condition1.H3K4me3.bam', package='crupR'),
#'           system.file('extdata', 'Condition1.H3K27ac.bam', package='crupR'),
#'           system.file('extdata', 'Condition2.H3K4me1.bam', package='crupR'),
#'           system.file('extdata', 'Condition2.H3K4me3.bam', package='crupR'),
#'           system.file('extdata', 'Condition2.H3K27ac.bam', package='crupR'))
#'
#'inputs <- c(rep(system.file('extdata', 'Condition1.Input.bam', package='crupR'), 3), 
#'            rep(system.file('extdata', 'Condition2.Input.bam', package='crupR'),3))
#'                                                        
#'metaData <- data.frame(HM = rep(c('H3K4me1','H3K4me3','H3K27ac'),2),
#'                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
#'                       bamFile = files, inputFile = inputs)
#'                                             
#' normalize(metaData = metaData, condition = 1, replicate = 1,
#'     genome = 'mm10', sequencing = 'paired')
#' 
#' #Example for a customized genome:
#' genome = GenomeInfoDb::Seqinfo( seqnames=c('chr3', 'chr4', 'chrM'),
#'                                 seqlengths=c(1000, 2000, 500))
#'
#' @export
#' @importFrom bamsignals bamProfile
#' @importFrom Rsamtools scanBamHeader BamFile
#' @importFrom stats knots
#' @importFrom GenomicRanges seqinfo mcols tileGenome width
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevels getChromInfoFromUCSC
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom S4Vectors metadata
#' @importFrom BiocParallel SerialParam

normalize <- function(metaData, condition, replicate, genome, mapq = 10, sequencing,
    input.free = FALSE, chroms = NULL, BPPARAM = BiocParallel::SerialParam()) {
    start_time <- Sys.time()
    message("\n")

    m <- metaData[which(metaData$condition == condition & metaData$replicate == replicate),
        ]
    if (nrow(m) != 3) {
        if (nrow(m) == 0)
            stop(paste0("The chosen combination of condition and replicate is not valid.\n There are no files for condition ",
                condition, " replicate ", replicate)) else if (nrow(m) > 3)
            stop(paste0("There are too many bam files for condition ", condition,
                " replicate ", replicate)) else stop(paste0("There are not enough bam files for condition ", condition,
            " replicate ", replicate))
    }

    for (path in m$bamFile) check_file(path)
    if (!input.free)
        for (path in m$inputFile) check_file(path)

    hm <- c()
    for (i in seq_along(hm_values)) hm <- c(hm, normalizePath(as.vector(m[which(m$HM ==
        hm_values[i]), ]$bamFile)[1]))
    if (input.free == FALSE)
        for (i in seq_along(hm_values)) hm <- c(hm, normalizePath(as.vector(m[which(m$HM ==
            hm_values[i]), ]$inputFile)[1]))

    bamHM <- hm[seq_len(3)]
    if (input.free == TRUE) {
        bamInp <- NULL
    } else {
        bamInp <- unique(hm[4:length(hm)])
    }

    ### EDIT: BSgenome objects have been replaced by using GenomeInfoDb to get
    ### the Seqinfo The Bsgenome packages have been removed from the
    ### Description EDIT END ####
    if (!is(genome, "Seqinfo") && !(genome %in% genome_values)) {
        {
            stop(paste0("Your genome is neither one of ", paste0(genome_values, collapse = ","),
                " nor is it a valid Seqinfo object."))
        }
    } else {
        genome <- GenomeInfoDb::getChromInfoFromUCSC(genome, as.Seqinfo = TRUE)
    }


    ################################################################## 
    # get and  adjust binned genome
    ################################################################## 

    message("Prepare the binned genome ...\n")  #EDIT: replaced cat() with message() 

    gr <- get_binned_genome(genome, chr = chroms)
    bf <- Rsamtools::BamFile(bamHM[1])
    si <- GenomicRanges::seqinfo(bf)
    if (!("chr1" %in% names(Rsamtools::scanBamHeader(bamHM[1])[[1]]$targets)))
        GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(si)[1]

    ################################################################## 
    # get counts from ChIP-seq experiments
    ################################################################## 

    message("Get summarized counts from ChIP-seq experiments ...\n")  #EDIT: replaced cat() with message()

    names <- c(hm_values, paste0("Input_", hm_values))
    if (input.free)
        names <- hm_values
    if (length(unique(bamInp)) == 1) {
        bamInp <- bamInp[1]
        names <- c(hm_values, "Input_All")
    }

    ### EDIT: replaced parallel::mcmapply with BiocParallel::bpmapply EDIT END
    ### ###
    if (sequencing == "paired") {
        counts <- BiocParallel::bpmapply(bamsignals::bamProfile, bampath = c(bamHM,
            bamInp), MoreArgs = list(gr = gr, binsize = 100, mapqual = mapq, ss = FALSE,
            paired.end = "midpoint", filteredFlag = 1024, verbose = FALSE), SIMPLIFY = FALSE,
            BPPARAM = BPPARAM)

    } else if (sequencing == "single") {
        counts <- BiocParallel::bpmapply(bamsignals::bamProfile, bampath = c(bamHM,
            bamInp), MoreArgs = list(gr = gr, binsize = 100, mapqual = mapq, ss = FALSE,
            shift = 100, filteredFlag = 1024, verbose = FALSE), SIMPLIFY = FALSE,
            BPPARAM = BPPARAM)
    } else stop("Sequencing parameter is not valid.\n Choose one of:", paste(sequencing_values,
        collapse = ","))

    for (i in seq_along(counts)) counts[[i]] <- unlist(as.list(counts[[i]]))
    names(counts) <- names

    ################################################################## 
    # get counts from ChIP-seq experiments
    ##################################################################
    if (!input.free) {
        message("Normalize histone modifications by Input ...\n")  #EDIT: replaced cat() with message()

        if ("Input_All" %in% names(counts))
            countsNorm <- lapply(hm_values, function(x) log2((counts[[x]] + 1)/(counts[[paste0("Input_All")]] +
                1))) else countsNorm <- lapply(hm_values, function(x) log2((counts[[x]] + 1)/(counts[[paste0("Input_",
            x)]] + 1)))
    } else countsNorm <- counts

    ################################################################## 
    # create data matrix for all normalized ChIP-seq experiments
    ##################################################################
    
    message("Create summarized data matrix ...\n")  #EDIT: replaced cat() with message()

    GenomicRanges::mcols(gr) <- matrix(unlist(countsNorm), ncol = 3, byrow = FALSE,
        dimnames = list(NULL, hm_values))
    GenomeInfoDb::seqlevels(gr) <- paste0("chr", gsub("chr|Chr", "", GenomeInfoDb::seqlevels(gr)))

    n <- GenomicRanges::mcols(gr)[, "H3K4me1"] + abs(min(GenomicRanges::mcols(gr)[,
        "H3K4me1"])) + 1
    d <- GenomicRanges::mcols(gr)[, "H3K4me3"] + abs(min(GenomicRanges::mcols(gr)[,
        "H3K4me3"])) + 1
    GenomicRanges::mcols(gr)[, "ratio"] <- log2(n/d)

    message(paste0("time: ", format(Sys.time() - start_time), "\n"))
    # EDIT: adjusted the output to a single GRanges object with meta data EDIT
    # END ####
    S4Vectors::metadata(gr) <- m
    return(gr)
}
