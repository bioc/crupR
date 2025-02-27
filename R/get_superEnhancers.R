# EDIT: Added more details to the man pages EDIT END
#' Find enhancer peaks and super enhancers
#' @description
#' An optional intermediate step to summarize the genome-wide enhancer predictions from the prior step (getEnhancers). First, enhancer peaks are detected.
#' Next, neighbouring enhancer peaks are further grouped together into super-enhancers, which we define as clusters of proximal enhancer regions.
#' 
#' @details
#' First, enhancer peaks are identified. All bins with an enhancer probability >=0.5 are sorted in a descending manner by their probabilities.
#' The bins are then extended by 5 bins up and downstream resulting in regions of size 1100bp. Overlapping regions are discarded, 
#' while keeping the region with the higher enhancer probability. The resulting list of enhancer peaks is further summarized 
#' into clusters by grouping peaks in close vicinity (default: max. 12.5kb up-/downstream distance) together. These clusters are supposed to reflect
#' super-enhancers.
#' 
#' @param data enhancer prediction values in a GRanges object (output of getEnhancers())
#' @param cutoff Cut-off for enhancer probabilities. Default is 0.5.
#' @param distance Maximum distance (in bp) for clustering. Default is 12500.
#' @param BPPARAM An object of class SerialParam that is used as input for the BiocParallel functions.
#' @return A list containing the enhancer prediction values as a GRanges object (like the input data),
#' the enhancer peak calls as a GRanges object (can be exported as a bedGraph) and the clusters of peaks 
#' (super-enhancers) in a GRanges object (can be exported as BED file)
#' 
#' @examples
#'
#' # first recreate the output of crupR:getPredictions
#' files <- c(system.file('extdata', 'Condition2.H3K4me1.bam', package='crupR'),
#'            system.file('extdata', 'Condition2.H3K4me3.bam', package='crupR'),
#'            system.file('extdata', 'Condition2.H3K27ac.bam', package='crupR'))
#' inputs <- rep(system.file('extdata', 'Condition2.Input.bam', package='crupR'), 3)
#' #create the metaData frame
#' metaData <- data.frame(HM = c('H3K4me1','H3K4me3','H3K27ac'),
#'                        condition = c(2,2,2), replicate = c(1,1,1),
#'                        bamFile = files, inputFile = inputs)
#' prediction <- readRDS(system.file('extdata', 'condition2_predictions.rds', package='crupR'))                                                                                           
#' S4Vectors::metadata(prediction) <- metaData
#' #run the function
#' se <- getSE(data = prediction)
#' @export
#' @importFrom GenomicRanges start reduce mcols findOverlaps width seqnames GRanges
#' @importFrom S4Vectors queryHits subjectHits metadata
#' @importFrom BiocParallel SerialParam


getSE <- function(data, cutoff = 0.5, distance = 12500, BPPARAM = BiocParallel::SerialParam()) {
    start_time <- Sys.time()
    message("\n")

    if (cutoff < 0 | cutoff > 1)
        stop(paste0(cutoff, " is not a valid cutoff. Please choose a cutoff between 0 and 1."))
    if (distance < 0)
        stop(paste0(distance, " is not a valid distance. Please choose a distance greater than 0."))
    # EDIT: Adjuested error messages
    if (!("prob" %in% colnames(GenomicRanges::mcols(data)))) {
        stop("Input needs to be the output of getEnhancers")
    }

    peaks <- GenomicRanges::GRanges()
    peaks <- c(peaks, get_enhancerPeaks(data, cutoff, BPPARAM))  #EDIT: use new BPPRAM 
    message(paste0(length(peaks), " single enhancer peak(s) found.\n"))  #EDIT: replaced cat() with message()

    cluster <- GenomicRanges::GRanges()
    cluster <- c(cluster, get_enhancerCluster(peaks, distance, BPPARAM))  #EDIT: use new BPPRAM 
    message(paste0(length(cluster), " enhancer cluster(s) found.\n"))  #EDIT: replaced cat() with message()

    message(paste0("time: ", format(Sys.time() - start_time), "\n"))  #EDIT: replaced cat() with message()
    # EDIT: adjusted the output to a list w/o meta data separately
    return(list(D = data, peaks = peaks, cluster = cluster))
}
