# EDIT: Added more details to the man pages EDIT END
#' Predicts the occurence of enhancers using the normalized HMs counts.
#' @description
#' This function gets the output of the prior normalization step as input and uses the three histone modifications to predict the probability of an active enhancer being present.
#' @details
#' First, the input-normalized histone modification counts are quantile-normalized. 
#' This way their distribution is matching to what the classifiers have been trained on.
#' Next, the two random forest classifiers go through every bin and use the 5 upstream and downstream flanking bins to make predictions. 
#' One classifiers predicts whether the current bin is an active region based on the surrounding histone modification patterns.
#' The other predicts whether the current bin is an enhancer or promoter, assuming it's already active. 
#' Last, the predicted probabilities of each classifier for every bin are mulitplied and their product forms the final binwise enhancer activity probability. 
#' 
#' @param data normalized ChIP-seq counts in a GRanges object (output of the normalize() step)
#' @param classifier The path of the classifier to use for the prediction. When set to NULL (default), the default classifier is used. 
#' Both classifiers are objects of class 'randomForest' as implemented by the randomForest package. For more information on this object, please check the documentation of randomForest.
#' @param all [LOGICAL] Whether to include the probabilities of the two individual random forests in the output. Default is FALSE.
#' @return GRanges object containing the enhancer probabilities for each 100bp bin
#' @examples
#' #first recreate the output of crupR::normalize (so skip this)
#' files <- c(system.file('extdata', 'Condition2.H3K4me1.bam', package='crupR'),
#'           system.file('extdata', 'Condition2.H3K4me3.bam', package='crupR'),
#'           system.file('extdata', 'Condition2.H3K27ac.bam', package='crupR'))
#' inputs <- rep(system.file('extdata', 'Condition2.Input.bam', package='crupR')) 
#' metaData <- data.frame(HM = c('H3K4me1','H3K4me3','H3K27ac'),
#'                     condition = c(2,2,2), replicate = c(1,1,1),
#'                     bamFile = files, inputFile = inputs)
#' norm <- readRDS(system.file('extdata', 'condition2_normalized.rds', package='crupR'))
#' S4Vectors::metadata(norm) <- metaData
#' #let's run the actual function
#' getEnhancers(data = norm)
#'
#' @export
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @import randomForest
#' @importFrom stats predict
#' @importFrom GenomicRanges mcols elementMetadata end reduce
#' @importFrom S4Vectors metadata
#' @importFrom BiocParallel SerialParam


getEnhancers <- function(data, classifier = NULL, all = FALSE) {
    start_time <- Sys.time()
    message("\n")

    if (!is.null(classifier) && (!dir.exists(classifier)))
        stop(paste0(classifier, " is not a valid directory"))
    if (all(colnames(GenomicRanges::mcols(data)) != c("H3K4me1", "H3K4me3", "H3K27ac",
        "ratio"))) {
        # EDIT: Adjusted error messages
        stop("Input must be an output of normalize() or contain counts for H3K4me1, H3K4me3, H3K27ac and the ratio of H3K4me1 and H3K4me3")
    }
    
    ################################################################## 
    # Read classifier files and ecdf
    ##################################################################

    message("Get classifier and empirical distribution function ...\n")  #EDIT: replaced cat() with message()

    if (is.null(classifier)) {
        classifierF1 <- system.file("extdata", "active_vs_inactive.rds", package = "crupR")
        classifierF2 <- system.file("extdata", "enhancer_vs_active_promoter.rds",
            package = "crupR")
    } else {
        classifierF1 <- file.path(classifier, "active_vs_inactive.rds")
        classifierF2 <- file.path(classifier, "enhancer_vs_active_promoter.rds")
    }

    check_file(classifierF1)
    check_file(classifierF2)
    classifier1 <- readRDS(classifierF1)
    classifier2 <- readRDS(classifierF2)

    feat1 <- unique(gsub("_.*", "", names(classifier1$forest$xlevels)))
    feat2 <- unique(gsub("_.*", "", names(classifier2$forest$xlevels)))
    featAll <- unique(c(feat1, feat2))

    ecdf_file <- system.file("extdata", "ecdf.rds", package = "crupR")
    check_file(ecdf_file)
    ecdf <- readRDS(ecdf_file)

    ################################################################## 
    # Quantile normalization
    ################################################################## 

    message(paste0("Quantile normalize counts for features ...\n"))  #EDIT: replaced cat() with message()

    # d <-data$D dNorm <- data$D
    d <- data
    dNorm <- data
    for (f in featAll) GenomicRanges::mcols(dNorm)[, f] <- preprocessCore::normalize.quantiles.use.target(matrix(GenomicRanges::mcols(dNorm)[,
        f]), get_targetQuantileNorm(ecdf[[f]]))

    ################################################################## 
    # Extend data matrix
    ##################################################################

    message("Extend data matrix ...\n")  #EDIT: replaced cat() with message()

    d_ext <- extend_dataMatrix(N = 5, df = data.frame(d), f = feat1)
    zero.idx <- which(rowSums(d_ext[, -c(seq_len(3))]) == 0)
    dNorm_ext <- extend_dataMatrix(N = 5, df = data.frame(dNorm), f = featAll, zero.idx = zero.idx)  #zero rows already removed
    # dNorm_ext <- dNorm_ext[-zero.idx,]

    ################################################################## 
    # make predictions
    ##################################################################
    #EDIT: replaced cat() with message()
    message("Predict enhancers for each 100 bp bin ...\n")  

    mid <- round(nrow(dNorm_ext)/2, 0)
    pred1 <- predict(classifier1, dNorm_ext[seq_len(mid), ], type = "prob")[, 2]
    pred1 <- c(pred1, predict(classifier1, dNorm_ext[(mid + 1):nrow(dNorm_ext), ],
        type = "prob")[, 2])
    pred2 <- predict(classifier2, dNorm_ext[seq_len(mid), ], type = "prob")[, 2]
    pred2 <- c(pred2, predict(classifier2, dNorm_ext[(mid + 1):nrow(dNorm_ext), ],
        type = "prob")[, 2])

    GenomicRanges::elementMetadata(d) <- NULL
    GenomicRanges::mcols(d)["prob"] <- 0
    GenomicRanges::mcols(d[-zero.idx])["prob"] <- pred1 * pred2
    if (all == TRUE) {
        GenomicRanges::mcols(d)["probA"] <- 0
        GenomicRanges::mcols(d)["probE"] <- 0
        GenomicRanges::mcols(d[-zero.idx])["probA"] <- pred1
        GenomicRanges::mcols(d[-zero.idx])["probE"] <- pred2
    }
    message(paste0("time: ", format(Sys.time() - start_time), "\n"))
    # EDIT: adjusted the output to a single GRanges object with meta data
    S4Vectors::metadata(d) <- S4Vectors::metadata(data)
    return(d)
}
