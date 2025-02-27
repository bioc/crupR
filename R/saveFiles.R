#' Function to save the different output objects of each step
#' @description
#' This functions provides an easy way to save the outcomes of every step in a suitable format. 
#' The type of format available depends on the outcome itself.
#' 
#' @param data An output object of getEnhancers(), getSE(), getDynamics() or getTargets()
#' @param modes Formats in which the GRanges object should be saved. Following modes are available:
#' for the output of getEnhancers(): 'bigWig' for a bigWig file, 'rds' for an .rds file
#' for the output of getDynamics(): 'beds' for saving each cluster in a seperate BED file
#' for the output of getTargets(): 'UCSC' for a UCSC interaction file
#' for the output of getSE(): 'bedGraph' for saving the peak calls in a 
#' bedGraph file, 'bed' for saving the clusters of peaks in a BED file
#' @param outdir Output directory in which the files should be saved
#' @param nearest Only relevant, if you want to save the output of enhancerTargets. 
#' Specifies if the output was produced by the nearest gene mode of the function or not. 
#' Default is FALSE.
#' @return True
#' @examples
#' \donttest{
#' # output directory
#' example_path <- file.path(tempdir(), 'crupR') #let's use a temporary direcotry for the outputs
#' dir.create(example_path) #create the directory
#' # recreate the output of getDynamics() to save
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
#' clusters <- readRDS(system.file('extdata', 'differential_enhancers.rds', 
#' package='crupR'))
#' S4Vectors::metadata(clusters) <- metaData
#' # save enhancer clusters as bed files
#' saveFiles(data = clusters, modes = 'beds', outdir = example_path)}
#' @export
#' @importFrom rtracklayer export
#' @importFrom GenomicRanges mcols seqnames
#' @importFrom utils write.table
#' @importFrom fs path_norm


saveFiles <- function(data, modes, outdir, nearest = FALSE) {
    possible.modes <- c("rds", "bigWig", "bed", "bedGraph", "beds", "UCSC")
    modes.pred <- c("bigWig", "bedGraph", "rds", "bed")
    if (!dir.exists(outdir))
        stop(paste0("Directory ", outdir, " doesn't exist!"))
    outdir <- paste0(fs::path_norm(outdir), "/")  #EDIT: added slash so files are saved correctly
    if (!all(modes %in% possible.modes)) {
        stop("One of the modes you chose is not supported by this function.\n Please choose from rds, bigWig, bed, bedGraph, beds and UCSC.")
    }
    # check if the mode and the data format are matching
    if ("UCSC" %in% modes) {
        if (is.list(data))
            stop("Only outputs of getTargets() can be saved in the mode UCSC.")
        # if not getTargets output stop()
        if (!("cluster" %in% colnames(GenomicRanges::mcols(data))) | !(any(grepl("GENE",
            colnames(mcols(data)))))) {
            stop("Only outputs of getTargets() can be saved in the mode UCSC.")
        }
    }
    if ("beds" %in% modes) {
        if (is.list(data))
            stop("Only outputs of getDynamics() can be saved in the mode beds.")
        # if not getDynamics output stop()
        if (!("cluster" %in% colnames(GenomicRanges::mcols(data))) | (any(grepl("GENE",
            colnames(mcols(data)))))) {
            stop("Only outputs of getDynamics() can be saved in the mode beds.")
        }

    }
    if (("bed" %in% modes) || ("bedGraph" %in% modes)) {
        # if not getSE output stop()
        if (!is.list(data)) {
            stop("Only outputs of getSE() can be saved in the modes bedGraph and/or bed.")
        } else {
            if ((length(data) < 3) || any(names(data) != c("D", "peaks", "cluster"))) {
                stop("Only outputs of getSE() can be saved in the modes bedGraph and/or bed.")
            }
        }
    }
    if (("rds" %in% modes) || ("bigWig" %in% modes)) {
        # if not getEnhancers or getSE output stop()
        if (is.list(data)) {
            # test if SE output
            if (any(names(data) != c("D", "peaks", "cluster"))) {
                stop("Only outputs of getEnhancers or getSE() can be saved in the modes rds and/or bigWig.")
            } else {
                if (!("prob" %in% colnames(GenomicRanges::mcols(data$D)))) {
                  stop("Only outputs of getEnhancers or getSE() can be saved in the modes rds and/or bigWig.")
                }
            }
        } else {
            # test if getEnhancers output
            if (!("prob" %in% colnames(GenomicRanges::mcols(data)))) {
                stop("Only outputs of getEnhancers or getSE() can be saved in the modes rds and/or bigWig.")
            }
        }
    }

    # for: enhancerPrediction
    if ("bedGraph" %in% modes) {
        out.bedgraph <- paste0(outdir, "singleEnh.bedGraph")
        rtracklayer::export(data$peaks, out.bedgraph)
    }
    if ("bed" %in% modes) {
        out.bed <- paste0(outdir, "clusterEnh.bed")
        rtracklayer::export(data$cluster, out.bed)
    }
    if ("bigWig" %in% modes) {
        # EDIT: data is now a GRanges and not a list anymore
        if (is.list(data)) {
            data.D <- data$D
        } else {
            data.D <- data
        }
        if (length(GenomicRanges::mcols(data.D)) > 1) {
            GenomicRanges::mcols(data.D) <- NULL
            metaCols <- GenomicRanges::mcols(data)
            # for the normal predictios
            out.bw <- paste0(outdir, "prediction.bw")
            data.D$score <- metaCols$prob
            rtracklayer::export(data.D, out.bw)
            # for probA
            out.bw <- paste0(outdir, "prediction_probA.bw")
            data.D$score <- metaCols$probA
            rtracklayer::export(data.D, out.bw)
            # for probE
            out.bw <- paste0(outdir, "prediction_probE.bw")
            data.D$score <- metaCols$probE
            rtracklayer::export(data.D, out.bw)
        } else {
            colnames(GenomicRanges::mcols(data.D)) <- c("score")
            out.bw <- paste0(outdir, "prediction.bw")
            rtracklayer::export(data.D, out.bw)
        }
    }
    if ("rds" %in% modes) {
        if (is.list(data)) {
            data.D <- data$D
        } else {
            data.D <- data
        }
        out.rds <- paste0(outdir, "prediction.rds")
        saveRDS(data.D, out.rds)  #EDIT: data is now a GRanges and not a list anymore
    }


    # for: enhancerDynamics
    if ("beds" %in% modes) {
        peaks <- data  #EDIT: data is now a GRanges and not a list anymore
        clusters <- unique(peaks$cluster)

        for (c in clusters) {
            if (!grepl("r", c)) {
                out.bed <- paste0(outdir, paste0("dynamicEnh__cluster_", c, ".bed"))
                write.table(data.frame(peaks)[which(peaks$cluster == c), GR_header_short],
                  file = out.bed, quote = FALSE, row.names = FALSE, col.names = FALSE,
                  sep = "\t")
            }
        }
    }

    # for: enhancerTargets
    if ("UCSC" %in% modes) {
        out.interaction <- ""
        units <- data  #EDIT: data is now a GRanges and not a list anymore
        header <- "track type=interact name=\"Dynamic Promoter-Enhancer Pairs\" description=\"Dynamic Promoter-Enhancer Pairs\" interactDirectional=true visibility=full"
        enhancer <- as.matrix(data.frame(units)[, c("start", "end")])
        promoter <- c()
        score <- 0
        if (nearest == FALSE) {
            out.interaction <- paste0(outdir, paste0("RegulatoryUnits.interaction"))
            promoter <- as.matrix(data.frame(GenomicRanges::mcols(units)$CORRELATED_GENE_PROMOTER_START,
                GenomicRanges::mcols(units)$CORRELATED_GENE_PROMOTER_END))
            score <- GenomicRanges::mcols(units)$CORRELATION
        } else {
            out.interaction <- paste0(outdir, paste0("RegulatoryUnitsNearestGene.interaction"))
            promoter <- as.matrix(data.frame(GenomicRanges::mcols(units)$NEAREST_GENE_PROMOTER_START,
                GenomicRanges::mcols(units)$NEAREST_GENE_PROMOTER_END))
            score <- GenomicRanges::mcols(units)$DISTANCE_TO_NEAREST
        }

        interaction <- cbind(as.character(GenomicRanges::seqnames(units)), apply(cbind(enhancer,
            promoter), 1, min), apply(cbind(enhancer, promoter), 1, max), rep(".",
            length(units)), rep(0, length(units)), score, rep(".", length(units)),
            rep("#7A67EE", length(units)), as.character(GenomicRanges::seqnames(units)),
            enhancer, rep(".", length(units)), rep(".", length(units)), as.character(GenomicRanges::seqnames(units)),
            promoter, rep(".", length(units)), rep(".", length(units)))

        write.table(header, file = out.interaction, quote = FALSE, col.names = FALSE,
            row.names = FALSE, sep = "\t")
        write.table(interaction, file = out.interaction, quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t", append = TRUE)
    }
    return(TRUE)
}
