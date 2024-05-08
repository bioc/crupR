# **crupR - Short Description**   

crupR is improved version of the enhancer detection pipeline CRUP ([C]ondition-specific [R]egulatory [U]nits [P]rediction). Its workflow consists of the same three main steps as the original pipeline (enhancer prediction, condition-specific enhancer detection, target gene detection) and an additional pre-preapring step (normalization). It uses different layers of epigenetic information in the form of ChIP-seq data from histone modification to provide a comprehensive list of regulatory units
consisting of dynamically changing enhancers and target genes.

#### *Installation*
Install crupR using following commad:

```
devtools::install_git("https://github.com/akbariomgba/crupR")
```


#### *Contact*

omgba@molgen.mpg.de

#### *Citation*

Since the crupR application note has not been submitted yet, you can cite this repository and the original paper if you use crupR:

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1860-7



## *crupR Overview*

*crupR* is a pipeline which aims to efficiently combine all its functions into a simple workflow. The main steps are:

* `normalize()`: Normalizes the ChIPseq counts
* `getPrediction()`: Uses the normalized ChIpseq counts to predict the occurence of active enhancers by applying a random forest based classifier
* `getDynamics()`: Uses the predicted probabilities to find dynamically changing, condition-specific enhancer clusters
* `getTargets(()`: Uses RNAseq experiments to find possible target genes for the dynamic enhancers.

Beyond that, *crupR* offers an additional function to summarize the single enhancers into super enhancers (`getSE()`), a function to visualise the condition-specific enhancer clusters (`plotSummary()`) and a function to export all the files produced by each *crupR* step as appropriate formats (i.e. bed, BigWig or bedGraph).

For a more detailed overview of *crupR* with more code examples, please check the vignette.

### *Getting started*

#### What will I need?

In order to run *crupR*, you need ChIPseq experiments for the histone modifications H3K27ac, H3K4me1 and H3K4me3. Additionally, it is recommended to also inclue an input experiment for those experiments, however if these aren't available for some reason, you can run *crupR* without them. All the experiments need to be already aligned and indexed.

#### Prepare the metaData file

After you prepared your ChIPseq files, you need to prepare the meta data frame. This data frame contains the paths to all the necessary bam files. It also contains further information about them, such as condition, replicate, histone modification and input file. 

Let's look at the metaData frame for the example files.
First, we'll need the paths to the bam files and input files.

```{r}
files <- c(system.file("extdata", "Condition1.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition1.H3K27ac.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me1.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K4me3.bam", package="crupR"),
           system.file("extdata", "Condition2.H3K27ac.bam", package="crupR"))
                     
inputs <- c(rep(system.file("extdata", "Condition1.Input.bam", package="crupR"), 3), rep(system.file("extdata", "Condition2.Input.bam", package="crupR"), 3))
```

Now we can build the meta data frame by adding the information about the hinstone modification, condition and replicate for each file. As there exists only one sample per condition, all of the files get replicate "1".

```{r}
metaData <- data.frame(HM = rep(c("H3K4me1","H3K4me3","H3K27ac"),2),
                       condition = c(1,1,1,2,2,2), replicate = c(1,1,1,1,1,1),
                       bamFile = files, inputFile = inputs)

metaData
```
After creating the meta data frame, you are ready to run *crupR*.


### *Run the crupR pipeline*

#### *Step 0: Normalize the ChIPseq counts*

*crupR* needs normalized ChIPseq counts for its prediction, so we offer an additional function to normalize the ChIPseq files beforehand. You need to normalize the ChIP-seq experiments for every replicate of every condition. In our example files two conditions with one replicate each are present. 
Thus, we must run `normalize()`twice. Please note, that *crupR* can only normalize BAM files for the genome builts mm9, mm10, hg19 and hg38 or for a costum genome (which needs to be a Seqinfo object).  

```{r}
normalized_1 <- crupR::normalize(metaData = metaData, condition = 1, replicate = 1, genome = "mm10", 
                                sequencing = "paired", C = 2)
normalized_2 <- crupR::normalize(metaData = metaData, condition = 2, replicate = 1, genome = "mm10", 
                                sequencing = "paired", C = 2)
```
If an input experiment is not availabe, *crupR* also offers and input-free mode (check the vignette or documentation of `normalize()` for more information).
The output of `normalize()`  should be a list of length 2 containing the meta data of the samples that were normalized and a GRanges object containing the normalized ChIPseq counts for the binned genome of your choice. 

#### *Step 1: Predict active enhancers with crupR*

Now, the normalized ChIP-seq counts can be used to predict the enhancer activity on the genome of your choice.

```{r}
prediction_1 <- crupR::getEnhancers(data = normalized_1, C = 2)
prediction_2 <- crupR::getEnhancers(data = normalized_2, C = 2)
```
The output of this step is a list containing the truncated meta data file containing the information about the respective condition and replicate and a GRanges file containing the the binned genome with the enhancer prediction values for each bin. 

#### *Step 1.5: Find enhancer peaks and super enhancers with crupR*

*crupR* offers the `getSE()` function as an additional step after the enhancer prediction. Enhancers that are in close proximity are summarized into super enhancers during this step.
`getSE()` uses the output of the last step (enhancerPrediction), meaning the list with the meta data and the GRanges object, as input. 
```{r}
se <- crupR::getSE(data = prediction_2, C = 2)
```
The output of this step is a list containing the same meta data file as the input, the same GRanges object that was produced by `getPrediction()`, a GRanges object containing the single enhancer peak calls (can be saved as a .bedGraph file) and a GRanges object containing the super enhancers or rather peak clusters (can be saved as a bed file).

This step is not necessary for the next ones and can be skipped if one isn't interested in the peak calls or peak clusters.

#### *Step 2: Find conditon-specific enhancer clusters*

In the second workflow step, the *crupR* function `getDynamics()` defines dynamic enhancer regions by applying a Kolmogorov-Smirnov test directly on the enhancer probabilities in a pairwise manner.

Before running the actual function, the output objects from the prior steps need to be put in one list.
```{r}
predictions <- list(prediction_1, prediction_2)
```
Now everything is ready for `enhancerDynamics()` to run.
```{r}
clusters <- crupR::getDynamics(data = predictions, C = 2)
```
The output of this step is a list containing the complete meta data file and the condition-specific clusters in a GRanges object.

##### *Visualization of the clusters using `plotSummary()`*

As trying to decode the patterns for each cluster can be a bit complicated and also take some time, *crupR* offers the visualization function `plotSummary()` for the detected enhancer clusters. The function shows the boxplots of the median probabilities of the enhancers. This is a simple way of checking which clusters are active in which conditions. 

```{r}
crupR::plotSummary(clusters)
```

#### *Step 3: Find target genes of the condition-specific enhancer clusters*

The last step in the pipeline is the identification of the target genes that get regulated by the condition-specific enhancer clusters found in the previous step.


##### *Gene expression counts*
For the function `getTargets()`the output of the step before is not sufficient. *crupR* requires additionally the gene expression counts for each condition and replicate all in one GRanges object. 
*crupR* offers such an expression file for the the example sets.

```{r}
expression <- readRDS(file = system.file("extdata", "expressions.rds", package="crupR"))
```

However, there is no in-built function that calculates the gene expression counts from RNAseq experiments. That means the user has to generate this file by themselves. Please make sure that the column names for the columns that contain the gene expression counts of each sample are the same as the column names of the columns that contain the enhancer probabilities of each sample, i.e. "cond1_1", "cond1_2", "cond2_1" and so on. 

##### *Running `getTargets()`*

Looking for the enhancer targets requires - besides the gene expression counts - a .bed file containing TADs for the genome of your choice. In case mm10 is used, *crupR* provides a suitable .bed file which is also used per default for mm10 data. In this case, simply running the code below is sufficient:

```{r}
targets <- crupR::getTargets(data=clusters, expr = expression, genome = "mm10", C = 2)
```

In case another genome is used or a different BED file is preferred, providing the path to the BED file is necessary.

#### *Exporting the files*
After running the *crupR* pipeline, you might want to save the resulting GRanges objects. *crupR* provides the function `saveFiles()`that exports the different files in suitable formats. 

In total, there are X possible formats: "bigWig", "bedGraph", "rds", "bed", "beds" and "UCSC". However, the format to export your file to depends on the *crupR* step:

*Output of `getEnhancers()`: The GRanges file that is produced in this step can be saved as bigWig file or an .rds file. Thus, the options "bigWig" and/or "rds" can be chosen.
*Output of `getSE()`: This step produces two new GRanges files: single enhancer peak calls () and peak clusters/super enhancers (). The single enhancer peak calls can be export as a bedGraph file and the peak clusters can be exported as a BED file.
*Output of `getDynamics()`: The GRanges file produced in this step can be exported as multiple BED files ("beds"). Each bed file contains the dynamic enhancer regions of one cluster.
*Output of `getTargets()`: The GRanges file produced in this step can be exported in (UCSC) interaction format ("UCSC").

If you want to export the files, you just have to use the output list of the respetive step as input and specify your desired format(s):

```{r}
out_dir <- ""
#save the GRanges object of the getEnhancers() step
saveFiles(data = prediction_1, modes = c("bigWig", "rds"), outdir = out_dir)
#save the GRanges object of the getSE() step
saveFiles(data = se, modes = c("bedGraph", "bed"), outdir = out_dir)
#save the GRanges object of the getDynamics() step
saveFiles(data = clusters, modes = "beds", outdir = out_dir)
#save the GRanges object of the getTargets() step
saveFiles(data = targets, modes = "UCSC", outdir = out_dir)
```
