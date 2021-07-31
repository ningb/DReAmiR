
# diffreg
miRNA targets differential regulation analysis
<!-- badges: start -->
<!-- badges: end -->

Rewiring of transcriptional regulatory networks has been implicated in many biological and pathological processes. However, most current methods for detecting rewiring events (differential network connectivity) are not optimized for miRNA-mediated gene regulation and fail to systematically examine predicted target genes in study designs with multiple groups. We developed a novel method to address the current. The method first estimates miRNA-gene expression correlations with Spatial Quantile Normalization to remove the mean-correlation relationship. Then, for each miRNA, genes are ranked by their correlation strength per group. Enrichment patterns of predicted target genes are compared using the Anderson-Darling test and significance levels are estimated via permutation. Finally, graph embedding or difference in enrichment score maximization is performed to prioritize group-specific target genes. 

## Installation

`diffreg` is still under active development.

``` r
library(devtools)
install_github("https://github.com/ningb/diffreg")
```

## Example

To run the main function:

``` r
library(diffreg)

# Read-in simulated data
target.mat <- sim.data[[1]]
group.label <- sim.data[[2]]
mrna.mat <- sim.data[[3]]
mirna.mat <- sim.data[[4]]

# Build correlation matrices by subtype
cor.list <- BuildCorList(
	mrna.mat    = mrna.mat,
	mirna.mat   = mirna.mat,
	group.label = group.label
	)

# Build target matrix	
target.mat <- BuildTargetMat(
	cor.list  = cor.list,
	Pair.df   = "TargetScan",
	gene.name = "hgnc_symbol"
	)

results <- Run_diffreg(
	cor.list          = cor.list,
	target.mat        = target.mat,
	min.target.number = 5,
	direction         = "Negative"
	)
```

