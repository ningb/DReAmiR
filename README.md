
# diffreg
miRNA targets differential regulation analysis
<!-- badges: start -->
<!-- badges: end -->

## Introduction
Rewiring of transcriptional regulatory networks has been implicated in many biological and pathological processes. However, most current methods for detecting rewiring events (differential network connectivity) are not optimized for miRNA-mediated gene regulation and fail to systematically examine predicted target genes in study designs with multiple groups. We developed a novel method to address the current. The method first estimates miRNA-gene expression correlations with Spatial Quantile Normalization to remove the mean-correlation relationship. Then, for each miRNA, genes are ranked by their correlation strength per group. Enrichment patterns of predicted target genes are compared using the Anderson-Darling test and significance levels are estimated via permutation. Finally, graph embedding or difference in enrichment score maximization is performed to prioritize group-specific target genes. 

<img src="media/Workflow.png" height="400px" class="center"/>

## Requirements
R (>= 4.0.0) is recommended.

## Installation

`diffreg` is still under active development.

``` r
library(devtools)
install_github("https://github.com/ningb/diffreg")
```

## Usage
Here is the basic workflow for `diffreg`. 
### Load package and read in toy data:
``` r
library(diffreg)

# Read-in simulated data
target.mat <- sim.data[[1]]
group.label <- sim.data[[2]]
mrna.mat <- sim.data[[3]]
mirna.mat <- sim.data[[4]]
```
### Build correlation matrices by subtype and extract target matrix
```r
cor.list <- BuildCorList(
	mrna.mat    = mrna.mat,
	mirna.mat   = mirna.mat,
	group.label = group.label
	)
# > lapply(cor.list, dim)
# $A
# [1] 1000   20

# $B
# [1] 1000   20

# $C
# [1] 1000   20

# Build target matrix	
target.mat <- BuildTargetMat(
	cor.list  = cor.list,
	Pair.df   = "TargetScan",
	gene.name = "hgnc_symbol"
	)
```
### Optional: run spatial quantile normalization on the correlation matrices
*Note: this is not going to work in this toy dat since it is too small*
```r
cor.list <- RunGroupSpQN(mrna.mat=mrna.mat, 
	mirna.mat=mirna.mat, 
	group.label=group.label, 
	cor.list=cor.list)
```
### Run the main function for edifferential regulation
```r
results <- Run_diffreg(
	cor.list          = cor.list,
	target.mat        = target.mat,
	min.target.number = 5,
	direction         = "Negative"
	)
# > results
#           ADval       Pval   A_Score   B_Score   C_Score
# miR_1   3.63230 8.1650e-03 0.7370506 0.6812772 0.4733412
# miR_2  10.92700 1.3873e-06 0.6731506 0.6729634 0.3102756
# miR_3  10.68900 1.9488e-06 0.6142691 0.6680506 0.2988116
# miR_4   8.22520 5.0877e-05 0.7027055 0.6946374 0.3628987
# miR_5   5.89900 6.9189e-04 0.6278100 0.7477924 0.4162563
# miR_6  13.09300 6.2400e-08 0.7722832 0.7524360 0.3041277
# miR_7   9.58270 9.5025e-06 0.7102737 0.7108067 0.3399027
# miR_8   4.73000 2.4806e-03 0.6721830 0.6941325 0.4255792
# miR_9   3.73600 7.3135e-03 0.6437979 0.6235318 0.4303972
# miR_10  7.30570 1.4629e-04 0.7087622 0.6415018 0.3662423
# miR_11 -1.23240 9.7939e-01 0.3690861 0.3762153 0.3422135
# miR_12  0.61642 2.0732e-01 0.2829414 0.4710669 0.4591127
# miR_13 -0.89984 8.5619e-01 0.2791646 0.3617085 0.3357241
# miR_14 -1.34350 9.9431e-01 0.4527785 0.4176795 0.4294675
# miR_15 -0.83689 8.2330e-01 0.3488505 0.4364934 0.3917992
# miR_16 -1.08150 9.3690e-01 0.4280941 0.3581289 0.4050904
# miR_17 -0.79527 8.0074e-01 0.4038234 0.4294814 0.3630465
# miR_18 -0.90107 8.5682e-01 0.3873749 0.4043080 0.4672643
# miR_19 -1.16770 9.6434e-01 0.4161954 0.4236637 0.3915652
# miR_20 -0.95749 8.8457e-01 0.4089479 0.3164258 0.4239376
```
### Plot an enrichment plot for an significant miRNA
```r
Plot_enrichment(cor.list=cor.list, 
	miR="miR_1", 
	target.mat=target.mat)
```
<img src="media/Diff_reg_example.png" height="400px" class="center"/>

###	Extract leading-edge gene and prioritize targets using either Graph Representation Learning or dES maximization