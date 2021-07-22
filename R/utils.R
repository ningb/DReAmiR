#' Generate a list of correlation matrices for later analysis
#' 
#' @param mrna.mat A matrix of mRNA expression (genes at rows and samples at columns)
#' @param mirna.mat A matrix of miRNA expression (genes at rows and samples at columns)
#' @param raw.mrna.mat A raw matrix of mRNA expression (genes at rows and samples at columns)
#' @param raw.mirna.mat A raw matrix of miRNA expression (genes at rows and samples at columns)
#' @param group.label A factor vector subtype label for the samples
#' @param cor.type A character string indicating type of correlation to do; options are "pearson", "kendall", "spearman"
#' @param do.libnorm A logic type indicating whether library size normalization (logCPM) is needed
#' @param do.spqn A logic type indicating whether mean-correlation relationship is needed to be removed
#' 
#' @export
BuildCorList <- function(mrna.mat, mirna.mat, group.label) {

	# Check dimension

	# Check raw vs. normalized matrices size

	# Check format

	# Form correlation matrix per group

	# Subset to a miRNA-mRNA correlation matrix

	# Attach name to the list

	names(cor.list) <- group.label
	
	return(cor.list)
}


#' Generate a list of correlation matrices for later analysis
BuildCorList <- function(mrna.mat, mirna.mat, group.label) {
	expr.mat <- rbind(mrna.mat, mirna.mat)
	cor.list <- lapply(levels(group.label), function(x) cor(t(expr.mat[, group.label == x])))
	cor.list <- lapply(cor.list, function(x) x[1:nrow(mrna.mat), (nrow(mrna.mat)+1):ncol(x)])
	names(cor.list) <- levels(group.label)
	return(cor.list)
}

#' Generate a list of correlation matrices for later analysis
#' 
#' @param mrna.mat A matrix of mRNA expression (genes at rows and samples at columns)
#' @param mirna.mat A matrix of miRNA expression (genes at rows and samples at columns)
#' @return A binary target matrix with mRNA on rows and miRNA on cols
#' 
#' @export
BuildTargetMat <- function(cor.list, Pair.df=c(), gene.name=c("hgnc_symbol", "ensembl_ID")) {}

#' Generate a list of LE genes to be prioritized
#' 
#' @param cor.list A list of correlation matrices with miRNA on columns and mRNA on rows
#' @param miR miRNA to plot
#' @param target.mat A binary matrix indicating target with miRNA on columns and mRNA on rows
#' @param direction A character string specifying the direction of enrichment, must be "Negative", "Positive" or "Both"
#' @param alpha Weight for running AD/KS test. Default is 0.1
#' @param with.weight Logic of whether to return weights as well
#' @return An list of leading edge genes for each group
#' 
#' @export
GetLE <- function(cor.list, target.mat, miR, direction=c("Negative", "Positive"), alpha=0.1, with.weight=TRUE) {

	if (!length(direction) == 1) {
		stop("There should be only one 'direction'.\n")
	}

	if (!direction %in% c("Negative", "Positive", "Both")) {
		stop("'direction' must be either Negative or Positive.\n")
	}

	target.set <- names(target.mat[,miR][target.mat[,miR]!=0])

	rank.list <- lapply(cor.list, function(x) sort(x[,miR]))

	ws.list <- lapply(rank.list, function(x) .CalcGseaStat(x, target.set, alpha))

	if (direction == "Negative") {
		le.gene <- lapply(ws.list, function(x) names(x[1:which(x == max(x))]))
	} else if (direction == "Positive") {
		le.gene <- lapply(ws.list, function(x) names(x[which(x == min(x)):length(x)]))
	} else {
		le.gene <- lapply(ws.list, function(x) c(names(x[1:which(x == max(x))]),
												names(x[which(x == min(x)):length(x)])))
	}

	le.gene <- lapply(le.gene, function(x) x[x %in% target.set])

	if (with.weight == TRUE) {
		le.gene <- sapply(names(le.gene), function(x) rank.list[[x]][le.gene[[x]]])
	}

	return(le.gene)
}

#' Generate a vector of colors to use
#' 
#' @param n Number of levels
#' @param set Pallete to choose
#' 
#' @import RColorBrewer
default.pal <- function(n, set) {
	getPalette <- colorRampPalette(brewer.pal(11, set))
	pal <- getPalette(n)
	return(pal)
}

