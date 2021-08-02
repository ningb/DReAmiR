#' Generate a list of correlation matrices for later analysis
#' 
#' @param mrna.mat A matrix of mRNA expression (genes at rows and samples at columns)
#' @param mirna.mat A matrix of miRNA expression (genes at rows and samples at columns)
#' @param group.label A factor vector subtype label for the samples
#' @param numCores An integer indicating the number of cores to use
#' 
#' @import parallel
#' @export
BuildCorList <- function(mrna.mat, mirna.mat, group.label, numCores=1) {

	# Check format of group.label: factor variable; if not, turn in to a factor
	if(class(group.label) != "factor") {
		group.label <- as.factor(group.label)
		warning("Changed class of group.label to factor.")
	} 

	# Check if names exist: sample names at columns of miRNA/mRNA matrices; miRNA and mRNA names at rows; names for group.label
	if(is.null(rownames(mrna.mat)) | is.null(rownames(mirna.mat))) {
		stop("Input misses gene or miRNA names")
	} 
	if(is.null(colnames(mrna.mat)) | is.null(colnames(mirna.mat)) | is.null(names(group.label))) {
		stop("Inputs misses sample names\n")
	} 

	# Select common samples
	common.samples <- Reduce(intersect, list(names(group.label), colnames(mrna.mat), colnames(mirna.mat)))

	# Check if duplicated sample name exists
	if(any(duplicated(common.samples))) {
		stop("Duplicated sample name exists. Please fix before running this function.\n")
	}

	# Subset data to only the overlapped samples
	if(length(common.samples) != length(group.label) | length(common.samples) != ncol(mrna.mat) | length(common.samples) != ncol(mirna.mat)) {
		warning("Detected unequal samples sizes. Only the overlapped sampels are used.\n")
		group.label <- group.label[common.samples]
		mrna.mat <- mrna.mat[, common.samples]
		mirna.mat <- mirna.mat[, common.samples]
	}

	# Check if samples sizes meet requirement
	if (min(table(group.label)) <= 30) {
		warning("One group has sample size smaller than 30. Correlation might not work well.\n")
	}

	# Form correlation matrix per group
	expr.mat <- rbind(mrna.mat, mirna.mat)
	cor.list <- mclapply(levels(group.label), 
						function(x) cor(t(expr.mat[, group.label == x])), 
						mc.cores = numCores)

	# Subset to a miRNA-mRNA correlation matrix
	cor.list <- mclapply(cor.list, 
						function(x) x[row.names(mrna.mat), row.names(mirna.mat)], 
						mc.cores = numCores)

	# Attach name to the list
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

