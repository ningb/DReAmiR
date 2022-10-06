#' Calculating random walking score
#' 
#' @param ranklist A named vector for all genes sorted by correlation
#' @param gs A character vector for the set of targets of interest
#' @param alpha A numeric variable for gsea parameter
#' @return A walking score vector
.CalcGseaStat <- function(ranklist, gs, alpha=NULL) {

	if (is.null(alpha)) {
		alpha <- 1
	}

	id.inside <- ifelse(names(ranklist) %in% gs, 1, 0)

	CDF.in <- cumsum((abs(ranklist)^alpha * id.inside)) /
	                  sum((abs(ranklist)^alpha * id.inside))

	CDF.notin <- cumsum(!id.inside) / sum(!id.inside)

	ws <- CDF.in - CDF.notin

	return(ws)
}

#' Run multiple-group differential enrichment using AD test
#' 
#' @param cor.list A list of correlation matrices with miRNA on columns and mRNA on rows
#' @param target.mat A binary matrix indicating target with miRNA on columns and mRNA on rows
#' @param min.target.number A integer indicating minimum number of targets requried for doing diffreg
#' @param direction A character string specifying the direction of enrichment, must be "Negative" or "Positive"
#' @param alpha Weight for running AD/KS test
#' @return A data frame with ES and test results for each miRNA
#' 
#' @importFrom kSamples ad.test
.diffreg_ad <- function(target.mat, cor.list, min.target.number=5, direction="Negative", alpha=NULL) {

	# Format the outputs first
	ngroup <- length(cor.list)
	ad.out.ncol <- ngroup + 3

	ad.out <- matrix(nrow=0, ncol=ad.out.ncol)

	for (i in colnames(target.mat)) {

		# Check if we have enough targets; if not, we will skip this miRNA
		if (colSums(target.mat)[i] <= min.target.number) {
			next

		} else {

			target.set <- names(target.mat[,i][target.mat[,i]!=0])

			rank.list <- lapply(cor.list, function(x) sort(x[,i]))

			pos.list <- lapply(rank.list, function(x) match(target.set, names(x)))

			ad.res <- ad.test(pos.list, method="asymptotic", dist=F)

			ws.list <- lapply(rank.list, function(x) .CalcGseaStat(x, target.set, alpha))

			es.max <- unlist(lapply(ws.list, max))
			es.min <- unlist(lapply(ws.list, min))

			# Check if enrichment direction is the same as input
			# Since we start ranking from the most negative side, positive ES score means negative enrichment
			es.direction <- sign(mean(es.max) + mean(es.min))
			test.direction <- ifelse(es.direction == -1, "Positive", "Negative")

			if (test.direction != direction) {
				warning(paste0("Enrichment patterns may not match with expectation for ", i, "\n", sep=""))
			}

			# Organize output		
			if (direction == "Negative") {
				es <- es.max
			} else {
				es <- es.min
			}			

			out <- c(i, ad.res$ad[1,2], ad.res$ad[1,3], es)
			ad.out <- rbind(ad.out, out)

		}
	}

	result <- as.data.frame(ad.out[, 2:ncol(ad.out)])
	result <- data.frame(sapply(result, function(x) as.numeric(x)))
	row.names(result) <- ad.out[,1]

	names(result) <- c("ADval", "Pval",
						paste(names(cor.list), "_Score", sep=""))

	return(result)

}

#' Run two-group differential enrichment using KS test
#' #' 
#' @param cor.list A list of correlation matrices with miRNA on columns and mRNA on rows
#' @param target.mat A binary matrix indicating target with miRNA on columns and mRNA on rows
#' @param min.target.number A integer indicating minimum number of targets requried for doing diffreg
#' @param direction A character string specifying the direction of enrichment, must be "Negative" or "Positive"
#' @param alpha Weight for running AD/KS test
#' @return A data frame with ES and test results for each miRNA
.diffreg_ks <- function(target.mat, cor.list, min.target.number=5, direction="Negative", alpha=NULL) {

	# Format the outputs first
	ks.out <- matrix(nrow=0, ncol=5)

	for (i in colnames(target.mat)) {

		# Check if we have enough targets; if not, we will skip this miRNA
		if (colSums(target.mat)[i] <= min.target.number) {
			next

		} else {

			target.set <- names(target.mat[,i][target.mat[,i]!=0])

			rank.list <- lapply(cor.list, function(x) sort(x[,i]))

			pos.list <- lapply(rank.list, function(x) match(target.set, names(x)))

			ks.res <- suppressWarnings(ks.test(pos.list[[1]], pos.list[[2]]))

			ws.list <- lapply(rank.list, function(x) .CalcGseaStat(x, target.set, alpha))

			es.max <- unlist(lapply(ws.list, max))
			es.min <- unlist(lapply(ws.list, min))

			# Check if enrichment direction is the same as input
			# Since we start ranking from the most negative side, positive ES score means negative enrichment
			es.direction <- sign(mean(es.max) + mean(es.min))
			test.direction <- ifelse(es.direction == -1, "Positive", "Negative")

			if (test.direction != direction) {
				warning(paste0("Enrichment patterns do not match with expectation for ", i, "\n", sep=""))
			}

			# Organize output		
			if (direction == "Negative") {
				es <- es.max
			} else {
				es <- es.min
			}			

			out <- c(i, ks.res$statistic, ks.res$p.value, es)
			ks.out <- rbind(ks.out, out)

		}
	}

	result <- as.data.frame(ks.out[, 2:ncol(ks.out)])
	result <- data.frame(sapply(result, function(x) as.numeric(x)))
	row.names(result) <- ks.out[,1]

	names(result) <- c("KSval", "Pval",
						paste(names(cor.list), "_Score", sep=""))

	return(result)

}

#' Given a list of correlation matrcies, return the SpQN normalized matrices
#' 
#' @param cor.list A list of correlation matrices with miRNA on columns and mRNA on rows
#' @param target.mat A binary matrix indicating target with miRNA on columns and mRNA on rows
#' @param min.target.number A integer indicating minimum number of targets requried for doing diffreg
#' @param direction A character string specifying the direction of enrichment, must be "Negative" or "Positive"
#' @param alpha Weight for running AD/KS test
#' @return A data frame with ES and test results for each miRNA
#' 
#' @importFrom kSamples ad.test
#' @export
#' 
Run_diffreg <- function(cor.list, target.mat, min.target.number=5, direction="Negative", alpha=1) {

	if (!length(direction) == 1) {
		stop("There should be only one 'direction'.\n")
	}

	if (!direction %in% c("Negative", "Positive")) {
		stop("'direction' must be either Negative or Positive.\n")
	}

	# Check whether target and correlation matrix match
	if (!(nrow(target.mat) == nrow(cor.list[[1]]) & ncol(target.mat) == ncol(cor.list[[1]]))) {
		stop("Dimensions of target matrix and correlation matrices do not match..\n")
	}
	
	# Check whether mRNA and miRNA names match between target.mat and cor.list
	if (!(identical(row.names(target.mat), row.names(cor.list[[1]])))) {
		stop("mRNA names in target matrix do not match with correlation matrices..\n")
	}

	if (!(identical(colnames(target.mat), colnames(cor.list[[1]])))) {
		stop("mRNA names in target matrix do not match with correlation matrices..\n")
	}

	# Check number of levels
	if (length(cor.list) < 2) {

	} else if (length(cor.list) == 2) {
		message("Comparing two-group diff-reg using KS test\n")
		result <- .diffreg_ks(target.mat, cor.list, min.target.number, direction, alpha)
	} else {
		message("Comparing multiple-group diff-reg using AD test\n")
		result <- .diffreg_ad(target.mat, cor.list, min.target.number, direction, alpha)
	}

	return(result)
}

#' Given a diff-reg results, do permutation on group label 
#' 
#' @param diffreg.res Differential regulation results from Run_diffreg()
#' @param mrna.mat A matrix of mRNA expression (genes at rows and samples at columns)
#' @param mirna.mat A matrix of miRNA expression (genes at rows and samples at columns)
#' @param raw.mrna.mat (optional) A raw matrix of mRNA expression (genes at rows and samples at columns)
#' @param raw.mirna.mat (optional) A raw matrix of miRNA expression (genes at rows and samples at columns)
#' @param do.spqn A logical variable indication whether SpQN should be done
#' @param target.mat A binary matrix indicating target with miRNA on columns and mRNA on rows
#' @param group.label A factor vector subtype label for the samples
#' @param do.parallel A logical variable indicating whether parallization will be performed
#' @param n.iter An integer indicating the number of interaions for n.iter
#' @param seed A integer for permutation seed
#' @return A data frame similar to output from Run_diffreg() with additional permutation significance values
#' 
#' @import progress
#' @importFrom kSamples ad.test
#' @export
#' 
Run_permutation <- function(diffreg.res, mrna.mat, mirna.mat, 
							raw.mrna.mat=NULL, raw.mirna.mat=NULL,
							do.spqn=TRUE,
							target.mat, group.label, 
							do.parallel = FALSE, 
							n.iter=1000, seed=1234) {

	# Check if raw matrix is there; if not the mrna.mat and mirna.mat has expression value needed for normalization
	if(is.null(raw.mrna.mat)){
		message("Raw gene expression missing. Using mrna.mat as raw.mrna.mat.\n")
		raw.mrna.mat <- mrna.mat
	}

	if(is.null(raw.mirna.mat)){
		message("Raw miRNA expression missing. Using mirna.mat as raw.mrna.mat.\n")
		raw.mirna.mat <- mirna.mat
	}

	# Do permutation test: shuffle each rank list
	perm.res <- NULL
	perm.res[row.names(diffreg.res)] <- list(NULL)

	# Setup parallization
	if (isTRUE(do.parallel)) {
		numCores <- detectCores()
	} else {
		numCores = 1
	}

	# Check number of levels
	if (length(levels(group.label)) == 2) {

		pb <- progress_bar$new(total = 100)
		for (iter in seq(1:n.iter)) {

			pb$tick()

			# Get Permutation cor mat
			perm.label <- sample(group.label)
			perm.cor.mat.list <- BuildCorList(mrna.mat=mrna.mat, 
												mirna.mat=mirna.mat, 
												group.label=perm.label,
												numCores=numCores)

			# Do miSpQN on permutation cor mat
			if (isTRUE(do.spqn)) {
				perm.spqn.cor.mat.list <- RunGroupSpQN(mrna.mat=raw.mrna.mat, 
													mirna.mat=raw.mirna.mat, 
													group.label=perm.label, 
													cor.list=perm.cor.mat.list,
													numCores=numCores)
			} else {
				perm.spqn.cor.mat.list <- perm.cor.mat.list
			}

			# Calculate ad on permutation cor mat and Append to perm.res list
			for (i in names(perm.res)) {
				target.set <- names(target.mat[,i][target.mat[,i]!=0])

				perm.rank.list <- lapply(perm.spqn.cor.mat.list, function(x) sort(x[,i]))

				perm.pos.list <- lapply(perm.rank.list, function(x) match(target.set, names(x)))

				this.ks <- ks.test(perm.pos.list[[1]], perm.pos.list[[2]])

				# Organize output		
				perm.res[[i]][iter] <- this.ks$statistic
				# print(i)
			}

			Sys.sleep(1 / 100)

		}

		diffreg.res$Perm.pval <- unlist(lapply(row.names(diffreg.res), function(x) 1 - sum(diffreg.res[x, ]$KSval > perm.res[[x]])/n.iter))
		diffreg.res$FDR <- p.adjust(diffreg.res$Perm.pval, "fdr")
		diffreg.res <- diffreg.res[order(diffreg.res$ADval, decreasing=TRUE), ]

	} else {

		pb <- progress_bar$new(total = 100)
		for (iter in seq(1:n.iter)) {

			pb$tick()

			# Get Permutation cor mat
			perm.label <- sample(group.label)
			perm.cor.mat.list <- BuildCorList(mrna.mat=mrna.mat, 
												mirna.mat=mirna.mat, 
												group.label=perm.label,
												numCores=numCores)

			# Do miSpQN on permutation cor mat
			if (isTRUE(do.spqn)) {
				perm.spqn.cor.mat.list <- RunGroupSpQN(mrna.mat=raw.mrna.mat, 
													mirna.mat=raw.mirna.mat, 
													group.label=perm.label, 
													cor.list=perm.cor.mat.list,
													numCores=numCores)
			} else {
				perm.spqn.cor.mat.list <- perm.cor.mat.list
			}

			# Calculate ad on permutation cor mat and Append to perm.res list
			for (i in names(perm.res)) {
				target.set <- names(target.mat[,i][target.mat[,i]!=0])

				perm.rank.list <- lapply(perm.spqn.cor.mat.list, function(x) sort(x[,i]))

				perm.pos.list <- lapply(perm.rank.list, function(x) match(target.set, names(x)))

				this.ad <- ad.test(perm.pos.list, method="asymptotic", dist=F)

				# Organize output		
				perm.res[[i]][iter] <- this.ad$ad[1,2]
				# print(i)
			}

			Sys.sleep(1 / 100)

		}

		diffreg.res$Perm.pval <- unlist(lapply(row.names(diffreg.res), function(x) 1 - sum(diffreg.res[x, ]$ADval > perm.res[[x]])/n.iter))
		diffreg.res$FDR <- p.adjust(diffreg.res$Perm.pval, "fdr")
		diffreg.res <- diffreg.res[order(diffreg.res$ADval, decreasing=TRUE), ]

	}
	return(diffreg.res)
}


#' Do post-hoc pair-wise ks tests for the significant diff-reg miRNAs
#' 
#' @param diffreg.res Differential regulation results from Run_diffreg()
#' @param cor.list A list of correlation matrices with miRNA on columns and mRNA on rows
#' @param target.mat A binary matrix indicating target with miRNA on columns and mRNA on rows
#' @param significance A character value for type of significance threshold
#' @param cutoff A numeric value for significance threshold
#' @param alpha Weight for running AD/KS test
#' @return A list with ks scores and p-values matrices for each pair-wise comparison
#' 
#' @export
#' 
Pairwise_postpoc <- function(diffreg.res, cor.list, target.mat, significance=c("Pval", "Perm.pval", "FDR"), cutoff=0.05, alpha=NULL) {
	if(is.null(significance)) {
		stop("Must select one significance value from: Pval, Perm.pval, FDR.\n")
	}
	
	if (!length(significance) == 1) {
		stop("There should be only one 'significance'.\n")
	}

	if (!significance %in% c("Pval", "Perm.pval", "FDR")) {
		stop("'significance' must be in: 'Pval', 'Perm.pval' or 'FDR'.\n")
	}

	if(!significance %in% names(diffreg.res)) {
		stop("Please check if permutation test has been run or use other significance values.\n")
	}

	diffreg.res.use <- diffreg.res[diffreg.res[,significance] <= cutoff, ]
	target.mat.use <- target.mat[, row.names(diffreg.res.use)]

	ks.mat <- p.mat <- matrix(nrow=length(cor.list), 
								ncol=length(cor.list),
								dimnames = list(names(cor.list), names(cor.list)))

	posthoc.res <- NULL
	posthoc.res[row.names(diffreg.res.use)] <- list(list(KS=ks.mat, Pval=p.mat))

	posthoc.res <- posthoc.p.res <- NULL
	posthoc.res[row.names(diffreg.res.use)] <- list(list(KS=ks.mat, Pval=p.mat))
	posthoc.p.res[row.names(diffreg.res.use)] <- list(p.mat)


	for (i in names(cor.list)) {
		for (j in names(cor.list)[which(names(cor.list)==i):length(cor.list)]) {
			if (i == j) {
				next
			} else {
				cor.list.this <- cor.list[c(i,j)]
				ks.res.this <- .diffreg_ks(target.mat.use, cor.list.this, alpha=alpha)

				for (k in row.names(ks.res.this)) {
					posthoc.res[[k]][["KS"]][i,j] <- ks.res.this[k, "KSval"]
					posthoc.res[[k]][["Pval"]][i,j] <- ks.res.this[k, "Pval"]
				}
			}
		}
	}
	return(posthoc.res)
}


#' Do post-hoc pair-wise ks tests for the significant diff-reg miRNAs
#' 
#' @param cor.list A list of correlation matrices with miRNA on columns and mRNA on rows
#' @param miR miRNA to plot
#' @param target.mat A binary matrix indicating target with miRNA on columns and mRNA on rows
#' @param target.set A vector of targets for plotting
#' @param direction A character string specifying the direction of enrichment, must be "Negative" or "Positive"
#' @param pal Pallette to use; default will be used if NULL
#' @param alpha Weight for running AD/KS test
#' @return An enrichment plot
#' 
#' @import ggplot2
#' @import ggrepel
#' @import grid
#' @import gridExtra
#' @export
Plot_enrichment <- function(cor.list, miR, target.mat=NULL, target.set=NULL, direction="Negative", pal=NULL, alpha=NULL) {

	if (!length(direction) == 1) {
		stop("There should be only one 'direction'.\n")
	}

	if (!direction %in% c("Negative", "Positive")) {
		stop("'direction' must be either Negative or Positive.\n")
	}

	if (is.null(target.mat) & is.null(target.set)) {
		stop("Must provide target information to plot.\n")
	}

	n.group <- length(cor.list)
	n.gene <- nrow(cor.list[[1]])

	if (is.null(pal)) {
		pal <- default.pal(n.group, "Set1")
		names(pal) <- names(cor.list)
	}

	if (is.null(target.set) & !is.null(target.mat)) {
		target.set <- names(target.mat[,miR][target.mat[,miR]!=0])
	}


	rank.list <- lapply(cor.list, function(x) sort(x[,miR]))

	pos.list <- lapply(rank.list, function(x) match(target.set, names(x)))

	ws.list <- lapply(rank.list, function(x) .CalcGseaStat(x, target.set, alpha))

	if (direction == "Negative") {
		es <- unlist(lapply(ws.list, max))
	} else {
		es <- unlist(lapply(ws.list, min))
	}

	es.this <- data.frame(
							Rank=rep(1:n.gene, n.group), 
							Gene=unlist(lapply(ws.list, function(x) names(x))),
							miRNA=rep(miR, n.gene*n.group),
							Cor=unlist(ws.list),
							Subtype=rep(names(cor.list), each=n.gene),
							ES=unlist(lapply(ws.list, function(x) rank(-x)))
		)

	es.this$Subtype <- factor(es.this$Subtype, levels = names(cor.list))

	es.this$ES <- ifelse(es.this$ES==1, 1, 0)
	es.this[es.this$ES==1, ]$ES <- as.numeric(es)

	es.this$Position <- ifelse(es.this$Gene %in% target.set, es.this$Rank, NA)

	ratio.values.1 <- n.gene/(max(es.this$Cor)-min(es.this$Cor))
	ratio.values.2 <- n.gene/length(levels(es.this$Subtype))

	p1 <- ggplot(es.this, aes(x = Rank, y = Cor, color = Subtype, fill=Subtype, label=round(ES, digits=3))) +
			geom_hline(yintercept = 0, linetype="dashed", color = "grey75") +
			geom_line() +
			geom_label_repel(data = subset(es.this, ES!=0), aes(color=Subtype, fill=Subtype), size=2, fontface = 'bold', color = 'white') +
			scale_color_manual(values = pal) +
			scale_fill_manual(values = pal) +
			coord_fixed(ratio.values.1 / 2) +
			ggtitle(miR) +
			xlab("") + ylab("Enrichment Score") +
			xlim(0, n.gene) +
			guides(fill=guide_legend(title="")) +
			theme_bw() + 
			theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.text.x=element_blank(),
				axis.ticks.x = element_blank(),
				plot.margin = unit(c(1,1,-0.5,1), "line"), 
				legend.position='none')

	p2 <- ggplot(es.this, aes(x=Position, y=Subtype, fill=Subtype)) +
			geom_tile() +
			scale_fill_manual(values = pal) +
			coord_fixed(ratio.values.2 / 6) +
			xlab("Position") + ylab("Group") + 
			xlim(0, n.gene) +
			guides(fill=guide_legend(title="")) +
			theme_ridges(grid = FALSE) + 
			theme_bw() +
			theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				plot.margin = unit(c(-0.5,1,1,1), "line"), 
				legend.position='bottom')

	grid.newpage()
	grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
}