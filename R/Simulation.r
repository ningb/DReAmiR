#' Generate Simulation Data
#' 
#' @param Ntargets An integer giving the number of targets per miRNA
#' @param Nvar.percent A numeric value giving the percentage of targets to be variable between groups
#' @param Nsample An integer indicating the number of samples per group
#' @return A list containing simulated data
#' 
#' @importFrom lqmm make.positive.definite
#' @importFrom MASS mvrnorm
#' @importFrom kSamples ad.test
Generate_simulation <- function(Ntargets, Nvar.percent, Nsample) {
	set.seed(1)
	#
	# Target matrices for three groups
	# --------------------------------------------------------------------------
	gene <- paste("gene_", seq(1000), sep="")
	mir <- paste("miR_", seq(20), sep="")

	# Assign genes to targets by miRNAs <- predicted targets
	# Randomly splitted the genes into n groups where n=length(mir)
	predict.target <- list()
	for (i in mir) {
		predict.target[[i]] <- sample(gene, Ntargets, replace=TRUE)
		# gene <- setdiff(gene, unlist(predict.target))
	}
	# predict.target[[i]] <- gene
	names(predict.target) <- mir

	# Reset
	gene <- paste("gene_", seq(1000), sep="")

	# Now subset to three groups
	# First five miRNAs will have different targets per group
	# Last five will have the same predicted targets
	target.A <- target.B <- target.C <- list()
	for (i in mir) {
		if (which(mir == i) <= 10) {
			target.A[[i]] <- sample(predict.target[[i]], round(Ntargets * (Nvar.percent / 100)))
			target.B[[i]] <- sample(predict.target[[i]], round(Ntargets * (Nvar.percent / 100)))
			target.C[[i]] <- sample(predict.target[[i]], round(Ntargets * (Nvar.percent / 100)))
		} else {
			target.A[[i]] <- target.B[[i]] <- target.C[[i]] <- predict.target[[i]]
		}
		# print(i)
	}

	#
	# Correlation matrices
	# Group-specific targets have stronger regulation
	# --------------------------------------------------------------------------
	all.names <- c(gene, mir)
	adj.A <- matrix(0, nrow=length(all.names), ncol=length(all.names), dimnames=list(all.names, all.names))
	adj.B <- matrix(0, nrow=length(all.names), ncol=length(all.names), dimnames=list(all.names, all.names))
	adj.C <- matrix(0, nrow=length(all.names), ncol=length(all.names), dimnames=list(all.names, all.names))

	for (i in mir) {
			adj.A[predict.target[[i]], i] <- rnorm(length(predict.target[[i]]), mean=-0.2, sd=0.3)
			adj.B[predict.target[[i]], i] <- rnorm(length(predict.target[[i]]), mean=-0.2, sd=0.3)
			adj.C[predict.target[[i]], i] <- rnorm(length(predict.target[[i]]), mean=-0.2, sd=0.3)
		if (which(mir == i) <= 10) {
			beta2 <- sample(1:15, 3)
			beta.a <- rbeta(length(target.A[[i]]), 20, beta2[1])
			adj.A[target.A[[i]], i] <- -0.2 + beta.a * (min(adj.A[,i]) + 0.2)
			beta.b <- rbeta(length(target.B[[i]]), 20, beta2[2])
			adj.B[target.B[[i]], i] <- -0.2 + beta.b * (min(adj.B[,i]) + 0.2)
		}
		# print(i)
	}

	#
	# Turn to expression matrix
	# --------------------------------------------------------------------------
	diag(adj.A) <- 1
	diag(adj.B) <- 1
	diag(adj.C) <- 1
	adj.A[lower.tri(adj.A)] = t(adj.A)[lower.tri(adj.A)]
	adj.B[lower.tri(adj.B)] = t(adj.B)[lower.tri(adj.B)]
	adj.C[lower.tri(adj.C)] = t(adj.C)[lower.tri(adj.C)]

	adj.A.pd <- make.positive.definite(adj.A)
	adj.B.pd <- make.positive.definite(adj.B)
	adj.C.pd <- make.positive.definite(adj.C)

	# Generate multivariate guassian networks
	expr.A <- t(mvrnorm(Nsample, mu=rep(0, length(all.names)), Sigma=adj.A.pd))
	expr.B <- t(mvrnorm(Nsample, mu=rep(0, length(all.names)), Sigma=adj.B.pd))
	expr.C <- t(mvrnorm(Nsample, mu=rep(0, length(all.names)), Sigma=adj.C.pd))

	colnames(expr.A) <- paste("Sample_A", seq(Nsample), sep="")
	colnames(expr.B) <- paste("Sample_B", seq(Nsample), sep="")
	colnames(expr.C) <- paste("Sample_C", seq(Nsample), sep="")

	mrna.mat <- cbind(expr.A[gene, ], expr.B[gene, ], expr.C[gene, ])
	mirna.mat <- cbind(expr.A[mir, ], expr.B[mir, ], expr.C[mir, ])

	#
	# Output Results
	# --------------------------------------------------------------------------
	group.label <- as.factor(rep(c("A", "B", "C"), each=Nsample))
	names(group.label) <- colnames(mrna.mat)

	sim.data <- list(predict.target = predict.target,
					group.label = group.label, 
					mrna.mat = mrna.mat, 
					mirna.mat = mirna.mat)

	return(sim.data)
}

# Example
# sim.data <- Generate_simulation(Ntargets=70, Nvar.percent=75, Nsample=100)