#' Given a list of correlation matrcies, return the SpQN normalized matrices
#' 
#' @param le.gene A list of correlation matrices with miRNA on columns and mRNA on rows
#' @param dim An integer indicating number embedded features returned from LINE
#' @param edge.sign An character for edge list weight; if negative, the reverse will be use
#' @param do.dim.reduction A logical variable, whether to return dimension reduction results
#' @return A matrix of graph embedding features for all the leading edge genes
#' 
#' @export
Run_LINE <- function(le.gene, dim=5, edge.sign=c("Negative", "Positive", "Both"), do.dim.reduction=TRUE) {

	if (!length(edge.sign) == 1) {
		stop("There should be only one 'edge.sign'.\n")
	}

	if (!edge.sign %in% c("Negative", "Positive", "Both")) {
		stop("'edge.sign' must be: Negative, Positive or Both.\n")
	}

	if (is.null(names(le.gene[[1]]))) {
		stop("le.gene should be the weights with gene as names to run graph embedding.\n")
	} else {
		all.targets <- unlist(lapply(le.gene, names))
	}

	# Build edgelist between subtype-miRNA and gene nodes
	# Using weight from cor.list
	target.edgelist <- data.frame(miRNA = rep(names(le.gene), times = unlist(lapply(le.gene, length))),
									Targets = all.targets,
									weight = unlist(le.gene))

	# Make the weight into the correct sign and normalize between 0 and 1

	if (edge.sign == "Negative") {
		target.edgelist$weight <- (-target.edgelist$weight + 1) / 2
	} else if (edge.sign == "Positive") {
		target.edgelist$weight <- (target.edgelist$weight + 1) / 2
	} else {
		target.edgelist$weight <- (abs(target.edgelist$weight) + 1) / 2
	}

	# Do bipartite network
	target.edgelist2 <- target.edgelist[,c(2,1,3)]
	names(target.edgelist2) <- names(target.edgelist)
	target.edgelist <- rbind(target.edgelist, target.edgelist2)
	names(target.edgelist) <- c("u", "v", "w")

	#
	# Run graph embedding
	# --------------------------------------------------------------------------
	reconstruct_df <- reconstruct(target.edgelist, max_depth = 2, max_k = 10)
	line_one_matrix <- line(reconstruct_df, dim = dim, order = 1, negative = 5, samples = 10, rho = 0.025)
	line_two_matrix <- line(reconstruct_df, dim = dim, order = 2, negative = 5, samples = 10, rho = 0.025)
	concatenate_matrix <- concatenate(line_one_matrix, line_two_matrix)
	normalize_matrix <- normalize(concatenate_matrix)

	# Check if doing PCA is needed (mostly for visualization)
	if (do.dim.reduction == TRUE) {

		count.degree <- table(all.targets)

		count_pca <- prcomp(scale(normalize_matrix))
		pc.mat <- count_pca$x
		embed.res <- data.frame(
						Dim1 = pc.mat[,1],
						Dim2 = pc.mat[,2],
						row.names = row.names(normalize_matrix))

		embed.res$Node.type <- ifelse(row.names(embed.res) %in% names(le.gene), row.names(embed.res), "Targets")
		embed.res$degree <- as.numeric(count.degree[row.names(embed.res)])

		le.mat <- matrix(nrow=nrow(embed.res), 
						ncol=length(le.gene), 
						dimnames=list(row.names(embed.res), 
						paste(names(le.gene), ".target", sep="")))

		for (i in 1:length(le.gene)) {
			le.mat[, i] <- ifelse(row.names(le.mat) %in% names(le.gene[[i]]), "1", "0")
		}

		embed.res <- cbind(embed.res, le.mat)

		output <- list(LINE=normalize_matrix, PCA=embed.res)

		return(output)
	} else {
		return(normalize_matrix)
	}
}

#' Try to find optimal cluster number fro fuzzy k meansand make plots for k selection
#' 
#' @param normalize_matrix A matrix from Run_LINE()
#' @param k.range A numeric vector for k numbers to try
#' @param do.plot A Logic variable indicating whether plots are shown
#' @return An enrichment plot
#' 
#' @import ggplot2
#' @importFrom e1071 cmeans
#' @importFrom cluster silhouette
#' @import grid
#' @import gridExtra
#' @export
selectK <- function(normalize_matrix, k.range=2:10, do.plot=TRUE) {

	output.stat <- sapply(k.range, function(k) {
		cl <- e1071::cmeans(scale(normalize_matrix), centers = k)
		ss <- silhouette(cl$cluster, dist(normalize_matrix))
		ave.ss <- ifelse(k == 1, 0, summary(ss)[["avg.width"]])
		c(cl$withinerror, ave.ss)
	})

	dat.plot <- data.frame(
						K             = k.range,
						tot_within_ss = output.stat[1,],
						ave.ss        = output.stat[2,]
					)

	if (do.plot) {

		p1 <- ggplot(data=dat.plot, aes(x=k.range, y=tot_within_ss)) + 
				geom_line(color="grey75") +
				geom_point(color="steelblue") +
				xlab("K") + ylab("Total Within-Cluster Variation") +
				scale_x_continuous(breaks = k.range) +
				theme_bw()
		p2 <- ggplot(data=dat.plot, aes(x=k.range, y=ave.ss)) + 
				geom_vline(xintercept = dat.plot[which(dat.plot$ave.ss == max(dat.plot$ave.ss)), ]$K, color="grey25") +
				geom_line(color="grey75") +
				geom_point(color="steelblue") +
				xlab("K") + ylab("Average Silhouette") +
				scale_x_continuous(breaks = k.range) +
				theme_bw()

		grid.newpage()
		grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
	}

	return(dat.plot)

}

#' Try to find optimal cluster number fro fuzzy k meansand make plots for k selection
#' 
#' @param normalize_matrix A matrix from Run_LINE()
#' @param select.k A numeric for k numbers to use
#' @param group.label A factor vector subtype label for the samples
#' @param p.cutoff A numeric variable indicating the p-value cutoff when assigning clusters to groups
#' @return A list in which GeneList shows targets genes to each group, Label stores the group assigned to each cluster
#' 
#' @importFrom e1071 cmeans
#' @export
AssignTargetCluster <- function(normalize_matrix, select.k=7, group.label, p.cutoff=0.05) {
	cluster.res <- cmeans(scale(normalize_matrix), centers = select.k)
	
	# Aissgn those can be done easily
	assigned.clusters <- cluster.res$cluster[group.label]

	# First test if all miR nodes are in the same clusters
	if (all(assigned.clusters==assigned.clusters[1])) {
		cluster.gene <- lapply(group.label, function(x) names(cluster.res$cluster)[cluster.res$cluster == assigned.clusters[1]])
		names(cluster.gene) <- group.label

		cluster.label <- lapply(as.character(seq(select.k)), function(x){
						if (x != assigned.clusters[1]) {
							label <- "Common"
						} else {
							label <- paste(group.label, collapse="+")
						}
						return(label)
						})
		names(cluster.label) <- as.character(seq(select.k))

	} else {

		# Clusters directly assigned by miRNA node membership or hasn't been clustered
		cluster.assigned <- as.character(unique(sort(assigned.clusters)))
		cluster.use <- names(table(cluster.res$cluster))
		cluster.tbd <- cluster.use[!cluster.use %in% cluster.assigned]

		# First assigned those that already labelled by miRNA position
		assign.res <- as.data.frame(cluster.res$membership[, cluster.assigned])
		assign.res$Cluster <- cluster.res$cluster
			
		# For the remaining usable clusters, calculate per miRNA-cluster
		# one-vs-others one sided KS, and label as one miRNA or multiple
		# or none
		ks.res <- matrix(nrow=length(cluster.tbd), ncol=length(cluster.assigned), dimnames=list(cluster.tbd, cluster.assigned))
		for (i in cluster.tbd) {
			use.res <- assign.res[assign.res$Cluster == i, ]
			use.res <- gather(use.res, Names, Value, cluster.assigned[1]:cluster.assigned[length(cluster.assigned)])
			
			this.ks <- c()
			for (j in cluster.assigned){
				this.res <- ks.test(use.res[use.res$Names == j, ]$Value, use.res[use.res$Names != j, ]$Value, alternative="less")
				this.ks <- c(this.ks, this.res$p.value)
			}

			ks.res[i, ] <- this.ks
		}
		ks.res <- ifelse(ks.res <= p.cutoff, 1, 0)

		# Assign the subtype labels to clusters
		assigned.dat <- as.matrix(data.frame(cluster=assigned.clusters, type=group.label))
		tbd.list <- lapply(cluster.tbd, function(x) colnames(ks.res)[ks.res[x, ] == 1])
		tbd.dat <- lapply(tbd.list, function(x) names(assigned.clusters)[as.character(assigned.clusters) %in% x])
		names(tbd.dat) <- row.names(ks.res)

		# Organize to binary matrix
		bin.assign.mat <- matrix(nrow=length(cluster.use), ncol=length(group.label), dimnames=list(cluster.use, group.label))
		bin.assign.mat[assigned.dat[, c(1,2)]] <- 1
		for (n in names(tbd.dat)) {
			bin.assign.mat[n, tbd.dat[[n]]] <- 1
		}

		# Output final cluster assignment as targets per subtype
		cluster.gene <- lapply(group.label, function(x) {
						type.vec <- bin.assign.mat[,x]
						cluster.pick <- names(type.vec)[!is.na(type.vec)]
						cluster.gene <- names(cluster.res$cluster)[cluster.res$cluster %in% as.numeric(cluster.pick)]
						cluster.gene[!cluster.gene %in% levels(group.label)]
						})
		names(cluster.gene) <- group.label
		
		cluster.label <- lapply(row.names(bin.assign.mat), function(x){
						cluster <- bin.assign.mat[x, ]
						label <- names(cluster)[!is.na(cluster)]
						if (length(label) == 0) {
							label <- "Common"
						} else if (length(label) == 1){
							label <- label
						} else {
							label <- paste(label, collapse="+")
						}
						return(label)
						})
		names(cluster.label) <- row.names(bin.assign.mat)

	}

	return(list(GeneList=cluster.gene, Label=cluster.label, Cluster=cluster.res$cluster))
}