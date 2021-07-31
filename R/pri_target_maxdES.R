#' Calculating background ES distribution for normalizing dES
#' 
#' @param gs A string vector for set of targets of interests. Could be an element from output of GetLE()
#' @param rk A numeric vector for all mRNAs ranked by correlation
#' @param direction A character for the direction of enrichment 
#' @param n.iter An integer for the number of background to generate. Default is 500
#' @param alpha A numeric value for the gsea param
#' @return A numeric vector of background normalized ES distribution
#' 
.Calc_NES_bkgd <- function(gs, rk, direction, n.iter, alpha) {

	if (direction == "Negative") {
		perm.es <- unlist(lapply(seq_len(n.iter), function(x) max(.CalcGseaStat(sample(rk), gs, alpha))))
	} else {
		perm.es <- unlist(lapply(seq_len(n.iter), function(x) min(.CalcGseaStat(sample(rk), gs, alpha))))
	}

	bkgd <- mean(perm.es)
	return(bkgd)

}

#' Remove each gene in the gs and see what is the changes in the dES
#' 
#' @param current.gs A character vector for a set of targets. One element from output of GetLE().
#' @param main.rk A numeric vector for all mRNAs ranked by correlation in the group of interest
#' @param ref.rk A numeric vector for all mRNAs ranked by correlation in the group to compare with
#' @param direction A character for the direction of enrichment 
#' @param n.iter An integer for the number of background to generate. Default is 500
#' @param alpha A numeric value for the gsea param
#' @return A numberic vector of the dES for removing each gene in gs
#' 
.Diff_by_Remove <- function(gs, main.rk, ref.rk, direction, n.iter, alpha) {

	if (direction == "Negative") {

		main.bkgd <- .Calc_NES_bkgd(gs[-1], main.rk, direction=direction, n.iter=n.iter, alpha=alpha)
		main.nes <- unlist(lapply(gs, function(x) max(.CalcGseaStat(main.rk, gs[gs != x], alpha)) / main.bkgd))

		ref.bkgd <- .Calc_NES_bkgd(gs[-1], ref.rk, direction=direction, n.iter=n.iter, alpha=alpha)
		ref.nes <- unlist(lapply(gs, function(x) max(.CalcGseaStat(ref.rk, gs[gs != x], alpha)) / ref.bkgd))

		remove.des <- main.nes - ref.nes

	} else {

		main.bkgd <- .Calc_NES_bkgd(gs[-1], main.rk, direction=direction, n.iter=n.iter, alpha=alpha)
		main.nes <- unlist(lapply(gs, function(x) min(.CalcGseaStat(main.rk, gs[gs != x], alpha)) / main.bkgd))

		ref.bkgd <- .Calc_NES_bkgd(gs[-1], ref.rk, direction=direction, n.iter=n.iter, alpha=alpha)
		ref.nes <- unlist(lapply(gs, function(x) min(.CalcGseaStat(ref.rk, gs[gs != x], alpha)) / ref.bkgd))

		remove.des <- ref.nes - main.nes
	}

	return(remove.des)

}

#' Do backward random search to find set of targets that maximize the normalized dES
#' 
#' @param current.gs A character vector for a set of targets. One element from output of GetLE().
#' @param cor.list A matrix from Run_LINE()
#' @param main A character for the name of the group to compared to. 
#' @param ref A character for the name of the group to compared to. 
#' @param direction A character for the direction of enrichment 
#' @param min.gs.size An integer indicating the minimum number of genes to used. Default is 20.
#' @param n.iter An integer indicating the number of background to generate for normlizing ES. Default is 50
#' @param alpha A numeric variable for gsea parameter
#' @return A character vector of target names
#' 
#' @examples
#' Find_Target_Max_dES(le.gene[[1]], cor.list, "A", "B", direction="Negative")
#' 
#' @import progress
#' @export
Find_Target_Max_dES <- function(gs, cor.list, main, ref, 
								direction=c("Negative", "Positive"), 
								min.gs.size=20, n.iter=500, alpha=0.1) {

	# Check if main and ref are in cor.list
	if(!main %in% names(cor.list) | !ref %in% names(cor.list)) {
		stop("The groups are not in cor.list.\n")
	}
	
	# Check if direction is correct
	if (!length(direction) == 1) {
		stop("There should be only one 'direction'.\n")
	}

	if (!direction %in% c("Negative", "Positive")) {
		stop("'direction' must be either Negative or Positive.\n")
	}

	# Check if the gs is in right format; if not, will use the names()
	if(!is.null(names(gs))) {
		gs <- names(gs)
	}

	# Check if the min.gs.size has already been reached
	if (length(gs) <= min.gs.size) {
		stop("The min.gs.size has already been reached. No maximization done.\n")
	}

	# Format gs and cor.list
	rank.list <- lapply(cor.list, function(x) sort(x[,miR]))

	ws.list <- lapply(rank.list, function(x) .CalcGseaStat(x, gs, alpha))

	# Calculate NES and dES from the full gs
	if (direction == "Negative") {
		es <- unlist(lapply(ws.list, max))
		
		main.nes <- es[main] / .Calc_NES_bkgd(gs, rank.list[[main]], direction=direction, n.iter=n.iter, alpha=alpha)
		ref.nes <- es[ref] / .Calc_NES_bkgd(gs, rank.list[[ref]], direction=direction, n.iter=n.iter, alpha=alpha)

		baseline.des <- main.nes - ref.nes
	} else {
		es <- unlist(lapply(ws.list, min))
		
		main.nes <- es[main] / .Calc_NES_bkgd(gs, rank.list[[main]], direction=direction, n.iter=n.iter, alpha=alpha)
		ref.nes <- es[ref] / .Calc_NES_bkgd(gs, rank.list[[ref]], direction=direction, n.iter=n.iter, alpha=alpha)

		baseline.des <- ref.nes - main.nes
	}

	# Start doing random backward search
	genes.des <- c()
	current.gs <- gs

	pb <- progress_bar$new(total = 100)
	for (i in seq_len(length(gs) - min.gs.size)) {
		  
		pb$tick()

		# Start from all, Calculate new remove.des
		new.es.res <- .Diff_by_Remove(gs=current.gs, main.rk=rank.list[[main]], ref.rk=rank.list[[ref]], direction=direction, n.iter=n.iter, alpha=alpha)
		
		# Decide which one to remove
		this.remove <- current.gs[which(new.es.res == max(new.es.res))]

		# Pick a ranodm one if ties appear
		this.remove <- ifelse(length(this.remove) > 1, this.remove[sample(1:length(this.remove), 1)], this.remove)

		# Update gene set by deleting the target that gives the largest change
		current.gs <- current.gs[current.gs != this.remove]
		
		# Update dES score
		new.des <- max(new.es.res)

		# Record the gene removed and the score at the current step
		genes.des <- setNames(c(genes.des, new.des), c(names(genes.des), this.remove))

		# print(paste("Iteration ", i, " removes ", this.remove, " new dES=", new.des, sep=""))

		# Once we reach the desired number of genes, add the remainig genes to the genes.des and exit loop
		if (length(current.gs) <= min.gs.size) {
			print("Reached minimum gs size.")
			remaining.des <- setNames(rep(0, length(current.gs)), current.gs)
			genes.des <- c(genes.des, remaining.des)
			break
		}

		Sys.sleep(1 / 100)

	}

	# Find the gene removing which gives the largest dES and select the remainings
	final.gs <- names(genes.des[(which(genes.des == max(genes.des))+1):length(genes.des)])

	# Compare with the full one to see if we actually have increase
	if (max(genes.des) <= baseline.des) {
		warnings("Can't find better gs for increasing dES. Full one should be use.")
		final.gs <- gs
	}

	return(final.gs)

}