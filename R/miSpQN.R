
# ****************************************************************************
# *           Try to modify SPQN to accomendate miRNA-mRNA network           *
# *                   1) Binning miRNA and mRNA separately                   *
# *              2) Examine most negative correlation as signal              *
# *                   3) Use non-square correlation matrix                   *
# *                    4) Normalize to top right reference                   *
# ****************************************************************************

# ==========================================================================
# Examine mean-correlation relationship
# ==========================================================================

.get_grp_loc <- function(cor_mat, ngrp=10){
    ngene <- nrow(cor_mat)
    idx <- seq_len(ngene)
    grp_label <- cut(idx, ngrp)
    grp_loc1 <- split(idx, grp_label)

    nmir <- ncol(cor_mat)
    idx <- seq_len(nmir)
    grp_label <- cut(idx, ngrp)
    grp_loc2 <- split(idx, grp_label)

    return(list(grp_loc1, grp_loc2))
}

#' Calculate IQR and Q1 per bin x bin; Calculate expression mean per bin
#' 
#' @param cor_mat A correlation matrix for a group with mRNA on the rows and miRNA on the columns
#' @param ave_exp1 A numeric vecotr with same length as the nrow(cor_mat) for average mRNA expressions
#' @param ave_exp2 A numeric vecotr with same length as the ncol(cor_mat) for average miRNA expressions
#' @return a list for IQR/Q1/average expression for mRNA and miRNA
#' 
#' @examples
#' IQR_list <- get_condition_exp(cor_mat, ave_exp1, ave_exp2)
#' 
#' @export
get_condition_exp <- function(cor_mat, ave_exp1, ave_exp2){
    idx1 <- order(ave_exp1)
    idx2 <- order(ave_exp2)

    cor_mat <- cor_mat[idx1, idx2]
    ave_exp1 <- ave_exp1[idx1]
    ave_exp2 <- ave_exp2[idx2]

    grp_loc <- .get_grp_loc(cor_mat)
    IQR_cor_mat <- array(dim=c(10,10))
    Q1_cor_mat <- array(dim=c(10,10))
    grp_mean1 <- grp_mean2 <- array(dim=10)

    for(i in seq_len(10)) {
        ave_exp_grp1 <- ave_exp1[grp_loc[[1]][[i]]]
        grp_mean1[i] <-  mean(ave_exp_grp1)
        for(j in seq_len(10)) {
            cor_tmp <- cor_mat[grp_loc[[1]][[i]],grp_loc[[2]][[j]]]
            IQR_cor_mat[i,j] <- IQR(cor_tmp)
            Q1_cor_mat[i,j] <- quantile(cor_tmp)[2]
            if (i == 10) {
                ave_exp_grp2 <- ave_exp2[grp_loc[[2]][[j]]]
                grp_mean2[j] <-  mean(ave_exp_grp2)
            }
        }
    }

    return(list(IQR_cor_mat=IQR_cor_mat,
                Q1_cor_mat=Q1_cor_mat,
                grp_mean1=grp_mean1,
                grp_mean2=grp_mean2))
}

#' Plot heatmap shwoing the Q1 and IQR of correlation coefficients distribution per bin x bin
#' 
#' @examples 
#' plot_IQR_condition_exp_Corplot(IQR_list)
#' 
#' @import corrplot
#' @import viridis
#' @export
plot_IQR_condition_exp_Corplot <- function(IQR_list) {
    IQR_cor_mat <- IQR_list$IQR_cor_mat
    Q1_cor_mat <- IQR_list$Q1_cor_mat
    grp_mean1 <- IQR_list$grp_mean1
    grp_mean2 <- IQR_list$grp_mean2

    row.names(IQR_cor_mat) <- row.names(Q1_cor_mat) <- as.character(round(grp_mean1, 3))
    colnames(IQR_cor_mat) <- colnames(Q1_cor_mat) <- as.character(round(grp_mean2, 3))

    IQR_cor_mat <- IQR_cor_mat[nrow(IQR_cor_mat):1, ]
    Q1_cor_mat <- Q1_cor_mat[nrow(Q1_cor_mat):1, ]

    c_colors <- viridis(15, option="A")

    corrplot(IQR_cor_mat, is.corr=FALSE, method="color",
        tl.col = "black", tl.srt = 45,
        pch.col = "white",
        col = c_colors)

    corrplot(Q1_cor_mat, is.corr=FALSE, method="color",
        tl.col = "black", tl.srt = 45,
        pch.col = "white",
        col = c_colors)
}

#' Plot density plot showing signal (most negative correlation) vs. background per bin
#' 
#' @examples 
#' plot_signal_condition_exp(cor_mat, ave_exp1, ave_exp2, signal=0.001)
#' 
#' @import ggplot2
#' @import ggridges
#' @export
plot_signal_condition_exp <- function(cor_mat, ave_exp1, ave_exp2, signal=0.001) { 
    exp_idx1 <- order(ave_exp1)
    exp_idx2 <- order(ave_exp2)

    cor_mat <- cor_mat[exp_idx1, exp_idx2]

    ncor <- nrow(cor_mat) * ncol(cor_mat)
    nsignal <- ncor*signal
    nsignal_grp <- round(nsignal/100)
    
    # Plot by mRNA bins (ave_exp1, rows) or miRNA bins (ave_exp2, cols)
    grp_locs <- .get_grp_loc(cor_mat)

    for (type in seq_len(length(grp_locs))) {
        list_cor_sig <- grp_sig <- numeric(10*nsignal_grp)
        nbackground_group <- lapply(grp_locs[[type]], length)
        nbackground_group <- unlist(nbackground_group)
        if (type == 1) {
            nbackground_group <- nbackground_group * ncol(cor_mat)
        } else {
            nbackground_group <- nbackground_group * nrow(cor_mat)            
        }
        list_cor_back_cumulate <- cumsum(nbackground_group)
        list_cor_back <- grp_back <- numeric(list_cor_back_cumulate[10])

        for(ngrp in seq_len(10)) {
            loc_back <- grp_locs[[type]][[ngrp]]
            if (type == 1) {
                cor_back <- cor_mat[loc_back, ]
            } else {
                cor_back <- cor_mat[, loc_back]
            }
            order_cor_ori_grp <- order(cor_back)
            idx <- seq_len(nsignal_grp)
            id_signal <- order_cor_ori_grp[idx]
            cor_signal_ori <- cor_back[id_signal]

            idx1 <- (ngrp-1)*nsignal_grp+1
            idx2 <- ngrp*nsignal_grp
            idx <- c(idx1 : idx2)
            list_cor_sig[idx] <- cor_signal_ori
            grp_sig[idx] <- ngrp
          
            if(ngrp==1){
                idx <- seq_len(list_cor_back_cumulate[ngrp])
                list_cor_back[idx] <- cor_back
                grp_back[idx] <- ngrp
            }else{
                idx1 <- list_cor_back_cumulate[ngrp-1]+1
                idx2 <- list_cor_back_cumulate[ngrp]
                idx <- c(idx1:idx2)
                list_cor_back[idx] <- cor_back
                grp_back[idx] <- ngrp }                                                                                                       
        }

        df_sig <- data.frame(correlation=list_cor_sig, bin=grp_sig, group="signal")
        df_back <- data.frame(correlation=list_cor_back, 
                              bin=grp_back, group="background")
        df_merge <- rbind(df_sig,df_back)
        df_merge$bin_group  <-  paste(df_merge$bin, df_merge$group)

        title <- ifelse(type == 1, "mRNA Bins", "miRNA Bins")

        p <- ggplot(df_merge, aes_string(y = "bin")) +
            geom_density_ridges(
            aes_string(x = "correlation", fill = "bin_group"),
            alpha = .8, color = "white") +
            labs(x = "correlation", y = "bin")  +
            scale_y_discrete(limits=seq_len(10)) +
            scale_fill_cyclical(
                labels = c("background","signal"),
                values = c("grey75", "firebrick"),
                name = "group", guide = "legend") +
            labs(title = title) +
            theme_ridges(grid = FALSE)
        print(p)
    }   
}

# ==========================================================================
# Remove mean-correlation relationship
# ==========================================================================

##### Asssign inner bins
## In each running group, get a inner group
.get_grps <- function(cor_mat, ngrp, size_grp, type="row"){
    if (type == "row") {
        ngene <- nrow(cor_mat)
    } else {
        ngene <- ncol(cor_mat)
    }
    idx <- seq_len(ngene-size_grp+1)
    grp_label <- cut(idx, ngrp-1)
    grp_loc1 <- split(idx, grp_label)

    if(size_grp-length(grp_loc1[[1]])<5){
        idx <- seq_len(ngene)
        grp_label <- cut(idx, ngrp)
        grp_loc <- split(idx, grp_label)
    }else{
        grp_loc <- vector(mode = "list", length = ngrp)
        idx <- seq_len(ngene-size_grp+1)
        grp_label <- cut(idx, ngrp-1)
        grp_loc0 <- split(idx, grp_label)
        for(i in seq_len(ngrp-1)){
            id1 <- grp_loc0[[i]][1]
            id2 <- grp_loc0[[i]][1]+size_grp-1
            grp_loc[[i]] <- c(id1:id2)
        }
        id1 <- ngene-size_grp+1
        grp_loc[[ngrp]] <- c(id1:ngene)
    }
    grp_loc
}

##### Asssign inner bins
## In each running group, get a inner group
.get_grps_inner <- function(grp_loc){
    grp_loc_inner <- list()
    ngrp <- length(grp_loc)
    size_bin <- length(grp_loc[[1]])
    ngene <- max(grp_loc[[ngrp]])
    
    width_tmp <- grp_loc[[2]][1] - grp_loc[[1]][1]
    length_inner_1 <- round(size_bin/2 + width_tmp/2)
    grp_loc_inner[[1]] <- seq_len(length_inner_1)
    
    for(i in 2:(ngrp-1)){
        id1 <- grp_loc[[i+1]][1]
        id2 <- grp_loc[[i]][1]
        width_tmp <- id1 - id2
        tail <- tail(grp_loc_inner[[i-1]], 1)
        grp_loc_inner[[i]] <- (tail + 1):(tail + width_tmp)
    }
    tail <- tail(grp_loc_inner[[ngrp-1]],1)
    grp_loc_inner[[ngrp]] <- c((tail+1) : ngene)
    
    grp_loc_inner
}

##### Get rank for each running bin
.get_bin_rank <- function(cor_obs, grp_loc, group_loc_adj, cor_ref){
    rank_bin <- rank_bin_pre <- array(dim=dim(cor_obs))
    ngrp <- length(grp_loc[[1]])
    
    l_cor_tmp_ref <- length(cor_ref)
    
    for(i in seq_len(ngrp)){
        for(j in seq_len(ngrp)){
            cor_bin_tmp <- cor_obs[grp_loc[[1]][[i]],grp_loc[[2]][[j]]]
            rank_bin_tmp <- array(dim=dim(cor_bin_tmp))
            l_cor_tmp <- length(cor_bin_tmp)
            
            idx <- seq_len(l_cor_tmp)
            rank_bin_tmp[idx] <- rank(cor_bin_tmp[idx])
            
            ## Scale the rank of each bin to same scale as rank(cor_ref),
            ## After scaling, rank_bin could contain non-integars
            rank_bin_pre[grp_loc[[1]][[i]], grp_loc[[2]][[j]]] <- rank_bin_tmp
            rank_replace <- rank_bin_pre[group_loc_adj[[1]][[i]], group_loc_adj[[2]][[j]]]
            # scale the rank to range [0,1]
            rank_replace <- (rank_replace-1) / l_cor_tmp
            # scale the rank to the same scale as the referent bin
            rank_replace <- 1 + (rank_replace * (l_cor_tmp_ref-1))
            # change the ranks that larger than the size of referent bin, to be the size of referent bin
            rank_replace[which(rank_replace>l_cor_tmp_ref)] <- l_cor_tmp_ref
            
            rank_bin[group_loc_adj[[1]][[i]],group_loc_adj[[2]][[j]]] <- rank_replace
        }
    }
    
    rank_bin
}

##### Transform rank to cor_est
.est_cor <- function(rank_bin, cor_ref){
    cor_adj <- array(dim=dim(rank_bin))
    
    cor_ref_sorted <- sort(cor_ref)
    
    ## Find to nearest integars to rank_bin, since rank_bin could contain non-integars
    rank_bin <- array(rank_bin)
    rank_bin1 <- rank_bin %/% 1

    ## Assign weights to each rank according to the distance to rank_bin
    rank_bin_w2 <- (rank_bin - rank_bin1)

    rank_bin2 <- rank_bin1+1
    rank_bin_w1 <- 1-rank_bin_w2

    ## Find the correlations in the cor_ref corresponding to the two nearest ranks to rank_bin
    ## Estimate the correlation using weighted average based on the distance to rank_bin
    length_cor_ref <- length(cor_ref_sorted)
    rank_bin2[which(rank_bin2>length_cor_ref)] <- length_cor_ref
    cor_adj[,] <- rank_bin_w1*cor_ref_sorted[rank_bin1] + rank_bin_w2*cor_ref_sorted[rank_bin2]

    cor_adj
}

#' Main function to remove mean-correlation relationship between miRNA and mRNA
#'
#' @param cor_mat A correlation matrix for a group with mRNA on the rows and miRNA on the columns
#' @param ave_exp1 A numeric vecotr with same length as the nrow(cor_mat) for average mRNA expressions
#' @param ave_exp2 A numeric vecotr with same length as the ncol(cor_mat) for average miRNA expressions
#' @param ngrp An integer for number of inner bins. Default is 20
#' @param size_grp A numeric vector of length 2 for number of mRNA and miRNA per outer group. Default c(1000, 300)
#' @param ref_grp An integer for the index of bin to use as reference. Default is 18
#' @return A normalized correlation matrix with the same dimension as cor_mat
#' 
#' @import corrplot
#' @import viridis
#' @export
normalize_correlation <- function(cor_mat, ave_exp1, ave_exp2, ngrp=20, size_grp=c(1000,300), ref_grp=18){
    stopifnot(length(ngrp) == 1, length(size_grp) == 2, length(ref_grp) == 1)
    stopifnot(ref_grp <= ngrp)
    stopifnot(0 < ngrp, 0 < ref_grp)
    stopifnot(is.matrix(cor_mat))
    stopifnot(nrow(cor_mat) == length(ave_exp1), ncol(cor_mat) == length(ave_exp2))

    idx1 <- order(ave_exp1)
    idx2 <- order(ave_exp2)
    original.names <- dimnames(cor_mat)
    cor_mat <- cor_mat[idx1, idx2]

    ngene <- length(ave_exp1)
    nmir <- length(ave_exp2)
    ref_vec1 <- seq_len(ngene)
    ref_vec_rearranged1 <- ref_vec1[idx1]
    ref_vec2 <- seq_len(nmir)
    ref_vec_rearranged2 <- ref_vec2[idx2]

    group_loc1 <- .get_grps(cor_mat, ngrp, size_grp[1], type="row")
    group_loc2 <- .get_grps(cor_mat, ngrp, size_grp[2], type="col")
    group_loc <- list(group_loc1, group_loc2)

    ## Asssign inner bins
    group_loc_adj1 <- .get_grps_inner(group_loc1)
    group_loc_adj2 <- .get_grps_inner(group_loc2)
    group_loc_adj <- list(group_loc_adj1, group_loc_adj2)

    ## Get rank for each running bin
    cor_ref <- cor_mat[group_loc1[[ref_grp]], group_loc2[[ref_grp]]]
    rank_bin <- .get_bin_rank(cor_mat, group_loc, group_loc_adj, cor_ref)

    ## Transform rank to cor_adj
    cor_est <- .est_cor(rank_bin, cor_ref) # matrix

    idx1 <- order(ref_vec_rearranged1)
    idx2 <- order(ref_vec_rearranged2)
    cor_est <- cor_est[idx1, idx2]
    dimnames(cor_est) <- original.names

    return(cor_est)
}

# ==========================================================================
# Function for running in DReAmiR
# ==========================================================================

.GetAveExpbySubtype <- function(expr.mat, group.label) {
	group.label <- as.factor(group.label)
	ave.list <- lapply(levels(group.label), function(x) rowMeans(expr.mat[, group.label == x]))
	names(ave.list) <- levels(group.label)
	return(ave.list)
}

#' Given a list of correlation matrcies, return the SpQN normalized matrices
#' 
#' @param mrna.mat An expression matrix for mRNAs
#' @param mirna.mat An expression matrix for miRNAs
#' @param group.label A factor object for subtype/group
#' @param cor.list A list object from BuildCorList()
#' @param numCores An integer indicating the number of cores to use
#' @param ... Arguments passed to normalize_correlation()
#' @return A list of normalized correlation matrices
#' 
#' @import parallel
#' @export
RunGroupSpQN <- function(mrna.mat, mirna.mat, group.label, cor.list, numCores=1, ...) {

	# Check names
	if (!identical(levels(group.label), names(cor.list))) {
		stop("Levels of cor.list do not match with group labels\n")
	}

	# Check dimensions between correlation matrix and expression matrix
	if (nrow(mrna.mat) != nrow(cor.list[[1]]) | nrow(mirna.mat) != ncol(cor.list[[1]])) {
		stop("Dimensions of correlation matrix do not match with expression matrices\n")
	}

	# Check dimensions between expression matrices
	if (ncol(mrna.mat) != ncol(mirna.mat)) {
		stop("Dimensions of expression matrices do not match\n")
	}

	mrna_ave_exp <- .GetAveExpbySubtype(expr.mat=mrna.mat, group.label=group.label)

	mirna_ave_exp <- .GetAveExpbySubtype(expr.mat=mirna.mat, group.label=group.label)

	perm.spqn.cor.mat.list <- mclapply(levels(as.factor(group.label)), function(x) 
														normalize_correlation(cor.list[[x]], 
														ave_exp1=mrna_ave_exp[[x]], 
														ave_exp2=mirna_ave_exp[[x]], 
														...),
														mc.cores = numCores)

	return(perm.spqn.cor.mat.list)
}