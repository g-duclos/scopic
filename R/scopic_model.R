#' Function for fitting topic models to scRNA-Seq gene counts
#'

#' This function fits topic models to scRNA-Seq gene counts
#' @param counts.mat Matrix of scRNA-Seq gene counts (columns are cells, rows are genes)
#' @param type Options include "genes" (gene set topics) or "cells" (cell cluster topics)
#' @param n.topics Number of topics (transcriptional states)
#' @param n.starts Number of random starts (will select best model)
#' @param seeds Random seeds (number of seeds must equal n.starts)
#' @param verbose TRUE to print status updates; FALSE to silence status updates
#' @return Topic model of scRNA-Seq gene counts
#' @export
#' @examples
#' scopic_model(counts.mat = counts.mat, type = c("genes", "cells"), n.topics = 3, n.starts = 5, seeds = length(n.starts), verbose = TRUE)
#

scopic_model <- function(
	counts.mat = counts.mat,
	type = c("genes", "cells"),
	n.topics = 3,
	n.starts = 5,
	seeds = 1:length(n.starts),
	verbose = TRUE) {

	#
	if (type == "genes") {
		type.out = "Gene Set"
	} else if (type == "cells") {
		type.out = "Cell Cluster"
	}

	# Parameters for LDA
	lda.param <- list(estimate.alpha = TRUE, estimate.beta = TRUE,
		verbose = 100, save = 0, keep = 0,
		seed = seeds, nstart = n.starts, best = TRUE)

	# Run LDA
	if (verbose == TRUE) {cat(paste("Build Model - ", type.out, " Topics: ", n.topics, sep=""), "\n")}
	if (type == "genes") {
		result.lda <- topicmodels::LDA(t(counts.mat), k=n.states, method="VEM", control=lda.param)
		} else if (type == "cells") {
			result.lda <- topicmodels::LDA(counts.mat, k=n.states, method="VEM", control=lda.param)
			}
	if (verbose == TRUE) {cat(paste("Complete Model - ", type.out, " Topics: ", n.topics, sep=""), "\n")}
	
	#
	return(result.lda)
	#
}


