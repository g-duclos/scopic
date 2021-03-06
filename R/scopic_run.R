#' Function for fitting topic models to scRNA-Seq gene counts and assigning genes to gene sets OR cells to cell clusters
#'

#' This function fits topic models to scRNA-Seq gene counts and assigns genes to gene sets OR cells to cell clusters
#' @param counts.mat Matrix of scRNA-Seq gene counts (columns are cells, rows are genes)
#' @param type Options include "genes" (gene set topics) or "cells" (cell cluster topics)
#' @param n.topics Number of topics (transcriptional states)
#' @param n.starts Number of random starts (will select best model)
#' @param seeds Random seeds (number of seeds must equal n.starts)
#' @param max.q Max adjusted p-value (FDR q) considered for cell/gene-topic association (< 0.05 is recommended)
#' @param min.coef Minimum regression coefficient considered for cell/gene-topic association (> 1 is recommended)
#' @param min.spec Minimum topic specificity value considered for cell/gene-topic association (value between 0 and 1)
#' @param min.sim min.spec Minimum topic similarity value considered for cell/gene-topic association (value between 0 and 1)
#' @param model.return FALSE to return only results of model evaluation; TRUE to return results of model evaluation AND model output (topics & terms)
#' @param verbose TRUE to print status updates; FALSE to silence status updates
#' @return Topic model analysis of scRNA-Seq gene counts
#' @export
#' @examples
#' scopic_run(counts.mat = counts.mat, type = c("genes", "cells"), n.topics = 3, n.starts = 5, seeds = 1:length(n.starts), max.q = 0.05, min.coef = 1, min.spec = 0.1, min.sim = 0, model.return = FALSE, verbose = TRUE)
#

scopic_run <- function(
	counts.mat = counts.mat,
	type = c("genes", "cells"),
	n.topics = 3,
	n.starts = 5,
	seeds = 1:length(n.starts),
	max.q = 0.05,
	min.coef = 1,
	min.spec = 0.1,
	min.sim = 0,
	model.return = FALSE,
	verbose = FALSE) {

	#
	result.lda <- scopic_model(
		counts.mat = counts.mat,
		type = type,
		n.topics = n.topics,
		n.starts = n.starts,
		seeds = seeds,
		verbose = verbose)

	# Extract States Matrix from LDA output
	states <- topicmodels::posterior(result.lda)$topics

	# Negative Binomial Regression
	result.stats <- data.frame(row.names=colnames(counts.mat))
	result.stats[unlist(lapply(1:ncol(states), paste, c("p", "FDR q", "Coef", "Spec", "Sim")))] <- NA

	# Maximum number of model iterations
	max.it <- 100
	cat("Evaluate State-Feature Association \n")

	if (type == "genes") {
		result.stats <- data.frame(row.names=rownames(counts.mat))
		result.stats[unlist(lapply(1:ncol(states), paste, c("p", "FDR q", "Coef", "Spec", "Sim")))] <- NA
		for (k in 1:ncol(states)) {
			cat("State", k, sep=" ", "\n")
			for (i in 1:nrow(counts.mat)) {
			    if ((i %% 100) == 0) cat("State", k, i, sep=" ", "\r")
			    fit <- try(glm.nb(counts.mat[i,] ~ states[,k], maxit=max.it), silent=TRUE)
			    if (!inherits(fit, "try-error")) {
			        result.stats[i, paste(k, "p")] <- summary(fit)$coefficients[2,4]
			        result.stats[i, paste(k, "Coef")] <- summary(fit)$coefficients[2,1]
			    	}
				}
			}
		} else if (type == "cells") {
			result.stats <- data.frame(row.names=colnames(counts.mat))
			result.stats[unlist(lapply(1:ncol(states), paste, c("p", "FDR q", "Coef", "Spec", "Sim")))] <- NA
			for (k in 1:ncol(states)) {
				cat("State", k, sep=" ", "\n")
				for (i in 1:ncol(counts.mat)) {
		    		if ((i %% 100) == 0) cat("State", k, i, sep=" ", "\r")
		    		fit <- try(glm.nb(counts.mat[,i] ~ states[,k], maxit=max.it), silent=TRUE)
		    		if (!inherits(fit, "try-error")) {
		        		result.stats[i, paste(k, "p")] <- summary(fit)$coefficients[2,4]
		        		result.stats[i, paste(k, "Coef")] <- summary(fit)$coefficients[2,1]
		    		}
				}
			}
		}

	cat("Evaluation Complete \n")

	# FDR correction
	for (prefix in sub(" p$", "", grep(" p$", colnames(result.stats), value=TRUE))) {
	    result.stats[[paste(prefix, "FDR q")]] <- p.adjust(result.stats[[paste(prefix, "p")]], method="fdr")
	}

	# Calculate State-Specificity and add to results matrix
	result.spec <- t(sweep(topicmodels::posterior(result.lda)$terms, 2, apply(topicmodels::posterior(result.lda)$terms, 2, sum), "/"))
	for (k in 1:ncol(states)) {
		result.stats[, paste(k, "Spec")] <- result.spec[,k]
	}
	cat("Calculate State-Specificity \n")

	# Calculate relative expression values: counts divided by total counts per cell
	rel.exp.mat <- sweep(counts.mat, MARGIN=2, STATS=colSums(counts.mat), FUN="/")

	# Calculate State-Similarity and add to results matrix
	if (type == "genes") {
		for (i in 1:nrow(rel.exp.mat)) {
			for (k in 1:ncol(states)) {
				result.stats[i, paste(k, "Sim")] <- lsa::cosine(rel.exp.mat[i,], states[,k])
				}
			}
		} else if (type == "cells") {
			for (i in 1:ncol(rel.exp.mat)) {
				for (k in 1:ncol(states)) {
					result.stats[i, paste(k, "Sim")] <- lsa::cosine(rel.exp.mat[,i], states[,k])
				}
			}
		}
	cat("Calculate State-Similarity \n")

	# Feature Selection
	result.stats["Assignment"] <- NA

	# Extract FDR q-values and model coefficients from statistical analysis
	for (i in 1:nrow(result.stats)) {
		target.q <- result.stats[i,][grep("FDR",names(result.stats[i,]),value=T)]
		target.coef <- result.stats[i,][grep("Coef",names(result.stats[i,]),value=T)]
		target.spec <- result.stats[i,][grep("Spec",names(result.stats[i,]),value=T)]
		target.sim <- result.stats[i,][grep("Sim",names(result.stats[i,]),value=T)]

		target.sig <- Reduce(intersect, list(
			gsub("[^0-9]", "", names(target.q)[which(target.q < max.q)]),
			gsub("[^0-9]", "", names(target.coef)[which(target.coef > min.coef)]),
			gsub("[^0-9]", "", names(target.spec)[which(target.spec > min.spec)]),
			gsub("[^0-9]", "", names(target.sim)[which(target.sim > min.sim)])
			))

		#
		if (length(target.sig) == 0) {
			result.stats[i, "Assignment"] <- NA
			} else if (length(target.sig) == 1) {
				result.stats[i, "Assignment"] <- target.sig
				} else {
					result.stats[i, "Assignment"] <- gsub("[^0-9]", "", names(result.stats[i,][paste(target.sig, "FDR q")])[which(result.stats[i,][paste(target.sig, "FDR q")] == min(result.stats[i,][paste(target.sig, "FDR q")]) )] )
				}
			}

	cat("State Assignment Complete \n")	
	#
	if (model.return == TRUE) {
		#
		result.list <- list()
		result.list["stats"] <- list(result.stats)
		result.list["topics"] <- list(states)
		result.list["terms"] <- list(topicmodels::posterior(result.lda)$terms)
		#
		return(result.list)
		#
		} else {
			#
			return(result.stats)
			#
			}
}

