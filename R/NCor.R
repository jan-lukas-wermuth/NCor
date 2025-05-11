#' Proper Correlations for Nominal Random Variables
#'
#' `NCor()` computes proper correlation coefficients for nominal random variables together with optional confidence intervals and independence tests.
#'
#' @param X an n x 1 character vector or a contingency table with at least one nominal random variable.
#' @param Y NULL (default) or a n x 1 character or numeric vector.
#' @param alpha confidence level for the returned confidence interval. FALSE yields the coefficient without confidence intervals.
#' @param digits only applicable if a contingency table of probabilities is supplied. Each probability gets multiplied by 10^digits to obtain an artificial sample/population size.
#' @param nominal determines whether the rows ("r"), the columns ("c") or the rows and the columns ("rc") of the resulting contingency table should be treated as nominal. If X and Y are supplied separately, X refers to the rows and Y to the columns.
#' @param method determines the coefficient to be used. Default is Goodman-Kruskal's gamma \insertCite{Goodman1954}{NCor}.
#'
#' @return The value of the coefficient.
#' @import Rdpack
#' @import foreach
#' @import rstatix
#' @import arrangements
#' @import doParallel
#' @importFrom magrittr "%>%"
#' @import class
#' @export
#'
#' @references \insertRef{Goodman1954}{NCor}
#'
#' @examples
#' tab <- matrix(c(38, 0, 75, 16, 0, 43, 86, 84, 60), ncol = 3)
#' colnames(tab) <- c("Afghanistan", "Somalia", "Iran")
#' rownames(tab) <- c("Islam", "Judaism", "Christianity")
#' NCor(tab)
NCor <- function(X, Y = NULL, alpha = 0.1, digits = 5, nominal = "rc", method = "gamma"){
  if (is.null(Y)){
    ContTable <- X
  } else ContTable <- table(X, Y)
  if (any(ContTable %% 1 != 0)){
    round(ContTable, digits = digits)
    ContTable <- 10 ^ digits * ContTable
  }
  # Start cluster for parallel computing
  cl <- parallel::makeCluster(2, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl)) # Need to stop the parallel computing

  if (nominal == "r"){ # case 1: Rows are nominal
    dim_r <- nrow(ContTable)
    rows <- factorial(dim_r)
    results <- foreach::foreach(iperm_r = arrangements::ipermutations(dim_r, dim_r), i = 1:rows, .combine = 'c') %dopar% {
      rownames(ContTable) <- iperm_r
      colnames(ContTable) <- 1:ncol(ContTable)
      cases <- rstatix::counts_to_cases(ContTable)
      round(RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), method = method)$Gamma, digits = 10)
    }
    iperm_r <- arrangements::ipermutations(dim_r, dim_r)
    permutations <- iperm_r$collect()
    permutations_index <- which(results == max(results))[1]
    if (permutations_index[1] == 1){
      final_perm_row <- iperm_r$getnext(permutations_index[1])
    } else final_perm_row <- iperm_r$getnext(permutations_index[1])[permutations_index[1],]
    rownames(ContTable) <- final_perm_row
    colnames(ContTable) <- 1:ncol(ContTable)
    cases <- rstatix::counts_to_cases(ContTable)
    finalres <- RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), method = method)
    return(list(finalres, allres <- results))
  }
  if (nominal == "c"){ # case 2: Cols are nominal
    dim_c <- ncol(ContTable)
    cols <- factorial(dim_c)
    results <- foreach::foreach(iperm_c = arrangements::ipermutations(dim_c, dim_c), i = 1:cols, .combine = 'c') %dopar% {
      colnames(ContTable) <- iperm_c
      rownames(ContTable) <- 1:nrow(ContTable)
      cases <- rstatix::counts_to_cases(ContTable)
      round(RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), alpha = alpha, method = method)$Gamma, digits = 10)
    }
    iperm_c <- arrangements::ipermutations(dim_c, dim_c)
    permutations <- iperm_c$collect()
    permutations_index <- which(results == max(results))[1]
    if (permutations_index[2] == 1){
      final_perm_col <- iperm_c$getnext(permutations_index[2])
    } else final_perm_col <- iperm_c$getnext(permutations_index[2])[permutations_index[2],]
    colnames(ContTable) <- final_perm_col
    rownames(ContTable) <- 1:nrow(ContTable)
    cases <- rstatix::counts_to_cases(ContTable)
    finalres <- RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), alpha = alpha, method = method)
    return(list(finalres, allres <- results))
  }
  if (nominal == "rc"){ # case 3: Rows and cols are nominal
    if (ncol(ContTable) < nrow(ContTable)){ContTable <- t(ContTable)} # Flip longer side to columns to exploit parallel computing
    dim_r <- nrow(ContTable)
    dim_c <- ncol(ContTable)
    cols <- factorial(dim_c)
    rows <- factorial(dim_r)
    iperm_r <- arrangements::ipermutations(dim_r, dim_r)
    results <- matrix(NA, nrow = rows, ncol = cols)
    for (i in 1:rows) {
      rownames(ContTable) <- iperm_r$getnext()
      results[i,] <- foreach(iperm_c = arrangements::ipermutations(dim_c, dim_c), i = 1:cols, .combine = 'c') %dopar% {
        colnames(ContTable) <- iperm_c
        cases <- rstatix::counts_to_cases(ContTable)
        round(RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), method = method)$Gamma, digits = 10)
      }
    }
    iperm_r <- arrangements::ipermutations(dim_r, dim_r)
    iperm_c <- arrangements::ipermutations(dim_c, dim_c)
    permutations_row <- iperm_r$collect()
    permutations_col <- iperm_c$collect()
    permutations_index <- which(results == max(results), arr.ind = TRUE)[1,]
    if (permutations_index[1] == 1){
      final_perm_row <- iperm_r$getnext(permutations_index[1])
    } else final_perm_row <- iperm_r$getnext(permutations_index[1])[permutations_index[1],]
    if (permutations_index[2] == 1){
      final_perm_col <- iperm_c$getnext(permutations_index[2])
    } else final_perm_col <- iperm_c$getnext(permutations_index[2])[permutations_index[2],]
    rownames(ContTable) <- final_perm_row
    colnames(ContTable) <- final_perm_col
    cases <- rstatix::counts_to_cases(ContTable)
    finalres <- RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), alpha = alpha, method = method)
    return(list(finalres, allres <- results))
  }
}

