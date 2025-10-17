#' Proper Correlations for Nominal Random Variables \insertCite{wermuth2025nominal}{NCor}
#'
#' `NCor()` computes Goodman-Kruskal's gamma \insertCite{Goodman1954}{NCor} for nominal random variables together with an optional confidence interval and a p-value for the independence test.
#'
#' @param X an n x 1 character vector or a contingency table with at least one nominal random variable.
#' @param Y NULL (default) or a n x 1 character or numeric vector.
#' @param alpha a numeric value specifying the significance level. The confidence level will be 1 - alpha.
#' @param digits only applicable if a contingency table of probabilities is supplied. Each probability gets multiplied by 10^digits to obtain an artificial sample/population size.
#' @param nominal determines whether the rows ("r") or the rows and the columns ("rc") of the resulting contingency table should be treated as nominal. If X and Y are supplied separately, X refers to the rows and Y to the columns.
#' @param CIs Boolean variable determining whether confidence intervals shall be reported or not.
#' @param Test Boolean variable determining whether the p-value for the null hypothesis of independence between X and Y shall be reported or not.
#'
#' @return The value of the coefficient, its confidence interval and the p-value for the independence test.
#' @useDynLib NCor, .registration = TRUE
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom magrittr %>%
#' @importFrom rstatix counts_to_cases
#' @importFrom parallel detectCores
#' @export
#'
#'
#' @references
#'  - \insertRef{Goodman1954}{NCor}
#'  - \insertRef{wermuth2025nominal}{NCor}
#'
#' @examples
#' tab <- matrix(c(38, 0, 75, 16, 0, 43, 86, 84, 60), ncol = 3)
#' colnames(tab) <- c("Afghanistan", "Somalia", "Iran")
#' rownames(tab) <- c("Islam", "Judaism", "Christianity")
#' NCor(tab)
NCor <- function(X, Y = NULL, alpha = 0.1, digits = 5, nominal = "rc", CIs = FALSE, Test = FALSE){
  if (nominal == "r"){ # case 1: Rows (X) are nominal
    if (is.null(Y)){stop("Please insert valid values for Y!")}
    if (length(X) != length(Y)){stop("Please provide X and Y vectors of equal length!")}
    n <- length(X)
    X_fac <- factor(X)
    X_vals <- levels(X_fac)
    X_int <- as.integer(X_fac)
    K <- length(X_vals)
    Y <- as.numeric(Y)
    results <- max_gamma_nominal_cont_dp_parallel(X_int, Y, K, tol = 1e-12, num_threads = detectCores() - 1)
    estim <- results$value
    order <- X_vals[results$order]
    if (isTRUE(Test)){
      rows <- factorial(K)
      CT <- table(X_int, Y)
      RcppParallel::setThreadOptions(numThreads = detectCores() - 1)
      out <- build_A1_SigmaU1_lehmer(CT)
      Sigma1 <- tcrossprod(out$A1 %*% out$Sigma_U1, out$A1)
      p_val <- 1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(sqrt(n)*estim, rows), mean = rep(0, rows), sigma = Sigma1)[1]
      if (isTRUE(CIs)){
        X_num <- as.numeric(factor(X, levels = order))
        finalres_CIs <- RCor::RCor(X_num, Y, alpha = alpha, method = "gamma", Fisher = FALSE)[1:3]
        finalres <- list(Gamma = finalres_CIs[[1]], CI_lower = max(finalres_CIs[2], 0), CI_upper = min(finalres_CIs[3], 1), PValue = p_val, Order = order)
      } else if (isFALSE(CIs)){
        finalres <- list(Gamma = estim, Order = order, PValue = p_val)
      } else {stop("Please insert either TRUE or FALSE for the variable CIs!")}
    } else if (isFALSE(Test)){
      if (isTRUE(CIs)){
        X_num <- as.numeric(factor(X, levels = order))
        finalres_CIs <- RCor::RCor(X_num, Y, alpha = alpha, method = "gamma", Fisher = FALSE)[1:3]
        finalres <- list(Gamma = finalres_CIs[[1]], CI_lower = max(finalres_CIs[2], 0), CI_upper = min(finalres_CIs[3], 1), Order = order)
      }
      else if (isFALSE(CIs)){
        finalres <- list(Gamma = estim, Order = order)
      } else {stop("Please insert either TRUE or FALSE for the variable CIs!")}
    } else {stop("Please insert either TRUE or FALSE for the variable Test!")}
    return(finalres)
  }
  if (nominal == "rc"){ # case 2: Rows and cols are nominal
    if (is.null(Y)){
      ContTable <- X
    } else ContTable <- table(X, Y)
    n <- sum(ContTable)
    if (any(ContTable %% 1 != 0)){
      round(ContTable, digits = digits)
      ContTable <- 10 ^ digits * ContTable
    }
    if (nrow(ContTable) < ncol(ContTable)){ContTable <- t(ContTable)} # Flip longer side to columns to exploit parallel computing
    results <- fill_matrix_max_tbb_lehmer_cpp(ContTable/n)
    estim <- results$value
    table <- results$table
    dim_r <- nrow(ContTable)
    dim_c <- ncol(ContTable)
    rownames(table) <- 1:dim_r
    colnames(table) <- 1:dim_c
    row_perm <- results$row_perm
    col_perm <- results$col_perm
    K <- factorial(dim_r)
    L <- factorial(dim_c)
    if (isTRUE(Test)){
      RcppParallel::setThreadOptions(numThreads = detectCores() - 1)
      out <- build_A2_SigmaU2_lehmer(ContTable)
      Sigma2 <- tcrossprod(out$A2 %*% out$Sigma_U2, out$A2)
      p_val <- 1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(sqrt(n)*estim, K * L), mean = rep(0, K * L), sigma = Sigma2)[1]
      if (isTRUE(CIs)){
        cases <- rstatix::counts_to_cases(table*n)
        finalres_CIs <- RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), alpha = alpha, method = "gamma", Fisher = FALSE)[1:3]
        finalres <- list(Gamma = finalres_CIs[[1]], CI_lower = max(finalres_CIs[2], 0), CI_upper = min(finalres_CIs[3], 1), PValue = p_val, Table = table, Row_perm = row_perm, Col_perm = col_perm)
      } else if (isFALSE(CIs)){
        finalres <- list(Gamma = estim, PValue = p_val, Table = table, Row_perm = row_perm, Col_perm = col_perm)
      } else {stop("Please insert either TRUE or FALSE for the variable CIs!")}
    } else if (isFALSE(Test)){
      if (isTRUE(CIs)){
        cases <- rstatix::counts_to_cases(table*n)
        finalres_CIs <- RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), alpha = alpha, method = "gamma", Fisher = FALSE)[1:3]
        finalres <- list(Gamma = finalres_CIs[[1]], CI_lower = max(finalres_CIs[2], 0), CI_upper = min(finalres_CIs[3], 1), Table = table, Row_perm = row_perm, Col_perm = col_perm)
      } else if (isFALSE(CIs)){
        finalres <- list(Gamma = estim, Table = table, Row_perm = row_perm, Col_perm = col_perm)
      } else {stop("Please insert either TRUE or FALSE for the variable CIs!")}
    } else {stop("Please insert either TRUE or FALSE for the variable Test!")}
    return(finalres)
  }
}
