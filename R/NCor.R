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
#' @import Rdpack
#' @import foreach
#' @import rstatix
#' @import arrangements
#' @import RCor
#' @import doParallel
#' @importFrom magrittr "%>%"
#' @import class
#' @export
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
  if (is.null(Y)){
    ContTable <- X
  } else ContTable <- table(X, Y)
  if (any(ContTable %% 1 != 0)){
    round(ContTable, digits = digits)
    ContTable <- 10 ^ digits * ContTable
  }
  n <- sum(ContTable)

  if (nominal == "r"){ # case 1: Rows are nominal
    # Start cluster for parallel computing
    cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl)) # Need to stop the parallel computing
    dim_r <- nrow(ContTable)
    rows <- factorial(dim_r)
    results <- foreach::foreach(iperm_r = arrangements::ipermutations(dim_r, dim_r), i = 1:rows, .combine = 'c') %dopar% {
      rownames(ContTable) <- iperm_r
      cases <- rstatix::counts_to_cases(ContTable)
      gamma_info <- DescTools:::.DoCount(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])))
      (gamma_info$C - gamma_info$D) / (gamma_info$C + gamma_info$D)
    }
    if (isTRUE(Test)){
      iperm_r <- arrangements::ipermutations(dim_r, dim_r)
      permutations <- iperm_r$collect()
      A1 <- matrix(0, nrow = rows, ncol = 2 * rows)
      Sigma_U1 <- matrix(NA, nrow = 2 * rows, ncol = 2 * rows)
      for (i in 1:rows) {
        rownames(ContTable) <- permutations[i,]
        cases <- rstatix::counts_to_cases(ContTable)
        X <- as.numeric(as.vector.factor(cases[,1]))
        Y <- as.numeric(as.vector.factor(cases[,2]))
        tau_info <- DescTools:::.DoCount(X, Y)
        tau <- (tau_info$C - tau_info$D) / choose(n, 2)
        X_TieProb <- sum((table(X)/length(X))^2)
        Y_TieProb <- sum((table(Y)/length(Y))^2)
        XY_TieProb <- sum((table(X, Y)/length(X))^2)
        tie_prob <- X_TieProb + Y_TieProb - XY_TieProb
        A1[i,2*i-1] <- 1 / (1 - tie_prob)
        A1[i,2*i] <- tau / (1 - tie_prob)^2
        G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
        G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
        G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
        x_eq <- Vectorize(function(x_val) mean(X == x_val))
        y_eq <- Vectorize(function(y_val) mean(Y == y_val))
        x_eq_y_eq <- Vectorize(function(x_val, y_val) mean(X == x_val & Y == y_val))
        # Calculate Marc's variance estimator
        G_XYXY <- G_XY(X, Y)
        G_XX <- G_X(X)
        G_YY <- G_Y(Y)
        x_eqX <- x_eq(X)
        y_eqY <- y_eq(Y)
        x_eq_y_eqXY <- x_eq_y_eq(X, Y)
        Sigma_U1[2*i-1,2*i-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau)^2)
        Sigma_U1[2*i,2*i] <- 4 * mean((x_eqX + y_eqY - x_eq_y_eqXY - tie_prob)^2)
        Sigma_U1[2*i,2*i-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (x_eqX + y_eqY - x_eq_y_eqXY - tie_prob))
        if (i < rows){
          for (j in (i+1):rows) {
            rownames(ContTable) <- permutations[j,]
            cases <- rstatix::counts_to_cases(ContTable)
            X_j <- as.numeric(as.vector.factor(cases[,1]))
            Y_j <- as.numeric(as.vector.factor(cases[,2]))
            tau_info_j <- DescTools:::.DoCount(X_j, Y_j)
            tau_j <- (tau_info_j$C - tau_info_j$D) / choose(n, 2)
            X_TieProb_j <- sum((table(X_j)/length(X_j))^2)
            Y_TieProb_j <- sum((table(Y_j)/length(Y_j))^2)
            XY_TieProb_j <- sum((table(X_j, Y_j)/length(X_j))^2)
            tie_prob_j <- X_TieProb_j + Y_TieProb_j - XY_TieProb_j
            G_XY_j <- Vectorize(function(x_val, y_val) (mean(X_j <= x_val & Y_j <= y_val) + mean(X_j <= x_val & Y_j < y_val) + mean(X_j < x_val & Y_j <= y_val) + mean(X_j < x_val & Y_j < y_val)) / 4)
            G_X_j <- Vectorize(function(x_val) (mean(X_j < x_val) + mean(X_j <= x_val)) / 2)
            G_Y_j <- Vectorize(function(y_val) (mean(Y_j < y_val) + mean(Y_j <= y_val)) / 2)
            x_eq_j <- Vectorize(function(x_val) mean(X_j == x_val))
            y_eq_j <- Vectorize(function(y_val) mean(Y_j == y_val))
            x_eq_y_eq_j <- Vectorize(function(x_val, y_val) mean(X_j == x_val & Y_j == y_val))
            # Calculate Marc's variance estimator
            G_XYXY_j <- G_XY_j(X_j, Y_j)
            G_XX_j <- G_X_j(X_j)
            G_YY_j <- G_Y_j(Y_j)
            x_eqX_j <- x_eq_j(X_j)
            y_eqY_j <- y_eq_j(Y_j)
            x_eq_y_eqXY_j <- x_eq_y_eq_j(X_j, Y_j)
            Sigma_U1[2*j,2*i-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (x_eqX_j + y_eqY_j - x_eq_y_eqXY_j - tie_prob_j))
            Sigma_U1[2*j-1,2*i-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (4 * G_XYXY_j - 2 * (G_XX_j + G_YY_j) + 1 - tau_j))
            Sigma_U1[2*j,2*i] <- 4 * mean((x_eqX + y_eqY - x_eq_y_eqXY - tie_prob) * (x_eqX_j + y_eqY_j - x_eq_y_eqXY_j - tie_prob_j))
            Sigma_U1[2*j-1,2*i] <- 4 * mean((x_eqX + y_eqY - x_eq_y_eqXY - tie_prob) * (4 * G_XYXY_j - 2 * (G_XX_j + G_YY_j) + 1 - tau_j))
          }
        }
      }
      Sigma_U1 <- Matrix::forceSymmetric(Sigma_U1, uplo = "L")
      Sigma1 <- as.matrix(A1 %*% Sigma_U1 %*% t(A1))
      p_val <- 1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(sqrt(n)*max(results), rows), mean = rep(0, rows), sigma = Sigma1)[1]
      if (isTRUE(CIs)){
        permutations_index <- which(results == max(results))[1]
        if (permutations_index == 1){
          final_perm_row <- iperm_r$getnext(permutations_index)
        } else {final_perm_row <- iperm_r$getnext(permutations_index)[permutations_index,]}
        rownames(ContTable) <- final_perm_row
        cases <- rstatix::counts_to_cases(ContTable)
        finalres_CIs <- RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), alpha = alpha, method = "gamma", Fisher = FALSE)[1:3]
        finalres_CIs <- dplyr::tribble(~Gamma, ~CI_lower, ~CI_upper,
                                       finalres_CIs[[1]], max(finalres_CIs[2], 0), min(finalres_CIs[3], 1))
        finalres <- finalres_CIs %>% mutate(PValue = p_val)
      } else if (isFALSE(CIs)){
        finalres <- dplyr::tribble(~Gamma, ~PValue,
                                   max(results), p_val)
      } else {stop("Please insert either TRUE or FALSE for the variable CIs!")}
    } else if (isFALSE(Test)){
      if (isTRUE(CIs)){
        iperm_r <- arrangements::ipermutations(dim_r, dim_r)
        permutations <- iperm_r$collect()
        permutations_index <- which(results == max(results))[1]
        if (permutations_index == 1){
          final_perm_row <- iperm_r$getnext(permutations_index)
        } else {final_perm_row <- iperm_r$getnext(permutations_index)[permutations_index,]}
        rownames(ContTable) <- final_perm_row
        cases <- rstatix::counts_to_cases(ContTable)
        finalres <- RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), alpha = alpha, method = "gamma", Fisher = FALSE)[1:3]
        finalres <- dplyr::tribble(~Gamma, ~CI_lower, ~CI_upper,
                                   finalres[[1]], max(finalres[2], 0), min(finalres[3], 1))
      }
      else if (isFALSE(CIs)){
        finalres <- max(results)
      } else {stop("Please insert either TRUE or FALSE for the variable CIs!")}
    } else {stop("Please insert either TRUE or FALSE for the variable Test!")}
    return(list(finalres, allres <- results))
  }
  if (nominal == "rc"){ # case 2: Rows and cols are nominal
    if (ncol(ContTable) < nrow(ContTable)){ContTable <- t(ContTable)} # Flip longer side to columns to exploit parallel computing
    dim_r <- nrow(ContTable)
    dim_c <- ncol(ContTable)
    cols <- factorial(dim_c)
    rows <- factorial(dim_r)
    iperm_r <- arrangements::ipermutations(dim_r, dim_r)$collect()
    iperm_c <- arrangements::ipermutations(dim_c, dim_c)$collect()
    results <- matrix(NA, nrow = rows, ncol = cols)
    for (i in 1:rows) {
      for (j in 1:cols) {
        rownames(ContTable) <- iperm_r[i,]
        colnames(ContTable) <- iperm_c[j,]
        ContTable_2 <- ContTable[order(as.integer(rownames(ContTable))),order(as.integer(colnames(ContTable))),drop = FALSE] / n
        A1 <- sign(row(diag(dim_r)) - col(diag(dim_r)))
        B1 <- sign(row(diag(dim_c)) - col(diag(dim_c)))
        A2 <- 1 - diag(dim_r)
        B2 <- 1 - diag(dim_c)
        Q <- sum(ContTable_2 * (A1 %*% ContTable_2 %*% t(B1)))   # = trace(t(H2) %*% A %*% H1 %*% t(B))
        N <- sum(ContTable_2 * (A2 %*% ContTable_2 %*% t(B2)))
        results[i,j] <- Q/N
      }
    }
    if (isTRUE(Test)){
      iperm_r <- arrangements::ipermutations(dim_r, dim_r)
      iperm_c <- arrangements::ipermutations(dim_c, dim_c)
      permutations_row <- iperm_r$collect()
      permutations_col <- iperm_c$collect()
      A2 <- matrix(0, nrow = rows * cols, ncol = 2 * rows * cols)
      Sigma_U2 <- matrix(NA, nrow = 2 * rows * cols, ncol = 2 * rows * cols)
      for (i in 1:rows){
        k <- (i - 1) * cols
        for (j in 1:cols) {
          rownames(ContTable) <- permutations_row[i,]
          colnames(ContTable) <- permutations_col[j,]
          cases <- rstatix::counts_to_cases(ContTable)
          X <- as.numeric(as.vector.factor(cases[,1]))
          Y <- as.numeric(as.vector.factor(cases[,2]))
          tau_info <- DescTools:::.DoCount(X, Y)
          tau <- (tau_info$C - tau_info$D) / choose(n, 2)
          X_TieProb <- sum((table(X)/length(X))^2)
          Y_TieProb <- sum((table(Y)/length(Y))^2)
          XY_TieProb <- sum((table(X, Y)/length(X))^2)
          tie_prob <- X_TieProb + Y_TieProb - XY_TieProb
          A2[j+k,2*(j+k)-1] <- 1 / (1 - tie_prob)
          A2[j+k,2*(j+k)] <- tau / (1 - tie_prob)^2
          G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
          G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
          G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
          x_eq <- Vectorize(function(x_val) mean(X == x_val))
          y_eq <- Vectorize(function(y_val) mean(Y == y_val))
          x_eq_y_eq <- Vectorize(function(x_val, y_val) mean(X == x_val & Y == y_val))
          # Calculate Marc's variance estimator
          G_XYXY <- G_XY(X, Y)
          G_XX <- G_X(X)
          G_YY <- G_Y(Y)
          x_eqX <- x_eq(X)
          y_eqY <- y_eq(Y)
          x_eq_y_eqXY <- x_eq_y_eq(X, Y)
          Sigma_U2[2*(j+k)-1,2*(j+k)-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau)^2)
          Sigma_U2[2*(j+k),2*(j+k)] <- 4 * mean((x_eqX + y_eqY - x_eq_y_eqXY - tie_prob)^2)
          Sigma_U2[2*(j+k),2*(j+k)-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (x_eqX + y_eqY - x_eq_y_eqXY - tie_prob))
          # Calculate off-block-diag elements of Sigma_U2 that use a different row permutation
          if (i < rows){
            p <- 0
            for (m in (i+1):rows){
              for (o in 1:cols){
                rownames(ContTable) <- permutations_row[m,]
                colnames(ContTable) <- permutations_col[o,]
                cases <- rstatix::counts_to_cases(ContTable)
                X_mo <- as.numeric(as.vector.factor(cases[,1]))
                Y_mo <- as.numeric(as.vector.factor(cases[,2]))
                tau_info_mo <- DescTools:::.DoCount(X_mo, Y_mo)
                tau_mo <- (tau_info_mo$C - tau_info_mo$D) / choose(n, 2)
                X_TieProb_mo <- sum((table(X_mo)/length(X_mo))^2)
                Y_TieProb_mo <- sum((table(Y_mo)/length(Y_mo))^2)
                XY_TieProb_mo <- sum((table(X_mo, Y_mo)/length(X_mo))^2)
                tie_prob_mo <- X_TieProb_mo + Y_TieProb_mo - XY_TieProb_mo
                G_XY_mo <- Vectorize(function(x_val, y_val) (mean(X_mo <= x_val & Y_mo <= y_val) + mean(X_mo <= x_val & Y_mo < y_val) + mean(X_mo < x_val & Y_mo <= y_val) + mean(X_mo < x_val & Y_mo < y_val)) / 4)
                G_X_mo <- Vectorize(function(x_val) (mean(X_mo < x_val) + mean(X_mo <= x_val)) / 2)
                G_Y_mo <- Vectorize(function(y_val) (mean(Y_mo < y_val) + mean(Y_mo <= y_val)) / 2)
                x_eq_mo <- Vectorize(function(x_val) mean(X_mo == x_val))
                y_eq_mo <- Vectorize(function(y_val) mean(Y_mo == y_val))
                x_eq_y_eq_mo <- Vectorize(function(x_val, y_val) mean(X_mo == x_val & Y_mo == y_val))
                # Calculate Marc's variance estimator
                G_XYXY_mo <- G_XY_mo(X_mo, Y_mo)
                G_XX_mo <- G_X_mo(X_mo)
                G_YY_mo <- G_Y_mo(Y_mo)
                x_eqX_mo <- x_eq_mo(X_mo)
                y_eqY_mo <- y_eq_mo(Y_mo)
                x_eq_y_eqXY_mo <- x_eq_y_eq_mo(X_mo, Y_mo)
                Sigma_U2[2*(k+cols+o+p),2*(j+k)-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (x_eqX_mo + y_eqY_mo - x_eq_y_eqXY_mo - tie_prob_mo))
                Sigma_U2[2*(k+cols+o+p)-1,2*(j+k)-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (4 * G_XYXY_mo - 2 * (G_XX_mo + G_YY_mo) + 1 - tau_mo))
                Sigma_U2[2*(k+cols+o+p),2*(j+k)] <- 4 * mean((x_eqX + y_eqY - x_eq_y_eqXY - tie_prob) * (x_eqX_mo + y_eqY_mo - x_eq_y_eqXY_mo - tie_prob_mo))
                Sigma_U2[2*(k+cols+o+p)-1,2*(j+k)] <- 4 * mean((x_eqX + y_eqY - x_eq_y_eqXY - tie_prob) * (4 * G_XYXY_mo - 2 * (G_XX_mo + G_YY_mo) + 1 - tau_mo))
              }
              p <- (m - i) * cols
            }
          }
          if (j < cols){
            # Calculate off-block-diag elements of Sigma_U2 that use the same row permutation
            for (l in (j+1):cols){
              rownames(ContTable) <- permutations_row[i,]
              colnames(ContTable) <- permutations_col[l,]
              cases <- rstatix::counts_to_cases(ContTable)
              X_l <- as.numeric(as.vector.factor(cases[,1]))
              Y_l <- as.numeric(as.vector.factor(cases[,2]))
              tau_info_l <- DescTools:::.DoCount(X_l, Y_l)
              tau_l <- (tau_info_l$C - tau_info_l$D) / choose(n, 2)
              X_TieProb_l <- sum((table(X_l)/length(X_l))^2)
              Y_TieProb_l <- sum((table(Y_l)/length(Y_l))^2)
              XY_TieProb_l <- sum((table(X_l, Y_l)/length(X_l))^2)
              tie_prob_l <- X_TieProb_l + Y_TieProb_l - XY_TieProb_l
              G_XY_l <- Vectorize(function(x_val, y_val) (mean(X_l <= x_val & Y_l <= y_val) + mean(X_l <= x_val & Y_l < y_val) + mean(X_l < x_val & Y_l <= y_val) + mean(X_l < x_val & Y_l < y_val)) / 4)
              G_X_l <- Vectorize(function(x_val) (mean(X_l < x_val) + mean(X_l <= x_val)) / 2)
              G_Y_l <- Vectorize(function(y_val) (mean(Y_l < y_val) + mean(Y_l <= y_val)) / 2)
              x_eq_l <- Vectorize(function(x_val) mean(X_l == x_val))
              y_eq_l <- Vectorize(function(y_val) mean(Y_l == y_val))
              x_eq_y_eq_l <- Vectorize(function(x_val, y_val) mean(X_l == x_val & Y_l == y_val))
              # Calculate Marc's variance estimator
              G_XYXY_l <- G_XY_l(X_l, Y_l)
              G_XX_l <- G_X_l(X_l)
              G_YY_l <- G_Y_l(Y_l)
              x_eqX_l <- x_eq_l(X_l)
              y_eqY_l <- y_eq_l(Y_l)
              x_eq_y_eqXY_l <- x_eq_y_eq_l(X_l, Y_l)
              Sigma_U2[2*(k+l),2*(j+k)-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (x_eqX_l + y_eqY_l - x_eq_y_eqXY_l - tie_prob_l))
              Sigma_U2[2*(k+l)-1,2*(j+k)-1] <- 4 * mean((4 * G_XYXY - 2 * (G_XX + G_YY) + 1 - tau) * (4 * G_XYXY_l - 2 * (G_XX_l + G_YY_l) + 1 - tau_l))
              Sigma_U2[2*(k+l),2*(j+k)] <- 4 * mean((x_eqX + y_eqY - x_eq_y_eqXY - tie_prob) * (x_eqX_l + y_eqY_l - x_eq_y_eqXY_l - tie_prob_l))
              Sigma_U2[2*(k+l)-1,2*(j+k)] <- 4 * mean((x_eqX + y_eqY - x_eq_y_eqXY - tie_prob) * (4 * G_XYXY_l - 2 * (G_XX_l + G_YY_l) + 1 - tau_l))
            }
          }
        }
      }
      Sigma_U2 <- Matrix::forceSymmetric(Sigma_U2, uplo = "L")
      Sigma2 <- as.matrix(A2 %*% Sigma_U2 %*% t(A2))
      p_val <- 1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(sqrt(n)*max(results), rows * cols), mean = rep(0, rows * cols), sigma = Sigma2)[1]
      if (isTRUE(CIs)){
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
        finalres_CIs <- RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), alpha = alpha, method = "gamma", Fisher = FALSE)[1:3]
        finalres_CIs <- dplyr::tribble(~Gamma, ~CI_lower, ~CI_upper,
                                       finalres_CIs[[1]], max(finalres_CIs[2], 0), min(finalres_CIs[3], 1))
        finalres <- finalres_CIs %>% mutate(PValue = p_val)
      } else if (isFALSE(CIs)){
        finalres <- dplyr::tribble(~Gamma, ~PValue,
                                   max(results), p_val)
      } else {stop("Please insert either TRUE or FALSE for the variable CIs!")}
    } else if (isFALSE(Test)){
      if (isTRUE(CIs)){
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
        finalres <- RCor::RCor(as.numeric(as.vector.factor(cases[,1])), as.numeric(as.vector.factor(cases[,2])), alpha = alpha, method = "gamma", Fisher = FALSE)[1:3]
        finalres <- dplyr::tribble(~Gamma, ~CI_lower, ~CI_upper,
                                   finalres[[1]], max(finalres[2], 0), min(finalres[3], 1))
      } else if (isFALSE(CIs)){
        finalres <- max(results)
      } else {stop("Please insert either TRUE or FALSE for the variable CIs!")}
    } else {stop("Please insert either TRUE or FALSE for the variable Test!")}
    return(list(finalres, allres <- results))
  }
}
