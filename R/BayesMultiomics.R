#' BayesMultiomics.CV
#'
#' Run first and second stage models in *Hao Xue, Sounak Chakraborty & Tanujit Dey (2024), Bayesian Shrinkage Models for Integration and Analysis of Multiplatform High-dimensional Genomics Data*.
#' In the first stage, Expectation–Maximization Variable Selection (EMVS) is used to learn the regulating mechanism between epigenomics (e.g., gene methylation) and gene expression while considering functional gene annotations.
#' Genes are grouped based on the regulating mechanism learned in the first stage using the function, EMVS.
#' Then, a group-wise penalty is applied to select genes significantly associated with clinical outcomes while incorporating clinical features using Bayesian variable selection model with extended NEG priors in *Veronika Rockova & Edward I. George (2013): EMVS: The EM Approach to Bayesian Variable Selection, Journal of the American Statistical Association*.
#' K-fold Cross validation is performed to evaluate model performances. The coefficients from all models are averaged.
#'
#'
#' @param M DNA Methylation matrix, a matrix with rows corresponding to samples and columns corresponding to probes/sites.
#' @param G Gene expression level, a matrix with rows corresponding to samples and columns corresponding to genes.
#' @param grouping Gene grouping, a vector of group membership of each gene. This could be given by the user or obatined by performing DAVID Functional Classification of genes. See vignette("gene_grouping").
#' @param nu0 Parameter `nu0` for spike-and-slab Gaussian mixture prior on \eqn{\beta}
#' @param nu1 parameter `nu1` for spike-and-slab Gaussian mixture prior on \eqn{\beta}
#' @param lambda For the prior on \eqn{\sigma2}, an inverse gamma prior \eqn{\pi(\sigma2 | \gamma) = IG(\nu/2, \nu\lambda/2)}
#' @param a Parameter \eqn{a} for the beta prior \eqn{\pi(\theta) \propto \theta a-1(1-\theta)b-1}, \eqn{a, b > 0}
#' @param b Parameter \eqn{b} for the beta prior \eqn{\pi(\theta) \propto \theta a-1(1-\theta)b-1}, \eqn{a, b > 0}
#' @param EMVS_I Maximum number of iterations of EMVS
#' @param EMVS_thresh Threshold for convergence criterion for EMVS
#' @param Y clinical outcome
#' @param G gene expression level
#' @param C clinical features
#' @param a0 size sparsity
#' @param gstr prior, two options: "scale" or "1"
#' @param NEG_I number of maximum iterations for NEG_em
#' @param NEG_thresh convergence criterion fpr NEG_em
#' @param .mpmath (optional) function depends on mpmath package from python. pointer for mpmath package
#' @param n_fold number of folds in K-fold cross validation
#' @param random_seed random seed for splitting folds for cross validation
#' @param lower_R2 lower limit for R2 threshold, default = 0.2
#' @param upper_R2 upper limit for R2 threshold, default = 0.8
#' @param Delta Delta
#' @param transform_M transform methylation matrix `M`. Options: "L" = linear transformation, "Q" = quadratic transformation; "QS" = cubic spline transformation
#'
#'
#' @details
#' `gstr` will take two options "scale" or "1." if `gstr` == "scale" then g = 1/N^2 where N = number of genes
#'
#' @examples
#' # params
#' M <- GBM_data2$M
#' G <- GBM_data2$G
#' # fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
#' fp <- system.file("eg_fx_classification.txt", package="BayesMultiomics")
#' gene_grouping <- DAVIDGeneGrouping(eg_gene_symbols, fp)
#' Y <- GBM_data2$Y
#' C <- GBM_data2$C
#' a0 <- 0.1
#' gstr <- "scale"
#' Delta <- GBM_data2$Delta
#'
#' # run
#' multiomics_cv <- BayesMultiomics.CV(
#'   M=M,
#'   G=G,
#'   grouping=gene_grouping,
#'   Y=Y,
#'   C=C,
#'   a0=a0,
#'   gstr=gstr,
#'   Delta=Delta,
#'   n_fold=10,
#'   random_seed = 123,
#'   EMVS_I = 10,
#'   NEG_I=20,
#'   EMVS_thresh = 0.0001,
#'   NEG_thresh = 0.0001
#' )
#'
#' # plot
#' plot(multiomics_cv)
#'
#' @export
BayesMultiomics.CV <- function(
    M, G, grouping, nu0=0.5, nu1=10^3, nu=1, lambda=1, a=1, b=1, EMVS_I=10, EMVS_thresh=0.0001,
    Y, C, Delta, a0, gstr, n_fold=10, random_seed=NULL, NEG_I=10, NEG_thresh=0.0001,
    lower_R2=0.2, upper_R2=0.8, transform_M="L", .mpmath=setup_mpmath()
){

  # within MultiOmics, run EMVS
  EMVS_result <- .EMVS(M,
                       G,
                       grouping=grouping,
                       nu0 = nu0,
                       nu1 = nu1,
                       lambda = lambda,
                       a = a,
                       b = b,
                       I=EMVS_I,
                       thresh=EMVS_thresh,
                       transform_M = transform_M)

  # run 2nd stage
  N <- nrow(G)
  K <- ncol(G)
  J <- ncol(M)
  L <- ncol(C)
  R2 <- EMVS_result$R2

  # build Z matrix
  Zmatrix <- Zmat_builder(R2, G, lower_R2, upper_R2)

  # E-M loop step
  # run cross validation
  if(is.null(random_seed)) set.seed(random_seed)
  # folds <- caret::createFolds(1:N, k=n_fold)
  folds <- .folds(1:N, num_fold=n_fold)

  # initialize result objects
  selected_biomarkers <- c()
  res_cols <- c("n_selected_vars", "r2_train", "r2_test",
                "cindex_train", "cindex_test", "mse_train", "mse_test",
                colnames(C))

  res_table <- data.frame(matrix(0, n_fold, length(res_cols))) |>
    setNames(res_cols)

  G_coeffs <- data.frame(matrix(0, n_fold, ncol(G))) |>
    setNames(colnames(G))

  message("\n2nd stage modeling with a0 = ", a0, "and g = ", gstr)
  message("Running cross validation...")
  pb <- txtProgressBar(min=0, max=n_fold, initial=0, style=3)

  # loop
  for(jjj in 1:n_fold){
    test_indx <- folds[[jjj]]
    train_delta <- Delta[-test_indx]
    test_delta <- Delta[test_indx]

    NEG_output <- .NEG_em(Y=Y[-test_indx],
                          G=as.matrix(G[-test_indx,]),
                          C=as.matrix(C[-test_indx,]),
                          a0=a0,
                          gstr=gstr,
                          Zmatrix=Zmatrix,
                          I=NEG_I,
                          thresh=NEG_thresh,
                          .mpmath=.mpmath)

    # save result
    NEG_output$beta <- data.frame(NEG_output$beta)
    colnames(NEG_output$beta) <- c(colnames(G), colnames(C))
    estbeta <- NEG_output$beta[NEG_output$k, ]
    selected_biomarkers <- unique(c(selected_biomarkers,
                                    names(estbeta)[abs(estbeta) > 1e-5 & !names(estbeta) %in% colnames(C)]))
    # evaluate cross validation result
    pred_train <- cbind(G[-test_indx,], C[-test_indx,]) %*% as.numeric(estbeta)
    pred_test <- cbind(G[test_indx,], C[test_indx,]) %*% as.numeric(estbeta)
    res_table[jjj, "n_selected_vars"] <- sum(abs(estbeta) > 1e-5 & !names(estbeta) %in% colnames(C))
    res_table[jjj, "r2_train"] <- .rsq(pred_train, Y[-test_indx])
    res_table[jjj, "r2_test"] <- .rsq(pred_test, Y[test_indx])
    res_table[jjj, "cindex_train"] <- .cindx(pred_train, Y[-test_indx])
    res_table[jjj, "cindex_test"] <- .cindx(pred_test, Y[test_indx])
    res_table[jjj, "mse_train"] <- mean((pred_train - Y[-test_indx])^2)
    res_table[jjj, "mse_test"] <- mean((pred_test - Y[test_indx])^2)
    res_table[jjj, colnames(C)] <- estbeta[colnames(C)]
    G_coeffs[jjj, colnames(G)] <- estbeta[1, colnames(G)]

    setTxtProgressBar(pb, jjj)
  }

  # take the average
  performance_df <- data.frame(prior_a=a0,
                               prior_g=gstr,
                               as.list(colMeans(res_table)))
  performance_df2 <- performance_df[, ! colnames(performance_df) %in% colnames(C)] |>
    t() |>
    `colnames<-`(paste0("average of ", n_fold, "-fold"))

  # assign gene groups from loading matrix
  selected_Z <- Zmatrix[selected_biomarkers, -1] |>
    data.frame() |>
    setNames(c("M", "joint", "nonM"))

  M_group <- colnames(selected_Z) |>
    lapply(function(cc) {
      selected_Z[cc] <- ifelse(selected_Z[[cc]] == 1, cc, "")
    }) |>
    data.frame() |>
    apply(1, paste0, collapse="")

  # genes output
  g_beta <- G_coeffs[selected_biomarkers] |>
    colMeans() |>
    as.matrix()

  genes_table <- data.frame(
    gene_name = selected_biomarkers,
    coeff = g_beta,
    methylation_group = M_group,
    functional_group = grouping[names(grouping) %in% selected_biomarkers]
  ) |>
    `rownames<-`(NULL)

  # clinical var output
  c_beta <- res_table[colnames(C)] |>
    colMeans() |>
    as.matrix()

  clinical_table <- data.frame(
    clinical_variable = colnames(C),
    coeff = c_beta
  ) |>
    `rownames<-`(NULL)

  # save result
  final_result <- list(
    performance = performance_df2,
    `Clinical Variables` = clinical_table,
    Genes = genes_table,
    cv_result = res_table |> transform(fold = 1:n_fold)
  )

  cat("\n")
  print(final_result[c("performance", "Clinical Variables")])

  # assign class
  class(final_result) <- "BayesMultiomics.CV"

  # return
  final_result
}

#' plot.BayesMultiomics.CV
#'
#' Summary box plot of BayesMultiomics.CV model object
#'
#' @param BayesMultiomics.CV_mod BayesMultiomics.CV model object
#' @export
plot.BayesMultiomics.CV <- function(BayesMultiomics.CV_mod){
  res_table <- BayesMultiomics.CV_mod$cv_result
  unwanted_cols <- c("n_selected_vars", "r2_train", "r2_test",
                     "cindex_train", "cindex_test", "mse_train", "mse_test", "fold")
  res_table[, setdiff(names(res_table), unwanted_cols)] |>
    utils::stack() |>
    setNames(c("beta", "clinical_x")) |>
    ggplot2::ggplot(ggplot2::aes(x=beta)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~clinical_x, scales = "free_x", ncol=1) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::ggtitle(paste0("a0 = ", BayesMultiomics.CV_mod$performance[[1]], ", ",
                            "g = ", BayesMultiomics.CV_mod$performance[[2]]))

}

# plot(multiomics_cv)

# print.multiOmics
#
# print method for multiOmics model object
#
# @param multiOmics_mod multiOmics model object
# @export
# print.multiOmics <- function(multiOmics_mod){
#   print(multiOmics_mod[c("performance", "coeffs")])
#   print(multiOmics_mod$selected_genes)
# }

# coef.multiOmics
#
# coef method for multiOmics model object to obtain clinical feature coefficients
#
# @param multiOmics_mod
# @export
# coef.multiOmics <- function(multiOmics_mod){
#   coeffs <- multiOmics_mod[["coeffs"]][, "beta"]
#   # return
#   coeffs
# }

#' BayesMultiomicsPriorSensitivity
#'
#' @description
#' Run `BayesMultiomics` for given sets of parameters `g` and `a`
#'
#' @param gstr_vec a vector of `g` values, two options: "scale" or 1
#' @param a0_vec a vector of `a` values
#' @param ... Additional parameters for `BayesMultiomics`
#'
#' @examples
#' # params
#' M <- GBM_data2$M
#' G <- GBM_data2$G
#' # fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
#' fp <- system.file("eg_fx_classification.txt", package="BayesMultiomics")
#' gene_grouping <- DAVIDGeneGrouping(eg_gene_symbols, fp)
#' Y <- GBM_data2$Y
#' C <- GBM_data2$C
#' Delta <- GBM_data2$Delta
#'
#' a0 <- c(0.1, 1, 10, 50)
#' g <- list("scale", 1)
#'
#'
#' # run
#' multiomics_sensitivity <- BayesMultiomicsPriorSensitivity(
#'   M=M,
#'   G=G,
#'   grouping=gene_grouping,
#'   Y=Y,
#'   C=C,
#'   a0_vec=a0,
#'   gstr_vec=g,
#'   Delta=Delta,
#'   n_fold=10,
#'   random_seed = 123,
#'   EMVS_I = 10,
#'   NEG_I=20,
#'   EMVS_thresh = 0.0001,
#'   NEG_thresh = 0.0001
#' )
#'
#' @export
BayesMultiomicsPriorSensitivity <- function(
    gstr_vec, a0_vec, M, G, grouping, nu0=0.5, nu1=10^3, nu=1, lambda=1, a=1, b=1, EMVS_I=10, EMVS_thresh=0.0001,
    Y, C, Delta, n_fold=10, random_seed=NULL, NEG_I=10, NEG_thresh=0.0001,
    lower_R2=0.2, upper_R2=0.8, transform_M = "L", .mpmath=setup_mpmath()
){
  # within MultiOmics, run EMVS
  EMVS_result <- .EMVS(M,
                       G,
                       grouping=grouping,
                       nu0 = nu0,
                       nu1 = nu1,
                       lambda = lambda,
                       a = a,
                       b = b,
                       I=EMVS_I,
                       thresh=EMVS_thresh,
                       transform_M = transform_M)

  # run 2nd stage
  N <- nrow(G)
  K <- ncol(G)
  J <- ncol(M)
  L <- ncol(C)
  R2 <- EMVS_result$R2

  # build Z matrix
  Zmatrix <- .Zmat_builder(R2, G, lower_R2 = lower_R2, upper_R2 = upper_R2)

  # E-M loop step
  # loop params
  # a0 <- c(0.1, 1, 10, 50)
  # g <- list("scale", 1)
  combos <- expand.grid(list(a0=a0_vec, gstr=gstr_vec))
  combos <- within(combos, {
    nms <- paste0("a0=",a0," gstr=",gstr)
  })

  result_list <- Map(
    function(aa, gg){
      # run cross validation
      if(is.null(random_seed)) set.seed(random_seed)
      # folds <- caret::createFolds(1:N, k=n_fold)
      folds <- .folds(1:N, num_fold=n_fold)

      # initialize result objects
      selected_biomarkers <- c()
      res_cols <- c("n_selected_vars", "r2_train", "r2_test",
                    "cindex_train", "cindex_test", "mse_train", "mse_test",
                    colnames(C))

      res_table <- data.frame(matrix(0, n_fold, length(res_cols))) |>
        setNames(res_cols)

      G_coeffs <- data.frame(matrix(0, n_fold, ncol(G))) |>
        setNames(colnames(G))

      message("\n2nd stage modeling with a0 = ", aa, "and g = ", gg)
      message("Running cross validation...")
      pb <- txtProgressBar(min=0, max=n_fold, initial=0, style=3)

      # loop
      for(jjj in 1:n_fold){
        test_indx <- folds[[jjj]]
        train_delta <- Delta[-test_indx]
        test_delta <- Delta[test_indx]

        NEG_output <- .NEG_em(Y=Y[-test_indx],
                              G=as.matrix(G[-test_indx,]),
                              C=as.matrix(C[-test_indx,]),
                              a0=aa,
                              gstr=gg,
                              Zmatrix=Zmatrix,
                              I=NEG_I,
                              thresh=NEG_thresh,
                              .mpmath=.mpmath)

        # save result
        NEG_output$beta <- data.frame(NEG_output$beta)
        colnames(NEG_output$beta) <- c(colnames(G), colnames(C))
        estbeta <- NEG_output$beta[NEG_output$k, ]
        selected_biomarkers <- unique(c(selected_biomarkers,
                                        names(estbeta)[abs(estbeta) > 1e-5 & !names(estbeta) %in% colnames(C)]))
        # evaluate cross validation result
        pred_train <- cbind(G[-test_indx,], C[-test_indx,]) %*% as.numeric(estbeta)
        pred_test <- cbind(G[test_indx,], C[test_indx,]) %*% as.numeric(estbeta)
        res_table[jjj, "n_selected_vars"] <- sum(abs(estbeta) > 1e-5 & !names(estbeta) %in% colnames(C))
        res_table[jjj, "r2_train"] <- .rsq(pred_train, Y[-test_indx])
        res_table[jjj, "r2_test"] <- .rsq(pred_test, Y[test_indx])
        res_table[jjj, "cindex_train"] <- .cindx(pred_train, Y[-test_indx])
        res_table[jjj, "cindex_test"] <- .cindx(pred_test, Y[test_indx])
        res_table[jjj, "mse_train"] <- mean((pred_train - Y[-test_indx])^2)
        res_table[jjj, "mse_test"] <- mean((pred_test - Y[test_indx])^2)
        res_table[jjj, colnames(C)] <- estbeta[colnames(C)]
        setTxtProgressBar(pb, jjj)
      }

      # take the average
      performance_df <- data.frame(prior_a=aa,
                                   prior_g=gg,
                                   as.list(colMeans(res_table)))
      performance_df2 <- performance_df[, ! colnames(performance_df) %in% colnames(C)] |>
        t() |>
        `colnames<-`(paste0("average of ", n_fold, "-fold"))

      # assign gene groups from loading matrix
      selected_Z <- Zmatrix[selected_biomarkers, -1] |>
        data.frame() |>
        setNames(c("M", "joint", "nonM"))

      M_group <- colnames(selected_Z) |>
        lapply(function(cc) {
          selected_Z[cc] <- ifelse(selected_Z[[cc]] == 1, cc, "")
        }) |>
        data.frame() |>
        apply(1, paste0, collapse="")

      # genes output
      g_beta <- G_coeffs[selected_biomarkers] |>
        colMeans() |>
        as.matrix()

      genes_table <- data.frame(
        gene_name = selected_biomarkers,
        coeff = g_beta,
        methylation_group = M_group,
        functional_group = grouping[names(grouping) %in% selected_biomarkers]
      ) |>
        `rownames<-`(NULL)

      # clinical var output
      c_beta <- res_table[colnames(C)] |>
        colMeans() |>
        as.matrix()

      clinical_table <- data.frame(
        clinical_variable = colnames(C),
        coeff = c_beta
      ) |>
        `rownames<-`(NULL)

      # save result
      final_result <- list(
        performance = performance_df2,
        `Clinical Variables` = clinical_table,
        Genes = genes_table,
        cv_result = res_table |> transform(fold = 1:n_fold)
      )

      # assign class
      class(final_result) <- "BayesMultiomics.CV"

      # return
      final_result

    },
    combos$a0, combos$gstr
  ) |>
    setNames(combos$nms)



  # assign class
  class(result_list) <- "BayesMultiomicsPriorSensitivity"

  # return
  result_list
}

#' summary.BayesMultiomicsPriorSensitivity
#'
#' summary method for BayesMultiomicsPriorSensitivity object
#'
#' @param BayesMultiomicsPriorSensitivity_obj BayesMultiomicsPriorSensitivity object
#' @export
summary.BayesMultiomicsPriorSensitivity <- function(BayesMultiomicsPriorSensitivity_obj){
  summary_df <- lapply(BayesMultiomicsPriorSensitivity_obj, function(x){
    pfm <- x[["performance"]] |> t() |> as.data.frame()
    cf <- x[["Clinical Variables"]] |>
      `rownames<-`(x[["Clinical Variables"]][["clinical_variable"]])
    cf <- cf[["coeff"]] |> t() |> as.data.frame()

    # return
    cbind(pfm, cf)
  }) |>
    do.call(rbind, args=_) |>
    `rownames<-`(NULL)

  # return
  summary_df
}

#' plot.BayesMultiomicsPriorSensitivity
#'
#' plot function for BayesMultiomicsPriorSensitivity object
#'
#' @param BayesMultiomicsPriorSensitivity_obj BayesMultiomicsPriorSensitivity object
#' @export
plot.BayesMultiomicsPriorSensitivity <- function(BayesMultiomicsPriorSensitivity_obj){
  lapply(BayesMultiomicsPriorSensitivity_obj, function(x) {
    unwanted_cols <- c("n_selected_vars", "r2_train", "r2_test",
                       "cindex_train", "cindex_test", "mse_train", "mse_test", "fold")

    res_table <- x$cv_result

    res_table <- res_table[, setdiff(names(res_table), unwanted_cols)] |>
      utils::stack() |>
      setNames(c("beta", "clinical_x"))
    # res_table$a0 <- x$performance[[1]]
    # res_table$g <- x$performance[[2]]
    res_table$a_g <- paste0("a=",x$performance[[1]],",",
                            "g=",x$performance[[2]])
    res_table
  }) |>
    do.call(rbind, args=_) |>
    ggplot2::ggplot(ggplot2::aes(y=beta, x=a_g)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~clinical_x, scales = "free_y", ncol=1) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x =ggplot2::element_text(angle=15))
}






#' BayesMultiomics
#'
#' Run first and second stage models in *Hao Xue, Sounak Chakraborty & Tanujit Dey (2024), Bayesian Shrinkage Models for Integration and Analysis of Multiplatform High-dimensional Genomics Data*.
#' In the first stage, Expectation–Maximization Variable Selection (EMVS) is used to learn the regulating mechanism between epigenomics (e.g., gene methylation) and gene expression while considering functional gene annotations.
#' Genes are grouped based on the regulating mechanism learned in the first stage using the function, EMVS.
#' Then, a group-wise penalty is applied to select genes significantly associated with clinical outcomes while incorporating clinical features using Bayesian variable selection model with extended NEG priors in *Veronika Rockova & Edward I. George (2013): EMVS: The EM Approach to Bayesian Variable Selection, Journal of the American Statistical Association*.
#'
#'
#' @param M DNA Methylation matrix, a matrix with rows corresponding to samples and columns corresponding to probes/sites.
#' @param G Gene expression level, a matrix with rows corresponding to samples and columns corresponding to genes.
#' @param grouping Gene grouping, a vector of group membership of each gene. This could be given by the user or obatined by performing DAVID Functional Classification of genes. See vignette("gene_grouping").
#' @param nu0 Parameter `nu0` for spike-and-slab Gaussian mixture prior on \eqn{\beta}
#' @param nu1 parameter `nu1` for spike-and-slab Gaussian mixture prior on \eqn{\beta}
#' @param lambda For the prior on \eqn{\sigma2}, an inverse gamma prior \eqn{\pi(\sigma2 | \gamma) = IG(\nu/2, \nu\lambda/2)}
#' @param a Parameter \eqn{a} for the beta prior \eqn{\pi(\theta) \propto \theta a-1(1-\theta)b-1}, \eqn{a, b > 0}
#' @param b Parameter \eqn{b} for the beta prior \eqn{\pi(\theta) \propto \theta a-1(1-\theta)b-1}, \eqn{a, b > 0}
#' @param EMVS_I Maximum number of iterations of EMVS
#' @param EMVS_thresh Threshold for convergence criterion for EMVS
#' @param Y clinical outcome
#' @param G gene expression level
#' @param C clinical features
#' @param a0 size sparsity
#' @param gstr prior, two options: "scale" or "1"
#' @param NEG_I number of maximum iterations for NEG_em
#' @param NEG_thresh convergence criterion fpr NEG_em
#' @param .mpmath (Optional) function depends on mpmath package from python. pointer for mpmath package. Default: setup_mpmath()
#' @param lower_R2 lower limit for R2 threshold, default = 0.2
#' @param upper_R2 upper limit for R2 threshold, default = 0.8
#' @param transform_M transform methylation matrix `M`. Options: "L" = linear transformation, "Q" = quadratic transformation; "QS" = cubic spline transformation
#'
#'
#' @details
#' `gstr` will take two options "scale" or "1." if `gstr` == "scale" then g = 1/N^2 where N = number of genes
#'
#' @examples
#' # params
#' M <- GBM_data2$M
#' G <- GBM_data2$G
#' fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
#' # fp <- system.file("eg_fx_classification.txt", package="BayesMultiomics")
#' gene_grouping <- DAVIDGeneGrouping(eg_gene_symbols, fp)
#' Y <- GBM_data2$Y
#' C <- GBM_data2$C
#' a0 <- 0.1
#' gstr <- "scale"
#' mpmath <- setup_mpmath()
#'
#' # run
#' multiomics <- BayesMultiomics(
#'   M=M,
#'   G=G,
#'   grouping=gene_grouping,
#'   Y=Y,
#'   C=C,
#'   a0=a0,
#'   gstr=gstr,
#'   EMVS_I = 10,
#'   NEG_I= 10,
#'   EMVS_thresh = 0.0001,
#'   NEG_thresh = 0.0001,
#'   transform_M = "L",
#'   lower_R2 = 0.2,
#'   upper_R2 = 0.8,
#'   .mpmath = mpmath
#' )
#'
#' @export
BayesMultiomics <- function(
    M, G, grouping, nu0=0.5, nu1=10^3, nu=1, lambda=1, a=1, b=1, EMVS_I=10, EMVS_thresh=0.0001,
    transform_M="L",
    Y, C,
    a0, gstr,
    NEG_I=10, NEG_thresh=0.0001,
    lower_R2=0.2, upper_R2=0.8,
    .mpmath = setup_mpmath()
){

  # 1st stage
  EMVS_result <- .EMVS(M,
                       G,
                       grouping=grouping,
                       nu0 = nu0,
                       nu1 = nu1,
                       lambda = lambda,
                       a = a,
                       b = b,
                       I=EMVS_I,
                       thresh=EMVS_thresh,
                       transform_M = transform_M)


  R2 <- EMVS_result$R2

  # build Z matrix
  Zmatrix <- Zmat_builder(R2, G, lower_R2, upper_R2)

  message("\n2nd stage modeling with a0 = ", a0, "and g = ", gstr)

  # 2nd stage
  NEG_output <- .NEG_em(Y=Y,
                        G=G,
                        C=C,
                        a0=a0,
                        gstr=gstr,
                        Zmatrix=Zmatrix,
                        I=NEG_I,
                        thresh=NEG_thresh,
                        .mpmath=.mpmath)

  # estimated coefficients
  NEG_output$beta <- data.frame(NEG_output$beta)
  colnames(NEG_output$beta) <- c(colnames(G), colnames(C))
  estbeta <- NEG_output$beta[NEG_output$k, ]

  # format output
  genes <- estbeta[, setdiff(names(estbeta), colnames(C))]
  selected_genes <- genes[, abs(genes) > 1e-5] # the threshold should be parameterized?
  selected_Z <- Zmatrix[names(selected_genes), -1] |>
    data.frame() |>
    setNames(c("M", "joint", "nonM"))

  # combine M effect group columns into one
  M_group <- colnames(selected_Z) |>
    lapply(function(cc) {
      selected_Z[cc] <- ifelse(selected_Z[[cc]] == 1, cc, "")
    }) |>
    data.frame() |>
    apply(1, paste0, collapse="")

  genes_table <- data.frame(
    gene_name = names(selected_genes),
    coeff = t(selected_genes)[,],
    methylation_group = M_group,
    functional_group = grouping[names(grouping) %in% names(selected_genes)]
  ) |>
    `rownames<-`(NULL)

  clinical_table <- data.frame(
    clinical_variable = colnames(C),
    coeff = t(estbeta[, colnames(C)])[,]
  ) |>
    `rownames<-`(NULL)

  # return
  list(
    Genes = genes_table,
    `Clinical Variables` = clinical_table
  )
}
