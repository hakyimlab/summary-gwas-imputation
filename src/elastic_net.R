suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))
suppressMessages(library(dplyr))

generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

calc_R2 <- function(y, y_pred) {
  tss <- sum(y**2)
  rss <- sum((y - y_pred)**2)
  1 - rss/tss
}

calc_corr <- function(y, y_pred) {
  sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))
}


nested_cv_elastic_net_perf <- function(x, y, n_samples, n_train_test_folds, n_k_folds, alpha, observation_weights, penalty_factor) {
  # Gets performance estimates for k-fold cross-validated elastic-net models.
  # Splits data into n_train_test_folds disjoint folds, roughly equal in size,
  # and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is
  # cross validated. Then get performance measures for how the model predicts on the hold-out
  # fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis
  # is there is no correlation between prediction and observed.
  #
  # The mean and standard deviation of R^2 over all folds is then reported, and the p-values
  # are combined using Fisher's method.
    message(n_train_test_folds)
  R2_folds <- rep(0, n_train_test_folds)
  corr_folds <- rep(0, n_train_test_folds)
  zscore_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  converged_folds <- rep(0, n_train_test_folds)
  # Outer-loop split into training and test set
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    x_train <- x[train_idxs, ]
    y_train <- y[train_idxs]
    x_test <- x[test_idxs, ]
    y_test <- y[test_idxs]
    observation_weights_ <- observation_weights[train_idxs]
    # Inner-loop - split up training set for cross-validation to choose lambda.
    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
    res <- tryCatch({
      # Fit model with training data.
      fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids, weights = observation_weights_, penalty.factor = penalty_factor)
      best_lam_ind <- which.min(fit$cvm)
      status <-  if (fit$nzero[best_lam_ind] > 0) { "valid" } else { "null_model" }
      # Predict test data using model that had minimal mean-squared error in cross validation.
      list(y_pred=predict(fit, x_test, s = 'lambda.min'), status=status)
    },
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) {
          message(cond);
          list(y_pred=rep(mean(y_train), length(y_test)), status="error")
    })
    y_pred = res$y_pred
    converged_folds[test_fold] <- res$status == "valid"
    R2_folds[test_fold] <- calc_R2(y_test, y_pred)
    # Get p-value for correlation test between predicted y and actual y.
    # If there was no model, y_pred will have var=0, so cor.test will yield NA.
    # In that case, give a random number from uniform distribution, which is what would
    # usually happen under the null.
    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
  }
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  rho_avg <- mean(corr_folds)
  rho_se <- sd(corr_folds)
  rho_avg_squared <- rho_avg**2
  # Stouffer's method for combining z scores.
  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  list(R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval, converged=sum(converged_folds))
}

matrixify_ <- function(x) {
  x <- data.frame(x)
  column_labels <- colnames(x)
  row_labels <- rownames(x)
  # Convert cis_gt to a matrix for glmnet
  x <- matrix(as.matrix(x), ncol=ncol(x)) # R is such a bad language.
  colnames(x) <- column_labels
  rownames(x) <- row_labels
  x
}

set_seed <- function(seed = NA) {
    seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
    set.seed(seed)
    seed
}

train_elastic_net <- function(y, x, n_train_test_folds=5, n_k_folds=10, alpha=0.5, observation_weights=NULL, penalty_factor=NULL,  matrixify=FALSE) {
    if (matrixify) {
     x <- matrixify_(x)
     y <- as.double(unlist(data.frame(y)[1]))
    }
    if (is.null(observation_weights)) {
        observation_weights = rep(1, nrow(x))
    }
    if (is.null(penalty_factor)) {
        penalty_factor = rep(1, ncol(x))
    }

    perf_measures <- nested_cv_elastic_net_perf(x, y, length(y), n_train_test_folds, n_k_folds, alpha, observation_weights, penalty_factor)

    R2_avg <- perf_measures$R2_avg
    R2_sd <- perf_measures$R2_sd
    pval_est <- perf_measures$pval_est
    rho_avg <- perf_measures$rho_avg
    rho_se <- perf_measures$rho_se
    rho_zscore <- perf_measures$rho_zscore
    rho_avg_squared <- perf_measures$rho_avg_squared
    zscore_pval <- perf_measures$zscore_pval
    cv_converged <- perf_measures$converged

    cv_fold_ids <- generate_fold_ids(length(y), n_k_folds)
    fit <- tryCatch(
      cv.glmnet(x, y, nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids, keep = TRUE, weights = observation_weights, penalty.factor = penalty_factor),
      error = function(cond) {message('Error'); message(geterrmessage()); list()}
    )

    if (length(fit) > 0) {
        cv_R2_folds <- rep(0, n_k_folds)
        cv_corr_folds <- rep(0, n_k_folds)
        cv_zscore_folds <- rep(0, n_k_folds)
        cv_pval_folds <- rep(0, n_k_folds)
        best_lam_ind <- which.min(fit$cvm)
        for (j in 1:n_k_folds) {
          fold_idxs <- which(cv_fold_ids == j)
          adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
          cv_R2_folds[j] <- calc_R2(y[fold_idxs], adj_expr_fold_pred)
          cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, y[fold_idxs]), 0)
          cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(y[fold_idxs]) - 3) # Fisher transformation
          cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, y[fold_idxs])$p.value, runif(1))
        }
        cv_R2_avg <- mean(cv_R2_folds)
        cv_R2_sd <- sd(cv_R2_folds)
        y_pred <- predict(fit, as.matrix(x), s = 'lambda.min')
        training_R2 <- calc_R2(y, y_pred)

        cv_rho_avg <- mean(cv_corr_folds)
        cv_rho_se <- sd(cv_corr_folds)
        cv_rho_avg_squared <- cv_rho_avg**2
        # Stouffer's method for combining z scores.
        cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_k_folds)
        cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
        cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_k_folds, lower.tail = F)

        if (fit$nzero[best_lam_ind] > 0) {

          weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
          weights_n <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
          #weighted_f <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
          w <-  data.frame(weight=weights %>% as.double, feature=weights_n, stringsAsFactors=FALSE)
          model_summary <- data.frame(alpha=alpha, n_features=ncol(x), n_features_in_model=fit$nzero[best_lam_ind], lambda_min_mse=fit$lambda[best_lam_ind],
                            test_R2_avg=R2_avg, test_R2_sd=R2_sd, cv_R2_avg=cv_R2_avg, cv_R2_sd=cv_R2_sd, in_sample_R2=training_R2,
                            nested_cv_fisher_pval=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=rho_zscore, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval,
                            cv_rho_avg=cv_rho_avg, cv_rho_se=cv_rho_se, cv_rho_avg_squared=cv_rho_avg_squared, cv_zscore_est=cv_zscore_est, cv_zscore_pval=cv_zscore_pval, cv_pval_est=cv_pval_est,
                            cv_converged=cv_converged)
        } else {
          w <-  data.frame(weight = numeric(0), feature=character(0), stringsAsFactors=FALSE)
          model_summary <- data.frame(alpha=alpha, n_features=ncol(x), n_features_in_model=0, lambda_min_mse=fit$lambda[best_lam_ind],
                            test_R2_avg=R2_avg, test_R2_sd=R2_sd, cv_R2_avg=cv_R2_avg, cv_R2_sd=cv_R2_sd, in_sample_R2=training_R2,
                            nested_cv_fisher_pval=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=rho_zscore, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval,
                            cv_rho_avg=cv_rho_avg, cv_rho_se=cv_rho_se, cv_rho_avg_squared=cv_rho_avg_squared, cv_zscore_est=cv_zscore_est, cv_zscore_pval=cv_zscore_pval, cv_pval_est=cv_pval_est,
                            cv_converged=cv_converged)
        }
      } else {
        w <-  data.frame(weight = numeric(0), feature=character(0), stringsAsFactors=FALSE)
        model_summary <- list(alpha=alpha, n_features=col(x), n_features_in_model=0,lambda_min_mse=NA,
                        test_R2=R2_avg, test_R2_sd=R2_sd, cv_R2_avg=NA, cv_R2_sd=NA, in_sample_R2=NA,
                        nested_cv_fisher_pval=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=rho_zscore, rho_avg_squared=rho_avg_squared,
                        zscore_pval=zscore_pval,
                        cv_rho_avg=NA, cv_rho_se=NA, cv_rho_avg_squared=NA, cv_score_est=NA, cv_score_pval=NA, cv_pval_est=NA, cv_converged=NA)
      }
    return(list(weights=w, summary=model_summary))
}