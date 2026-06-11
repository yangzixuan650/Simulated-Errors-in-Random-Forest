# ============================================================
# Utility functions
# ============================================================




clamp_prob <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}


binary_to_numeric01 <- function(y, positive_class = NULL) {
  if (is.factor(y)) {
    if (length(levels(y)) != 2) stop("y must be binary.")
    lv <- levels(y)
    if (is.null(positive_class)) {
      positive_class <- lv[2]
    }
    if (!positive_class %in% lv) {
      stop("positive_class not found in y levels.")
    }
    y_num <- as.integer(y == positive_class)
    return(list(y_num = y_num, positive_class = positive_class, levels = lv))
  }
  
  if (!all(y %in% c(0, 1))) {
    stop("y must be a 2-level factor or numeric 0/1.")
  }
  return(list(y_num = as.integer(y), positive_class = 1, levels = c(0, 1)))
}


brier_score <- function(y, p) {
  y_info <- binary_to_numeric01(y)
  mean((y_info$y_num - p)^2)
}


log_loss <- function(y, p, eps = 1e-6) {
  y_info <- binary_to_numeric01(y)
  p <- clamp_prob(p, eps = eps)
  -mean(y_info$y_num * log(p) + (1 - y_info$y_num) * log(1 - p))
}


auc_binary <- function(y, p) {
  y_info <- binary_to_numeric01(y)
  y_num <- y_info$y_num
  
  n1 <- sum(y_num == 1)
  n0 <- sum(y_num == 0)
  
  if (n1 == 0 || n0 == 0) {
    warning("AUC is undefined because only one class is present.")
    return(NA_real_)
  }
  
  r <- rank(p, ties.method = "average")
  auc <- (sum(r[y_num == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  auc
}


confusion_binary <- function(y, p, threshold = 0.5, positive_class = NULL) {
  y_info <- binary_to_numeric01(y, positive_class = positive_class)
  y_num <- y_info$y_num
  
  pred <- ifelse(p >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & y_num == 1)
  tn <- sum(pred == 0 & y_num == 0)
  fp <- sum(pred == 1 & y_num == 0)
  fn <- sum(pred == 0 & y_num == 1)
  
  out <- matrix(c(tn, fp, fn, tp), nrow = 2, byrow = TRUE)
  rownames(out) <- c("Observed_0", "Observed_1")
  colnames(out) <- c("Pred_0", "Pred_1")
  out
}


metric_summary_binary <- function(y, p, threshold = 0.5, positive_class = NULL) {
  y_info <- binary_to_numeric01(y, positive_class = positive_class)
  y_num <- y_info$y_num
  pred <- ifelse(p >= threshold, 1, 0)
  
  tp <- sum(pred == 1 & y_num == 1)
  tn <- sum(pred == 0 & y_num == 0)
  fp <- sum(pred == 1 & y_num == 0)
  fn <- sum(pred == 0 & y_num == 1)
  
  sensitivity <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  precision   <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  npv         <- if ((tn + fn) == 0) NA_real_ else tn / (tn + fn)
  accuracy    <- mean(pred == y_num)
  f1          <- if (is.na(precision) || is.na(sensitivity) || (precision + sensitivity) == 0) {
    NA_real_
  } else {
    2 * precision * sensitivity / (precision + sensitivity)
  }
  
  c(
    auc = auc_binary(y, p),
    brier = brier_score(y, p),
    logloss = log_loss(y, p),
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    npv = npv,
    f1 = f1
  )
}

f <- function(FFMC.star, param){
  for(x0 in 0:101){
    # if FFMC.star is within a (x0-0.5, x0+0.5) neighborhood of x0:
    if((FFMC.star - x0 > -0.5) & (FFMC.star - x0 < 0.5)){
      a <- param[x0 + 1, 1]
      b <- param[x0 + 1, 2]
      mu <- param[x0 + 1, 3]
      sd <- param[x0 + 1, 4]
    }
  }
  return(list(a = a, b = b, mu = mu, sd = sd))
}



# ============================================================
# Main training function
# ============================================================

noiseAvgRF <- function(formula,
                       data,
                       error.variables,
                       measurement.error,
                       B = 100,
                       zeta = 1,
                       ntree = 500,
                       mtry = NULL,
                       nodesize = 1,
                       maxnodes = NULL,
                       replace = TRUE,
                       sampsize = NULL,
                       balance_classes = FALSE,
                       importance = FALSE,
                       keep.forest = TRUE,
                       positive_class = NULL,
                       seed = NULL,
                       verbose = TRUE) {
  
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required.")
  }
  
  if (!inherits(formula, "formula")) {
    stop("formula must be a formula.")
  }
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }
  if (length(error.variables) == 0) {
    stop("error.variables must contain at least one variable.")
  }
  if (!all(error.variables %in% names(data))) {
    missing_vars <- error.variables[!error.variables %in% names(data)]
    stop("These error.variables are missing from data: ",
         paste(missing_vars, collapse = ", "))
  }
  if (zeta < 0) {
    stop("zeta must be >= 0.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  mf <- model.frame(formula, data = data, na.action = na.fail)
  y <- model.response(mf)
  
  if (!is.factor(y) || length(levels(y)) != 2) {
    stop("This version currently supports binary classification only. y must be a 2-level factor.")
  }
  
  if (is.null(positive_class)) {
    positive_class <- levels(y)[2]
  }
  if (!positive_class %in% levels(y)) {
    stop("positive_class must be one of the factor levels of y.")
  }
  
  y_name <- names(mf)[1]
  x_terms <- attr(stats::terms(formula), "term.labels")
  
  p_total <- ncol(model.matrix(formula, data = data)) - 1L
  if (is.null(mtry)) {
    mtry <- max(floor(sqrt(max(p_total, 1))), 1)
  }
  
  #a <- measurement.error$a
  #b <- measurement.error$b
  #mean <- measurement.error$mean
  #sd <- measurement.error$sd
  
  #extra_a <- a
  #extra_b <- b
  #extra_mean <- mean
  #extra_sd <- sqrt(zeta) * sd
  
  forests <- vector("list", B)
  
  # balanced class sampling inside each RF, if requested
  rf_sampsize <- sampsize
  if (balance_classes) {
    tab <- table(y)
    min_n <- min(tab)
    rf_sampsize <- rep(min_n, length(tab))
    names(rf_sampsize) <- names(tab)
  }
  
  if (verbose) {
    message("Training ", B, " noise-augmented random forests ...")
  }
  
  for (b in seq_len(B)) {
    dat_b <- data
    
    for (v in error.variables) {
      if (!is.numeric(dat_b[[v]])) {
        stop("Currently, error.variables must be numeric. Problem variable: ", v)
      }
      
      #noise_b <- rtruncnorm(
      #  n = nrow(dat_b),
      #  a = extra_a,
      #  b = extra_b,
      #  mean = extra_mean,
      #  sd = extra_sd
      #)
      
      res <- sapply(dat_b[[v]], FUN = f, param = measurement.error)
      
      dat_b <- dat_b %>%
        mutate(Ui = rtruncnorm(nrow(dat_b), a = res[1,], b = res[2,], 
                               mean = res[3,], sd = res[4,]))
      
      dat_b[[v]] <- dat_b[[v]] - dat_b$Ui
    }
    
    rf_args <- list(
      formula = formula,
      data = dat_b,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      maxnodes = maxnodes,
      replace = replace,
      importance = importance,
      keep.forest = keep.forest
    )
    
    if (!is.null(rf_sampsize)) {
      rf_args$sampsize <- rf_sampsize
    }
    
    forests[[b]] <- do.call(randomForest::randomForest, rf_args)
    
    if (verbose && (b %% max(1, floor(B / 10)) == 0 || b == B)) {
      message("  finished ", b, "/", B)
    }
  }
  
  out <- list(
    call = match.call(),
    formula = formula,
    response_name = y_name,
    predictor_names = x_terms,
    error.variables = error.variables,
    measurement.error = measurement.error,
    zeta = zeta,
    B = B,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    maxnodes = maxnodes,
    replace = replace,
    sampsize = rf_sampsize,
    balance_classes = balance_classes,
    importance = importance,
    keep.forest = keep.forest,
    positive_class = positive_class,
    class_levels = levels(y),
    forests = forests,
    train_data = data,
    calibration = NULL
  )
  
  class(out) <- "noiseAvgRF"
  out
}


# ============================================================
# Raw prediction helper
# ============================================================

predict_noiseAvgRF_raw_prob <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    newdata <- object$train_data
  }
  
  class_levels <- object$class_levels
  pos_class <- object$positive_class
  B <- length(object$forests)
  
  prob_array <- array(
    0,
    dim = c(nrow(newdata), length(class_levels), B),
    dimnames = list(NULL, class_levels, NULL)
  )
  
  for (b in seq_len(B)) {
    prob_b <- stats::predict(object$forests[[b]], newdata = newdata, type = "prob")
    
    prob_b2 <- matrix(0, nrow = nrow(newdata), ncol = length(class_levels))
    colnames(prob_b2) <- class_levels
    prob_b2[, colnames(prob_b)] <- prob_b
    
    prob_array[, , b] <- prob_b2
  }
  
  avg_prob <- apply(prob_array, c(1, 2), mean)
  if (is.vector(avg_prob)) {
    avg_prob <- matrix(avg_prob, ncol = length(class_levels))
    colnames(avg_prob) <- class_levels
  }
  
  list(
    prob = avg_prob,
    prob_positive = avg_prob[, pos_class],
    individual_prob = prob_array
  )
}


# ============================================================
# Calibration
# ============================================================

calibrate_noiseAvgRF <- function(object,
                                 newdata,
                                 y,
                                 method = c("platt", "none"),
                                 eps = 1e-6) {
  method <- match.arg(method)
  
  if (!inherits(object, "noiseAvgRF")) {
    stop("object must be a noiseAvgRF object.")
  }
  
  if (method == "none") {
    object$calibration <- list(
      method = "none",
      model = NULL,
      eps = eps
    )
    return(object)
  }
  
  y_info <- binary_to_numeric01(y, positive_class = object$positive_class)
  y_num <- y_info$y_num
  
  raw <- predict_noiseAvgRF_raw_prob(object, newdata = newdata)
  p_raw <- clamp_prob(raw$prob_positive, eps = eps)
  
  df_cal <- data.frame(
    y = y_num,
    score = qlogis(p_raw)
  )
  
  calib_fit <- stats::glm(y ~ score, data = df_cal, family = stats::binomial())
  
  object$calibration <- list(
    method = "platt",
    model = calib_fit,
    eps = eps
  )
  
  object
}


# ============================================================
# Main prediction method
# ============================================================

predict.noiseAvgRF <- function(object,
                               newdata = NULL,
                               type = c("prob", "raw_prob", "class", "link", "all"),
                               threshold = 0.5,
                               calibrated = TRUE,
                               ...) {
  type <- match.arg(type)
  
  raw <- predict_noiseAvgRF_raw_prob(object, newdata = newdata)
  p_raw <- raw$prob_positive
  class_levels <- object$class_levels
  pos_class <- object$positive_class
  neg_class <- setdiff(class_levels, pos_class)
  
  p_use <- p_raw
  link_use <- qlogis(clamp_prob(p_raw))
  
  if (calibrated &&
      !is.null(object$calibration) &&
      identical(object$calibration$method, "platt")) {
    
    eps <- object$calibration$eps
    df_new <- data.frame(score = qlogis(clamp_prob(p_raw, eps = eps)))
    link_use <- stats::predict(object$calibration$model, newdata = df_new, type = "link")
    p_use <- stats::plogis(link_use)
  }
  
  prob_mat <- cbind(1 - p_use, p_use)
  colnames(prob_mat) <- c(as.character(neg_class), as.character(pos_class))
  
  pred_class <- ifelse(p_use >= threshold, as.character(pos_class), as.character(neg_class))
  pred_class <- factor(pred_class, levels = class_levels)
  
  if (type == "raw_prob") {
    out <- cbind(1 - p_raw, p_raw)
    colnames(out) <- c(as.character(neg_class), as.character(pos_class))
    return(out)
  }
  
  if (type == "prob") {
    return(prob_mat)
  }
  
  if (type == "link") {
    return(link_use)
  }
  
  if (type == "class") {
    return(pred_class)
  }
  
  list(
    prob = prob_mat,
    raw_prob = cbind(1 - p_raw, p_raw),
    class = pred_class,
    individual_prob = raw$individual_prob
  )
}


# ============================================================
# Evaluation
# ============================================================

evaluate_noiseAvgRF <- function(object,
                                newdata,
                                y,
                                threshold = 0.5,
                                calibrated = TRUE) {
  p <- predict(object, newdata = newdata, type = "prob",
               threshold = threshold, calibrated = calibrated)[, as.character(object$positive_class)]
  
  metrics <- metric_summary_binary(
    y = y,
    p = p,
    threshold = threshold,
    positive_class = object$positive_class
  )
  
  cm <- confusion_binary(
    y = y,
    p = p,
    threshold = threshold,
    positive_class = object$positive_class
  )
  
  list(
    metrics = metrics,
    confusion_matrix = cm
  )
}


# ============================================================
# Print / summary methods
# ============================================================

print.noiseAvgRF <- function(x, ...) {
  cat("noiseAvgRF object\n")
  cat("  Formula            :", deparse(x$formula), "\n")
  cat("  Error variables    :", paste(x$error.variables, collapse = ", "), "\n")
  cat("  B                  :", x$B, "\n")
  cat("  zeta               :", x$zeta, "\n")
  cat("  ntree per RF       :", x$ntree, "\n")
  cat("  mtry               :", x$mtry, "\n")
  cat("  Balanced sampling  :", x$balance_classes, "\n")
  cat("  Positive class     :", as.character(x$positive_class), "\n")
  cat("  Calibration        :",
      if (is.null(x$calibration)) "not fitted" else x$calibration$method, "\n")
  invisible(x)
}


summary.noiseAvgRF <- function(object, ...) {
  print(object)
  
  if (object$importance) {
    imp_list <- lapply(object$forests, function(f) randomForest::importance(f))
    cat("\nVariable importance was requested during fitting.\n")
    cat("You can inspect individual forests via object$forests.\n")
  }
  
  invisible(object)
}




