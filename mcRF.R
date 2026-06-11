

mcRF <- function(formula,
                 data,
                 B = 100,
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
  
  forests <- vector("list", B)
  
  rf_sampsize <- sampsize
  if (balance_classes) {
    tab <- table(y)
    min_n <- min(tab)
    rf_sampsize <- rep(min_n, length(tab))
    names(rf_sampsize) <- names(tab)
  }
  
  if (!is.null(rf_sampsize)) {
    if (any(!is.finite(rf_sampsize)) || any(rf_sampsize <= 0)) {
      stop("sampsize must be NULL or positive finite integer(s).")
    }
  }
  
  if (verbose) {
    message("Training ", B, " Monte Carlo baseline random forests ...")
  }
  
  for (b in seq_len(B)) {
    rf_args <- list(
      formula = formula,
      data = data,
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
  
  class(out) <- "mcRF"
  out
}


calibrate.mcRF <- function(object,
                           newdata,
                           y,
                           eps = 1e-6) {
  if (!inherits(object, "mcRF")) {
    stop("object must be an mcRF object.")
  }
  
  if (is.factor(y)) {
    if (length(levels(y)) != 2) stop("y must be binary.")
    y_num <- as.integer(y == levels(y)[2])
  } else {
    if (!all(y %in% c(0, 1))) stop("y must be coded as 0/1 or a 2-level factor.")
    y_num <- y
  }
  
  p_raw <- predict(object, newdata = newdata, type = "prob")[, 2]
  p_raw <- pmin(pmax(p_raw, eps), 1 - eps)
  
  df_cal <- data.frame(
    y = y_num,
    score = qlogis(p_raw)
  )
  
  calib_fit <- glm(y ~ score, data = df_cal, family = binomial())
  
  object$calibration <- list(
    method = "platt",
    model = calib_fit,
    eps = eps
  )
  
  object
}


predict.mcRF <- function(object,
                         newdata = NULL,
                         type = c("prob", "class", "all"),
                         threshold = 0.5,
                         ...) {
  type <- match.arg(type)
  
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
  
  pred_class <- factor(
    class_levels[max.col(avg_prob, ties.method = "first")],
    levels = class_levels
  )
  
  if (type == "prob") return(avg_prob)
  if (type == "class") return(pred_class)
  
  list(
    prob = avg_prob,
    class = pred_class,
    individual_prob = prob_array
  )
}






predict.calibrated.mcRF <- function(object,
                                    newdata,
                                    type = c("prob", "class", "link"),
                                    threshold = 0.5) {
  type <- match.arg(type)
  
  if (is.null(object$calibration) || object$calibration$method != "platt") {
    stop("mcRF object has no Platt calibration fitted.")
  }
  
  p_raw <- predict(object, newdata = newdata, type = "prob")[, 2]
  p_raw <- pmin(pmax(p_raw, object$calibration$eps), 1 - object$calibration$eps)
  
  df_new <- data.frame(score = qlogis(p_raw))
  lp <- predict(object$calibration$model, newdata = df_new, type = "link")
  p_cal <- plogis(lp)
  
  if (type == "prob") {
    return(cbind(`0` = 1 - p_cal, `1` = p_cal))
  }
  if (type == "link") {
    return(lp)
  }
  
  factor(ifelse(p_cal >= threshold, 1, 0), levels = c(0, 1))
}





