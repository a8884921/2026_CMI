# Scoring functions. Each reproduces the original prediction path from the model parameters alone
# (tests/test_equivalence.R).

load_assets <- function(ref_dir) {
  rd <- function(f) readRDS(file.path(ref_dir, f))
  list(ig = rd("ig_gbm.rds"), mtp = rd("mtp_gbm.rds"), cmi = rd("cmi_gbm.rds"),
       ppi = rd("ccat_ppi.rds")$adjacency, cnv = rd("cnv_reference.rds"),
       ann = rd("annotation_snapshot.rds"))
}

# ---- caret preProcess (center, scale, then median imputation, in caret's order) ------------------
apply_preprocess <- function(x, pre) {
  x <- as.data.frame(x, check.names = FALSE)
  if (length(pre$center)) x[, pre$center] <- sweep(x[, pre$center, drop = FALSE], 2, pre$mean[pre$center], "-")
  if (length(pre$scale))  x[, pre$scale]  <- sweep(x[, pre$scale,  drop = FALSE], 2, pre$std[pre$scale], "/")
  if (length(pre$median_impute) && any(!stats::complete.cases(x))) {
    for (v in pre$median_impute) { miss <- is.na(x[[v]]); if (any(miss)) x[miss, v] <- pre$median[[v]] }
  }
  x
}

# gbm_pred (gbm 2.x C routine) in plain R: F = initF, then F += value of the reached leaf, tree by
# tree; continuous split x < code -> left, NA -> missing node. Bit-identical to gbm::predict.gbm.
gbm_link <- function(model, X) {
  X <- as.matrix(X); storage.mode(X) <- "double"; n <- nrow(X); f <- rep(model$init_f, n)
  for (tr in model$trees) {
    node <- rep.int(1L, n)
    repeat {
      sv <- tr$split_var[node]; act <- which(sv > 0L)
      if (!length(act)) break
      nd <- node[act]; v <- X[cbind(act, sv[act])]
      node[act] <- ifelse(is.na(v), tr$missing[nd], ifelse(v < tr$split_code[nd], tr$left[nd], tr$right[nd]))
    }
    f <- f + tr$split_code[node]
  }
  f
}
gbm_tumor_prob <- function(asset, x) {
  x <- apply_preprocess(x[, asset$features, drop = FALSE], asset$preprocess)
  p <- 1 / (1 + exp(-gbm_link(asset$model, x)))
  probs <- cbind(p, 1 - p); colnames(probs) <- asset$obs_levels   # caret bernoulli prob()
  as.numeric(probs[, asset$positive])
}

# ---- CMI: gradient boosting model evaluated from its exported trees ------------------------------
cmi_eta <- function(features, model) {
  X <- as.matrix(features[, model$feature_names, drop = FALSE]); storage.mode(X) <- "double"
  n <- nrow(X); eta <- rep(0, n)
  for (tree in model$trees) {
    node <- rep(1L, n)
    repeat {
      active <- node > 0L
      if (!any(active)) break
      idx <- which(active)
      v <- X[cbind(idx, tree$f[node[idx]])]
      missing <- is.na(v)
      go_left <- !missing & v < tree$thr[node[idx]]
      if (any(missing)) go_left[missing] <- tree$na_left[node[idx[missing]]]
      nxt <- integer(length(idx))
      nxt[go_left] <- tree$left[node[idx[go_left]]]
      nxt[!go_left] <- tree$right[node[idx[!go_left]]]
      node[idx] <- nxt
      leaf_idx <- idx[tree$leaf[node[idx]]]
      if (length(leaf_idx)) { eta[leaf_idx] <- eta[leaf_idx] + tree$pred[node[leaf_idx]]; node[leaf_idx] <- 0L }
    }
  }
  eta + model$offset
}
cmi_predict <- function(features, model) {
  if (any(!is.finite(as.matrix(features[, model$feature_names, drop = FALSE])))) stop("All five indices must be finite before CMI scoring.")
  1 / (1 + exp(-cmi_eta(features, model)))
}

# ---- CMI-riskscore: portable random survival forest -----------------------------------------------
# Split rule of randomForestSRC (getDaughterPolaritySimpleNonFactor): x <= split value -> left.
# Prediction: mean over the 800 trees of the reached leaf's mortality.
rsf_predict <- function(asset, X) {
  X <- as.matrix(X[, asset$xvar_names, drop = FALSE]); storage.mode(X) <- "double"
  if (any(!is.finite(X))) stop("REO features must be finite after imputation.")
  n <- nrow(X); acc <- numeric(n)
  for (b in seq_len(asset$ntree)) {
    node <- rep.int(asset$root[b], n)
    repeat {
      internal <- which(asset$split_var[node] != 0L)
      if (!length(internal)) break
      nd <- node[internal]
      go_left <- X[cbind(internal, asset$split_var[nd])] <= asset$split_value[nd]
      node[internal] <- ifelse(go_left, asset$left[nd], asset$right[nd])
    }
    acc <- acc + asset$leaf_value[node]
    if (b %% 100L == 0L) invisible(gc(verbose = FALSE))
  }
  acc / asset$ntree
}
