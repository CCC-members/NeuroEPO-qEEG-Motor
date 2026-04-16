# =============================================================================
# ExportLocfdrCcaWeights.R
# =============================================================================
# Part 1 of 3: locfdr + CCA weights → NIfTI brain maps → BrainNet renders
#
# Purpose:
#   Apply local false discovery rate (locfdr) to the raw CCA canonical
#   weight vector (X-block; qEEG features). Identifies significant weights
#   at a given FDR threshold and exports all quantities (raw, normalized,
#   FDR p-values, significance mask, and scaled versions) to a MATLAB .mat
#   file for downstream frequency-band reshaping and NIfTI generation.
#
# Input:
#   cca.out -- whitening::scca() output object containing:
#     $WX           -- canonical weight matrix (typically 1 x n_features)
#     $lambda       -- canonical correlation(s)
#
# Output:
#   cca_wx_locfdr_componentN.mat -- MATLAB file containing:
#     wx_raw, wx_normal, wx_locfdr_fdr, wx_sig_mask, wx_sig_weights,
#     wx_sig_weights_lambda, plus metadata (component, fdr_threshold,
#     n_voxels, n_freqs, lambda_component)
#
# Workflow:
#   1. ExportLocfdrCcaWeights.R (this file) → cca_wx_locfdr_component0_1.mat
#   2. CreateNiisLocfdrFigureMaxMin.m → NIfTI files
#   3. PlotBrainNetLocfdrDisplayBatch4ViewsBipolar.m → PNG renders
#
# =============================================================================

export_locfdr_cca_weights <- function(cca.out,
                                      component = 1,
                                      fdr_threshold = 0.2,
                                      dims = c(3244, 49),
                                      out_file = sprintf("cca_wx_locfdr_component%d.mat", component)) {
  if (!requireNamespace("locfdr", quietly = TRUE)) {
    install.packages("locfdr", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    install.packages("R.matlab", repos = "https://cloud.r-project.org")
  }

  if (is.null(cca.out$WX)) {
    stop("cca.out$WX is missing.")
  }
  if (component < 1 || component > nrow(cca.out$WX)) {
    stop(sprintf("component must be between 1 and %d.", nrow(cca.out$WX)))
  }

  wx_raw <- as.numeric(cca.out$WX[component, ])
  expected_length <- prod(dims)
  if (length(wx_raw) != expected_length) {
    stop(sprintf(
      "Expected %d weights for dims %s, but found %d.",
      expected_length,
      paste(dims, collapse = " x "),
      length(wx_raw)
    ))
  }

  wx_normal <- as.numeric(scale(wx_raw))
  fit <- locfdr::locfdr(wx_normal, plot = 1, nulltype = 1)

  if (is.null(fit$fdr) || length(fit$fdr) != length(wx_raw)) {
    stop("locfdr did not return an fdr vector with the expected length.")
  }

  sig_mask <- fit$fdr < fdr_threshold
  wx_sig_weights <- wx_raw
  wx_sig_weights[!sig_mask] <- 0

  lambda_component <- NA_real_
  wx_sig_weights_lambda <- wx_sig_weights
  if (!is.null(cca.out$lambda) && length(cca.out$lambda) >= component) {
    lambda_component <- as.numeric(cca.out$lambda[component])
    wx_sig_weights_lambda <- wx_sig_weights * lambda_component
  }

  R.matlab::writeMat(
    con = out_file,
    wx_raw = wx_raw,
    wx_normal = wx_normal,
    wx_locfdr_fdr = as.numeric(fit$fdr),
    wx_sig_mask = as.numeric(sig_mask),
    wx_sig_weights = wx_sig_weights,
    wx_sig_weights_lambda = wx_sig_weights_lambda,
    component = as.numeric(component),
    fdr_threshold = as.numeric(fdr_threshold),
    n_voxels = as.numeric(dims[1]),
    n_freqs = as.numeric(dims[2]),
    lambda_component = lambda_component
  )

  message(sprintf(
    "Saved %s with %d significant voxel-frequency entries at FDR < %.3f.",
    normalizePath(out_file, winslash = "/", mustWork = FALSE),
    sum(sig_mask),
    fdr_threshold
  ))

  invisible(list(
    out_file = normalizePath(out_file, winslash = "/", mustWork = FALSE),
    n_significant = sum(sig_mask),
    fdr_threshold = fdr_threshold
  ))
}

# ---------------------------------------------------------------------------
# Example: Run after whitening::scca() has produced cca.out
# ---------------------------------------------------------------------------
# source("ExportLocfdrCcaWeights.R")
export_locfdr_cca_weights(
  cca.out,
  component = 1,           # First canonical component
  fdr_threshold = 0.1,     # FDR < 0.1 for statistical significance
  out_file = "cca_wx_locfdr_component0_1.mat"
)