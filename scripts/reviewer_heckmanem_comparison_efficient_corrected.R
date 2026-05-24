#!/usr/bin/env Rscript
## =============================================================================
## reviewer_heckmanem_comparison_efficient_corrected.R
##
## Corrected reviewer benchmark comparison for the ESE/ESN/EST Heckman paper.
##
## Models compared
##   NH  : classical normal Heckman model, fitted with HeckmanEM
##   HT  : symmetric Student-t Heckman model, fitted with HeckmanEM
##   CN  : contaminated-normal Heckman model, fitted with HeckmanEM
##   ESN : proposed hidden-truncation skew-normal Heckman model
##   EST : proposed hidden-truncation skew-t Heckman model
##
## Main corrections in this version
##   1. AIC/BIC are always recomputed on a common likelihood scale using the full
##      parameter count of each model.
##   2. Raw HeckmanEM AIC/BIC are used only diagnostically or to recover a
##      missing log-likelihood; they are never used as final criteria.
##   3. Numeric columns are coerced robustly before aggregation to avoid character
##      AIC/BIC arithmetic errors.
##   4. Finite objective values are separated from regular convergence. Numerical
##      penalties are finite by construction, but never count as valid fits.
##   5. Bivariate probability evaluations try deterministic TVPACK first and then
##      fall back to Genz--Bretz.
##   6. Student-t degrees of freedom are continuous. No rounding of nu is used.
##   7. EST starts use non-oracle generic degrees-of-freedom starts.
##   8. Mean-IC winners are reported only among models with regular convergence
##      rate at least CONV_THRESHOLD_FOR_IC.
##   9. Model-selection frequencies are computed replication by replication;
##      non-converged fits are counted as not selected.
##  10. Tables, diagnostics, metadata, and ggplot figures are written to the
##      report folder.
##
## Typical Deucalion run through SLURM
##   export PROFILE=final
##   export REPS_PER_SCEN=50
##   export N_PER_DATASET=1000
##   export MASTER_SEED=20260601
##   export ROOT_DIR=/projects/F202500010HPCVLABUMINHO/cecilia/ese_heckman/results/reviewer_heckmanem_comparison_B_final_seed20260601
##   export RESUME_EXISTING=FALSE
##   Rscript scripts/reviewer_heckmanem_comparison_efficient_corrected.R
##
## Report-only use
##   REPORT_ONLY=TRUE ROOT_DIR=<existing-root> Rscript scripts/reviewer_heckmanem_comparison_efficient_corrected.R
## =============================================================================

options(stringsAsFactors = FALSE)

## -----------------------------------------------------------------------------
## 0. User configuration
## -----------------------------------------------------------------------------

SCRIPT_VERSION <- "2026-05-18-efficient-v4-B-final-corrected"

PROFILE <- Sys.getenv("PROFILE", unset = "quick")
MASTER_SEED <- as.integer(Sys.getenv("MASTER_SEED", unset = "20260601"))

ROOT_DIR <- Sys.getenv(
  "ROOT_DIR",
  unset = paste0("reviewer_heckmanem_comparison_", PROFILE, "_seed", MASTER_SEED)
)
RESULTS_DIR <- file.path(ROOT_DIR, "results")
REPORT_DIR  <- file.path(ROOT_DIR, "report")

REPORT_ONLY <- as.logical(Sys.getenv("REPORT_ONLY", unset = "FALSE"))
RESUME_EXISTING <- as.logical(Sys.getenv("RESUME_EXISTING", unset = "TRUE"))
AUTO_INSTALL <- as.logical(Sys.getenv("AUTO_INSTALL", unset = "FALSE"))

profile_defaults <- function(profile) {
  switch(
    profile,
    quick = list(
      reps = 2L, n = 300L, workers = 4L,
      cdf_maxpts = 3000L, cdf_abseps = 1e-5, cdf_releps = 1e-5,
      stage1_maxit = 60L, stage2_maxit = 150L,
      stage2_keep = 1L, timeout_esn = 600L, timeout_est = 900L
    ),
    main = list(
      reps = 10L, n = 1000L, workers = 16L,
      cdf_maxpts = 10000L, cdf_abseps = 1e-5, cdf_releps = 1e-5,
      stage1_maxit = 120L, stage2_maxit = 400L,
      stage2_keep = 2L, timeout_esn = 1800L, timeout_est = 2700L
    ),
    final = list(
      reps = 50L, n = 1000L, workers = 32L,
      cdf_maxpts = 20000L, cdf_abseps = 5e-6, cdf_releps = 5e-6,
      stage1_maxit = 150L, stage2_maxit = 600L,
      stage2_keep = 2L, timeout_esn = 3600L, timeout_est = 5400L
    ),
    stop("Unknown PROFILE. Use quick, main, or final.")
  )
}

def <- profile_defaults(PROFILE)

REPS_PER_SCEN <- as.integer(Sys.getenv("REPS_PER_SCEN", unset = as.character(def$reps)))
N_PER_DATASET <- as.integer(Sys.getenv("N_PER_DATASET", unset = as.character(def$n)))
MAX_WORKERS <- as.integer(Sys.getenv(
  "MAX_WORKERS",
  unset = Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(def$workers))
))

FIT_MODELS <- strsplit(Sys.getenv("FIT_MODELS", unset = "NH,HT,CN,ESN,EST"), ",")[[1]]
FIT_MODELS <- trimws(toupper(FIT_MODELS))

SCENARIO_FILTER <- strsplit(Sys.getenv("SCENARIO_FILTER", unset = ""), ",")[[1]]
SCENARIO_FILTER <- trimws(SCENARIO_FILTER)
SCENARIO_FILTER <- SCENARIO_FILTER[nzchar(SCENARIO_FILTER)]
FIT_MODELS <- FIT_MODELS[FIT_MODELS %in% c("NH", "HT", "CN", "ESN", "EST")]
if (!length(FIT_MODELS)) stop("FIT_MODELS is empty after validation.")

## A model is eligible for best-by-mean IC only when regular convergence is high.
CONV_THRESHOLD_FOR_IC <- as.numeric(Sys.getenv("CONV_THRESHOLD_FOR_IC", unset = "0.80"))

## DGP settings for reviewer benchmark.
BETA_SEL <- c(-0.70, 0.15, 1.00)
BETA_OUT <- c(5.00, 1.50, 0.80)
SIGMA_S <- 1.00
SIGMA_Y <- 3.00
RHO <- 0.40
ETA_S <- -0.70
ETA_Y <-  0.70
TAU_BAR <- 0.50
NU_HEAVY <- 5.00
CONTAM_EPS <- 0.08
CONTAM_KAPPA <- 16.00

## HeckmanEM settings.
HECKMANEM_ITER_MAX <- as.integer(Sys.getenv("HECKMANEM_ITER_MAX", unset = "500"))
HECKMANEM_ERROR <- as.numeric(Sys.getenv("HECKMANEM_ERROR", unset = "1e-5"))
HECKMANEM_IM <- as.logical(Sys.getenv("HECKMANEM_IM", unset = "FALSE"))
HECKMANEM_VERBOSE <- as.logical(Sys.getenv("HECKMANEM_VERBOSE", unset = "FALSE"))
HECKMANEM_NU_T_INIT <- as.numeric(Sys.getenv("HECKMANEM_NU_T_INIT", unset = "6"))
HECKMANEM_NU_CN_INIT <- as.numeric(strsplit(Sys.getenv("HECKMANEM_NU_CN_INIT", unset = "0.05,0.20"), ",")[[1]])

## Numerical CDF and optimizer settings.
CDF_MAXPTS <- as.integer(Sys.getenv("CDF_MAXPTS", unset = as.character(def$cdf_maxpts)))
CDF_ABSEPS <- as.numeric(Sys.getenv("CDF_ABSEPS", unset = as.character(def$cdf_abseps)))
CDF_RELEPS <- as.numeric(Sys.getenv("CDF_RELEPS", unset = as.character(def$cdf_releps)))
STAGE1_MAXIT <- as.integer(Sys.getenv("STAGE1_MAXIT", unset = as.character(def$stage1_maxit)))
STAGE2_MAXIT <- as.integer(Sys.getenv("STAGE2_MAXIT", unset = as.character(def$stage2_maxit)))
STAGE2_KEEP <- as.integer(Sys.getenv("STAGE2_KEEP", unset = as.character(def$stage2_keep)))
TIMEOUT_ESN <- as.integer(Sys.getenv("TIMEOUT_ESN", unset = as.character(def$timeout_esn)))
TIMEOUT_EST <- as.integer(Sys.getenv("TIMEOUT_EST", unset = as.character(def$timeout_est)))
COMPUTE_HESSIAN <- as.logical(Sys.getenv("COMPUTE_HESSIAN", unset = "FALSE"))

## Bounds in eta_j = L_ETA tanh(raw_j), tau_bar = L_TAU tanh(raw_tau).
L_ETA <- as.numeric(Sys.getenv("L_ETA", unset = "1.5"))
L_TAU <- as.numeric(Sys.getenv("L_TAU", unset = "1.5"))

## Numerical penalties and validity threshold.
NLL_PENALTY <- as.numeric(Sys.getenv("NLL_PENALTY", unset = "1e12"))
NLL_VALID_MAX <- as.numeric(Sys.getenv("NLL_VALID_MAX", unset = as.character(NLL_PENALTY / 10)))

## Generic non-oracle starts for EST degrees of freedom.
EST_NU_STARTS <- as.numeric(strsplit(Sys.getenv("EST_NU_STARTS", unset = "4,10,20"), ",")[[1]])
EST_NU_STARTS <- EST_NU_STARTS[is.finite(EST_NU_STARTS) & EST_NU_STARTS > 2]
if (!length(EST_NU_STARTS)) EST_NU_STARTS <- c(4, 10, 20)

## For mean time summaries. proc.time()[3] is seconds, so do not auto-divide.
TIME_DIVISOR <- as.numeric(Sys.getenv("TIME_DIVISOR", unset = "1"))

## Threading hygiene: one BLAS/LAPACK thread per worker.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

## -----------------------------------------------------------------------------
## 1. Packages
## -----------------------------------------------------------------------------

CORE_PKGS <- c(
  "HeckmanEM", "mvtnorm", "parallel", "dplyr", "tidyr", "readr", "purrr",
  "stringr", "tibble", "ggplot2", "fs", "knitr"
)
OPT_PKGS <- c("kableExtra", "svglite", "R.utils")

ensure_pkg <- function(pkg, critical = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (AUTO_INSTALL) {
      message(sprintf("[setup] installing '%s' from CRAN ...", pkg))
      tryCatch(
        utils::install.packages(
          pkg,
          repos = "https://cloud.r-project.org",
          dependencies = TRUE,
          lib = Sys.getenv("R_LIBS_USER", unset = .libPaths()[1])
        ),
        error = function(e) message(sprintf("[setup] install error for '%s': %s", pkg, conditionMessage(e)))
      )
    }
  }
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok && critical) stop(sprintf("Required package '%s' is not available", pkg))
  ok
}

invisible(lapply(CORE_PKGS, ensure_pkg, critical = TRUE))
invisible(lapply(OPT_PKGS, ensure_pkg, critical = FALSE))

suppressPackageStartupMessages({
  for (p in CORE_PKGS) library(p, character.only = TRUE)
})

HAS_KABLE_EXTRA <- requireNamespace("kableExtra", quietly = TRUE)
HAS_SVGLITE <- requireNamespace("svglite", quietly = TRUE)
HAS_RUTILS <- requireNamespace("R.utils", quietly = TRUE)

## -----------------------------------------------------------------------------
## 2. Output folders and general helpers
## -----------------------------------------------------------------------------

DIR_TABLES <- file.path(REPORT_DIR, "tables")
DIR_FIGS   <- file.path(REPORT_DIR, "figures")
DIR_META   <- file.path(REPORT_DIR, "metadata")
DIR_OBJECT <- file.path(REPORT_DIR, "objects")

for (d in c(RESULTS_DIR, REPORT_DIR, DIR_TABLES, DIR_FIGS, DIR_META, DIR_OBJECT)) {
  fs::dir_create(d, recurse = TRUE)
}

clamp01 <- function(x, eps = 1e-12) pmin(pmax(x, eps), 1 - eps)
safe_log <- function(x) log(pmax(x, .Machine$double.xmin))

as_numeric_safely <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x)
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  sqrt(mean(x^2))
}

write_csv_atomic <- function(df, path) {
  fs::dir_create(fs::path_dir(path), recurse = TRUE)
  tmp <- paste0(path, ".tmp")
  readr::write_csv(df, tmp)
  if (file.exists(path)) unlink(path)
  ok <- file.rename(tmp, path)
  if (!ok) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

write_latex_table <- function(df, path, caption, label = NULL, digits = 4) {
  fs::dir_create(fs::path_dir(path), recurse = TRUE)
  label_clean <- if (!is.null(label)) sub("^tab:", "", label) else NULL
  df2 <- df %>% dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, digits)))
  tab <- knitr::kable(df2, format = "latex", booktabs = TRUE,
                      caption = caption, label = label_clean, escape = FALSE)
  if (HAS_KABLE_EXTRA) {
    tab <- tab %>% kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))
  }
  writeLines(tab, path)
  invisible(path)
}

theme_bench <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "plain"),
      strip.text = ggplot2::element_text(face = "plain"),
      axis.title = ggplot2::element_text(face = "plain")
    )
}

save_plot_all <- function(p, file_base, width = 7, height = 5) {
  fs::dir_create(fs::path_dir(file_base), recurse = TRUE)
  ggplot2::ggsave(paste0(file_base, ".pdf"), p, width = width, height = height,
                  device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(file_base, ".png"), p, width = width, height = height, dpi = 300)
  if (HAS_SVGLITE) {
    svglite::svglite(paste0(file_base, ".svg"), width = width, height = height)
    print(p)
    grDevices::dev.off()
  }
  invisible(file_base)
}

scenario_seed <- function(name) as.integer(MASTER_SEED + sum(utf8ToInt(name)))

make_rep_seeds <- function(n_reps, seed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed))
  as.integer(sample.int(1e9, size = n_reps, replace = FALSE))
}

is_valid_nll_value <- function(x) is.finite(x) & x < NLL_VALID_MAX

is_abnormal_message <- function(message) {
  stringr::str_detect(
    dplyr::coalesce(as.character(message), ""),
    stringr::regex("ERROR|ABNORMAL|LNSRCH|LINESEARCH|FAIL|FAILED|TIMEOUT", ignore_case = TRUE)
  )
}

regular_convergence_flag <- function(success, convergence, message = "", nll_value = NA_real_) {
  success_bool <- success %in% TRUE | tolower(as.character(success)) %in% c("true", "1", "yes", "y")
  conv_num <- suppressWarnings(as.integer(convergence))
  success_bool &
    is.finite(conv_num) & conv_num == 0L &
    !is_abnormal_message(message) &
    (is.na(nll_value) | is_valid_nll_value(nll_value))
}

with_timeout <- function(expr, timeout, label = "operation") {
  if (HAS_RUTILS && is.finite(timeout) && timeout > 0) {
    return(R.utils::withTimeout(expr, timeout = timeout, onTimeout = "error"))
  }
  expr
}

## -----------------------------------------------------------------------------
## 3. Matrix/distribution helpers and bivariate probability evaluations
## -----------------------------------------------------------------------------

Omega_mat <- function(sigmaS, sigmaY, rho) {
  matrix(c(
    sigmaS^2,              rho * sigmaS * sigmaY,
    rho * sigmaS * sigmaY, sigmaY^2
  ), nrow = 2, byrow = TRUE)
}

lambda_eff_from_eta <- function(Omega, eta) {
  Omega <- (Omega + t(Omega)) / 2
  R <- chol(Omega)
  drop(backsolve(R, eta))
}

df_mvt <- function(nu) {
  if (!is.finite(nu)) return(Inf)
  max(as.numeric(nu), 1e-6)
}

make_cdf_algorithms_2d <- function() {
  tvp <- tryCatch(
    mvtnorm::TVPACK(maxpts = CDF_MAXPTS, abseps = CDF_ABSEPS, releps = CDF_RELEPS),
    error = function(e) NULL
  )
  gb <- mvtnorm::GenzBretz(maxpts = CDF_MAXPTS, abseps = CDF_ABSEPS, releps = CDF_RELEPS)
  Filter(Negate(is.null), list(tvp, gb))
}

pmvnorm2_safe <- function(lower, upper, rho) {
  rho <- max(min(rho, 0.999999), -0.999999)
  corr <- matrix(c(1, rho, rho, 1), 2, 2)
  for (alg in make_cdf_algorithms_2d()) {
    out <- tryCatch(
      as.numeric(mvtnorm::pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = alg)),
      error = function(e) NA_real_
    )
    if (is.finite(out)) return(clamp01(out))
  }
  NA_real_
}

pmvt2_safe <- function(lower, upper, rho, nu) {
  rho <- max(min(rho, 0.999999), -0.999999)
  if (!is.finite(nu)) return(pmvnorm2_safe(lower = lower, upper = upper, rho = rho))
  corr <- matrix(c(1, rho, rho, 1), 2, 2)
  for (alg in make_cdf_algorithms_2d()) {
    out <- tryCatch(
      as.numeric(mvtnorm::pmvt(lower = lower, upper = upper, corr = corr,
                               df = df_mvt(nu), algorithm = alg)),
      error = function(e) NA_real_
    )
    if (is.finite(out)) return(clamp01(out))
  }
  NA_real_
}

PE_gauss <- function(tau_bar, eta) {
  clamp01(stats::pnorm(tau_bar / sqrt(1 + sum(eta^2))))
}

PE_t <- function(tau_bar, eta, nu) {
  if (!is.finite(nu)) return(PE_gauss(tau_bar, eta))
  clamp01(stats::pt(tau_bar / sqrt(1 + sum(eta^2)), df = df_mvt(nu)))
}

## -----------------------------------------------------------------------------
## 4. ESN/EST likelihood blocks, identity transformation G = I
## -----------------------------------------------------------------------------

pS_le0_ESN <- function(muS, sigmaY, rho, etaS, etaY, tau_bar) {
  if (!all(is.finite(c(muS, sigmaY, rho, etaS, etaY, tau_bar))) || sigmaY <= 0) return(NA_real_)
  Omega <- Omega_mat(1, sigmaY, rho)
  lam <- tryCatch(lambda_eff_from_eta(Omega, c(etaS, etaY)), error = function(e) c(NA_real_, NA_real_))
  if (anyNA(lam)) return(NA_real_)

  lamS <- lam[1]
  lamY <- lam[2]
  varW <- 1 + etaS^2 + etaY^2
  sdW <- sqrt(max(varW, 1e-12))
  covSW <- -lamS - lamY * rho * sigmaY
  rSW <- covSW / sdW
  if (!is.finite(rSW)) return(NA_real_)

  PE <- PE_gauss(tau_bar, c(etaS, etaY))
  if (!is.finite(PE) || PE <= 0) return(NA_real_)

  num <- pmvnorm2_safe(lower = c(-Inf, -Inf), upper = c(-muS, tau_bar / sdW), rho = rSW)
  if (!is.finite(num)) return(NA_real_)
  clamp01(num / PE)
}

joint_y_spos_ESN <- function(y, muS, muY, sigmaY, rho, etaS, etaY, tau_bar) {
  if (!all(is.finite(c(y, muS, muY, sigmaY, rho, etaS, etaY, tau_bar))) || sigmaY <= 0) return(NA_real_)
  Omega <- Omega_mat(1, sigmaY, rho)
  lam <- tryCatch(lambda_eff_from_eta(Omega, c(etaS, etaY)), error = function(e) c(NA_real_, NA_real_))
  if (anyNA(lam)) return(NA_real_)

  lamS <- lam[1]
  lamY <- lam[2]

  meanS_y <- muS + rho / sigmaY * (y - muY)
  varS_cond <- pmax(1 - rho^2, 1e-12)
  sdS <- sqrt(varS_cond)

  varW <- 1 + etaS^2 + etaY^2
  covYW <- -lamS * rho * sigmaY - lamY * sigmaY^2
  covSW <- -lamS - lamY * rho * sigmaY
  varW_cond <- pmax(varW - covYW^2 / sigmaY^2, 1e-12)
  sdW <- sqrt(varW_cond)
  covSW_cond <- covSW - (rho * sigmaY) * covYW / sigmaY^2
  r_cond <- covSW_cond / (sdS * sdW)
  if (!is.finite(r_cond)) return(NA_real_)

  meanW_y <- -(lamY + lamS * rho / sigmaY) * (y - muY)
  hS <- (0 - meanS_y) / sdS
  hW <- (tau_bar - meanW_y) / sdW
  if (!all(is.finite(c(hS, hW)))) return(NA_real_)

  fY <- stats::dnorm((y - muY) / sigmaY) / sigmaY
  PE <- PE_gauss(tau_bar, c(etaS, etaY))
  if (!is.finite(PE) || PE <= 0) return(NA_real_)

  p_sel_cond <- pmvnorm2_safe(lower = c(hS, -Inf), upper = c(Inf, hW), rho = r_cond)
  if (!is.finite(p_sel_cond)) return(NA_real_)
  pmax(fY * p_sel_cond / PE, .Machine$double.xmin)
}

pS_le0_EST <- function(muS, sigmaY, rho, etaS, etaY, tau_bar, nu) {
  if (!all(is.finite(c(muS, sigmaY, rho, etaS, etaY, tau_bar, nu))) || sigmaY <= 0) return(NA_real_)
  Omega <- Omega_mat(1, sigmaY, rho)
  lam <- tryCatch(lambda_eff_from_eta(Omega, c(etaS, etaY)), error = function(e) c(NA_real_, NA_real_))
  if (anyNA(lam)) return(NA_real_)

  lamS <- lam[1]
  lamY <- lam[2]
  varW <- 1 + etaS^2 + etaY^2
  sdW <- sqrt(max(varW, 1e-12))
  covSW <- -lamS - lamY * rho * sigmaY
  rSW <- covSW / sdW
  if (!is.finite(rSW)) return(NA_real_)

  PE <- PE_t(tau_bar, c(etaS, etaY), nu)
  if (!is.finite(PE) || PE <= 0) return(NA_real_)
  num <- pmvt2_safe(lower = c(-Inf, -Inf), upper = c(-muS, tau_bar / sdW), rho = rSW, nu = nu)
  if (!is.finite(num)) return(NA_real_)
  clamp01(num / PE)
}

joint_y_spos_EST <- function(y, muS, muY, sigmaY, rho, etaS, etaY, tau_bar, nu) {
  if (!all(is.finite(c(y, muS, muY, sigmaY, rho, etaS, etaY, tau_bar, nu))) || sigmaY <= 0) return(NA_real_)
  Omega <- Omega_mat(1, sigmaY, rho)
  lam <- tryCatch(lambda_eff_from_eta(Omega, c(etaS, etaY)), error = function(e) c(NA_real_, NA_real_))
  if (anyNA(lam)) return(NA_real_)

  lamS <- lam[1]
  lamY <- lam[2]
  nu_c <- df_mvt(nu)
  zY <- (y - muY) / sigmaY
  qy <- zY^2
  infl <- if (is.finite(nu_c)) (nu_c + qy) / (nu_c + 1) else 1

  meanS_y <- muS + rho / sigmaY * (y - muY)
  varS_cond <- pmax(1 - rho^2, 1e-12)

  varW <- 1 + etaS^2 + etaY^2
  covYW <- -lamS * rho * sigmaY - lamY * sigmaY^2
  covSW <- -lamS - lamY * rho * sigmaY
  varW_cond <- pmax(varW - covYW^2 / sigmaY^2, 1e-12)
  covSW_cond <- covSW - (rho * sigmaY) * covYW / sigmaY^2

  sdS <- sqrt(varS_cond * infl)
  sdW <- sqrt(varW_cond * infl)
  r_cond <- covSW_cond / sqrt(varS_cond * varW_cond)
  if (!is.finite(r_cond)) return(NA_real_)

  meanW_y <- -(lamY + lamS * rho / sigmaY) * (y - muY)
  hS <- (0 - meanS_y) / sdS
  hW <- (tau_bar - meanW_y) / sdW
  if (!all(is.finite(c(hS, hW)))) return(NA_real_)

  fY <- if (is.finite(nu_c)) stats::dt(zY, df = nu_c) / sigmaY else stats::dnorm(zY) / sigmaY
  PE <- PE_t(tau_bar, c(etaS, etaY), nu_c)
  if (!is.finite(PE) || PE <= 0) return(NA_real_)

  p_sel_cond <- pmvt2_safe(lower = c(hS, -Inf), upper = c(Inf, hW), rho = r_cond, nu = nu_c + 1)
  if (!is.finite(p_sel_cond)) return(NA_real_)
  pmax(fY * p_sel_cond / PE, .Machine$double.xmin)
}

## -----------------------------------------------------------------------------
## 5. Data generation
## -----------------------------------------------------------------------------

generate_base_data <- function(n, betaSel, betaOut, sigmaS, sigmaY, rho,
                               family = c("NH", "HT", "CN"), nu = Inf,
                               contam_eps = CONTAM_EPS, contam_kappa = CONTAM_KAPPA,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  family <- match.arg(family)

  X1 <- rnorm(n, 0, 3)
  X2 <- rbinom(n, 1, 0.5)
  Z1 <- rnorm(n, 0, 1)
  Z2 <- rpois(n, 2)

  Xsel <- cbind(1, X1, X2)
  Xout <- cbind(1, Z1, Z2)
  Omega <- Omega_mat(sigmaS, sigmaY, rho)

  if (family == "NH") {
    E <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Omega)
  } else if (family == "HT") {
    Z <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Omega)
    mix <- sqrt(df_mvt(nu) / rchisq(n, df = df_mvt(nu)))
    E <- Z * mix
  } else {
    comp <- rbinom(n, 1, contam_eps)
    scale <- ifelse(comp == 1L, sqrt(contam_kappa), 1)
    Z <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Omega)
    E <- Z * scale
  }

  muS <- drop(Xsel %*% betaSel)
  muY <- drop(Xout %*% betaOut)
  S_star <- muS + E[, 1]
  Y_star <- muY + E[, 2]
  U <- as.integer(S_star > 0)

  data.frame(U = U, Yobs = ifelse(U == 1L, Y_star, NA_real_),
             X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2,
             S_star = S_star, Y_star = Y_star)
}

generate_ESE_data <- function(n, betaSel, betaOut, sigmaS, sigmaY, rho,
                              etaS, etaY, tau_bar, family = c("ESN", "EST"),
                              nu = Inf, seed = NULL, max_tries_per_obs = 100000L) {
  if (!is.null(seed)) set.seed(seed)
  family <- match.arg(family)

  X1 <- rnorm(n, 0, 3)
  X2 <- rbinom(n, 1, 0.5)
  Z1 <- rnorm(n, 0, 1)
  Z2 <- rpois(n, 2)

  Xsel <- cbind(1, X1, X2)
  Xout <- cbind(1, Z1, Z2)
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam <- lambda_eff_from_eta(Omega, c(etaS, etaY))

  S_star <- numeric(n)
  Y_star <- numeric(n)
  accepted <- 0L
  total_tries <- 0L

  for (i in seq_len(n)) {
    mu_i <- c(sum(Xsel[i, ] * betaSel), sum(Xout[i, ] * betaOut))
    tries <- 0L
    repeat {
      if (family == "EST") {
        nu_c <- df_mvt(nu)
        scale_mix <- sqrt(nu_c / rchisq(1, df = nu_c))
        z_xy <- as.numeric(mvtnorm::rmvnorm(1, mean = c(0, 0), sigma = Omega))
        z_t <- rnorm(1)
        x <- mu_i + scale_mix * z_xy
        t0 <- scale_mix * z_t
      } else {
        x <- as.numeric(mvtnorm::rmvnorm(1, mean = mu_i, sigma = Omega))
        t0 <- rnorm(1)
      }
      L <- t0 - sum(lam * (x - mu_i))
      tries <- tries + 1L
      if (L <= tau_bar) {
        S_star[i] <- x[1]
        Y_star[i] <- x[2]
        accepted <- accepted + 1L
        total_tries <- total_tries + tries
        break
      }
      if (tries >= max_tries_per_obs) stop("Hidden-truncation sampler reached max_tries_per_obs")
    }
  }

  U <- as.integer(S_star > 0)
  df <- data.frame(U = U, Yobs = ifelse(U == 1L, Y_star, NA_real_),
                   X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2,
                   S_star = S_star, Y_star = Y_star)
  attr(df, "acceptance_rate") <- if (total_tries > 0) accepted / total_tries else NA_real_
  df
}

## -----------------------------------------------------------------------------
## 6. Scenarios
## -----------------------------------------------------------------------------

define_scenarios <- function() {
  common <- list(n = N_PER_DATASET, reps = REPS_PER_SCEN,
                 betaSel = BETA_SEL, betaOut = BETA_OUT,
                 sigmaS = SIGMA_S, sigmaY = SIGMA_Y, rho = RHO,
                 fit_models = FIT_MODELS)

  list(
    nearly_gaussian = c(common, list(label = "nearly Gaussian", dgp_family = "NH",
                                     etaS = 0, etaY = 0, tau_bar = 0, nu = Inf,
                                     contam_eps = 0, contam_kappa = 1)),
    symmetric_heavy_tails = c(common, list(label = "symmetric heavy tails", dgp_family = "HT",
                                           etaS = 0, etaY = 0, tau_bar = 0, nu = NU_HEAVY,
                                           contam_eps = 0, contam_kappa = 1)),
    contaminated_normal = c(common, list(label = "contaminated normal", dgp_family = "CN",
                                         etaS = 0, etaY = 0, tau_bar = 0, nu = Inf,
                                         contam_eps = CONTAM_EPS, contam_kappa = CONTAM_KAPPA)),
    hidden_truncation_skewness = c(common, list(label = "hidden-truncation skewness", dgp_family = "ESN",
                                                etaS = ETA_S, etaY = ETA_Y, tau_bar = TAU_BAR, nu = Inf,
                                                contam_eps = 0, contam_kappa = 1)),
    skewness_heavy_tails = c(common, list(label = "skewness and heavy tails", dgp_family = "EST",
                                          etaS = ETA_S, etaY = ETA_Y, tau_bar = TAU_BAR, nu = NU_HEAVY,
                                          contam_eps = 0, contam_kappa = 1))
  )
}

generate_scenario_data <- function(scen, seed = NULL) {
  if (scen$dgp_family %in% c("NH", "HT", "CN")) {
    return(generate_base_data(n = scen$n, betaSel = scen$betaSel, betaOut = scen$betaOut,
                              sigmaS = scen$sigmaS, sigmaY = scen$sigmaY, rho = scen$rho,
                              family = scen$dgp_family, nu = scen$nu,
                              contam_eps = scen$contam_eps, contam_kappa = scen$contam_kappa,
                              seed = seed))
  }
  generate_ESE_data(n = scen$n, betaSel = scen$betaSel, betaOut = scen$betaOut,
                    sigmaS = scen$sigmaS, sigmaY = scen$sigmaY, rho = scen$rho,
                    etaS = scen$etaS, etaY = scen$etaY, tau_bar = scen$tau_bar,
                    family = scen$dgp_family, nu = scen$nu, seed = seed)
}

## -----------------------------------------------------------------------------
## 7. Model matrices
## -----------------------------------------------------------------------------

model_matrices <- function(data) {
  Xsel <- model.matrix(~ X1 + X2, data = data)
  Xout <- model.matrix(~ Z1 + Z2, data = data)
  colnames(Xsel) <- c("intercept", "X1", "X2")
  colnames(Xout) <- c("intercept", "Z1", "Z2")
  U <- as.integer(data$U)
  Y <- as.numeric(data$Yobs)
  if (nrow(Xsel) != length(U)) stop("selection design matrix row mismatch")
  if (nrow(Xout) != length(U)) stop("outcome design matrix row mismatch")
  if (!all(U %in% c(0, 1))) stop("U must be binary")
  list(Xsel = Xsel, Xout = Xout, U = U, Y = Y)
}

## -----------------------------------------------------------------------------
## 8. Parameter counts and HeckmanEM wrappers
## -----------------------------------------------------------------------------

npar_model <- function(model, nB, nG) {
  switch(toupper(model),
         NH  = nB + nG + 2,
         HT  = nB + nG + 2,  # nu is supplied to HeckmanEM and is not recovered as an estimated parameter
         CN  = nB + nG + 4,
         ESN = nB + nG + 5,
         EST = nB + nG + 6,
         stop("Unknown model: ", model))
}

flatten_numeric <- function(x, prefix = "") {
  out <- numeric(0)
  if (is.numeric(x) || is.integer(x) || is.logical(x)) {
    vals <- as.numeric(x)
    nms <- names(x)
    if (is.null(nms)) nms <- seq_along(vals)
    names(vals) <- paste0(prefix, nms)
    return(vals)
  }
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- paste0("field", seq_along(x))
    for (nm in nms) out <- c(out, flatten_numeric(x[[nm]], paste0(prefix, nm, ".")))
  }
  out
}

get_first_named_vector <- function(obj, candidates, len = NULL) {
  for (nm in candidates) {
    if (!is.null(obj[[nm]]) && is.numeric(obj[[nm]])) {
      v <- as.numeric(obj[[nm]])
      if (is.null(len) || length(v) == len) return(v)
    }
  }
  NULL
}

get_first_scalar <- function(flat, patterns) {
  if (!length(flat)) return(NA_real_)
  nm <- names(flat)
  for (pat in patterns) {
    idx <- grep(pat, nm, ignore.case = TRUE)
    if (length(idx)) {
      val <- flat[idx[1]]
      if (is.finite(val)) return(as.numeric(val))
    }
  }
  NA_real_
}

recompute_ic_rows <- function(df) {
  if (!all(c("fit_family", "logLik", "n") %in% names(df))) return(df)

  ## This reviewer benchmark uses three selection coefficients and three outcome
  ## coefficients: intercept plus two covariates in each equation.
  nB_current <- 3L
  nG_current <- 3L

  df %>%
    dplyr::mutate(
      fit_family = toupper(as.character(fit_family)),
      k_params = vapply(
        fit_family,
        function(m) {
          tryCatch(as.numeric(npar_model(m, nB_current, nG_current)),
                   error = function(e) NA_real_)
        },
        numeric(1)
      ),
      AIC = ifelse(
        is.finite(logLik) & is.finite(k_params),
        -2 * logLik + 2 * k_params,
        NA_real_
      ),
      BIC = ifelse(
        is.finite(logLik) & is.finite(k_params) & is.finite(n),
        -2 * logLik + log(n) * k_params,
        NA_real_
      ),
      AICc = ifelse(
        is.finite(AIC) & is.finite(k_params) & is.finite(n) & n > k_params + 1,
        AIC + (2 * k_params * (k_params + 1)) / (n - k_params - 1),
        NA_real_
      )
    )
}


extract_criteria <- function(obj, model, n, nB, nG) {
  k <- npar_model(model, nB, nG)
  cr <- tryCatch(HeckmanEM::HeckmanEM.criteria(obj), error = function(e) NULL)
  flat_cr <- if (!is.null(cr)) flatten_numeric(cr) else numeric(0)
  flat_obj <- flatten_numeric(obj)
  flat_all <- c(flat_cr, flat_obj)

  raw_AIC <- get_first_scalar(flat_cr, c("^AIC$", "\\.AIC$", "AIC"))
  raw_BIC <- get_first_scalar(flat_cr, c("^BIC$", "\\.BIC$", "BIC"))

  logLik_method <- tryCatch(as.numeric(stats::logLik(obj))[1], error = function(e) NA_real_)
  logLik <- logLik_method
  source <- "stats::logLik"

  if (!is.finite(logLik)) {
    ll_field <- get_first_scalar(flat_all, c("^logLik$", "\\.logLik$", "logLik", "loglik", "LogLik", "log_lik", "ll", "LL"))
    if (is.finite(ll_field)) {
      logLik <- ll_field
      source <- "object_logLik_field"
    }
  }

  k_pkg <- NA_real_
  if (is.finite(raw_AIC) && is.finite(raw_BIC) && n > 1 && abs(log(n) - 2) > 1e-8) {
    k_pkg <- (raw_BIC - raw_AIC) / (log(n) - 2)
  }

  if (!is.finite(logLik) && is.finite(raw_AIC) && is.finite(k_pkg) && k_pkg >= 0) {
    logLik <- (2 * k_pkg - raw_AIC) / 2
    source <- "raw_AIC_BIC_implied_logLik"
  }

  AIC <- if (is.finite(logLik)) -2 * logLik + 2 * k else NA_real_
  BIC <- if (is.finite(logLik)) -2 * logLik + k * log(n) else NA_real_
  AICc <- if (is.finite(AIC) && n > k + 1) AIC + (2 * k * (k + 1)) / (n - k - 1) else NA_real_

  list(logLik = logLik, AIC = AIC, BIC = BIC, AICc = AICc, k = k,
       raw_AIC = raw_AIC, raw_BIC = raw_BIC, k_package_implied = k_pkg,
       criteria_source = source)
}

extract_heckmanem_convergence <- function(obj, crit) {
  flat <- flatten_numeric(obj)
  nms <- names(flat)
  em_iter <- get_first_scalar(flat, c("^iter$", "\\.iter$", "iterations", "niter", "n\\.iter", "numiter"))

  conv <- NA_integer_
  source <- "not_reported"

  idx_bool <- grep("converged|success", nms, ignore.case = TRUE)
  if (length(idx_bool)) {
    val <- flat[idx_bool[1]]
    if (is.finite(val)) {
      conv <- ifelse(val != 0, 0L, 1L)
      source <- nms[idx_bool[1]]
    }
  }

  if (is.na(conv)) {
    idx_code <- grep("convergence|convcode|conv\\.code|status", nms, ignore.case = TRUE)
    if (length(idx_code)) {
      val <- flat[idx_code[1]]
      if (is.finite(val)) {
        conv <- ifelse(val == 0, 0L, as.integer(abs(val)))
        source <- nms[idx_code[1]]
      }
    }
  }

  if (is.finite(em_iter) && em_iter >= HECKMANEM_ITER_MAX) {
    conv <- 1L
    source <- "reached_iter_max"
  }

  if (is.na(conv)) {
    conv <- ifelse(is.finite(crit$logLik), 0L, 99L)
    source <- ifelse(is.finite(crit$logLik), "finite_logLik_no_explicit_code", "missing_logLik")
  }

  list(convergence = as.integer(conv), conv_source = source, em_iterations = em_iter)
}

extract_heckmanem_estimates <- function(obj, model, nB, nG) {
  flat <- flatten_numeric(obj)
  betaOut <- get_first_named_vector(obj, c("beta", "Beta", "beta_hat", "Beta_hat"), nG)
  betaSel <- get_first_named_vector(obj, c("gamma", "Gamma", "gamma_hat", "Gamma_hat"), nB)

  cf <- tryCatch(stats::coef(obj), error = function(e) NULL)
  if (is.null(betaOut) && !is.null(cf) && !is.null(names(cf))) {
    idx <- grep("beta|outcome", names(cf), ignore.case = TRUE)
    if (length(idx) >= nG) betaOut <- as.numeric(cf[idx[seq_len(nG)]])
  }
  if (is.null(betaSel) && !is.null(cf) && !is.null(names(cf))) {
    idx <- grep("gamma|selection|select", names(cf), ignore.case = TRUE)
    if (length(idx) >= nB) betaSel <- as.numeric(cf[idx[seq_len(nB)]])
  }

  if (is.null(betaOut)) betaOut <- rep(NA_real_, nG)
  if (is.null(betaSel)) betaSel <- rep(NA_real_, nB)

  sigma2 <- get_first_scalar(flat, c("^sigma2$", "\\.sigma2$", "sigma2", "sigma_sq", "sigma\\.2"))
  sigma  <- get_first_scalar(flat, c("^sigma$", "\\.sigma$", "sigma"))
  sigmaY <- if (is.finite(sigma2) && sigma2 >= 0) sqrt(sigma2) else if (is.finite(sigma)) sigma else NA_real_
  rho <- get_first_scalar(flat, c("^rho$", "\\.rho$", "rho"))
  nu_hat <- if (toupper(model) == "HT") get_first_scalar(flat, c("^nu$", "\\.nu$", "df", "dof")) else NA_real_

  list(betaSel = betaSel, betaOut = betaOut, sigmaY = sigmaY, rho = rho,
       etaS = NA_real_, etaY = NA_real_, tau_bar = NA_real_, alpha = NA_real_,
       nu = nu_hat, raw_fields = paste(names(flat)[seq_len(min(30, length(flat)))], collapse = ";"))
}

empty_est <- function(nB, nG) {
  list(betaSel = rep(NA_real_, nB), betaOut = rep(NA_real_, nG), sigmaY = NA_real_, rho = NA_real_,
       etaS = NA_real_, etaY = NA_real_, tau_bar = NA_real_, alpha = NA_real_, nu = NA_real_, raw_fields = "")
}

fit_HeckmanEM_model <- function(model, data) {
  model <- toupper(model)
  mm <- model_matrices(data)
  y_na <- mm$Y
  y_zero <- y_na
  y_zero[is.na(y_zero)] <- 0
  cc <- mm$U
  nB <- ncol(mm$Xsel)
  nG <- ncol(mm$Xout)

  family_arg <- switch(model, NH = "Normal", HT = "T", CN = "CN")
  nu_arg <- switch(model, NH = 4, HT = HECKMANEM_NU_T_INIT, CN = HECKMANEM_NU_CN_INIT)

  attempts <- list(y_na, y_zero)
  last_error <- NULL
  obj <- NULL

  for (yy in attempts) {
    obj <- tryCatch(
      HeckmanEM::HeckmanEM(
        y = yy, x = mm$Xout, w = mm$Xsel, cc = cc,
        nu = nu_arg, family = family_arg,
        error = HECKMANEM_ERROR, iter.max = HECKMANEM_ITER_MAX,
        im = HECKMANEM_IM, criteria = TRUE, verbose = HECKMANEM_VERBOSE
      ),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(obj)) break
  }

  if (is.null(obj)) {
    k <- npar_model(model, nB, nG)
    return(list(model = model, success = FALSE, convergence = 99L,
                message = paste("ERROR:", last_error), logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
                AICc = NA_real_, k = k, regular_convergence = FALSE, finite_fit = FALSE,
                est = empty_est(nB, nG), criteria_source = "fit_failed", conv_source = "fit_failed",
                em_iterations = NA_real_, raw_package_AIC = NA_real_, raw_package_BIC = NA_real_,
                k_package_implied = NA_real_, start_id = NA_integer_, stage = "HeckmanEM",
                nll_value = NA_real_, all_stage1_values = ""))
  }

  crit <- extract_criteria(obj, model, length(cc), nB, nG)
  conv <- extract_heckmanem_convergence(obj, crit)
  est <- extract_heckmanem_estimates(obj, model, nB, nG)
  nll_value <- if (is.finite(crit$logLik)) -crit$logLik else NA_real_
  finite_fit <- is.finite(crit$logLik) & (is.na(nll_value) | is_valid_nll_value(nll_value))
  reg <- regular_convergence_flag(finite_fit, conv$convergence, "", nll_value)

  list(model = model, success = finite_fit, convergence = conv$convergence, message = "",
       logLik = crit$logLik, AIC = crit$AIC, BIC = crit$BIC, AICc = crit$AICc,
       k = crit$k, regular_convergence = reg, finite_fit = finite_fit,
       est = est, criteria_source = crit$criteria_source, conv_source = conv$conv_source,
       em_iterations = conv$em_iterations, raw_package_AIC = crit$raw_AIC,
       raw_package_BIC = crit$raw_BIC, k_package_implied = crit$k_package_implied,
       start_id = NA_integer_, stage = "HeckmanEM", nll_value = nll_value,
       all_stage1_values = "")
}

## -----------------------------------------------------------------------------
## 9. Proposed ESN/EST fitting
## -----------------------------------------------------------------------------

initial_values <- function(Xsel, Xout, U, Y, family) {
  nB <- ncol(Xsel)
  nG <- ncol(Xout)

  betaS0 <- rep(0, nB)
  datS <- data.frame(U = U, Xsel[, -1, drop = FALSE])
  names(datS) <- c("U", paste0("x", seq_len(nB - 1)))
  fitS <- try(stats::glm(U ~ ., data = datS, family = binomial(link = "probit")), silent = TRUE)
  if (!inherits(fitS, "try-error")) {
    cs <- stats::coef(fitS)
    betaS0[seq_along(cs)] <- ifelse(is.finite(cs), cs, 0)
  }

  betaY0 <- rep(0, nG)
  y_obs <- Y[U == 1L & !is.na(Y)]
  Xout_obs <- Xout[U == 1L & !is.na(Y), , drop = FALSE]
  sigmaY0 <- stats::sd(y_obs)

  if (length(y_obs) > nG) {
    fitY <- try(stats::lm.fit(Xout_obs, y_obs), silent = TRUE)
    if (!inherits(fitY, "try-error")) {
      betaY0 <- ifelse(is.finite(fitY$coefficients), fitY$coefficients, 0)
      res <- y_obs - drop(Xout_obs %*% betaY0)
      sigmaY0 <- stats::sd(res)
    }
  }
  if (!is.finite(sigmaY0) || sigmaY0 <= 1e-4) sigmaY0 <- 1

  base <- c(betaS0, betaY0, logSigmaY = log(sigmaY0), rho_t = atanh(0.1),
            etaS_t = 0, etaY_t = 0, tau_bar_t = 0)
  if (family == "ESN") return(base)
  if (family == "EST") return(c(base, log_nu_minus2 = log(8)))
  stop("unknown family")
}

make_start_candidates <- function(model, Xsel, Xout, U, Y, external_starts = NULL) {
  model <- toupper(model)
  base <- initial_values(Xsel, Xout, U, Y, model)

  perturb <- rbind(
    c( 0.00,  0.00,  0.00),
    c(-0.20,  0.20,  0.00),
    c( 0.20, -0.20,  0.00),
    c(-0.50,  0.50,  0.00),
    c( 0.50, -0.50,  0.00),
    c(-0.50,  0.50,  0.35),
    c(-0.50,  0.50, -0.35),
    c( 0.50, -0.50,  0.35),
    c( 0.50, -0.50, -0.35),
    c(-0.80,  0.80,  0.35),
    c( 0.80, -0.80,  0.35),
    c( 0.00,  0.00,  0.35),
    c( 0.00,  0.00, -0.35)
  )

  starts <- lapply(seq_len(nrow(perturb)), function(j) {
    st <- base
    st["etaS_t"] <- perturb[j, 1]
    st["etaY_t"] <- perturb[j, 2]
    st["tau_bar_t"] <- perturb[j, 3]
    if ("log_nu_minus2" %in% names(st)) st["log_nu_minus2"] <- log(8)
    st
  })

  if (!is.null(external_starts) && is.numeric(external_starts)) {
    st <- base
    len <- min(length(st), length(external_starts))
    st[seq_len(len)] <- external_starts[seq_len(len)]
    starts[[length(starts) + 1L]] <- st
  }

  if (model == "EST") {
    starts <- unlist(
      lapply(starts, function(st) {
        lapply(EST_NU_STARTS, function(nu0) {
          st2 <- st
          st2["log_nu_minus2"] <- log(max(nu0 - 2, 0.05))
          st2
        })
      }),
      recursive = FALSE
    )
  }

  starts
}

negloglik_ESN <- function(par, Xsel, Xout, U, Y, nB, nG) {
  betaS <- par[1:nB]
  betaY <- par[(nB + 1):(nB + nG)]
  sigmaY <- exp(par[nB + nG + 1])
  rho <- tanh(par[nB + nG + 2])
  etaS <- L_ETA * tanh(par[nB + nG + 3])
  etaY <- L_ETA * tanh(par[nB + nG + 4])
  tau_bar <- L_TAU * tanh(par[nB + nG + 5])

  if (!is.finite(sigmaY) || sigmaY <= 0) return(NLL_PENALTY)
  ll <- 0
  for (i in seq_along(U)) {
    muS <- sum(Xsel[i, ] * betaS)
    muY <- sum(Xout[i, ] * betaY)
    if (U[i] == 0L) {
      p0 <- pS_le0_ESN(muS, sigmaY, rho, etaS, etaY, tau_bar)
      if (!is.finite(p0)) return(NLL_PENALTY)
      ll <- ll + safe_log(p0)
    } else {
      if (is.na(Y[i])) return(NLL_PENALTY)
      joint <- joint_y_spos_ESN(Y[i], muS, muY, sigmaY, rho, etaS, etaY, tau_bar)
      if (!is.finite(joint)) return(NLL_PENALTY)
      ll <- ll + safe_log(joint)
    }
  }
  nll <- -ll
  if (!is.finite(nll)) NLL_PENALTY else nll
}

negloglik_EST <- function(par, Xsel, Xout, U, Y, nB, nG) {
  betaS <- par[1:nB]
  betaY <- par[(nB + 1):(nB + nG)]
  sigmaY <- exp(par[nB + nG + 1])
  rho <- tanh(par[nB + nG + 2])
  etaS <- L_ETA * tanh(par[nB + nG + 3])
  etaY <- L_ETA * tanh(par[nB + nG + 4])
  tau_bar <- L_TAU * tanh(par[nB + nG + 5])
  nu <- 2 + exp(par[nB + nG + 6])

  if (!is.finite(sigmaY) || sigmaY <= 0 || !is.finite(nu)) return(NLL_PENALTY)
  ll <- 0
  for (i in seq_along(U)) {
    muS <- sum(Xsel[i, ] * betaS)
    muY <- sum(Xout[i, ] * betaY)
    if (U[i] == 0L) {
      p0 <- pS_le0_EST(muS, sigmaY, rho, etaS, etaY, tau_bar, nu)
      if (!is.finite(p0)) return(NLL_PENALTY)
      ll <- ll + safe_log(p0)
    } else {
      if (is.na(Y[i])) return(NLL_PENALTY)
      joint <- joint_y_spos_EST(Y[i], muS, muY, sigmaY, rho, etaS, etaY, tau_bar, nu)
      if (!is.finite(joint)) return(NLL_PENALTY)
      ll <- ll + safe_log(joint)
    }
  }
  nll <- -ll
  if (!is.finite(nll)) NLL_PENALTY else nll
}

parameter_names <- function(model, nB, nG) {
  base <- c(paste0("betaSel_", seq_len(nB)), paste0("betaOut_", seq_len(nG)),
            "logSigmaY", "rho_t", "etaS_t", "etaY_t", "tau_bar_t")
  if (model == "ESN") return(base)
  if (model == "EST") return(c(base, "log_nu_minus2"))
  stop("unknown model")
}

bounds_for_model <- function(model, pnames) {
  lower <- rep(-Inf, length(pnames))
  upper <- rep( Inf, length(pnames))
  names(lower) <- names(upper) <- pnames
  lower["logSigmaY"] <- log(1e-2)
  upper["logSigmaY"] <- log(1e2)
  lower["rho_t"] <- atanh(-0.995)
  upper["rho_t"] <- atanh( 0.995)
  lower[c("etaS_t", "etaY_t", "tau_bar_t")] <- -3
  upper[c("etaS_t", "etaY_t", "tau_bar_t")] <-  3
  if (model == "EST") {
    lower["log_nu_minus2"] <- log(0.05)
    upper["log_nu_minus2"] <- log(98)
  }
  list(lower = lower, upper = upper)
}

decode_proposed <- function(model, par, nB, nG) {
  out <- list(
    betaSel = par[1:nB],
    betaOut = par[(nB + 1):(nB + nG)],
    sigmaY = exp(par[nB + nG + 1]),
    rho = tanh(par[nB + nG + 2]),
    etaS = L_ETA * tanh(par[nB + nG + 3]),
    etaY = L_ETA * tanh(par[nB + nG + 4]),
    tau_bar = L_TAU * tanh(par[nB + nG + 5]),
    alpha = L_TAU * tanh(par[nB + nG + 5]),
    nu = NA_real_,
    raw_fields = ""
  )
  if (model == "EST") out$nu <- 2 + exp(par[nB + nG + 6])
  out
}

fit_proposed_model <- function(model, data, external_starts = NULL) {
  model <- toupper(model)
  mm <- model_matrices(data)
  Xsel <- mm$Xsel
  Xout <- mm$Xout
  U <- mm$U
  Y <- mm$Y
  nB <- ncol(Xsel)
  nG <- ncol(Xout)
  pnames <- parameter_names(model, nB, nG)
  bnd <- bounds_for_model(model, pnames)
  nll <- if (model == "ESN") negloglik_ESN else negloglik_EST
  timeout <- if (model == "ESN") TIMEOUT_ESN else TIMEOUT_EST

  starts <- make_start_candidates(model, Xsel, Xout, U, Y, external_starts = external_starts)
  clean_start <- function(st) {
    names(st) <- pnames
    st <- st[pnames]
    st[!is.finite(st)] <- 0
    pmin(pmax(st, bnd$lower), bnd$upper)
  }
  starts <- lapply(starts, clean_start)

  run_optim <- function(st, maxit, stage, start_id) {
    tryCatch(
      with_timeout(
        stats::optim(par = st, fn = nll, Xsel = Xsel, Xout = Xout, U = U, Y = Y,
                     nB = nB, nG = nG, method = "L-BFGS-B",
                     lower = bnd$lower, upper = bnd$upper,
                     control = list(maxit = maxit, factr = 1e7, pgtol = 1e-6),
                     hessian = COMPUTE_HESSIAN),
        timeout = timeout,
        label = paste(model, stage, start_id)
      ),
      error = function(e) list(par = st, value = Inf, convergence = 99L,
                               message = paste("ERROR:", conditionMessage(e)))
    )
  }

  stage1 <- lapply(seq_along(starts), function(j) {
    f <- run_optim(starts[[j]], STAGE1_MAXIT, "stage1", j)
    f$stage <- "stage1"
    f$start_id <- j
    f
  })
  vals1 <- vapply(stage1, function(f) f$value, numeric(1))
  valid1 <- is_valid_nll_value(vals1)

  if (!any(valid1)) {
    return(list(model = model, success = FALSE, convergence = 99L,
                message = "all stage-1 starts failed", logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
                AICc = NA_real_, k = npar_model(model, nB, nG), regular_convergence = FALSE,
                finite_fit = FALSE, est = empty_est(nB, nG), start_id = NA_integer_, stage = "none",
                nll_value = NA_real_, all_stage1_values = paste(vals1, collapse = ";"),
                criteria_source = "proposed_likelihood", conv_source = "stage1_failed",
                em_iterations = NA_real_, raw_package_AIC = NA_real_, raw_package_BIC = NA_real_,
                k_package_implied = NA_real_))
  }

  keep_ids <- which(valid1)
  keep_ids <- keep_ids[order(vals1[keep_ids])]
  keep_ids <- keep_ids[seq_len(min(STAGE2_KEEP, length(keep_ids)))]

  stage2 <- lapply(keep_ids, function(j) {
    f <- run_optim(stage1[[j]]$par, STAGE2_MAXIT, "stage2", j)
    f$stage <- "stage2"
    f$start_id <- j
    f
  })

  fits <- c(stage2, stage1)
  vals <- vapply(fits, function(f) f$value, numeric(1))
  conv <- vapply(fits, function(f) if (is.null(f$convergence)) 99L else as.integer(f$convergence), integer(1))
  msg <- vapply(fits, function(f) if (is.null(f$message)) "" else as.character(f$message), character(1))
  valid <- is_valid_nll_value(vals)
  regular <- valid & regular_convergence_flag(TRUE, conv, msg, vals)

  if (any(regular)) {
    best_pos <- which(regular)[which.min(vals[regular])]
  } else if (any(valid)) {
    best_pos <- which(valid)[which.min(vals[valid])]
  } else {
    return(list(model = model, success = FALSE, convergence = 99L,
                message = "all starts invalid after stage 2", logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
                AICc = NA_real_, k = npar_model(model, nB, nG), regular_convergence = FALSE,
                finite_fit = FALSE, est = empty_est(nB, nG), start_id = NA_integer_, stage = "none",
                nll_value = NA_real_, all_stage1_values = paste(vals1, collapse = ";"),
                criteria_source = "proposed_likelihood", conv_source = "all_invalid",
                em_iterations = NA_real_, raw_package_AIC = NA_real_, raw_package_BIC = NA_real_,
                k_package_implied = NA_real_))
  }

  best <- fits[[best_pos]]
  nll_best <- tryCatch(nll(best$par, Xsel, Xout, U, Y, nB, nG), error = function(e) NLL_PENALTY)
  logLik <- if (is_valid_nll_value(nll_best)) -nll_best else NA_real_
  k <- npar_model(model, nB, nG)
  conv_code <- if (is.null(best$convergence)) 99L else as.integer(best$convergence)
  msg_best <- if (is.null(best$message)) "" else as.character(best$message)
  reg <- is_valid_nll_value(nll_best) & regular_convergence_flag(TRUE, conv_code, msg_best, nll_best)

  list(model = model, success = is.finite(logLik), convergence = conv_code, message = msg_best,
       logLik = logLik, AIC = if (is.finite(logLik)) -2 * logLik + 2 * k else NA_real_,
       BIC = if (is.finite(logLik)) -2 * logLik + k * log(length(U)) else NA_real_,
       AICc = if (is.finite(logLik) && length(U) > k + 1) -2 * logLik + 2 * k + (2 * k * (k + 1)) / (length(U) - k - 1) else NA_real_,
       k = k, regular_convergence = reg, finite_fit = is.finite(logLik),
       est = decode_proposed(model, best$par, nB, nG), start_id = best$start_id,
       stage = best$stage, nll_value = nll_best,
       all_stage1_values = paste(vals1, collapse = ";"),
       criteria_source = "proposed_likelihood", conv_source = best$stage,
       em_iterations = NA_real_, raw_package_AIC = NA_real_, raw_package_BIC = NA_real_,
       k_package_implied = NA_real_)
}

fit_one_model <- function(model, data, external_starts = NULL) {
  model <- toupper(model)
  if (model %in% c("NH", "HT", "CN")) return(fit_HeckmanEM_model(model, data))
  if (model %in% c("ESN", "EST")) return(fit_proposed_model(model, data, external_starts = external_starts))
  stop("Unknown model: ", model)
}

## -----------------------------------------------------------------------------
## 10. Metrics and row-level output
## -----------------------------------------------------------------------------

truth_for_scenario <- function(scen) {
  list(betaSel = scen$betaSel / scen$sigmaS,
       betaOut = scen$betaOut,
       sigmaY = scen$sigmaY,
       rho = scen$rho,
       etaS = ifelse(scen$dgp_family %in% c("ESN", "EST"), scen$etaS, NA_real_),
       etaY = ifelse(scen$dgp_family %in% c("ESN", "EST"), scen$etaY, NA_real_),
       tau_bar = ifelse(scen$dgp_family %in% c("ESN", "EST"), scen$tau_bar, NA_real_),
       alpha = ifelse(scen$dgp_family %in% c("ESN", "EST"), scen$tau_bar, NA_real_),
       nu = ifelse(scen$dgp_family %in% c("HT", "EST"), scen$nu, NA_real_))
}

metrics_row <- function(scen_name, scen, rep_id, data, fit, elapsed_sec) {
  tr <- truth_for_scenario(scen)
  est <- fit$est
  nB <- length(tr$betaSel)
  nG <- length(tr$betaOut)
  acc <- attr(data, "acceptance_rate")
  if (is.null(acc) || length(acc) == 0L || !is.finite(acc)) acc <- NA_real_

  nll_value <- ifelse(is.null(fit$nll_value), NA_real_, fit$nll_value)
  finite_fit <- isTRUE(fit$finite_fit) & is.finite(fit$logLik) & (is.na(nll_value) | is_valid_nll_value(nll_value))
  reg <- isTRUE(fit$regular_convergence) & finite_fit

  row <- list(
    scenario = scen_name,
    scenario_label = scen$label,
    dgp_family = scen$dgp_family,
    fit_family = fit$model,
    rep = rep_id,
    n = nrow(data),
    n_selected = sum(data$U == 1L),
    selection_rate = mean(data$U == 1L),
    hidden_acceptance_rate = acc,
    success = isTRUE(fit$success),
    finite_fit = finite_fit,
    regular_convergence = reg,
    convergence = ifelse(is.null(fit$convergence), 99L, fit$convergence),
    message = ifelse(is.null(fit$message), "", fit$message),
    logLik = fit$logLik,
    AIC = fit$AIC,
    BIC = fit$BIC,
    AICc = fit$AICc,
    k_params = fit$k,
    seconds = elapsed_sec,
    em_iterations = ifelse(is.null(fit$em_iterations), NA_real_, fit$em_iterations),
    raw_package_AIC = ifelse(is.null(fit$raw_package_AIC), NA_real_, fit$raw_package_AIC),
    raw_package_BIC = ifelse(is.null(fit$raw_package_BIC), NA_real_, fit$raw_package_BIC),
    k_package_implied = ifelse(is.null(fit$k_package_implied), NA_real_, fit$k_package_implied),
    start_id = ifelse(is.null(fit$start_id), NA_integer_, fit$start_id),
    stage = ifelse(is.null(fit$stage), "", fit$stage),
    all_stage1_values = ifelse(is.null(fit$all_stage1_values), "", fit$all_stage1_values),
    criteria_source = ifelse(is.null(fit$criteria_source), "", fit$criteria_source),
    conv_source = ifelse(is.null(fit$conv_source), "", fit$conv_source),
    sigmaY_true = tr$sigmaY,
    sigmaY_hat = est$sigmaY,
    sigmaY_bias = est$sigmaY - tr$sigmaY,
    rho_true = tr$rho,
    rho_hat = est$rho,
    rho_bias = est$rho - tr$rho,
    etaS_true = tr$etaS,
    etaS_hat = est$etaS,
    etaS_bias = est$etaS - tr$etaS,
    etaY_true = tr$etaY,
    etaY_hat = est$etaY,
    etaY_bias = est$etaY - tr$etaY,
    tau_bar_true = tr$tau_bar,
    tau_bar_hat = est$tau_bar,
    tau_bar_bias = est$tau_bar - tr$tau_bar,
    alpha_legacy_true = tr$alpha,
    alpha_legacy_hat = est$alpha,
    nu_true = tr$nu,
    nu_hat = est$nu,
    nu_bias = est$nu - tr$nu
  )

  for (j in seq_len(nB)) {
    row[[paste0("betaSel_", j, "_true")]] <- tr$betaSel[j]
    row[[paste0("betaSel_", j, "_hat")]] <- est$betaSel[j]
    row[[paste0("betaSel_", j, "_bias")]] <- est$betaSel[j] - tr$betaSel[j]
  }
  for (j in seq_len(nG)) {
    row[[paste0("betaOut_", j, "_true")]] <- tr$betaOut[j]
    row[[paste0("betaOut_", j, "_hat")]] <- est$betaOut[j]
    row[[paste0("betaOut_", j, "_bias")]] <- est$betaOut[j] - tr$betaOut[j]
  }

  row$BIC_minus_AIC <- row$BIC - row$AIC
  row$expected_BIC_minus_AIC <- row$k_params * (log(row$n) - 2)
  row$IC_gap_error <- row$BIC_minus_AIC - row$expected_BIC_minus_AIC

  as.data.frame(row, check.names = FALSE)
}

error_row <- function(scen_name, scen, rep_id, fit_family = NA_character_, elapsed_sec, message) {
  tr <- truth_for_scenario(scen)
  row <- list(
    scenario = scen_name,
    scenario_label = scen$label,
    dgp_family = scen$dgp_family,
    fit_family = fit_family,
    rep = rep_id,
    n = scen$n,
    success = FALSE,
    finite_fit = FALSE,
    regular_convergence = FALSE,
    convergence = 99L,
    message = paste("ERROR:", message),
    n_selected = NA_real_,
    selection_rate = NA_real_,
    hidden_acceptance_rate = NA_real_,
    logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, AICc = NA_real_,
    k_params = NA_real_, seconds = elapsed_sec,
    sigmaY_true = tr$sigmaY, rho_true = tr$rho,
    etaS_true = tr$etaS, etaY_true = tr$etaY,
    tau_bar_true = tr$tau_bar, alpha_legacy_true = tr$alpha,
    nu_true = tr$nu
  )
  as.data.frame(row, check.names = FALSE)
}

## -----------------------------------------------------------------------------
## 11. Data cleaning and summaries
## -----------------------------------------------------------------------------

ensure_benchmark_numeric_columns <- function(df) {
  numeric_cols <- c(
    "rep", "n", "n_selected", "selection_rate", "hidden_acceptance_rate",
    "logLik", "AIC", "BIC", "AICc", "k_params", "seconds", "em_iterations",
    "raw_package_AIC", "raw_package_BIC", "k_package_implied",
    "sigmaY_true", "sigmaY_hat", "sigmaY_bias",
    "rho_true", "rho_hat", "rho_bias",
    "etaS_true", "etaS_hat", "etaS_bias",
    "etaY_true", "etaY_hat", "etaY_bias",
    "tau_bar_true", "tau_bar_hat", "tau_bar_bias",
    "alpha_legacy_true", "alpha_legacy_hat",
    "nu_true", "nu_hat", "nu_bias",
    "BIC_minus_AIC", "expected_BIC_minus_AIC", "IC_gap_error",
    paste0("betaSel_", 1:20, "_true"),
    paste0("betaSel_", 1:20, "_hat"),
    paste0("betaSel_", 1:20, "_bias"),
    paste0("betaOut_", 1:20, "_true"),
    paste0("betaOut_", 1:20, "_hat"),
    paste0("betaOut_", 1:20, "_bias")
  )

  for (cc in numeric_cols) {
    if (!cc %in% names(df)) df[[cc]] <- NA_real_
    df[[cc]] <- as_numeric_safely(df[[cc]])
  }
  if (!"fit_family" %in% names(df)) df$fit_family <- NA_character_
  if (!"message" %in% names(df)) df$message <- ""
  if (!"success" %in% names(df)) df$success <- FALSE
  if (!"finite_fit" %in% names(df)) df$finite_fit <- FALSE
  if (!"regular_convergence" %in% names(df)) df$regular_convergence <- FALSE
  if (!"convergence" %in% names(df)) df$convergence <- 99L
  if (!"stage" %in% names(df)) df$stage <- ""

  df$fit_family <- toupper(as.character(df$fit_family))
  df$dgp_family <- toupper(as.character(df$dgp_family))
  df$convergence <- suppressWarnings(as.integer(df$convergence))
  df
}

regularize_replications <- function(reps) {
  reps <- ensure_benchmark_numeric_columns(reps)
  reps %>%
    dplyr::mutate(
      success_chr = tolower(as.character(success)),
      success_bool = success %in% TRUE | success_chr %in% c("true", "1", "yes", "y"),
      finite_fit = as.logical(finite_fit) & is.finite(logLik) & is_valid_nll_value(-logLik),
      regular_convergence = as.logical(regular_convergence) & finite_fit &
        is.finite(convergence) & convergence == 0L & !is_abnormal_message(message),
      converged_regular = regular_convergence,
      AIC_correct = ifelse(is.finite(logLik) & is.finite(k_params), -2 * logLik + 2 * k_params, NA_real_),
      BIC_correct = ifelse(is.finite(logLik) & is.finite(k_params), -2 * logLik + log(n) * k_params, NA_real_),
      AICc_correct = ifelse(is.finite(AIC_correct) & n > k_params + 1,
                            AIC_correct + (2 * k_params * (k_params + 1)) / (n - k_params - 1), NA_real_),
      raw_BIC_minus_AIC = BIC - AIC,
      corrected_BIC_minus_AIC = BIC_correct - AIC_correct,
      expected_BIC_minus_AIC = k_params * (log(n) - 2),
      IC_gap_error = corrected_BIC_minus_AIC - expected_BIC_minus_AIC,
      k_implied_by_raw_BIC = raw_BIC_minus_AIC / (log(n) - 2),
      raw_ic_has_wrong_k = is.finite(k_implied_by_raw_BIC) & is.finite(k_params) & abs(k_implied_by_raw_BIC - k_params) > 0.25,
      seconds_correct = seconds / TIME_DIVISOR,
      etaS_near_boundary = abs(etaS_hat) > 0.95 * L_ETA,
      etaY_near_boundary = abs(etaY_hat) > 0.95 * L_ETA,
      tau_near_boundary = abs(tau_bar_hat) > 0.95 * L_TAU
    )
}

summarize_model_comparison <- function(reps) {
  reps <- recompute_ic_rows(reps)
  reps <- regularize_replications(reps)

  summary <- reps %>%
    dplyr::filter(!is.na(fit_family), fit_family %in% c("NH", "HT", "CN", "ESN", "EST")) %>%
    dplyr::group_by(scenario, scenario_label, dgp_family, fit_family, n) %>%
    dplyr::summarise(
      reps = dplyr::n(),
      finite_fit_rate = safe_mean(as.numeric(finite_fit)),
      conv_rate = safe_mean(as.numeric(converged_regular)),
      n_converged = sum(converged_regular, na.rm = TRUE),
      abnormal_rate = safe_mean(as.numeric(is_abnormal_message(message))),
      stage1_failed_rate = safe_mean(as.numeric(stage == "none" | stringr::str_detect(message, "stage-1"))),
      eligible_for_mean_IC = FALSE,
      mean_selection_rate = safe_mean(selection_rate),
      mean_logLik = safe_mean(ifelse(converged_regular, logLik, NA_real_)),
      mean_AIC = safe_mean(ifelse(converged_regular, AIC_correct, NA_real_)),
      mean_BIC = safe_mean(ifelse(converged_regular, BIC_correct, NA_real_)),
      mean_AICc = safe_mean(ifelse(converged_regular, AICc_correct, NA_real_)),
      mean_k_params = safe_mean(ifelse(converged_regular, k_params, NA_real_)),
      mean_k_package_implied = safe_mean(ifelse(finite_fit, k_package_implied, NA_real_)),
      rmse_betaOut_1 = safe_rmse(ifelse(converged_regular, betaOut_1_bias, NA_real_)),
      rmse_betaOut_2 = safe_rmse(ifelse(converged_regular, betaOut_2_bias, NA_real_)),
      rmse_betaOut_3 = safe_rmse(ifelse(converged_regular, betaOut_3_bias, NA_real_)),
      rmse_sigmaY = safe_rmse(ifelse(converged_regular, sigmaY_bias, NA_real_)),
      rmse_rho = safe_rmse(ifelse(converged_regular, rho_bias, NA_real_)),
      rmse_etaS = safe_rmse(ifelse(converged_regular, etaS_bias, NA_real_)),
      rmse_etaY = safe_rmse(ifelse(converged_regular, etaY_bias, NA_real_)),
      rmse_tau_bar = safe_rmse(ifelse(converged_regular, tau_bar_bias, NA_real_)),
      rmse_nu = safe_rmse(ifelse(converged_regular, nu_bias, NA_real_)),
      bias_sigmaY = safe_mean(ifelse(converged_regular, sigmaY_bias, NA_real_)),
      bias_rho = safe_mean(ifelse(converged_regular, rho_bias, NA_real_)),
      mean_etaS_hat = safe_mean(ifelse(converged_regular, etaS_hat, NA_real_)),
      mean_etaY_hat = safe_mean(ifelse(converged_regular, etaY_hat, NA_real_)),
      mean_tau_bar_hat = safe_mean(ifelse(converged_regular, tau_bar_hat, NA_real_)),
      mean_nu_hat = safe_mean(ifelse(converged_regular, nu_hat, NA_real_)),
      share_etaS_near_boundary = safe_mean(ifelse(converged_regular, as.numeric(etaS_near_boundary), NA_real_)),
      share_etaY_near_boundary = safe_mean(ifelse(converged_regular, as.numeric(etaY_near_boundary), NA_real_)),
      share_tau_near_boundary = safe_mean(ifelse(converged_regular, as.numeric(tau_near_boundary), NA_real_)),
      mean_seconds = safe_mean(seconds_correct),
      sd_seconds = safe_sd(seconds_correct),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      eligible_for_mean_IC = conv_rate >= CONV_THRESHOLD_FOR_IC &
        n_converged > 0 & is.finite(mean_AIC) & is.finite(mean_BIC)
    )

  best <- best_model_table(summary)
  summary %>%
    dplyr::left_join(best, by = c("scenario", "scenario_label", "dgp_family"))
}

best_model_table <- function(summary) {
  summary %>%
    dplyr::group_by(scenario, scenario_label, dgp_family) %>%
    dplyr::summarise(
      best_AIC = {
        idx <- which(eligible_for_mean_IC & is.finite(mean_AIC))
        if (length(idx)) fit_family[idx[which.min(mean_AIC[idx])]] else NA_character_
      },
      best_AIC_conv_rate = {
        idx <- which(eligible_for_mean_IC & is.finite(mean_AIC))
        if (length(idx)) conv_rate[idx[which.min(mean_AIC[idx])]] else NA_real_
      },
      best_BIC = {
        idx <- which(eligible_for_mean_IC & is.finite(mean_BIC))
        if (length(idx)) fit_family[idx[which.min(mean_BIC[idx])]] else NA_character_
      },
      best_BIC_conv_rate = {
        idx <- which(eligible_for_mean_IC & is.finite(mean_BIC))
        if (length(idx)) conv_rate[idx[which.min(mean_BIC[idx])]] else NA_real_
      },
      .groups = "drop"
    )
}

safe_min_model <- function(df, criterion_col) {
  if (!criterion_col %in% names(df)) return(NA_character_)
  df2 <- df %>% dplyr::filter(converged_regular, is.finite(.data[[criterion_col]]))
  if (!nrow(df2)) return(NA_character_)
  df2 %>% dplyr::arrange(.data[[criterion_col]], fit_family) %>% dplyr::slice(1) %>% dplyr::pull(fit_family)
}

model_selection_frequencies <- function(reps) {
  reps <- recompute_ic_rows(reps)
  reps <- regularize_replications(reps)
  rep_counts <- reps %>%
    dplyr::group_by(scenario, scenario_label, dgp_family) %>%
    dplyr::summarise(total_reps = dplyr::n_distinct(rep), .groups = "drop")

  best_aic <- reps %>%
    dplyr::group_by(scenario, scenario_label, dgp_family, rep) %>%
    dplyr::group_modify(~ tibble::tibble(selected_model = safe_min_model(.x, "AIC_correct"))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(criterion = "AIC")

  best_bic <- reps %>%
    dplyr::group_by(scenario, scenario_label, dgp_family, rep) %>%
    dplyr::group_modify(~ tibble::tibble(selected_model = safe_min_model(.x, "BIC_correct"))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(criterion = "BIC")

  long <- dplyr::bind_rows(best_aic, best_bic) %>%
    dplyr::mutate(selected_model = ifelse(is.na(selected_model), "no_regular_fit", selected_model))

  freq_long <- long %>%
    dplyr::count(scenario, scenario_label, dgp_family, criterion, selected_model, name = "wins") %>%
    dplyr::right_join(
      rep_counts %>% tidyr::crossing(criterion = c("AIC", "BIC"), selected_model = c("NH", "HT", "CN", "ESN", "EST", "no_regular_fit")),
      by = c("scenario", "scenario_label", "dgp_family", "criterion", "selected_model")
    ) %>%
    dplyr::mutate(wins = tidyr::replace_na(wins, 0L), frequency = wins / total_reps) %>%
    dplyr::arrange(scenario, criterion, selected_model)

  freq_wide <- freq_long %>%
    dplyr::select(scenario, scenario_label, dgp_family, criterion, selected_model, frequency) %>%
    tidyr::pivot_wider(names_from = selected_model, values_from = frequency, values_fill = 0) %>%
    dplyr::arrange(scenario, criterion)

  list(long = freq_long, wide = freq_wide, best_per_rep = long)
}

ic_gap_check <- function(reps) {
  reps <- recompute_ic_rows(reps)
  reps <- regularize_replications(reps)
  reps %>%
    dplyr::filter(is.finite(AIC_correct), is.finite(BIC_correct), is.finite(k_params), is.finite(n)) %>%
    dplyr::group_by(scenario, scenario_label, fit_family) %>%
    dplyr::summarise(
      rows = dplyr::n(),
      mean_k_correct = safe_mean(k_params),
      mean_observed_gap = safe_mean(BIC_correct - AIC_correct),
      mean_expected_gap = safe_mean(k_params * (log(n) - 2)),
      max_abs_gap_error = {
        vals <- abs((BIC_correct - AIC_correct) - k_params * (log(n) - 2))
        vals <- vals[is.finite(vals)]
        if (!length(vals)) NA_real_ else max(vals)
      },
      .groups = "drop"
    )
}

bic_parameter_count_diagnostics <- function(reps) {
  reps <- regularize_replications(reps)
  reps %>%
    dplyr::filter(!is.na(fit_family)) %>%
    dplyr::group_by(fit_family) %>%
    dplyr::summarise(
      mean_k_correct = safe_mean(k_params),
      mean_k_implied_by_raw_BIC = safe_mean(k_implied_by_raw_BIC),
      mean_raw_BIC_minus_AIC = safe_mean(raw_BIC_minus_AIC),
      mean_correct_BIC_minus_AIC = safe_mean(corrected_BIC_minus_AIC),
      share_raw_ic_has_wrong_k = safe_mean(as.numeric(raw_ic_has_wrong_k)),
      .groups = "drop"
    )
}

convergence_summary <- function(reps) {
  reps <- regularize_replications(reps)
  reps %>%
    dplyr::group_by(scenario, scenario_label, dgp_family, fit_family) %>%
    dplyr::summarise(
      reps = dplyr::n(),
      finite_fit_rate = safe_mean(as.numeric(finite_fit)),
      conv_rate = safe_mean(as.numeric(converged_regular)),
      abnormal_or_error_rate = safe_mean(as.numeric(is_abnormal_message(message))),
      stage1_failed_rate = safe_mean(as.numeric(stage == "none" | stringr::str_detect(message, "stage-1"))),
      mean_seconds = safe_mean(seconds_correct),
      sd_seconds = safe_sd(seconds_correct),
      .groups = "drop"
    )
}

boundary_check <- function(reps) {
  reps <- regularize_replications(reps)
  reps %>%
    dplyr::filter(converged_regular, fit_family %in% c("ESN", "EST")) %>%
    dplyr::group_by(scenario, scenario_label, dgp_family, fit_family) %>%
    dplyr::summarise(
      rows = dplyr::n(),
      share_etaS_near_boundary = safe_mean(as.numeric(etaS_near_boundary)),
      share_etaY_near_boundary = safe_mean(as.numeric(etaY_near_boundary)),
      share_tau_near_boundary = safe_mean(as.numeric(tau_near_boundary)),
      .groups = "drop"
    )
}

## -----------------------------------------------------------------------------
## 12. Running scenarios
## -----------------------------------------------------------------------------

run_one_rep <- function(scen_name, scen, rep_id, seed) {
  t0 <- proc.time()[3]
  tryCatch({
    data <- generate_scenario_data(scen, seed = seed)
    rows <- list()
    external_starts <- NULL

    for (model in scen$fit_models) {
      t_model <- proc.time()[3]
      fit <- tryCatch(
        fit_one_model(model, data, external_starts = external_starts),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        rows[[length(rows) + 1L]] <- error_row(scen_name, scen, rep_id, fit_family = model,
                                                elapsed_sec = proc.time()[3] - t_model,
                                                message = "model fit failed without returned object")
      } else {
        rows[[length(rows) + 1L]] <- metrics_row(scen_name, scen, rep_id, data, fit,
                                                  elapsed_sec = proc.time()[3] - t_model)
        if (model == "NH" && isTRUE(fit$finite_fit) && is.numeric(fit$est$betaSel)) {
          ## Keep a simple external start vector for proposed models if wanted.
          ## Proposed models still use data-based multi-starts; this is only an
          ## additional non-oracle start derived from the same dataset.
          mm <- model_matrices(data)
          nB <- ncol(mm$Xsel)
          nG <- ncol(mm$Xout)
          external_starts <- c(fit$est$betaSel, fit$est$betaOut, log(max(fit$est$sigmaY, 1e-4)), atanh(max(min(fit$est$rho, 0.99), -0.99)))
          if (length(external_starts) < nB + nG + 2) external_starts <- NULL
        }
      }
    }
    dplyr::bind_rows(rows)
  }, error = function(e) {
    error_row(scen_name, scen, rep_id, fit_family = NA_character_,
              elapsed_sec = proc.time()[3] - t0, message = conditionMessage(e))
  })
}

run_scenario <- function(scen_name, scen, n_workers) {
  n_workers <- min(n_workers, scen$reps)
  scen_dir <- file.path(RESULTS_DIR, scen_name)
  rep_dir <- file.path(scen_dir, "per_replication")

  if (!RESUME_EXISTING && dir.exists(rep_dir)) fs::dir_delete(rep_dir)
  fs::dir_create(rep_dir, recurse = TRUE)

  message("--------------------------------------------------------------------")
  message("scenario: ", scen_name, " | dgp: ", scen$dgp_family, " | n=", scen$n, " | reps=", scen$reps)
  message("models: ", paste(scen$fit_models, collapse = ", "))
  message("workers: ", n_workers)

  rep_ids <- seq_len(scen$reps)
  rep_seeds <- make_rep_seeds(scen$reps, scenario_seed(scen_name))

  worker_fun <- function(rep_id) {
    out_path <- file.path(rep_dir, sprintf("rep_%05d.csv", rep_id))
    if (isTRUE(RESUME_EXISTING) && file.exists(out_path)) return(TRUE)
    row <- run_one_rep(scen_name, scen, rep_id, rep_seeds[rep_id])
    write_csv_atomic(row, out_path)
    TRUE
  }

  if (n_workers <= 1L) {
    invisible(lapply(rep_ids, worker_fun))
  } else {
    cl <- parallel::makeCluster(n_workers, type = "PSOCK", setup_timeout = 300)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(HeckmanEM)
        library(mvtnorm)
        library(dplyr)
        library(tidyr)
        library(readr)
        library(purrr)
        library(tibble)
        library(fs)
        library(stringr)
      })
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1"
      )
      TRUE
    })
    parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    invisible(parallel::parLapplyLB(cl, rep_ids, worker_fun))
  }

  csvs <- fs::dir_ls(rep_dir, regexp = "rep_[0-9]+\\.csv$", type = "file")
  reps <- purrr::map_dfr(csvs, ~ suppressMessages(readr::read_csv(.x, show_col_types = FALSE)))
  reps <- regularize_replications(reps)
  write_csv_atomic(reps, file.path(scen_dir, "replicates_all.csv"))

  summary <- summarize_model_comparison(reps)
  write_csv_atomic(summary, file.path(scen_dir, "summary.csv"))

  errors <- reps %>% dplyr::filter(!converged_regular)
  if (nrow(errors)) write_csv_atomic(errors, file.path(scen_dir, "diagnostic_nonregular.csv"))

  message("completed: ", scen_name)
  invisible(reps)
}

## -----------------------------------------------------------------------------
## 13. Report generation
## -----------------------------------------------------------------------------

read_all_replicates <- function(root = RESULTS_DIR) {
  files <- fs::dir_ls(root, recurse = TRUE, type = "file", regexp = "replicates_all\\.csv$")
  if (!length(files)) stop("no replicates_all.csv files found under ", root)
  purrr::map_dfr(files, ~ suppressMessages(readr::read_csv(.x, show_col_types = FALSE)))
}

make_tables <- function(reps) {
  reps <- regularize_replications(reps)
  summary <- summarize_model_comparison(reps)
  best <- best_model_table(summary)
  freq <- model_selection_frequencies(reps)
  gap <- ic_gap_check(reps)
  bicdiag <- bic_parameter_count_diagnostics(reps)
  conv <- convergence_summary(reps)
  boundary <- boundary_check(reps)

  write_csv_atomic(reps, file.path(REPORT_DIR, "benchmark_replications_all_corrected.csv"))
  write_csv_atomic(summary, file.path(DIR_TABLES, "benchmark_summary_converged_fits.csv"))
  write_csv_atomic(best, file.path(DIR_TABLES, "best_model_by_mean_corrected.csv"))
  write_csv_atomic(freq$long, file.path(DIR_TABLES, "model_selection_frequencies_long.csv"))
  write_csv_atomic(freq$wide, file.path(DIR_TABLES, "model_selection_frequencies_wide.csv"))
  write_csv_atomic(freq$best_per_rep, file.path(DIR_TABLES, "best_model_by_replication.csv"))
  write_csv_atomic(gap, file.path(DIR_TABLES, "ic_gap_check_corrected.csv"))
  write_csv_atomic(bicdiag, file.path(DIR_TABLES, "bic_parameter_count_diagnostics.csv"))
  write_csv_atomic(conv, file.path(DIR_TABLES, "convergence_summary.csv"))
  write_csv_atomic(boundary, file.path(DIR_TABLES, "boundary_check.csv"))

  compact <- summary %>%
    dplyr::select(scenario_label, dgp_family, fit_family,
                  n, reps, finite_fit_rate, conv_rate, n_converged, eligible_for_mean_IC,
                  mean_selection_rate, mean_logLik, mean_AIC, mean_BIC,
                  rmse_betaOut_1, rmse_betaOut_2, rmse_betaOut_3,
                  rmse_sigmaY, rmse_rho,
                  mean_etaS_hat, mean_etaY_hat, mean_tau_bar_hat, mean_nu_hat,
                  mean_seconds, best_AIC, best_AIC_conv_rate, best_BIC, best_BIC_conv_rate) %>%
    dplyr::arrange(scenario_label, fit_family)
  write_csv_atomic(compact, file.path(DIR_TABLES, "compact_article_table.csv"))

  write_latex_table(
    compact,
    file.path(DIR_TABLES, "compact_article_table.tex"),
    caption = paste0("Benchmark comparison of Heckman-type specifications. Mean likelihood, AIC, BIC and RMSE summaries are computed over regularly converged fits. Models are eligible for mean-IC ranking only when their regular convergence rate is at least ", CONV_THRESHOLD_FOR_IC, "."),
    label = "tab:benchmark-comparison",
    digits = 4
  )
  write_latex_table(
    best,
    file.path(DIR_TABLES, "best_model_by_mean_ic.tex"),
    caption = "Best-performing models by corrected mean information criterion among eligible regularly converged fits.",
    label = "tab:benchmark-best-mean",
    digits = 4
  )
  write_latex_table(
    freq$wide,
    file.path(DIR_TABLES, "model_selection_frequencies.tex"),
    caption = "Model-selection frequencies by replication. Non-converged fits are counted as not selected.",
    label = "tab:benchmark-selection-frequency",
    digits = 4
  )
  write_latex_table(
    gap,
    file.path(DIR_TABLES, "ic_gap_check_corrected.tex"),
    caption = "Information-criterion gap diagnostic for the benchmark comparison.",
    label = "tab:benchmark-ic-gap-check",
    digits = 6
  )
  write_latex_table(
    conv,
    file.path(DIR_TABLES, "convergence_summary.tex"),
    caption = "Convergence diagnostics for the benchmark comparison.",
    label = "tab:benchmark-convergence",
    digits = 4
  )
  if (nrow(boundary)) {
    write_latex_table(
      boundary,
      file.path(DIR_TABLES, "boundary_check.tex"),
      caption = "Boundary diagnostics for ESN/EST hidden-truncation parameters in the benchmark comparison.",
      label = "tab:benchmark-boundary-check",
      digits = 4
    )
  }

  invisible(list(reps = reps, summary = summary, best = best, freq = freq, gap = gap, conv = conv, boundary = boundary))
}

make_figures <- function(tabs) {
  summary <- tabs$summary
  freq <- tabs$freq$wide

  if (nrow(summary)) {
    p_conv <- summary %>%
      ggplot2::ggplot(ggplot2::aes(x = fit_family, y = conv_rate)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::facet_wrap(~ scenario_label) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = "model", y = "regular convergence rate") +
      theme_bench()
    save_plot_all(p_conv, file.path(DIR_FIGS, "convergence_by_scenario"), width = 9, height = 5)

    p_aic <- summary %>%
      dplyr::filter(is.finite(mean_AIC)) %>%
      ggplot2::ggplot(ggplot2::aes(x = fit_family, y = mean_AIC, shape = eligible_for_mean_IC)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::facet_wrap(~ scenario_label, scales = "free_y") +
      ggplot2::labs(x = "model", y = "mean AIC", shape = "eligible") +
      theme_bench()
    save_plot_all(p_aic, file.path(DIR_FIGS, "mean_aic_by_scenario"), width = 9, height = 6)

    rmse_long <- summary %>%
      dplyr::select(scenario_label, fit_family, rmse_betaOut_1, rmse_betaOut_2, rmse_betaOut_3, rmse_sigmaY, rmse_rho) %>%
      tidyr::pivot_longer(cols = dplyr::starts_with("rmse_"), names_to = "parameter", values_to = "rmse") %>%
      dplyr::mutate(parameter = dplyr::recode(parameter,
        rmse_betaOut_1 = "outcome intercept",
        rmse_betaOut_2 = "outcome slope 1",
        rmse_betaOut_3 = "outcome slope 2",
        rmse_sigmaY = "sigma_y",
        rmse_rho = "rho"
      )) %>%
      dplyr::filter(is.finite(rmse))
    if (nrow(rmse_long)) {
      p_rmse <- rmse_long %>%
        ggplot2::ggplot(ggplot2::aes(x = fit_family, y = rmse)) +
        ggplot2::geom_point(size = 2.4) +
        ggplot2::facet_grid(parameter ~ scenario_label, scales = "free_y") +
        ggplot2::labs(x = "model", y = "RMSE") +
        theme_bench()
      save_plot_all(p_rmse, file.path(DIR_FIGS, "rmse_by_scenario"), width = 11, height = 7)
    }
  }

  if (nrow(freq)) {
    freq_long <- freq %>%
      tidyr::pivot_longer(cols = dplyr::any_of(c("NH", "HT", "CN", "ESN", "EST", "no_regular_fit")),
                          names_to = "model", values_to = "frequency") %>%
      dplyr::filter(model != "no_regular_fit")
    p_freq <- freq_long %>%
      ggplot2::ggplot(ggplot2::aes(x = model, y = frequency)) +
      ggplot2::geom_col() +
      ggplot2::facet_grid(criterion ~ scenario_label) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = "model", y = "selection frequency") +
      theme_bench()
    save_plot_all(p_freq, file.path(DIR_FIGS, "selection_frequencies"), width = 10, height = 5.5)
  }
}

build_report <- function() {
  message(">> reading replication outputs ...")
  reps <- read_all_replicates(RESULTS_DIR)
  reps <- regularize_replications(reps)
  write_csv_atomic(reps, file.path(REPORT_DIR, "all_replications.csv"))

  message(">> building tables ...")
  tabs <- make_tables(reps)

  message(">> building figures ...")
  make_figures(tabs)

  tabs
}

## -----------------------------------------------------------------------------
## 14. Metadata and main
## -----------------------------------------------------------------------------

scenario_definitions_table <- function(scenarios) {
  purrr::imap_dfr(scenarios, function(s, nm) {
    tibble::tibble(
      scenario = nm,
      scenario_label = s$label,
      dgp_family = s$dgp_family,
      n = s$n,
      reps = s$reps,
      betaSel = paste(s$betaSel, collapse = ";"),
      betaOut = paste(s$betaOut, collapse = ";"),
      sigmaS = s$sigmaS,
      sigmaY = s$sigmaY,
      rho = s$rho,
      etaS = s$etaS,
      etaY = s$etaY,
      tau_bar = s$tau_bar,
      alpha_legacy = s$tau_bar,
      nu = s$nu,
      contam_eps = s$contam_eps,
      contam_kappa = s$contam_kappa,
      fit_models = paste(s$fit_models, collapse = ",")
    )
  })
}

write_metadata <- function(scenarios, workers) {
  fs::dir_create(DIR_META, recurse = TRUE)
  cfg <- c(
    sprintf("date: %s", as.character(Sys.time())),
    sprintf("SCRIPT_VERSION: %s", SCRIPT_VERSION),
    sprintf("PROFILE: %s", PROFILE),
    sprintf("ROOT_DIR: %s", ROOT_DIR),
    sprintf("RESULTS_DIR: %s", RESULTS_DIR),
    sprintf("REPORT_DIR: %s", REPORT_DIR),
    sprintf("MASTER_SEED: %s", MASTER_SEED),
    sprintf("REPS_PER_SCEN: %s", REPS_PER_SCEN),
    sprintf("N_PER_DATASET: %s", N_PER_DATASET),
    sprintf("MAX_WORKERS: %s", MAX_WORKERS),
    sprintf("workers used: %s", workers),
    sprintf("FIT_MODELS: %s", paste(FIT_MODELS, collapse = ",")),
    sprintf("CONV_THRESHOLD_FOR_IC: %s", CONV_THRESHOLD_FOR_IC),
    sprintf("HECKMANEM_ITER_MAX: %s", HECKMANEM_ITER_MAX),
    sprintf("HECKMANEM_ERROR: %s", HECKMANEM_ERROR),
    sprintf("EST_NU_STARTS: %s", paste(EST_NU_STARTS, collapse = ",")),
    sprintf("CDF_MAXPTS: %s", CDF_MAXPTS),
    sprintf("CDF_ABSEPS: %s", CDF_ABSEPS),
    sprintf("CDF_RELEPS: %s", CDF_RELEPS),
    sprintf("STAGE1_MAXIT: %s", STAGE1_MAXIT),
    sprintf("STAGE2_MAXIT: %s", STAGE2_MAXIT),
    sprintf("STAGE2_KEEP: %s", STAGE2_KEEP),
    sprintf("TIMEOUT_ESN: %s", TIMEOUT_ESN),
    sprintf("TIMEOUT_EST: %s", TIMEOUT_EST),
    sprintf("NLL_PENALTY: %s", NLL_PENALTY),
    sprintf("NLL_VALID_MAX: %s", NLL_VALID_MAX),
    sprintf("L_ETA: %s", L_ETA),
    sprintf("L_TAU: %s", L_TAU),
    sprintf("SLURM_JOB_ID: %s", Sys.getenv("SLURM_JOB_ID", unset = "")),
    sprintf("SLURM_CPUS_PER_TASK: %s", Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")),
    "NOTE: tau_bar is the manuscript-facing normalized effective threshold; alpha_legacy is written only as an alias.",
    "NOTE: AIC/BIC are recomputed with full model parameter counts.",
    "NOTE: mean-IC rankings exclude models below CONV_THRESHOLD_FOR_IC."
  )
  writeLines(cfg, file.path(DIR_META, "run_configuration.txt"))
  writeLines(capture.output(sessionInfo()), file.path(DIR_META, "session_info.txt"))
  write_csv_atomic(scenario_definitions_table(scenarios), file.path(DIR_META, "scenario_definitions.csv"))
  invisible(NULL)
}

main <- function() {
  logical_cores <- parallel::detectCores(logical = TRUE)
  physical_cores <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
  workers <- min(MAX_WORKERS, ifelse(is.na(logical_cores), 1L, logical_cores), REPS_PER_SCEN)

  cat("====================================================================\n")
  cat("Efficient benchmark comparison: HeckmanEM + ESN/EST\n")
  cat("--------------------------------------------------------------------\n")
  cat(sprintf("profile      : %s\n", PROFILE))
  cat(sprintf("root folder  : %s\n", ROOT_DIR))
  cat(sprintf("workers used : %s\n", workers))
  cat(sprintf("n, reps      : %s, %s per scenario\n", N_PER_DATASET, REPS_PER_SCEN))
  cat(sprintf("logical cores : %s\n", ifelse(is.na(logical_cores), "unknown", logical_cores)))
  cat(sprintf("physical cores: %s\n", ifelse(is.na(physical_cores), "unknown", physical_cores)))
  cat("====================================================================\n")

  scenarios <- define_scenarios()

  if (length(SCENARIO_FILTER) > 0L) {
    missing_scenarios <- setdiff(SCENARIO_FILTER, names(scenarios))
    if (length(missing_scenarios) > 0L) {
      stop("Unknown SCENARIO_FILTER value(s): ", paste(missing_scenarios, collapse = ", "))
    }
    scenarios <- scenarios[SCENARIO_FILTER]
  }
  write_metadata(scenarios, workers)

  if (REPORT_ONLY) {
    build_report()
    return(invisible(NULL))
  }

  for (nm in names(scenarios)) {
    run_scenario(nm, scenarios[[nm]], workers)
  }

  build_report()

  cat("====================================================================\n")
  cat("completed benchmark comparison\n")
  cat(sprintf("results: %s\n", normalizePath(RESULTS_DIR, mustWork = FALSE)))
  cat(sprintf("report : %s\n", normalizePath(REPORT_DIR, mustWork = FALSE)))
  cat("====================================================================\n")
}

if (sys.nframe() == 0) {
  main()
}
