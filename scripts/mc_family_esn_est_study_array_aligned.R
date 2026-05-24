#!/usr/bin/env Rscript
## =============================================================================
## mc_family_esn_est_study_array_aligned.R
##
## Array-based Monte Carlo recovery study for the proposed ESN/EST Heckman family.
##
## Purpose
##   Reproduce the correctly specified Monte Carlo recovery study reported in
##   the paper: data are generated from ESN or EST errors and fitted under the
##   corresponding ESN or EST likelihood.
##
## Alignment with the article
##   - The hidden-truncation threshold reported in article tables is tau_bar,
##     the normalized effective threshold in the paper.
##   - In the canonical latent normalization used for simulation and fitting,
##       mu_T = 0, Sigma_XT = 0, sigma_{T|X} = 1,
##     so the raw threshold alpha equals tau_bar. The script keeps legacy
##     alpha_* aliases only for backward compatibility with older outputs.
##   - The canonical shape eta is reported on the identified probit scale.
##     If data are generated with sigmaS_true != 1, betaSel and eta are both
##     transformed to the probit-normalized scale before computing bias/RMSE.
##   - Student-t degrees of freedom are treated as continuous. No rounding of
##     nu is used in likelihood evaluations.
##   - Invalid likelihood evaluations are marked by a large penalty but are not
##     allowed to count as finite/regular fits.
##   - Bivariate normal/t probabilities try deterministic TVPACK first and fall
##     back to GenzBretz.
##   - EST uses non-oracle generic starts for nu; the default start is not set
##     equal to the typical DGP value nu=6.
##   - Convergence summaries distinguish finite returned fits from regular
##     optimizer convergence. Abnormal line-search termination is non-regular.
##
## Not in this script
##   - Benchmark comparisons against NH, Student-t Heckman, or contaminated-
##     normal Heckman models. Those reviewer-requested comparisons are handled
##     by the separate benchmark script.
##
## HPC safeguards
##   - Uses SLURM_CPUS_PER_TASK when available.
##   - Uses PSOCK clusters with setup_timeout = 300.
##   - Writes each replication before aggregation.
##   - Supports RUN_MODE=list, RUN_MODE=scenario, RUN_MODE=report, and RUN_MODE=all.
##   - In scenario mode, SLURM_ARRAY_TASK_ID selects one scenario.
##   - REPL_TIMEOUT prevents a single difficult replication from blocking an entire scenario.
##   - RESUME_EXISTING defaults to TRUE for safe resubmission on Deucalion.
##   - Bias/RMSE summaries are computed over regularly converged fits only.
## =============================================================================

options(stringsAsFactors = FALSE)

## -----------------------------------------------------------------------------
## 0. User settings
## -----------------------------------------------------------------------------

SCRIPT_VERSION <- "2026-05-19-array-aligned-v1"

RUN_MODE <- tolower(Sys.getenv("RUN_MODE", unset = "scenario"))
RUN_MODE_ALLOWED <- c("list", "scenario", "report", "all")
if (!RUN_MODE %in% RUN_MODE_ALLOWED) {
  stop("RUN_MODE must be one of: ", paste(RUN_MODE_ALLOWED, collapse = ", "))
}

SCENARIO_INDEX <- as.integer(Sys.getenv(
  "SCENARIO_INDEX",
  unset = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1")
))

RUN_TAG <- Sys.getenv("RUN_TAG", unset = "array_final100")
ARCH_TAG <- Sys.getenv("ARCH_TAG", unset = Sys.getenv("SLURM_JOB_PARTITION", unset = "local"))
MASTER_SEED <- as.integer(Sys.getenv("MASTER_SEED", unset = "20260611"))

ROOT_DIR <- Sys.getenv(
  "ROOT_DIR",
  unset = paste0("mc_family_esn_est_A_", RUN_TAG, "_", ARCH_TAG, "_seed", MASTER_SEED)
)
RESULTS_DIR <- file.path(ROOT_DIR, "results")
REPORT_DIR  <- file.path(ROOT_DIR, "report")

## On Deucalion this is set by --cpus-per-task. Locally, it defaults to 16.
MAX_WORKERS <- as.integer(Sys.getenv(
  "MAX_WORKERS",
  unset = Sys.getenv("SLURM_CPUS_PER_TASK", unset = "16")
))

RUN_GRID_STUDY   <- as.logical(Sys.getenv("RUN_GRID_STUDY",   unset = "TRUE"))
RUN_NSWEEP_STUDY <- as.logical(Sys.getenv("RUN_NSWEEP_STUDY", unset = "TRUE"))

REPS_GRID   <- as.integer(Sys.getenv("REPS_GRID",   unset = "100"))
REPS_NSWEEP <- as.integer(Sys.getenv("REPS_NSWEEP", unset = "100"))

N_GRID <- as.integer(Sys.getenv("N_GRID", unset = "1500"))
N_SWEEP <- as.integer(strsplit(Sys.getenv("N_SWEEP", unset = "500,1000,1500,2000,2500,3000"), ",")[[1]])

GRID_RHOS       <- as.numeric(strsplit(Sys.getenv("GRID_RHOS",       unset = "0.0,0.3,0.7"), ",")[[1]])
GRID_ETA_MAGS   <- as.numeric(strsplit(Sys.getenv("GRID_ETA_MAGS",   unset = "0.3,0.6,1.0"), ",")[[1]])
GRID_SEL_LEVELS <- strsplit(Sys.getenv("GRID_SEL_LEVELS", unset = "COMMON,MODERATE,RARE"), ",")[[1]]

## Backward compatibility: if GRID_TAU_BAR is not supplied, use GRID_ALPHA.
GRID_TAU_BAR <- as.numeric(Sys.getenv("GRID_TAU_BAR", unset = Sys.getenv("GRID_ALPHA", unset = "0.4")))
GRID_EST_NU  <- as.numeric(Sys.getenv("GRID_EST_NU",  unset = "6"))

## One direction (-,+) is used by default to keep the grid interpretable.
## Set INCLUDE_OPPOSITE_SIGN_DIRECTION=TRUE to also add (+,-).
INCLUDE_OPPOSITE_SIGN_DIRECTION <- as.logical(
  Sys.getenv("INCLUDE_OPPOSITE_SIGN_DIRECTION", unset = "FALSE")
)

## Data can be generated with sigmaS_true > 1. Estimation uses the probit
## normalization sigmaS = 1. True selection coefficients and canonical eta are
## reported on the corresponding identified scale.
SIGMA_S_TRUE <- as.numeric(Sys.getenv("SIGMA_S_TRUE", unset = "4.0"))
SIGMA_Y_TRUE <- as.numeric(Sys.getenv("SIGMA_Y_TRUE", unset = "3.0"))

BETA_SEL_BASE <- c(-5.0, 0.15, 1.00)
BETA_OUT_BASE <- c( 5.0, 1.50, 0.80)
## Selection-level labels are defined by selection prevalence/stringency.
## COMMON has the largest baseline selection rate; RARE has the smallest.
## Legacy LOW/MID/HIGH names are retained only for backward compatibility and
## are mapped according to selection rate, not stringency.
SEL_INTERCEPTS <- c(
  COMMON = -3.0, MODERATE = -5.0, RARE = -7.0,
  HIGH = -3.0, MID = -5.0, LOW = -7.0
)

## Shape/threshold bounds used in eta_j = L_ETA tanh(raw_j),
## tau_bar = L_TAU tanh(raw_tau).
L_ETA <- as.numeric(Sys.getenv("L_ETA", unset = "1.5"))
L_TAU <- as.numeric(Sys.getenv("L_TAU", unset = Sys.getenv("L_ALPHA", unset = "1.5")))
L_ALPHA <- L_TAU  ## legacy alias; alpha equals tau_bar under canonical normalization

## Numerical CDF settings.
CDF_MAXPTS <- as.integer(Sys.getenv("CDF_MAXPTS", unset = "100000"))
CDF_ABSEPS <- as.numeric(Sys.getenv("CDF_ABSEPS", unset = "1e-7"))
CDF_RELEPS <- as.numeric(Sys.getenv("CDF_RELEPS", unset = "1e-7"))

## Optimizer settings.
OPT_MAXIT <- as.integer(Sys.getenv("OPT_MAXIT", unset = "2000"))
COMPUTE_HESSIAN <- as.logical(Sys.getenv("COMPUTE_HESSIAN", unset = "FALSE"))

## Starts:
##   data_multistart: data-based starts plus non-oracle perturbations.
##   truth_warm: adds a truth-centered start; use only for internal diagnostics.
START_STRATEGY <- Sys.getenv("START_STRATEGY", unset = "data_multistart")

## Resume existing per-replication files if present. For array jobs this should
## normally be TRUE, because a failed or timed-out array task can then be safely
## resubmitted without overwriting finished replications.
RESUME_EXISTING <- as.logical(Sys.getenv("RESUME_EXISTING", unset = "TRUE"))

## Per-replication timeout in seconds. A timeout writes a non-regular row instead
## of blocking the whole scenario. Set to 0 to disable. R.utils is required.
REPL_TIMEOUT <- as.integer(Sys.getenv("REPL_TIMEOUT", unset = "43200"))

## Summaries are computed over regularly converged fits; this threshold flags
## scenarios that are stable enough for direct article interpretation.
CONV_THRESHOLD_FOR_REPORT <- as.numeric(Sys.getenv("CONV_THRESHOLD_FOR_REPORT", unset = "0.80"))

AUTO_INSTALL <- as.logical(Sys.getenv("AUTO_INSTALL", unset = "FALSE"))

## Avoid BLAS/LAPACK oversubscription inside workers.
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
  "mvtnorm", "parallel", "dplyr", "tidyr", "readr", "purrr", "stringr",
  "tibble", "ggplot2", "scales", "forcats", "fs", "knitr"
)

OPT_PKGS <- c("kableExtra", "svglite", "R.utils")

install_one <- function(pkg, critical = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE) && AUTO_INSTALL) {
    message(sprintf("[setup] installing '%s' from CRAN ...", pkg))
    ok <- TRUE
    tryCatch(
      install.packages(
        pkg,
        repos = "https://cloud.r-project.org",
        dependencies = TRUE,
        lib = Sys.getenv("R_LIBS_USER", unset = .libPaths()[1])
      ),
      error = function(e) {
        ok <<- FALSE
        message(sprintf("[setup] install error for '%s': %s", pkg, conditionMessage(e)))
      }
    )
    if (requireNamespace(pkg, quietly = TRUE)) return(ok)
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (critical) {
      stop(sprintf(
        "required package '%s' is not available. Install it before submitting the job, or rerun with AUTO_INSTALL=TRUE if installation is allowed.",
        pkg
      ))
    }
    message(sprintf("optional package '%s' is not available", pkg))
    return(FALSE)
  }
  TRUE
}

invisible(lapply(CORE_PKGS, install_one, critical = TRUE))
invisible(lapply(OPT_PKGS,  install_one, critical = FALSE))

suppressPackageStartupMessages({
  for (p in CORE_PKGS) library(p, character.only = TRUE)
})

HAS_KABLE_EXTRA <- requireNamespace("kableExtra", quietly = TRUE)
HAS_SVGLITE <- requireNamespace("svglite", quietly = TRUE)
HAS_RUTILS <- requireNamespace("R.utils", quietly = TRUE)

if (START_STRATEGY == "truth_warm") {
  warning(
    "START_STRATEGY='truth_warm' uses truth-centered starts and must not be used ",
    "for article-reported Monte Carlo results. Use only for internal diagnostics."
  )
}

## -----------------------------------------------------------------------------
## 2. Output folders and helpers
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

## Large numerical penalty for invalid negative log-likelihood evaluations.
## A penalized value is finite by construction, but it must not be treated as a
## valid fit or as regular convergence.
NLL_PENALTY <- 1e12
NLL_VALID_MAX <- NLL_PENALTY / 10

is_valid_nll_value <- function(x) {
  is.finite(x) & x < NLL_VALID_MAX
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x)
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  sqrt(mean(x^2))
}

write_csv_atomic <- function(df, path) {
  fs::dir_create(fs::path_dir(path), recurse = TRUE)
  tmp <- paste0(path, ".tmp")
  readr::write_csv(df, tmp)
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  invisible(path)
}

write_rds_atomic <- function(object, path) {
  fs::dir_create(fs::path_dir(path), recurse = TRUE)
  tmp <- paste0(path, ".tmp")
  saveRDS(object, tmp)
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  invisible(path)
}

with_timeout <- function(expr, timeout, label = "operation") {
  if (HAS_RUTILS && is.finite(timeout) && timeout > 0) {
    return(R.utils::withTimeout(expr, timeout = timeout, onTimeout = "error"))
  }
  expr
}


write_latex_table <- function(df, path, caption, label = NULL, digits = 4) {
  df2 <- df %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, digits)))

  label_clean <- if (!is.null(label)) sub("^tab:", "", label) else NULL

  if (HAS_KABLE_EXTRA) {
    tab <- knitr::kable(
      df2, format = "latex", booktabs = TRUE,
      caption = caption, label = label_clean, escape = FALSE
    ) %>%
      kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))
  } else {
    tab <- knitr::kable(
      df2, format = "latex", booktabs = TRUE,
      caption = caption, label = label_clean, escape = FALSE
    )
  }

  writeLines(tab, path)
  invisible(path)
}

theme_mc <- function() {
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
  ggplot2::ggsave(
    paste0(file_base, ".pdf"), p,
    width = width, height = height, device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(paste0(file_base, ".png"), p, width = width, height = height, dpi = 300)
  if (HAS_SVGLITE) {
    svglite::svglite(paste0(file_base, ".svg"), width = width, height = height)
    print(p)
    grDevices::dev.off()
  }
  invisible(file_base)
}

scenario_seed <- function(name) {
  as.integer(MASTER_SEED + sum(utf8ToInt(name)))
}

make_rep_seeds <- function(n_reps, seed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed))
  as.integer(sample.int(1e9, size = n_reps, replace = FALSE))
}

fmt_num <- function(x) {
  gsub("\\+", "p", gsub("-", "m", gsub("\\.", "", as.character(x))))
}

regular_convergence_flag_vec <- function(convergence, message = "") {
  conv_num <- suppressWarnings(as.integer(convergence))
  msg <- dplyr::coalesce(as.character(message), "")
  is.finite(conv_num) &
    conv_num == 0L &
    !stringr::str_detect(
      msg,
      stringr::regex("ERROR|ABNORMAL|LNSRCH|LINESEARCH|FAIL|FAILED", ignore_case = TRUE)
    )
}

## -----------------------------------------------------------------------------
## 3. Elliptical helpers and generator-aware CDF blocks
## -----------------------------------------------------------------------------

Omega_mat <- function(sigmaS, sigmaY, rho) {
  matrix(
    c(
      sigmaS^2,              rho * sigmaS * sigmaY,
      rho * sigmaS * sigmaY, sigmaY^2
    ),
    nrow = 2, byrow = TRUE
  )
}

## Omega = R'R, with R = chol(Omega). If eta is the canonical shape vector,
## lambda_eff solves R lambda_eff = eta.
lambda_eff_from_eta <- function(Omega, eta) {
  Omega <- (Omega + t(Omega)) / 2
  R <- chol(Omega)
  drop(backsolve(R, eta))
}

eta_from_lambda <- function(Omega, lambda) {
  Omega <- (Omega + t(Omega)) / 2
  R <- chol(Omega)
  drop(R %*% lambda)
}

normalize_eta_to_probit_scale <- function(sigmaS, sigmaY, rho, eta_raw) {
  ## If raw errors are e_raw = (eps_S, eps_Y), and the identified probit-scale
  ## errors are e_norm = (eps_S/sigmaS, eps_Y), then e_raw = A e_norm with
  ## A = diag(sigmaS, 1). The same hidden-truncation index satisfies
  ## lambda_raw' e_raw = lambda_norm' e_norm, where lambda_norm = A' lambda_raw.
  ## We then convert lambda_norm to the canonical eta for Omega_norm.
  Omega_raw <- Omega_mat(sigmaS, sigmaY, rho)
  lambda_raw <- lambda_eff_from_eta(Omega_raw, eta_raw)
  A <- diag(c(sigmaS, 1), nrow = 2)
  lambda_norm <- drop(t(A) %*% lambda_raw)
  Omega_norm <- Omega_mat(1, sigmaY, rho)
  eta_norm <- eta_from_lambda(Omega_norm, lambda_norm)
  eta_norm
}

## Continuous degrees of freedom. Do not round nu.
df_mvt <- function(nu) {
  if (!is.finite(nu)) return(1e6)
  max(1e-6, as.numeric(nu))
}

make_cdf_algorithms_2d <- function() {
  alg_tvp <- tryCatch(
    mvtnorm::TVPACK(
      maxpts = CDF_MAXPTS,
      abseps = CDF_ABSEPS,
      releps = CDF_RELEPS
    ),
    error = function(e) NULL
  )

  alg_gb <- mvtnorm::GenzBretz(
    maxpts = CDF_MAXPTS,
    abseps = CDF_ABSEPS,
    releps = CDF_RELEPS
  )

  Filter(Negate(is.null), list(alg_tvp, alg_gb))
}

pmvnorm2_safe <- function(lower, upper, rho) {
  rho <- max(min(rho, 0.999999), -0.999999)
  corr <- matrix(c(1, rho, rho, 1), 2, 2)

  for (alg in make_cdf_algorithms_2d()) {
    out <- tryCatch(
      as.numeric(mvtnorm::pmvnorm(
        lower = lower,
        upper = upper,
        corr = corr,
        algorithm = alg
      )),
      error = function(e) NA_real_
    )
    if (is.finite(out)) return(clamp01(out))
  }

  NA_real_
}

pmvt2_safe <- function(lower, upper, rho, nu) {
  rho <- max(min(rho, 0.999999), -0.999999)

  if (!is.finite(nu)) {
    return(pmvnorm2_safe(lower = lower, upper = upper, rho = rho))
  }

  corr <- matrix(c(1, rho, rho, 1), 2, 2)

  for (alg in make_cdf_algorithms_2d()) {
    out <- tryCatch(
      as.numeric(mvtnorm::pmvt(
        lower = lower,
        upper = upper,
        corr = corr,
        df = df_mvt(nu),
        algorithm = alg
      )),
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
## 4. ESN/EST likelihood blocks
## -----------------------------------------------------------------------------

pS_le0_ESN <- function(muS, sigmaY, rho, etaS, etaY, tau_bar) {
  if (!all(is.finite(c(muS, sigmaY, rho, etaS, etaY, tau_bar))) || sigmaY <= 0) {
    return(NA_real_)
  }

  Omega <- Omega_mat(1, sigmaY, rho)
  lam <- tryCatch(lambda_eff_from_eta(Omega, c(etaS, etaY)),
                  error = function(e) c(NA_real_, NA_real_))
  if (anyNA(lam)) return(NA_real_)

  lamS <- lam[1]
  lamY <- lam[2]

  varW <- 1 + etaS^2 + etaY^2
  sdW <- sqrt(max(varW, 1e-12))
  covSW <- -lamS - lamY * rho * sigmaY
  rSW <- covSW / sdW
  if (!is.finite(rSW)) return(NA_real_)

  hS <- -muS
  hW <- tau_bar / sdW
  PE <- PE_gauss(tau_bar, c(etaS, etaY))
  if (!is.finite(PE) || PE <= 0) return(NA_real_)

  num <- pmvnorm2_safe(lower = c(-Inf, -Inf), upper = c(hS, hW), rho = rSW)
  if (!is.finite(num)) return(NA_real_)

  clamp01(num / PE)
}

joint_y_spos_ESN <- function(y, muS, muY, sigmaY, rho, etaS, etaY, tau_bar) {
  if (!all(is.finite(c(y, muS, muY, sigmaY, rho, etaS, etaY, tau_bar))) || sigmaY <= 0) {
    return(NA_real_)
  }

  Omega <- Omega_mat(1, sigmaY, rho)
  lam <- tryCatch(lambda_eff_from_eta(Omega, c(etaS, etaY)),
                  error = function(e) c(NA_real_, NA_real_))
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
  if (!all(is.finite(c(muS, sigmaY, rho, etaS, etaY, tau_bar, nu))) || sigmaY <= 0) {
    return(NA_real_)
  }

  Omega <- Omega_mat(1, sigmaY, rho)
  lam <- tryCatch(lambda_eff_from_eta(Omega, c(etaS, etaY)),
                  error = function(e) c(NA_real_, NA_real_))
  if (anyNA(lam)) return(NA_real_)

  lamS <- lam[1]
  lamY <- lam[2]

  varW <- 1 + etaS^2 + etaY^2
  sdW <- sqrt(max(varW, 1e-12))
  covSW <- -lamS - lamY * rho * sigmaY
  rSW <- covSW / sdW
  if (!is.finite(rSW)) return(NA_real_)

  hS <- -muS
  hW <- tau_bar / sdW
  PE <- PE_t(tau_bar, c(etaS, etaY), nu)
  if (!is.finite(PE) || PE <= 0) return(NA_real_)

  num <- pmvt2_safe(lower = c(-Inf, -Inf), upper = c(hS, hW), rho = rSW, nu = nu)
  if (!is.finite(num)) return(NA_real_)

  clamp01(num / PE)
}

joint_y_spos_EST <- function(y, muS, muY, sigmaY, rho, etaS, etaY, tau_bar, nu) {
  if (!all(is.finite(c(y, muS, muY, sigmaY, rho, etaS, etaY, tau_bar, nu))) || sigmaY <= 0) {
    return(NA_real_)
  }

  Omega <- Omega_mat(1, sigmaY, rho)
  lam <- tryCatch(lambda_eff_from_eta(Omega, c(etaS, etaY)),
                  error = function(e) c(NA_real_, NA_real_))
  if (anyNA(lam)) return(NA_real_)

  lamS <- lam[1]
  lamY <- lam[2]

  zY <- (y - muY) / sigmaY
  qy <- zY^2
  infl <- (df_mvt(nu) + qy) / (df_mvt(nu) + 1)

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

  fY <- stats::dt(zY, df = df_mvt(nu)) / sigmaY
  PE <- PE_t(tau_bar, c(etaS, etaY), nu)
  if (!is.finite(PE) || PE <= 0) return(NA_real_)

  p_sel_cond <- pmvt2_safe(
    lower = c(hS, -Inf),
    upper = c(Inf, hW),
    rho = r_cond,
    nu = df_mvt(nu) + 1
  )
  if (!is.finite(p_sel_cond)) return(NA_real_)

  pmax(fY * p_sel_cond / PE, .Machine$double.xmin)
}

## -----------------------------------------------------------------------------
## 5. Data generation
## -----------------------------------------------------------------------------

generate_data_ESN <- function(scen, seed = NULL, max_tries_per_obs = 100000L) {
  if (!is.null(seed)) set.seed(seed)

  n <- scen$n
  X1 <- stats::rnorm(n, 0, 3)
  X2 <- stats::rbinom(n, 1, 0.5)
  Z1 <- stats::rnorm(n, 0, 1)
  Z2 <- stats::rpois(n, 2)

  Xsel <- cbind(1, X1, X2)
  Xout <- cbind(1, Z1, Z2)

  Omega <- Omega_mat(scen$sigmaS, scen$sigmaY, scen$rho)
  lam <- lambda_eff_from_eta(Omega, c(scen$etaS_raw, scen$etaY_raw))

  S_star <- numeric(n)
  Y_star <- numeric(n)
  accepted <- 0L
  total_tries <- 0L

  for (i in seq_len(n)) {
    mu_i <- c(sum(Xsel[i, ] * scen$betaSel), sum(Xout[i, ] * scen$betaOut))
    tries <- 0L

    repeat {
      x <- as.numeric(mvtnorm::rmvnorm(1, mean = mu_i, sigma = Omega))
      t0 <- stats::rnorm(1)
      L <- t0 - sum(lam * (x - mu_i))
      tries <- tries + 1L

      if (L <= scen$tau_bar) {
        S_star[i] <- x[1]
        Y_star[i] <- x[2]
        accepted <- accepted + 1L
        total_tries <- total_tries + tries
        break
      }

      if (tries >= max_tries_per_obs) {
        stop("ESN generation reached max_tries_per_obs")
      }
    }
  }

  U <- as.integer(S_star > 0)
  Yobs <- ifelse(U == 1L, Y_star, NA_real_)

  df <- data.frame(
    U = U, Yobs = Yobs,
    X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2,
    S_star = S_star, Y_star = Y_star
  )
  attr(df, "acceptance_rate") <- if (total_tries > 0) accepted / total_tries else NA_real_
  df
}

generate_data_EST <- function(scen, seed = NULL, max_tries_per_obs = 100000L) {
  if (!is.null(seed)) set.seed(seed)

  n <- scen$n
  X1 <- stats::rnorm(n, 0, 3)
  X2 <- stats::rbinom(n, 1, 0.5)
  Z1 <- stats::rnorm(n, 0, 1)
  Z2 <- stats::rpois(n, 2)

  Xsel <- cbind(1, X1, X2)
  Xout <- cbind(1, Z1, Z2)

  Omega <- Omega_mat(scen$sigmaS, scen$sigmaY, scen$rho)
  lam <- lambda_eff_from_eta(Omega, c(scen$etaS_raw, scen$etaY_raw))

  S_star <- numeric(n)
  Y_star <- numeric(n)
  accepted <- 0L
  total_tries <- 0L

  for (i in seq_len(n)) {
    mu_i <- c(sum(Xsel[i, ] * scen$betaSel), sum(Xout[i, ] * scen$betaOut))
    tries <- 0L

    repeat {
      scale_mix <- sqrt(scen$nu / stats::rchisq(1, df = df_mvt(scen$nu)))
      z_xy <- as.numeric(mvtnorm::rmvnorm(1, mean = c(0, 0), sigma = Omega))
      z_t <- stats::rnorm(1)
      x <- mu_i + scale_mix * z_xy
      t0 <- scale_mix * z_t
      L <- t0 - sum(lam * (x - mu_i))
      tries <- tries + 1L

      if (L <= scen$tau_bar) {
        S_star[i] <- x[1]
        Y_star[i] <- x[2]
        accepted <- accepted + 1L
        total_tries <- total_tries + tries
        break
      }

      if (tries >= max_tries_per_obs) {
        stop("EST generation reached max_tries_per_obs")
      }
    }
  }

  U <- as.integer(S_star > 0)
  Yobs <- ifelse(U == 1L, Y_star, NA_real_)

  df <- data.frame(
    U = U, Yobs = Yobs,
    X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2,
    S_star = S_star, Y_star = Y_star
  )
  attr(df, "acceptance_rate") <- if (total_tries > 0) accepted / total_tries else NA_real_
  df
}

generate_data <- function(scen, seed = NULL) {
  if (scen$family == "ESN") return(generate_data_ESN(scen, seed = seed))
  if (scen$family == "EST") return(generate_data_EST(scen, seed = seed))
  stop("unknown scenario family: ", scen$family)
}

## -----------------------------------------------------------------------------
## 6. Scenarios
## -----------------------------------------------------------------------------

make_scenario <- function(
  experiment, family, rho, etaS, etaY, tau_bar, sel_level, n, reps,
  nu = NA_real_, name_suffix = NULL
) {
  betaSel <- BETA_SEL_BASE
  betaSel[1] <- SEL_INTERCEPTS[[sel_level]]

  eta_ident <- normalize_eta_to_probit_scale(
    sigmaS = SIGMA_S_TRUE,
    sigmaY = SIGMA_Y_TRUE,
    rho = rho,
    eta_raw = c(etaS, etaY)
  )

  parts <- c(
    experiment,
    tolower(family),
    paste0("r", fmt_num(sprintf("%.2f", rho))),
    paste0("etaSraw", fmt_num(sprintf("%+.1f", etaS))),
    paste0("etaYraw", fmt_num(sprintf("%+.1f", etaY))),
    paste0("tau", fmt_num(sprintf("%.1f", tau_bar))),
    paste0("sel", sel_level)
  )

  if (family == "EST") {
    parts <- c(parts, paste0("nu", fmt_num(nu)))
  }

  if (!is.null(name_suffix) && length(name_suffix) > 0 && !is.na(name_suffix)) {
    parts <- c(parts, as.character(name_suffix))
  }

  name <- paste(parts, collapse = "_")

  list(
    name = name,
    experiment = experiment,
    family = family,
    n = as.integer(n),
    reps = as.integer(reps),
    betaSel = betaSel,
    betaOut = BETA_OUT_BASE,
    sigmaS = SIGMA_S_TRUE,
    sigmaY = SIGMA_Y_TRUE,
    rho = rho,
    etaS_raw = etaS,
    etaY_raw = etaY,
    etaS_identified = eta_ident[1],
    etaY_identified = eta_ident[2],
    tau_bar = tau_bar,
    alpha = tau_bar,  ## legacy alias, canonical alpha = tau_bar
    nu = if (family == "EST") nu else NA_real_,
    L_eta = L_ETA,
    L_tau = L_TAU,
    selection_level = sel_level
  )
}

define_grid_scenarios <- function() {
  sign_pairs <- list(c(-1, 1))
  if (INCLUDE_OPPOSITE_SIGN_DIRECTION) sign_pairs <- c(sign_pairs, list(c(1, -1)))

  out <- list()

  for (rho in GRID_RHOS) {
    for (eta_mag in GRID_ETA_MAGS) {
      for (sp in sign_pairs) {
        for (sel in GRID_SEL_LEVELS) {
          etaS <- sp[1] * eta_mag
          etaY <- sp[2] * eta_mag

          sc_esn <- make_scenario(
            experiment = "grid_n1500",
            family = "ESN",
            rho = rho, etaS = etaS, etaY = etaY,
            tau_bar = GRID_TAU_BAR, sel_level = sel,
            n = N_GRID, reps = REPS_GRID
          )
          out[[sc_esn$name]] <- sc_esn

          sc_est <- make_scenario(
            experiment = "grid_n1500",
            family = "EST",
            rho = rho, etaS = etaS, etaY = etaY,
            tau_bar = GRID_TAU_BAR, sel_level = sel,
            n = N_GRID, reps = REPS_GRID,
            nu = GRID_EST_NU
          )
          out[[sc_est$name]] <- sc_est
        }
      }
    }
  }

  out
}

define_nsweep_scenarios <- function() {
  out <- list()

  for (nval in N_SWEEP) {
    sc_esn <- make_scenario(
      experiment = "n_sweep",
      family = "ESN",
      rho = 0.0, etaS = -0.3, etaY = 0.3,
      tau_bar = GRID_TAU_BAR, sel_level = "COMMON",
      n = nval, reps = REPS_NSWEEP,
      name_suffix = paste0("n", nval)
    )
    out[[sc_esn$name]] <- sc_esn

    sc_est <- make_scenario(
      experiment = "n_sweep",
      family = "EST",
      rho = 0.0, etaS = -0.3, etaY = 0.3,
      tau_bar = GRID_TAU_BAR, sel_level = "COMMON",
      n = nval, reps = REPS_NSWEEP,
      nu = GRID_EST_NU,
      name_suffix = paste0("n", nval)
    )
    out[[sc_est$name]] <- sc_est
  }

  out
}

define_all_scenarios <- function() {
  scs <- list()
  if (RUN_GRID_STUDY)   scs <- c(scs, define_grid_scenarios())
  if (RUN_NSWEEP_STUDY) scs <- c(scs, define_nsweep_scenarios())
  scs
}

## -----------------------------------------------------------------------------
## 7. Model matrices, starts, likelihoods, and fitting
## -----------------------------------------------------------------------------

model_matrices <- function(data) {
  Xsel <- stats::model.matrix(~ X1 + X2, data = data, na.action = stats::na.pass)
  Xout <- stats::model.matrix(~ Z1 + Z2, data = data, na.action = stats::na.pass)
  U <- as.integer(data$U)
  Y <- as.numeric(data$Yobs)

  if (nrow(Xsel) != length(U)) stop("selection design matrix row mismatch")
  if (nrow(Xout) != length(U)) stop("outcome design matrix row mismatch")
  if (!all(U %in% c(0, 1))) stop("U must be binary")

  list(Xsel = Xsel, Xout = Xout, U = U, Y = Y)
}

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

  base <- c(
    betaS0,
    betaY0,
    logSigmaY = log(sigmaY0),
    rho_t = atanh(0.1),
    etaS_t = 0,
    etaY_t = 0,
    tau_bar_t = 0
  )

  if (family == "ESN") return(base)
  if (family == "EST") return(c(base, log_nu_minus2 = log(8)))
  stop("unknown family")
}

raw_from_bounded <- function(value, bound) {
  value <- max(min(value, bound * 0.999), -bound * 0.999)
  atanh(value / bound)
}

make_start_candidates <- function(scen, Xsel, Xout, U, Y) {
  base <- initial_values(Xsel, Xout, U, Y, scen$family)

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
    st["etaS_t"]    <- perturb[j, 1]
    st["etaY_t"]    <- perturb[j, 2]
    st["tau_bar_t"] <- perturb[j, 3]
    if ("log_nu_minus2" %in% names(st)) st["log_nu_minus2"] <- log(8)
    st
  })

  ## Non-oracle generic starts for Student-t degrees of freedom. These values
  ## avoid starting exactly at the common DGP value nu=6.
  if (scen$family == "EST") {
    nu_start_values <- c(4, 10, 20)
    starts <- unlist(
      lapply(starts, function(st) {
        lapply(nu_start_values, function(nu0) {
          st2 <- st
          st2["log_nu_minus2"] <- log(max(nu0 - 2, 0.05))
          st2
        })
      }),
      recursive = FALSE
    )
  }

  if (START_STRATEGY == "truth_warm") {
    st <- base
    st["etaS_t"]    <- raw_from_bounded(scen$etaS_identified, L_ETA)
    st["etaY_t"]    <- raw_from_bounded(scen$etaY_identified, L_ETA)
    st["tau_bar_t"] <- raw_from_bounded(scen$tau_bar, L_TAU)
    st["rho_t"]     <- atanh(max(min(scen$rho, 0.995), -0.995))
    if ("log_nu_minus2" %in% names(st)) {
      st["log_nu_minus2"] <- log(max(scen$nu - 2, 0.05))
    }
    starts[[length(starts) + 1L]] <- st
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

  -ll
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

  -ll
}

parameter_names <- function(family, nB, nG) {
  base <- c(
    paste0("betaSel_", seq_len(nB)),
    paste0("betaOut_", seq_len(nG)),
    "logSigmaY", "rho_t", "etaS_t", "etaY_t", "tau_bar_t"
  )
  if (family == "ESN") return(base)
  if (family == "EST") return(c(base, "log_nu_minus2"))
  stop("unknown family")
}

bounds_for_family <- function(family, pnames) {
  lower <- rep(-Inf, length(pnames))
  upper <- rep( Inf, length(pnames))
  names(lower) <- names(upper) <- pnames

  lower["logSigmaY"] <- log(1e-2)
  upper["logSigmaY"] <- log(1e2)

  lower["rho_t"] <- atanh(-0.995)
  upper["rho_t"] <- atanh( 0.995)

  lower[c("etaS_t", "etaY_t", "tau_bar_t")] <- -3
  upper[c("etaS_t", "etaY_t", "tau_bar_t")] <-  3

  if (family == "EST") {
    lower["log_nu_minus2"] <- log(0.05)
    upper["log_nu_minus2"] <- log(98)
  }

  list(lower = lower, upper = upper)
}

decode_estimates <- function(family, par, nB, nG) {
  out <- list(
    betaSel = par[1:nB],
    betaOut = par[(nB + 1):(nB + nG)],
    sigmaY = exp(par[nB + nG + 1]),
    rho = tanh(par[nB + nG + 2]),
    etaS = L_ETA * tanh(par[nB + nG + 3]),
    etaY = L_ETA * tanh(par[nB + nG + 4]),
    tau_bar = L_TAU * tanh(par[nB + nG + 5]),
    alpha = L_TAU * tanh(par[nB + nG + 5]),  ## legacy alias
    nu = NA_real_
  )

  if (family == "EST") out$nu <- 2 + exp(par[nB + nG + 6])
  out
}

npar_family <- function(family, nB, nG) {
  if (family == "ESN") return(nB + nG + 5)
  if (family == "EST") return(nB + nG + 6)
  stop("unknown family")
}

empty_fit_result <- function(scen, nB, nG, pnames, message) {
  list(
    par = rep(NA_real_, length(pnames)),
    est = list(
      betaSel = rep(NA_real_, nB),
      betaOut = rep(NA_real_, nG),
      sigmaY = NA_real_,
      rho = NA_real_,
      etaS = NA_real_,
      etaY = NA_real_,
      tau_bar = NA_real_,
      alpha = NA_real_,
      nu = NA_real_
    ),
    logLik = NA_real_,
    aic = NA_real_,
    bic = NA_real_,
    convergence = 99L,
    message = message,
    regular_convergence = FALSE,
    start_id = NA_integer_,
    start_values = NA_character_,
    mm = NULL
  )
}

fit_family <- function(scen, data) {
  mm <- model_matrices(data)
  Xsel <- mm$Xsel
  Xout <- mm$Xout
  U <- mm$U
  Y <- mm$Y
  nB <- ncol(Xsel)
  nG <- ncol(Xout)

  family <- scen$family
  pnames <- parameter_names(family, nB, nG)
  bnd <- bounds_for_family(family, pnames)
  nll <- if (family == "ESN") negloglik_ESN else negloglik_EST

  starts <- make_start_candidates(scen, Xsel, Xout, U, Y)

  fits <- lapply(starts, function(st) {
    names(st) <- pnames
    st[!is.finite(st)] <- 0
    st <- pmin(pmax(st, bnd$lower), bnd$upper)

    tryCatch(
      stats::optim(
        par = st,
        fn = nll,
        Xsel = Xsel, Xout = Xout, U = U, Y = Y,
        nB = nB, nG = nG,
        method = "L-BFGS-B",
        lower = bnd$lower,
        upper = bnd$upper,
        control = list(maxit = OPT_MAXIT, factr = 1e7, pgtol = 1e-6),
        hessian = COMPUTE_HESSIAN
      ),
      error = function(e) {
        list(par = st, value = Inf, convergence = 99L, message = conditionMessage(e))
      }
    )
  })

  values <- vapply(fits, function(f) f$value, numeric(1))
  conv <- vapply(fits, function(f) {
    if (is.null(f$convergence)) 99L else as.integer(f$convergence)
  }, integer(1))
  messages <- vapply(fits, function(f) {
    if (is.null(f$message)) "" else as.character(f$message)
  }, character(1))

  valid_value <- is_valid_nll_value(values)
  regular <- valid_value & regular_convergence_flag_vec(conv, messages)
  finite <- valid_value

  if (any(regular)) {
    best_id <- which(regular)[which.min(values[regular])]
  } else if (any(finite)) {
    best_id <- which(finite)[which.min(values[finite])]
  } else {
    return(empty_fit_result(scen, nB, nG, pnames, "all multi-start fits failed"))
  }

  best <- fits[[best_id]]
  logLik <- -nll(best$par, Xsel, Xout, U, Y, nB, nG)
  k <- npar_family(family, nB, nG)
  msg <- ifelse(is.null(best$message), "", best$message)
  conv_code <- ifelse(is.null(best$convergence), 99L, as.integer(best$convergence))
  reg <- regular_convergence_flag_vec(conv_code, msg) &
    is.finite(logLik) &
    is_valid_nll_value(-logLik)

  list(
    par = best$par,
    est = decode_estimates(family, best$par, nB, nG),
    logLik = logLik,
    aic = -2 * logLik + 2 * k,
    bic = -2 * logLik + k * log(length(U)),
    k = k,
    convergence = conv_code,
    message = msg,
    regular_convergence = reg,
    start_id = best_id,
    start_values = paste(round(values, 6), collapse = ";"),
    mm = mm
  )
}

## -----------------------------------------------------------------------------
## 8. Truth, metrics, and summaries
## -----------------------------------------------------------------------------

truth_for_scenario <- function(scen) {
  eta_ident <- c(scen$etaS_identified, scen$etaY_identified)
  eta_raw <- c(scen$etaS_raw, scen$etaY_raw)

  list(
    betaSel = scen$betaSel / scen$sigmaS,
    betaOut = scen$betaOut,
    sigmaY = scen$sigmaY,
    rho = scen$rho,
    etaS = eta_ident[1],
    etaY = eta_ident[2],
    etaS_raw = eta_raw[1],
    etaY_raw = eta_raw[2],
    tau_bar = scen$tau_bar,
    alpha = scen$tau_bar,
    nu = ifelse(scen$family == "EST", scen$nu, NA_real_),
    PE = if (scen$family == "ESN") {
      PE_gauss(scen$tau_bar, eta_raw)
    } else {
      PE_t(scen$tau_bar, eta_raw, scen$nu)
    }
  )
}

metrics_row <- function(scen, rep_id, data, fit, elapsed_sec) {
  tr <- truth_for_scenario(scen)
  est <- fit$est

  acc <- attr(data, "acceptance_rate")
  if (!is.finite(acc)) acc <- NA_real_

  finite_fit <- is.finite(fit$logLik) & is_valid_nll_value(-fit$logLik)
  regular_fit <- isTRUE(fit$regular_convergence) & finite_fit

  row <- list(
    scenario = scen$name,
    experiment = scen$experiment,
    family = scen$family,
    rep = rep_id,
    n = nrow(data),
    selection_level = scen$selection_level,
    n_selected = sum(data$U == 1L),
    selection_rate = mean(data$U == 1L),
    hidden_acceptance_rate = acc,
    PE_theoretical = tr$PE,
    acceptance_gap = acc - tr$PE,
    logLik = fit$logLik,
    AIC = fit$aic,
    BIC = fit$bic,
    k_params = fit$k,
    finite_fit = finite_fit,
    regular_convergence = regular_fit,
    convergence = fit$convergence,
    start_id = fit$start_id,
    seconds = elapsed_sec,
    sigmaS_true_raw = scen$sigmaS,
    sigmaY_true = tr$sigmaY,
    rho_true = tr$rho,
    etaS_raw_true = tr$etaS_raw,
    etaY_raw_true = tr$etaY_raw,
    etaS_true = tr$etaS,
    etaY_true = tr$etaY,
    tau_bar_true = tr$tau_bar,
    alpha_legacy_true = tr$alpha,
    nu_true = tr$nu,
    sigmaY_hat = est$sigmaY,
    rho_hat = est$rho,
    etaS_hat = est$etaS,
    etaY_hat = est$etaY,
    tau_bar_hat = est$tau_bar,
    alpha_legacy_hat = est$alpha,
    nu_hat = est$nu
  )

  for (j in seq_along(tr$betaSel)) {
    row[[paste0("betaSel_", j, "_true")]] <- tr$betaSel[j]
    row[[paste0("betaSel_", j, "_hat")]]  <- est$betaSel[j]
    row[[paste0("betaSel_", j, "_bias")]] <- est$betaSel[j] - tr$betaSel[j]
  }

  for (j in seq_along(tr$betaOut)) {
    row[[paste0("betaOut_", j, "_true")]] <- tr$betaOut[j]
    row[[paste0("betaOut_", j, "_hat")]]  <- est$betaOut[j]
    row[[paste0("betaOut_", j, "_bias")]] <- est$betaOut[j] - tr$betaOut[j]
  }

  row$sigmaY_bias <- est$sigmaY - tr$sigmaY
  row$rho_bias <- est$rho - tr$rho
  row$etaS_bias <- est$etaS - tr$etaS
  row$etaY_bias <- est$etaY - tr$etaY
  row$tau_bar_bias <- est$tau_bar - tr$tau_bar
  row$alpha_legacy_bias <- row$tau_bar_bias
  row$nu_bias <- ifelse(is.na(tr$nu), NA_real_, est$nu - tr$nu)
  row$error <- ifelse(regular_fit, "ok", "nonregular_or_failed")
  row$message <- fit$message
  row$start_values <- fit$start_values

  as.data.frame(row, check.names = FALSE)
}

error_row <- function(scen, rep_id, elapsed_sec, message) {
  tr <- truth_for_scenario(scen)
  data.frame(
    scenario = scen$name,
    experiment = scen$experiment,
    family = scen$family,
    rep = rep_id,
    n = scen$n,
    selection_level = scen$selection_level,
    finite_fit = FALSE,
    regular_convergence = FALSE,
    convergence = 99L,
    error = "fit_failed",
    message = message,
    seconds = elapsed_sec,
    PE_theoretical = tr$PE,
    tau_bar_true = tr$tau_bar,
    alpha_legacy_true = tr$tau_bar,
    stringsAsFactors = FALSE
  )
}

add_missing_summary_columns <- function(df) {
  needed <- c(
    "selection_rate", "hidden_acceptance_rate", "PE_theoretical", "acceptance_gap",
    "logLik", "AIC", "BIC", "k_params", "seconds", "convergence",
    "finite_fit", "regular_convergence", "message",
    "sigmaY_true", "sigmaY_hat", "sigmaY_bias",
    "rho_true", "rho_hat", "rho_bias",
    "etaS_true", "etaS_hat", "etaS_bias",
    "etaY_true", "etaY_hat", "etaY_bias",
    "etaS_raw_true", "etaY_raw_true",
    "tau_bar_true", "tau_bar_hat", "tau_bar_bias",
    "alpha_legacy_true", "alpha_legacy_hat", "alpha_legacy_bias",
    "nu_true", "nu_hat", "nu_bias",
    paste0("betaSel_", 1:3, "_true"),
    paste0("betaSel_", 1:3, "_hat"),
    paste0("betaSel_", 1:3, "_bias"),
    paste0("betaOut_", 1:3, "_true"),
    paste0("betaOut_", 1:3, "_hat"),
    paste0("betaOut_", 1:3, "_bias")
  )

  for (cc in needed) {
    if (!cc %in% names(df)) df[[cc]] <- NA_real_
  }

  df
}

ensure_convergence_columns <- function(df) {
  df <- add_missing_summary_columns(df)
  if (!"message" %in% names(df)) df$message <- ""
  df %>%
    dplyr::mutate(
      convergence = suppressWarnings(as.integer(convergence)),
      finite_fit = ifelse(
        is.na(finite_fit),
        is.finite(logLik) & is_valid_nll_value(-logLik),
        as.logical(finite_fit) & is.finite(logLik) & is_valid_nll_value(-logLik)
      ),
      regular_convergence = regular_convergence_flag_vec(convergence, message) & finite_fit,
      ok = regular_convergence
    )
}

summarize_by_scenario <- function(reps) {
  reps <- ensure_convergence_columns(reps)

  reps %>%
    dplyr::group_by(experiment, scenario, family, n, selection_level) %>%
    dplyr::summarise(
      reps = dplyr::n(),
      finite_fits = sum(finite_fit, na.rm = TRUE),
      converged = sum(ok, na.rm = TRUE),
      finite_fit_rate = safe_mean(as.numeric(finite_fit)),
      conv_rate = safe_mean(as.numeric(ok)),
      valid_for_report = conv_rate >= CONV_THRESHOLD_FOR_REPORT,

      mean_selection_rate = safe_mean(selection_rate),
      mean_hidden_acceptance = safe_mean(hidden_acceptance_rate),
      mean_PE_theoretical = safe_mean(PE_theoretical),
      mean_acceptance_gap = safe_mean(acceptance_gap),

      mean_logLik = safe_mean(ifelse(ok, logLik, NA_real_)),
      mean_AIC = safe_mean(ifelse(ok, AIC, NA_real_)),
      mean_BIC = safe_mean(ifelse(ok, BIC, NA_real_)),
      time_sec_avg = safe_mean(seconds),
      time_sec_sd = safe_sd(seconds),
      time_min_avg = safe_mean(seconds) / 60,

      rho_true = safe_mean(rho_true),
      etaS_raw_true = safe_mean(etaS_raw_true),
      etaY_raw_true = safe_mean(etaY_raw_true),
      etaS_true = safe_mean(etaS_true),
      etaY_true = safe_mean(etaY_true),
      tau_bar_true = safe_mean(tau_bar_true),
      alpha_legacy_true = safe_mean(alpha_legacy_true),
      nu_true = safe_mean(nu_true),

      mean_sigmaY_hat = safe_mean(ifelse(ok, sigmaY_hat, NA_real_)),
      mean_rho_hat = safe_mean(ifelse(ok, rho_hat, NA_real_)),
      mean_etaS_hat = safe_mean(ifelse(ok, etaS_hat, NA_real_)),
      mean_etaY_hat = safe_mean(ifelse(ok, etaY_hat, NA_real_)),
      mean_tau_bar_hat = safe_mean(ifelse(ok, tau_bar_hat, NA_real_)),
      mean_alpha_legacy_hat = safe_mean(ifelse(ok, alpha_legacy_hat, NA_real_)),
      mean_nu_hat = safe_mean(ifelse(ok, nu_hat, NA_real_)),

      bias_betaSel_1 = safe_mean(ifelse(ok, betaSel_1_bias, NA_real_)),
      rmse_betaSel_1 = safe_rmse(ifelse(ok, betaSel_1_bias, NA_real_)),
      bias_betaSel_2 = safe_mean(ifelse(ok, betaSel_2_bias, NA_real_)),
      rmse_betaSel_2 = safe_rmse(ifelse(ok, betaSel_2_bias, NA_real_)),
      bias_betaSel_3 = safe_mean(ifelse(ok, betaSel_3_bias, NA_real_)),
      rmse_betaSel_3 = safe_rmse(ifelse(ok, betaSel_3_bias, NA_real_)),

      bias_betaOut_1 = safe_mean(ifelse(ok, betaOut_1_bias, NA_real_)),
      rmse_betaOut_1 = safe_rmse(ifelse(ok, betaOut_1_bias, NA_real_)),
      bias_betaOut_2 = safe_mean(ifelse(ok, betaOut_2_bias, NA_real_)),
      rmse_betaOut_2 = safe_rmse(ifelse(ok, betaOut_2_bias, NA_real_)),
      bias_betaOut_3 = safe_mean(ifelse(ok, betaOut_3_bias, NA_real_)),
      rmse_betaOut_3 = safe_rmse(ifelse(ok, betaOut_3_bias, NA_real_)),

      bias_sigmaY = safe_mean(ifelse(ok, sigmaY_bias, NA_real_)),
      rmse_sigmaY = safe_rmse(ifelse(ok, sigmaY_bias, NA_real_)),
      bias_rho = safe_mean(ifelse(ok, rho_bias, NA_real_)),
      rmse_rho = safe_rmse(ifelse(ok, rho_bias, NA_real_)),
      bias_etaS = safe_mean(ifelse(ok, etaS_bias, NA_real_)),
      rmse_etaS = safe_rmse(ifelse(ok, etaS_bias, NA_real_)),
      bias_etaY = safe_mean(ifelse(ok, etaY_bias, NA_real_)),
      rmse_etaY = safe_rmse(ifelse(ok, etaY_bias, NA_real_)),
      bias_tau_bar = safe_mean(ifelse(ok, tau_bar_bias, NA_real_)),
      rmse_tau_bar = safe_rmse(ifelse(ok, tau_bar_bias, NA_real_)),
      bias_alpha_legacy = safe_mean(ifelse(ok, alpha_legacy_bias, NA_real_)),
      rmse_alpha_legacy = safe_rmse(ifelse(ok, alpha_legacy_bias, NA_real_)),
      bias_nu = safe_mean(ifelse(ok, nu_bias, NA_real_)),
      rmse_nu = safe_rmse(ifelse(ok, nu_bias, NA_real_)),
      .groups = "drop"
    )
}

summarize_by_family <- function(summary) {
  summary %>%
    dplyr::group_by(experiment, family) %>%
    dplyr::summarise(
      scenarios = dplyr::n(),
      mean_conv_rate = safe_mean(conv_rate),
      mean_finite_fit_rate = safe_mean(finite_fit_rate),
      mean_selection_rate = safe_mean(mean_selection_rate),
      mean_hidden_acceptance = safe_mean(mean_hidden_acceptance),
      mean_acceptance_gap = safe_mean(mean_acceptance_gap),
      mean_time_min = safe_mean(time_min_avg),
      mean_rmse_betaOut_1 = safe_mean(rmse_betaOut_1),
      mean_rmse_betaOut_2 = safe_mean(rmse_betaOut_2),
      mean_rmse_betaOut_3 = safe_mean(rmse_betaOut_3),
      mean_rmse_sigmaY = safe_mean(rmse_sigmaY),
      mean_rmse_rho = safe_mean(rmse_rho),
      mean_rmse_etaS = safe_mean(rmse_etaS),
      mean_rmse_etaY = safe_mean(rmse_etaY),
      mean_rmse_tau_bar = safe_mean(rmse_tau_bar),
      .groups = "drop"
    )
}

global_bias_rmse <- function(reps) {
  reps <- ensure_convergence_columns(reps)

  params <- c(
    "betaSel_1", "betaSel_2", "betaSel_3",
    "betaOut_1", "betaOut_2", "betaOut_3",
    "sigmaY", "rho", "etaS", "etaY", "tau_bar", "nu"
  )

  purrr::map_dfr(params, function(p) {
    bias_col <- paste0(p, "_bias")
    if (!bias_col %in% names(reps)) return(NULL)

    reps %>%
      dplyr::group_by(experiment, family) %>%
      dplyr::summarise(
        parameter = p,
        mean_bias = safe_mean(ifelse(ok, .data[[bias_col]], NA_real_)),
        RMSE = safe_rmse(ifelse(ok, .data[[bias_col]], NA_real_)),
        n_converged = sum(ok, na.rm = TRUE),
        .groups = "drop"
      )
  }) %>%
    dplyr::arrange(experiment, family, parameter)
}

ic_gap_check_recovery <- function(reps) {
  reps <- add_missing_summary_columns(reps)
  reps %>%
    dplyr::filter(is.finite(AIC), is.finite(BIC), is.finite(k_params), is.finite(n)) %>%
    dplyr::mutate(
      observed_gap = BIC - AIC,
      expected_gap = k_params * (log(n) - 2),
      gap_error = observed_gap - expected_gap
    ) %>%
    dplyr::group_by(experiment, family) %>%
    dplyr::summarise(
      rows = dplyr::n(),
      mean_k = safe_mean(k_params),
      mean_observed_gap = safe_mean(observed_gap),
      mean_expected_gap = safe_mean(expected_gap),
      max_abs_gap_error = max(abs(gap_error), na.rm = TRUE),
      .groups = "drop"
    )
}

convergence_summary <- function(reps) {
  reps <- ensure_convergence_columns(reps)
  reps %>%
    dplyr::group_by(experiment, family) %>%
    dplyr::summarise(
      reps = dplyr::n(),
      finite_fit_rate = safe_mean(as.numeric(finite_fit)),
      conv_rate = safe_mean(as.numeric(ok)),
      abnormal_message_rate = safe_mean(as.numeric(stringr::str_detect(
        dplyr::coalesce(as.character(message), ""),
        stringr::regex("ERROR|ABNORMAL|LNSRCH|LINESEARCH|FAIL|FAILED", ignore_case = TRUE)
      ))),
      mean_seconds = safe_mean(seconds),
      sd_seconds = safe_sd(seconds),
      .groups = "drop"
    )
}

boundary_check <- function(reps) {
  reps <- ensure_convergence_columns(reps)

  reps %>%
    dplyr::filter(ok) %>%
    dplyr::mutate(
      etaS_near_boundary = abs(etaS_hat) > 0.95 * L_ETA,
      etaY_near_boundary = abs(etaY_hat) > 0.95 * L_ETA,
      tau_near_boundary  = abs(tau_bar_hat) > 0.95 * L_TAU
    ) %>%
    dplyr::group_by(experiment, family) %>%
    dplyr::summarise(
      rows = dplyr::n(),
      share_etaS_near_boundary = safe_mean(as.numeric(etaS_near_boundary)),
      share_etaY_near_boundary = safe_mean(as.numeric(etaY_near_boundary)),
      share_tau_near_boundary  = safe_mean(as.numeric(tau_near_boundary)),
      .groups = "drop"
    )
}

## -----------------------------------------------------------------------------
## 9. Running scenarios
## -----------------------------------------------------------------------------

run_one_rep_body <- function(scen, rep_id, seed) {
  t0 <- proc.time()[3]
  data <- generate_data(scen, seed = seed)
  fit <- fit_family(scen, data)
  metrics_row(scen, rep_id, data, fit, elapsed_sec = proc.time()[3] - t0)
}

run_one_rep <- function(scen, rep_id, seed) {
  t0 <- proc.time()[3]

  tryCatch({
    with_timeout(
      run_one_rep_body(scen, rep_id, seed),
      timeout = REPL_TIMEOUT,
      label = paste0(scen$name, "_rep_", rep_id)
    )
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("timed out|timeout", msg, ignore.case = TRUE)) {
      msg <- paste0("TIMEOUT: ", msg)
    } else {
      msg <- paste0("ERROR: ", msg)
    }
    error_row(scen, rep_id, elapsed_sec = proc.time()[3] - t0, message = msg)
  })
}

run_scenario <- function(scen, n_workers) {
  n_workers <- min(n_workers, scen$reps)

  scen_dir <- file.path(RESULTS_DIR, scen$experiment, scen$name)
  rep_dir <- file.path(scen_dir, "per_replication")

  if (!RESUME_EXISTING && dir.exists(rep_dir)) {
    fs::dir_delete(rep_dir)
  }

  fs::dir_create(rep_dir, recurse = TRUE)

  message("--------------------------------------------------------------------")
  message("scenario: ", scen$name)
  message("experiment: ", scen$experiment)
  message("family: ", scen$family, " | n=", scen$n, " | reps=", scen$reps)
  message("rho=", scen$rho,
          " | etaS_raw=", scen$etaS_raw,
          " | etaY_raw=", scen$etaY_raw,
          " | etaS_ident=", signif(scen$etaS_identified, 4),
          " | etaY_ident=", signif(scen$etaY_identified, 4),
          " | tau_bar=", scen$tau_bar,
          " | nu=", ifelse(is.na(scen$nu), "NA", scen$nu))
  message("workers: ", n_workers)

  rep_ids <- seq_len(scen$reps)
  rep_seeds <- make_rep_seeds(scen$reps, scenario_seed(scen$name))

  worker_fun <- function(rep_id) {
    out_path <- file.path(rep_dir, sprintf("rep_%05d.csv", rep_id))
    if (isTRUE(RESUME_EXISTING) && file.exists(out_path)) return(TRUE)

    row <- run_one_rep(scen, rep_id, rep_seeds[rep_id])
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

    parallel::clusterExport(
      cl,
      varlist = c(
        "scen", "rep_dir", "rep_seeds", "RESUME_EXISTING", "REPL_TIMEOUT", "HAS_RUTILS",
        "CDF_MAXPTS", "CDF_ABSEPS", "CDF_RELEPS", "OPT_MAXIT", "COMPUTE_HESSIAN",
        "L_ETA", "L_TAU", "L_ALPHA", "START_STRATEGY", "CONV_THRESHOLD_FOR_REPORT",
        "NLL_PENALTY", "NLL_VALID_MAX",
        "clamp01", "safe_log", "is_valid_nll_value", "safe_mean", "safe_sd", "safe_rmse", "write_csv_atomic", "with_timeout",
        "regular_convergence_flag_vec",
        "Omega_mat", "lambda_eff_from_eta", "eta_from_lambda", "normalize_eta_to_probit_scale",
        "df_mvt", "make_cdf_algorithms_2d", "pmvnorm2_safe", "pmvt2_safe", "PE_gauss", "PE_t",
        "pS_le0_ESN", "joint_y_spos_ESN", "pS_le0_EST", "joint_y_spos_EST",
        "generate_data_ESN", "generate_data_EST", "generate_data",
        "model_matrices", "initial_values", "raw_from_bounded", "make_start_candidates",
        "negloglik_ESN", "negloglik_EST", "parameter_names", "bounds_for_family",
        "decode_estimates", "npar_family", "empty_fit_result", "fit_family",
        "truth_for_scenario", "metrics_row", "error_row",
        "run_one_rep_body", "run_one_rep", "worker_fun"
      ),
      envir = environment()
    )

    invisible(parallel::parLapplyLB(cl, rep_ids, worker_fun))
  }

  csvs <- fs::dir_ls(rep_dir, regexp = "rep_[0-9]+\\.csv$", type = "file")
  reps <- purrr::map_dfr(csvs, ~ suppressMessages(readr::read_csv(.x, show_col_types = FALSE)))
  write_csv_atomic(reps, file.path(scen_dir, "replicates_all.csv"))

  summary <- summarize_by_scenario(reps)
  write_csv_atomic(summary, file.path(scen_dir, "summary.csv"))

  errors <- ensure_convergence_columns(reps) %>%
    dplyr::filter(!ok)
  if (nrow(errors)) write_csv_atomic(errors, file.path(scen_dir, "diagnostic_nonregular.csv"))

  message("completed: ", scen$name)
  invisible(reps)
}

## -----------------------------------------------------------------------------
## 10. Report generation
## -----------------------------------------------------------------------------

read_all_replicates <- function(root = RESULTS_DIR) {
  rep_files <- fs::dir_ls(
    root,
    recurse = TRUE,
    type = "file",
    regexp = "per_replication[/\\\\]rep_[0-9]+\\.csv$"
  )

  if (length(rep_files)) {
    return(purrr::map_dfr(
      rep_files,
      ~ suppressMessages(readr::read_csv(.x, show_col_types = FALSE))
    ))
  }

  files <- fs::dir_ls(root, recurse = TRUE, type = "file", regexp = "replicates_all\\.csv$")
  if (!length(files)) stop("no replication files found under ", root)
  purrr::map_dfr(files, ~ suppressMessages(readr::read_csv(.x, show_col_types = FALSE)))
}

make_tables <- function(reps) {
  reps <- ensure_convergence_columns(reps)
  summary <- summarize_by_scenario(reps)
  family <- summarize_by_family(summary)
  global <- global_bias_rmse(reps)
  ic_gap <- ic_gap_check_recovery(reps)
  conv <- convergence_summary(reps)
  boundary <- boundary_check(reps)

  write_csv_atomic(summary, file.path(DIR_TABLES, "summary_by_scenario.csv"))
  write_csv_atomic(family,  file.path(DIR_TABLES, "family_aggregates.csv"))
  write_csv_atomic(global,  file.path(DIR_TABLES, "global_bias_rmse.csv"))
  write_csv_atomic(ic_gap, file.path(DIR_TABLES, "ic_gap_check_recovery.csv"))
  write_csv_atomic(conv,   file.path(DIR_TABLES, "convergence_summary.csv"))
  write_csv_atomic(boundary, file.path(DIR_TABLES, "boundary_check.csv"))

  write_rds_atomic(summary, file.path(DIR_OBJECT, "summary_by_scenario.rds"))
  write_rds_atomic(family,  file.path(DIR_OBJECT, "family_aggregates.rds"))
  write_rds_atomic(global,  file.path(DIR_OBJECT, "global_bias_rmse.rds"))
  write_rds_atomic(ic_gap, file.path(DIR_OBJECT, "ic_gap_check_recovery.rds"))
  write_rds_atomic(conv,   file.path(DIR_OBJECT, "convergence_summary.rds"))
  write_rds_atomic(boundary, file.path(DIR_OBJECT, "boundary_check.rds"))

  write_latex_table(
    ic_gap,
    file.path(DIR_TABLES, "ic_gap_check_recovery.tex"),
    caption = "Information-criterion diagnostic for the ESN/EST recovery study.",
    label = "tab:mc-ic-gap-check",
    digits = 6
  )

  write_latex_table(
    boundary,
    file.path(DIR_TABLES, "boundary_check.tex"),
    caption = "Boundary diagnostics for the bounded hidden-truncation parameters.",
    label = "tab:mc-boundary-check",
    digits = 4
  )

  grid_summary <- summary %>% dplyr::filter(experiment == "grid_n1500")
  nsweep_summary <- summary %>% dplyr::filter(experiment == "n_sweep")

  if (nrow(grid_summary)) {
    scenario_compact <- grid_summary %>%
      dplyr::select(
        scenario, family, n, rho_true,
        etaS_raw_true, etaY_raw_true, etaS_true, etaY_true,
        tau_bar_true, nu_true, selection_level,
        finite_fit_rate, conv_rate, valid_for_report,
        mean_selection_rate, mean_hidden_acceptance,
        mean_PE_theoretical, mean_acceptance_gap,
        mean_etaS_hat, mean_etaY_hat, mean_tau_bar_hat, mean_nu_hat,
        rmse_betaOut_1, rmse_betaOut_2, rmse_betaOut_3,
        rmse_sigmaY, rmse_rho, rmse_etaS, rmse_etaY, rmse_tau_bar, rmse_nu,
        time_min_avg
      ) %>%
      dplyr::arrange(family, rho_true, etaS_raw_true, selection_level)

    write_csv_atomic(scenario_compact, file.path(DIR_TABLES, "grid_scenario_compact.csv"))
    write_latex_table(
      scenario_compact,
      file.path(DIR_TABLES, "grid_scenario_compact.tex"),
      caption = "Scenario-level Monte Carlo summaries for the ESN/EST recovery grid. Shape parameters and tau_bar are reported on the identified probit-normalized scale.",
      label = "tab:mc-grid-scenarios",
      digits = 4
    )
  }

  if (nrow(nsweep_summary)) {
    nsweep_compact <- nsweep_summary %>%
      dplyr::select(
        family, n, finite_fit_rate, conv_rate, valid_for_report,
        mean_selection_rate, mean_acceptance_gap,
        mean_etaS_hat, mean_etaY_hat, mean_tau_bar_hat, mean_nu_hat,
        rmse_betaOut_1, rmse_betaOut_2, rmse_betaOut_3,
        rmse_sigmaY, rmse_rho, rmse_etaS, rmse_etaY, rmse_tau_bar, rmse_nu,
        time_min_avg
      ) %>%
      dplyr::arrange(family, n)

    write_csv_atomic(nsweep_compact, file.path(DIR_TABLES, "nsweep_summary.csv"))
    write_latex_table(
      nsweep_compact,
      file.path(DIR_TABLES, "nsweep_summary.tex"),
      caption = "Sample-size sweep for matched ESN and EST scenarios.",
      label = "tab:mc-nsweep",
      digits = 4
    )
  }

  write_latex_table(
    family,
    file.path(DIR_TABLES, "family_aggregates.tex"),
    caption = "Family-level Monte Carlo aggregates for the proposed ESN/EST models.",
    label = "tab:mc-family-aggregates",
    digits = 4
  )

  global_key <- global %>%
    dplyr::filter(parameter %in% c(
      "betaSel_1", "betaSel_2", "betaSel_3",
      "betaOut_1", "betaOut_2", "betaOut_3",
      "sigmaY", "rho", "etaS", "etaY", "tau_bar", "nu"
    ))

  write_latex_table(
    global_key,
    file.path(DIR_TABLES, "global_bias_rmse.tex"),
    caption = "Global bias and RMSE for selected parameters. The threshold parameter is the normalized effective tau_bar.",
    label = "tab:mc-global-bias-rmse",
    digits = 4
  )

  invisible(list(summary = summary, family = family, global = global, ic_gap = ic_gap, convergence = conv, boundary = boundary))
}

make_figures <- function(reps, summary) {
  reps <- ensure_convergence_columns(reps)

  params <- c("betaSel_2", "betaSel_3", "betaOut_2", "betaOut_3")
  long_bias <- purrr::map_dfr(params, function(p) {
    bc <- paste0(p, "_bias")
    if (!bc %in% names(reps)) return(NULL)
    reps %>%
      dplyr::filter(ok) %>%
      dplyr::transmute(experiment, family, parameter = p, bias = .data[[bc]]) %>%
      dplyr::filter(is.finite(bias))
  })

  if (nrow(long_bias)) {
    p_bias <- ggplot2::ggplot(long_bias, ggplot2::aes(x = parameter, y = bias, fill = family)) +
      ggplot2::geom_violin(trim = TRUE, alpha = 0.75) +
      ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
      ggplot2::facet_wrap(~ experiment, scales = "free_y") +
      ggplot2::labs(x = "parameter", y = "bias", fill = "family") +
      theme_mc()
    save_plot_all(p_bias, file.path(DIR_FIGS, "bias_distributions"), width = 8, height = 5.5)
  }

  acc_df <- summary %>%
    dplyr::filter(is.finite(mean_acceptance_gap))

  if (nrow(acc_df)) {
    p_acc <- ggplot2::ggplot(acc_df, ggplot2::aes(x = family, y = mean_acceptance_gap, shape = family)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
      ggplot2::geom_jitter(width = 0.08, height = 0, size = 2.3, alpha = 0.85) +
      ggplot2::facet_wrap(~ experiment, scales = "free_y") +
      ggplot2::labs(x = "family", y = "acceptance gap", shape = "family") +
      theme_mc()
    save_plot_all(p_acc, file.path(DIR_FIGS, "acceptance_gap"), width = 7, height = 5)
  }

  ct_df <- summary %>%
    dplyr::filter(is.finite(conv_rate), is.finite(time_min_avg))

  if (nrow(ct_df)) {
    p_ct <- ggplot2::ggplot(ct_df, ggplot2::aes(x = conv_rate, y = time_min_avg, shape = family)) +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::facet_wrap(~ experiment, scales = "free_y") +
      ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      ggplot2::labs(x = "convergence", y = "time (min)", shape = "family") +
      theme_mc()
    save_plot_all(p_ct, file.path(DIR_FIGS, "convergence_time"), width = 7, height = 5)
  }

  ns <- summary %>%
    dplyr::filter(experiment == "n_sweep") %>%
    dplyr::select(family, n, rmse_betaOut_1, rmse_betaOut_2, rmse_betaOut_3,
                  rmse_sigmaY, rmse_rho, rmse_etaS, rmse_etaY, rmse_tau_bar) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("rmse_"), names_to = "parameter", values_to = "rmse") %>%
    dplyr::filter(is.finite(rmse)) %>%
    dplyr::mutate(
      parameter = dplyr::recode(
        parameter,
        rmse_betaOut_1 = "outcome intercept",
        rmse_betaOut_2 = "outcome slope 1",
        rmse_betaOut_3 = "outcome slope 2",
        rmse_sigmaY = "sigma_y",
        rmse_rho = "rho",
        rmse_etaS = "eta_s",
        rmse_etaY = "eta_y",
        rmse_tau_bar = "tau_bar"
      )
    )

  if (nrow(ns)) {
    p_ns <- ggplot2::ggplot(ns, ggplot2::aes(x = n, y = rmse, linetype = family, shape = family)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2.2) +
      ggplot2::facet_wrap(~ parameter, scales = "free_y") +
      ggplot2::labs(x = "sample size", y = "rmse", linetype = "family", shape = "family") +
      theme_mc()
    save_plot_all(p_ns, file.path(DIR_FIGS, "nsweep_rmse"), width = 9, height = 6)
  }
}

build_report <- function() {
  message(">> reading replication outputs ...")
  reps <- read_all_replicates(RESULTS_DIR)
  reps <- ensure_convergence_columns(reps)
  write_csv_atomic(reps, file.path(REPORT_DIR, "all_replications.csv"))
  write_rds_atomic(reps, file.path(DIR_OBJECT, "all_replications.rds"))

  message(">> building tables ...")
  tabs <- make_tables(reps)

  message(">> building figures ...")
  make_figures(reps, tabs$summary)

  invisible(list(reps = reps, summary = tabs$summary, family = tabs$family, global = tabs$global))
}

## -----------------------------------------------------------------------------
## 11. Metadata
## -----------------------------------------------------------------------------

scenario_definitions_table <- function(scenarios) {
  purrr::map2_dfr(scenarios, seq_along(scenarios), function(s, idx) {
    tibble::tibble(
      scenario_index = idx,
      scenario = s$name,
      experiment = s$experiment,
      family = s$family,
      n = s$n,
      reps = s$reps,
      selection_level = s$selection_level,
      betaSel_raw = paste(s$betaSel, collapse = ";"),
      betaSel_identified = paste(s$betaSel / s$sigmaS, collapse = ";"),
      betaOut = paste(s$betaOut, collapse = ";"),
      sigmaS_raw = s$sigmaS,
      sigmaY = s$sigmaY,
      rho = s$rho,
      etaS_raw = s$etaS_raw,
      etaY_raw = s$etaY_raw,
      etaS_identified = s$etaS_identified,
      etaY_identified = s$etaY_identified,
      tau_bar = s$tau_bar,
      alpha_legacy = s$alpha,
      nu = s$nu,
      L_eta = s$L_eta,
      L_tau = s$L_tau
    )
  })
}

write_metadata <- function(scenarios) {
  fs::dir_create(DIR_META, recurse = TRUE)

  cfg <- c(
    sprintf("date: %s", as.character(Sys.time())),
    sprintf("SCRIPT_VERSION: %s", SCRIPT_VERSION),
    sprintf("RUN_MODE: %s", RUN_MODE),
    sprintf("SCENARIO_INDEX: %s", SCENARIO_INDEX),
    sprintf("RUN_TAG: %s", RUN_TAG),
    sprintf("ARCH_TAG: %s", ARCH_TAG),
    sprintf("ROOT_DIR: %s", ROOT_DIR),
    sprintf("RESULTS_DIR: %s", RESULTS_DIR),
    sprintf("REPORT_DIR: %s", REPORT_DIR),
    sprintf("MASTER_SEED: %s", MASTER_SEED),
    sprintf("MAX_WORKERS: %s", MAX_WORKERS),
    sprintf("SLURM_JOB_ID: %s", Sys.getenv("SLURM_JOB_ID", unset = "")),
    sprintf("SLURM_CPUS_PER_TASK: %s", Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")),
    sprintf("RUN_GRID_STUDY: %s", RUN_GRID_STUDY),
    sprintf("RUN_NSWEEP_STUDY: %s", RUN_NSWEEP_STUDY),
    sprintf("REPS_GRID: %s", REPS_GRID),
    sprintf("REPS_NSWEEP: %s", REPS_NSWEEP),
    sprintf("N_GRID: %s", N_GRID),
    sprintf("N_SWEEP: %s", paste(N_SWEEP, collapse = ",")),
    sprintf("GRID_RHOS: %s", paste(GRID_RHOS, collapse = ",")),
    sprintf("GRID_ETA_MAGS: %s", paste(GRID_ETA_MAGS, collapse = ",")),
    sprintf("GRID_SEL_LEVELS: %s", paste(GRID_SEL_LEVELS, collapse = ",")),
    sprintf("GRID_TAU_BAR: %s", GRID_TAU_BAR),
    sprintf("GRID_EST_NU: %s", GRID_EST_NU),
    sprintf("INCLUDE_OPPOSITE_SIGN_DIRECTION: %s", INCLUDE_OPPOSITE_SIGN_DIRECTION),
    sprintf("SIGMA_S_TRUE: %s", SIGMA_S_TRUE),
    sprintf("SIGMA_Y_TRUE: %s", SIGMA_Y_TRUE),
    sprintf("BETA_SEL_BASE: %s", paste(BETA_SEL_BASE, collapse = ",")),
    sprintf("BETA_OUT_BASE: %s", paste(BETA_OUT_BASE, collapse = ",")),
    sprintf("SEL_INTERCEPTS: %s", paste(names(SEL_INTERCEPTS), SEL_INTERCEPTS, sep = "=", collapse = ",")),
    sprintf("NLL_PENALTY: %s", NLL_PENALTY),
    sprintf("NLL_VALID_MAX: %s", NLL_VALID_MAX),
    sprintf("L_ETA: %s", L_ETA),
    sprintf("L_TAU: %s", L_TAU),
    sprintf("CDF_MAXPTS: %s", CDF_MAXPTS),
    sprintf("CDF_ABSEPS: %s", CDF_ABSEPS),
    sprintf("CDF_RELEPS: %s", CDF_RELEPS),
    sprintf("OPT_MAXIT: %s", OPT_MAXIT),
    sprintf("COMPUTE_HESSIAN: %s", COMPUTE_HESSIAN),
    sprintf("START_STRATEGY: %s", START_STRATEGY),
    sprintf("RESUME_EXISTING: %s", RESUME_EXISTING),
    sprintf("REPL_TIMEOUT: %s", REPL_TIMEOUT),
    sprintf("AUTO_INSTALL: %s", AUTO_INSTALL),
    sprintf("CONV_THRESHOLD_FOR_REPORT: %s", CONV_THRESHOLD_FOR_REPORT),
    "NOTE: selection labels COMMON/MODERATE/RARE denote prevalence, with COMMON having the largest baseline selection probability.",
    "NOTE: legacy LOW/MID/HIGH labels map to selection-rate levels LOW=-7, MID=-5, HIGH=-3 if explicitly requested.",
    "NOTE: alpha_legacy equals tau_bar under the canonical latent normalization used here.",
    "NOTE: eta_true in summaries is the identified probit-normalized canonical shape."
  )

  writeLines(cfg, file.path(DIR_META, "run_configuration.txt"))
  writeLines(capture.output(sessionInfo()), file.path(DIR_META, "session_info.txt"))

  scen_tab <- scenario_definitions_table(scenarios)
  write_csv_atomic(scen_tab, file.path(DIR_META, "scenario_definitions.csv"))
  write_csv_atomic(scen_tab, file.path(DIR_META, "scenario_index.csv"))
  write_rds_atomic(scen_tab, file.path(DIR_OBJECT, "scenario_definitions.rds"))

  invisible(NULL)
}

## -----------------------------------------------------------------------------
## 12. Main
## -----------------------------------------------------------------------------

pipeline_main <- function() {
  logical_cores <- parallel::detectCores(logical = TRUE)
  physical_cores <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
  workers <- min(MAX_WORKERS, ifelse(is.na(logical_cores), 1L, logical_cores))

  cat("====================================================================\n")
  cat("array-based Monte Carlo study for ESN/EST family\n")
  cat("--------------------------------------------------------------------\n")
  cat(sprintf("script       : %s\n", SCRIPT_VERSION))
  cat(sprintf("run mode     : %s\n", RUN_MODE))
  cat(sprintf("scenario idx : %s\n", SCENARIO_INDEX))
  cat(sprintf("logical cores : %s\n", ifelse(is.na(logical_cores), "unknown", logical_cores)))
  cat(sprintf("physical cores: %s\n", ifelse(is.na(physical_cores), "unknown", physical_cores)))
  cat(sprintf("workers used  : %d\n", workers))
  cat(sprintf("rep timeout   : %s sec\n", REPL_TIMEOUT))
  cat(sprintf("root folder   : %s\n", ROOT_DIR))
  cat(sprintf("start strategy: %s\n", START_STRATEGY))
  cat(sprintf("resume        : %s\n", RESUME_EXISTING))
  cat("====================================================================\n")

  scenarios <- define_all_scenarios()
  if (!length(scenarios)) stop("no scenarios were defined")

  ## Avoid repeated simultaneous metadata writes from all array tasks. Task 1
  ## writes metadata in scenario mode; list/report/all modes always write it.
  if (RUN_MODE != "scenario" || SCENARIO_INDEX == 1L) {
    write_metadata(scenarios)
  }

  if (RUN_MODE == "list") {
    scen_tab <- scenario_definitions_table(scenarios)
    print(scen_tab)
    cat(sprintf("\nnumber of scenarios: %d\n", length(scenarios)))
    cat(sprintf("scenario index file: %s\n", normalizePath(file.path(DIR_META, "scenario_index.csv"), mustWork = FALSE)))
    return(invisible(scen_tab))
  }

  if (RUN_MODE == "report") {
    out <- build_report()
    cat("====================================================================\n")
    cat("completed report generation\n")
    cat(sprintf("report : %s\n", normalizePath(REPORT_DIR, mustWork = FALSE)))
    cat("====================================================================\n")
    return(invisible(out))
  }

  if (RUN_MODE == "scenario") {
    if (!is.finite(SCENARIO_INDEX) || SCENARIO_INDEX < 1L || SCENARIO_INDEX > length(scenarios)) {
      stop("SCENARIO_INDEX/SLURM_ARRAY_TASK_ID must be between 1 and ", length(scenarios))
    }
    scen <- scenarios[[SCENARIO_INDEX]]
    run_scenario(scen, workers)
    cat("====================================================================\n")
    cat("completed scenario task\n")
    cat(sprintf("scenario index: %d / %d\n", SCENARIO_INDEX, length(scenarios)))
    cat(sprintf("scenario      : %s\n", scen$name))
    cat(sprintf("results       : %s\n", normalizePath(file.path(RESULTS_DIR, scen$experiment, scen$name), mustWork = FALSE)))
    cat("====================================================================\n")
    return(invisible(scen))
  }

  if (RUN_MODE == "all") {
    for (nm in names(scenarios)) {
      run_scenario(scenarios[[nm]], workers)
    }
    out <- build_report()
    cat("====================================================================\n")
    cat("completed all scenarios and report generation\n")
    cat(sprintf("results: %s\n", normalizePath(RESULTS_DIR, mustWork = FALSE)))
    cat(sprintf("report : %s\n", normalizePath(REPORT_DIR, mustWork = FALSE)))
    cat("====================================================================\n")
    return(invisible(out))
  }
}

if (sys.nframe() == 0) {
  pipeline_main()
}
