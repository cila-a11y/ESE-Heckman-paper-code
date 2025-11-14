#!/usr/bin/env Rscript

## =============================================================================
## mc_pipeline_all.R
## End-to-end Monte Carlo pipeline for ESE (ESN/EST) sample-selection models
## - Detects/prints CPU cores; tries to detect GPU and prints it
## - Runs selected scenarios; each scenario runs its reps in parallel
## - Writes per-rep CSVs atomically + per-scenario replicates_all.csv + summary.csv
## - After all scenarios, builds mc_report/ with article-ready tables & ggplots
##   (plots: ggplot2, no titles; all text in English and lowercase)
## =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------#
# 0) Auto-install and load (robust; per-package; with fallbacks)
# -----------------------------------------------------------------------------#
CORE_PKGS <- c(
  "mvtnorm","dplyr","tidyr","readr","stringr","purrr","tibble","ggplot2",
  "scales","forcats","glue","fs","parallel","knitr"
)
OPT_PKGS  <- c("kableExtra","svglite")  # optional: LaTeX & SVG

.install_one <- function(pkg, critical = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("[setup] installing '%s' ...", pkg))
    ok <- TRUE
    tryCatch(
      install.packages(pkg,
                       repos = "https://cloud.r-project.org",
                       dependencies = TRUE,
                       Ncpus = max(1L, getOption("Ncpus", 1L)),
                       type = if (.Platform$OS.type == "windows") "binary" else getOption("pkgType","both")
      ),
      error = function(e) { ok <<- FALSE; message(sprintf("[setup] install error for '%s': %s", pkg, conditionMessage(e))) }
    )
    if (!requireNamespace(pkg, quietly = TRUE)) {
      ok <- FALSE
      if (critical) stop(sprintf("[setup] required package '%s' is not available after install.", pkg))
      message(sprintf("[setup] optional package '%s' not available; features that depend on it will be skipped.", pkg))
    }
    invisible(ok)
  } else TRUE
}

invisible(lapply(CORE_PKGS, .install_one, critical = TRUE))
invisible(lapply(OPT_PKGS,  .install_one, critical = FALSE))

suppressPackageStartupMessages({
  for (p in CORE_PKGS) library(p, character.only = TRUE)
  if (requireNamespace("kableExtra", quietly = TRUE)) library(kableExtra)
})
svg_ok <- requireNamespace("svglite", quietly = TRUE)

# -----------------------------------------------------------------------------#
# 0b) Hardware info (CPU/GPU) — print banner
# -----------------------------------------------------------------------------#
.hw_cmd <- function(cmd, args = character()) {
  out <- tryCatch(system2(cmd, args = args, stdout = TRUE, stderr = FALSE), error = function(e) NULL)
  if (is.null(out) || length(out) == 0) character(0) else out
}

gpu_info <- function() {
  # Try NVIDIA
  nv <- .hw_cmd("nvidia-smi", c("--query-gpu=name,memory.total,driver_version",
                                "--format=csv,noheader"))
  if (length(nv)) {
    return(list(available=TRUE, vendor="NVIDIA", lines=nv))
  }
  # Try AMD ROCm
  ro <- .hw_cmd("rocm-smi", c("--showproductname"))
  if (length(ro) && any(grepl("GPU", ro, ignore.case=TRUE))) {
    return(list(available=TRUE, vendor="AMD", lines=ro))
  }
  # Fallback: list VGA controllers
  vga <- .hw_cmd("lspci", c("-nnk"))
  vga <- if (length(vga)) grep("VGA compatible controller", vga, value=TRUE, ignore.case=TRUE) else character(0)
  if (length(vga)) {
    return(list(available=FALSE, vendor="unknown", lines=vga))
  }
  list(available=FALSE, vendor=NA_character_, lines=character(0))
}

print_hw_banner <- function(use_cores) {
  cat("--------------------------------------------------------------------\n")
  # CPU info
  lscpu <- .hw_cmd("lscpu")
  cpu_logical  <- tryCatch(parallel::detectCores(logical = TRUE),  error=function(e) NA_integer_)
  cpu_physical <- tryCatch(parallel::detectCores(logical = FALSE), error=function(e) NA_integer_)
  cpu_model <- if (length(lscpu)) {
    ln <- lscpu[grepl("^Model name:", lscpu)]; if (length(ln)) sub("^Model name:\\s*", "", ln[1]) else NA_character_
  } else NA_character_
  cat(sprintf("machine cores (logical): %s\n", ifelse(is.na(cpu_logical), "unknown", cpu_logical)))
  cat(sprintf("machine cores (physical): %s\n", ifelse(is.na(cpu_physical), "unknown", cpu_physical)))
  cat(sprintf("workers per scenario (PSOCK/LB): %d\n", use_cores))
  if (!is.na(cpu_model)) cat(sprintf("cpu model: %s\n", cpu_model))
  
  # GPU info
  gi <- gpu_info()
  if (gi$available) {
    cat(sprintf("gpu detected (%s):\n", gi$vendor))
    for (ln in gi$lines) cat("  - ", ln, "\n", sep = "")
    if (gi$vendor == "NVIDIA" && Sys.info()[["sysname"]] != "Windows") {
      cat("note: mvtnorm/optim are cpu-bound. if you want to try nvblas (experimental, nvidia-only), launch r with:\n")
      cat("      LD_PRELOAD=libnvblas.so NVBLAS_CONFIG_FILE=/path/to/nvblas.conf Rscript mc_pipeline_all.R\n")
    }
  } else {
    cat("gpu: none detected or drivers not available. (this pipeline is cpu-bound.)\n")
    if (length(gi$lines)) for (ln in gi$lines) cat("  hint vga: ", ln, "\n", sep = "")
  }
  cat("--------------------------------------------------------------------\n")
  invisible(NULL)
}

# Avoid BLAS/LAPACK oversubscription inside parallel workers
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# -----------------------------------------------------------------------------#
# 1) Global configuration
# -----------------------------------------------------------------------------#
RESULTS_ROOT   <- "results_grid"
REPORT_ROOT    <- "mc_report"
MASTER_SEED    <- 13579
USE_CORES_MAX  <- 64L
REPS_PER_SCEN  <- 200L
N_PER_DATASET  <- 1500L  # used only in full-grid branch
RUN_FULL_GRID  <- TRUE
TEST_SUBSET_REGEX <- "^esn_r(000|030)_eS(m03|p03)_eY(p03|m03)_a04_sel(LOW|MID)$"

# Sentinel n-sweep (ONLY these two scenarios across N_SWEEP)
RUN_SENTINELS_NSWEEP <- TRUE
N_SWEEP <- c(500L, 1000L, 1500L, 2000L, 2500L, 3000L)
SENTINEL_ESN_ID <- "esn_r000_eSm03_eYp03_a04_selLOW"
SENTINEL_EST_ID <- "est_r000_eSm03_eYp03_a04_selLOW_nu6"  # use _nu10 if you prefer lighter tails

# -----------------------------------------------------------------------------#
# 2) Small numeric helpers and mappings
# -----------------------------------------------------------------------------#
.clamp01 <- function(x, eps = 1e-12) pmin(pmax(x, eps), 1 - eps)

Omega_mat <- function(sigmaS, sigmaY, rho) {
  matrix(c(
    sigmaS^2,              rho * sigmaS * sigmaY,
    rho * sigmaS * sigmaY,        sigmaY^2
  ), nrow = 2, byrow = TRUE)
}

# Ω = R'R (R upper chol). Ω^{-1/2} = R^{-T}. λ_eff = (Ω^{-1/2})^T η = R^{-1} η
lambda_eff_from_eta <- function(Omega, eta) {
  Omega <- (Omega + t(Omega))/2
  R <- chol(Omega)
  drop(backsolve(R, eta))
}

.df_mvt <- function(nu) if (!is.finite(nu)) 100000L else as.integer(max(1L, round(nu)))

# -----------------------------------------------------------------------------#
# 3) ESN / EST observed-probability building blocks
# -----------------------------------------------------------------------------#
# -- ESN (Gaussian core) --
PE_gauss_canonical <- function(alpha, eta) .clamp01(pnorm(alpha / sqrt(1 + sum(eta^2))))

pS_le0_ESE_gauss <- function(muS, sigmaS, muY, sigmaY, rho, etaS, etaY, alpha) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam   <- lambda_eff_from_eta(Omega, c(etaS, etaY))
  lamS <- lam[1]; lamY <- lam[2]
  varW <- 1 + sum(c(etaS, etaY)^2); sdW <- sqrt(varW)
  covSW <- - lamS * sigmaS^2 - lamY * rho * sigmaS * sigmaY
  rSW  <- max(min(covSW / (sigmaS * sdW), 0.999999), -0.999999)
  hS <- (0 - muS)/sigmaS; hW <- alpha/sdW
  PE <- PE_gauss_canonical(alpha, c(etaS, etaY))
  alg <- mvtnorm::GenzBretz(maxpts = 2e5, abseps = 1e-8, releps = 1e-8)
  num <- as.numeric(mvtnorm::pmvnorm(
    upper = c(hS, hW),
    corr  = matrix(c(1, rSW, rSW, 1), 2, 2),
    algorithm = alg
  ))
  .clamp01(num / PE)
}

int_bivar_spos_ESE_gauss <- function(y, muS, muY, sigmaS, sigmaY, rho, etaS, etaY, alpha) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam   <- lambda_eff_from_eta(Omega, c(etaS, etaY))
  lamS <- lam[1]; lamY <- lam[2]
  mean_S_y  <- muS + rho * sigmaS/sigmaY * (y - muY)
  varS_cond <- pmax(sigmaS^2 * (1 - rho^2), 1e-12); sdS <- sqrt(varS_cond)
  varW      <- 1 + sum(c(etaS, etaY)^2)
  covYW     <- - lamS * rho * sigmaS * sigmaY - lamY * sigmaY^2
  varW_cond <- pmax(varW - (covYW^2)/(sigmaY^2), 1e-12); sdW <- sqrt(varW_cond)
  covSW      <- - lamS * sigmaS^2 - lamY * rho * sigmaS * sigmaY
  covSW_cond <- covSW - (rho * sigmaS * sigmaY) * covYW / (sigmaY^2)
  r_cond     <- max(min(covSW_cond / (sdS * sdW), 0.999999), -0.999999)
  mean_W_y <- - (lamY + lamS * rho * sigmaS/sigmaY) * (y - muY)
  h1 <- (0 - mean_S_y)/sdS; h2 <- (alpha - mean_W_y)/sdW
  phiY <- dnorm((y - muY)/sigmaY)/sigmaY
  PE   <- PE_gauss_canonical(alpha, c(etaS, etaY))
  alg <- mvtnorm::GenzBretz(maxpts = 2e5, abseps = 1e-8, releps = 1e-8)
  p_sel_cond <- as.numeric(mvtnorm::pmvnorm(
    lower = c(h1, -Inf),
    upper = c(Inf, h2),
    corr  = matrix(c(1, r_cond, r_cond, 1), 2, 2),
    algorithm = alg
  ))
  joint <- phiY * p_sel_cond / PE
  pmax(joint, .Machine$double.xmin)
}

# -- EST (Student-t core) --
PE_t_canonical <- function(alpha, eta, nu) {
  if (!is.finite(nu)) return(.clamp01(pnorm(alpha / sqrt(1 + sum(eta^2)))))
  .clamp01(stats::pt(alpha / sqrt(1 + sum(eta^2)), df = nu))
}

pS_le0_ESE_t <- function(muS, sigmaS, muY, sigmaY, rho, etaS, etaY, alpha, nu) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam   <- lambda_eff_from_eta(Omega, c(etaS, etaY))
  lamS <- lam[1]; lamY <- lam[2]
  varW <- 1 + sum(c(etaS, etaY)^2); sdW <- sqrt(varW)
  covSW <- - lamS * sigmaS^2 - lamY * rho * sigmaS * sigmaY
  rSW  <- max(min(covSW / (sigmaS * sdW), 0.999999), -0.999999)
  hS <- (0 - muS)/sigmaS; hW <- alpha/sdW
  PE <- PE_t_canonical(alpha, c(etaS, etaY), nu)
  alg <- mvtnorm::GenzBretz(maxpts = 2e5, abseps = 1e-8, releps = 1e-8)
  num <- as.numeric(mvtnorm::pmvt(
    lower = c(-Inf, -Inf), upper = c(hS, hW),
    corr  = matrix(c(1, rSW, rSW, 1), 2, 2),
    df    = .df_mvt(nu), algorithm = alg
  ))
  .clamp01(num / PE)
}

int_bivar_spos_ESE_t <- function(y, muS, muY, sigmaS, sigmaY, rho, etaS, etaY, alpha, nu) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam   <- lambda_eff_from_eta(Omega, c(etaS, etaY))
  lamS <- lam[1]; lamY <- lam[2]
  mean_S_y  <- muS + rho * sigmaS/sigmaY * (y - muY)
  varS_cond <- pmax(sigmaS^2 * (1 - rho^2), 1e-12)
  varW      <- 1 + sum(c(etaS, etaY)^2)
  covYW     <- - lamS * rho * sigmaS * sigmaY - lamY * sigmaY^2
  covSW     <- - lamS * sigmaS^2 - lamY * rho * sigmaS * sigmaY
  varW_cond  <- pmax(varW - (covYW^2)/(sigmaY^2), 1e-12)
  covSW_cond <- covSW - (rho * sigmaS * sigmaY) * covYW / (sigmaY^2)
  qy   <- ((y - muY)/sigmaY)^2
  infl <- if (is.finite(nu)) (nu + qy)/(nu + 1) else 1
  sdS   <- sqrt(varS_cond * infl)
  sdW   <- sqrt(varW_cond * infl)
  r_cond <- max(min(covSW_cond / (sqrt(varS_cond) * sqrt(varW_cond)), 0.999999), -0.999999)
  mean_W_y <- - (lamY + lamS * rho * sigmaS/sigmaY) * (y - muY)
  h1 <- (0 - mean_S_y)/sdS; h2 <- (alpha - mean_W_y)/sdW
  fY <- if (is.finite(nu)) stats::dt((y - muY)/sigmaY, df = nu)/sigmaY else
    dnorm((y - muY)/sigmaY)/sigmaY
  PE <- PE_t_canonical(alpha, c(etaS, etaY), nu)
  alg <- mvtnorm::GenzBretz(maxpts = 2e5, abseps = 1e-8, releps = 1e-8)
  p_sel_cond <- as.numeric(mvtnorm::pmvt(
    lower = c(h1, -Inf), upper = c(Inf, h2),
    corr  = matrix(c(1, r_cond, r_cond, 1), 2, 2),
    df    = .df_mvt(if (is.finite(nu)) nu + 1 else 100000L),
    algorithm = alg
  ))
  joint <- fY * p_sel_cond / PE
  pmax(joint, .Machine$double.xmin)
}

# -----------------------------------------------------------------------------#
# 4) Data generation (ESE ESN/EST) and model matrices
# -----------------------------------------------------------------------------#
generate_data_ESE_ESN <- function(n, betaSel, betaOut, sigmaS, sigmaY, rho,
                                  etaS, etaY, alpha, seed = NULL, max_tries_per_obs = 1e5) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- rnorm(n, 0, 3); X2 <- rbinom(n, 1, 0.5)
  Z1 <- rnorm(n);       Z2 <- rpois(n, 2)
  Xsel <- cbind(1, X1, X2); Xout <- cbind(1, Z1, Z2)
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  eta <- c(etaS, etaY); lam_eff <- lambda_eff_from_eta(Omega, eta)
  S_star <- numeric(n); Y_star <- numeric(n)
  acc <- 0L; total_tries <- 0L
  for (i in seq_len(n)) {
    mu_i <- c(sum(Xsel[i,]*betaSel), sum(Xout[i,]*betaOut))
    tries <- 0L
    repeat {
      x  <- as.numeric(mvtnorm::rmvnorm(1, mean = mu_i, sigma = Omega))
      t0 <- rnorm(1)
      L  <- t0 - sum(lam_eff * (x - mu_i))
      tries <- tries + 1L
      if (L <= alpha) {
        S_star[i] <- x[1]; Y_star[i] <- x[2]
        acc <- acc + 1L; total_tries <- total_tries + tries
        break
      }
      if (tries >= max_tries_per_obs) stop("ESN: max_tries reached; tune ||eta||/alpha/Omega.")
    }
  }
  U <- as.integer(S_star > 0)
  Yobs <- ifelse(U==1L, Y_star, NA_real_)
  attr(U, "acceptance_rate") <- if (total_tries>0) acc/total_tries else NA_real_
  data.frame(U, Yobs, X1, X2, Z1, Z2, S_star = S_star, Y_star = Y_star)
}

generate_data_ESE_EST <- function(n, betaSel, betaOut, sigmaS, sigmaY, rho,
                                  etaS, etaY, alpha, nu, seed = NULL, max_tries_per_obs = 1e5) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- rnorm(n, 0, 3); X2 <- rbinom(n, 1, 0.5)
  Z1 <- rnorm(n);       Z2 <- rpois(n, 2)
  Xsel <- cbind(1, X1, X2); Xout <- cbind(1, Z1, Z2)
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  eta <- c(etaS, etaY); lam_eff <- lambda_eff_from_eta(Omega, eta)
  S_star <- numeric(n); Y_star <- numeric(n)
  acc <- 0L; total_tries <- 0L
  for (i in seq_len(n)) {
    mu_i <- c(sum(Xsel[i,]*betaSel), sum(Xout[i,]*betaOut))
    tries <- 0L
    repeat {
      s_mix <- if (is.finite(nu)) sqrt(nu / stats::rchisq(1, df = nu)) else 1
      z_xy  <- as.numeric(mvtnorm::rmvnorm(1, mean = c(0,0), sigma = Omega))
      z_t   <- rnorm(1)
      x  <- mu_i + s_mix * z_xy
      t0 <-          s_mix * z_t
      L  <- t0 - sum(lam_eff * (x - mu_i))
      tries <- tries + 1L
      if (L <= alpha) {
        S_star[i] <- x[1]; Y_star[i] <- x[2]
        acc <- acc + 1L; total_tries <- total_tries + tries
        break
      }
      if (tries >= max_tries_per_obs) stop("EST: max_tries reached; tune ||eta||/alpha/nu/Omega.")
    }
  }
  U <- as.integer(S_star > 0)
  Yobs <- ifelse(U==1L, Y_star, NA_real_)
  attr(U, "acceptance_rate") <- if (total_tries>0) acc/total_tries else NA_real_
  data.frame(U, Yobs, X1, X2, Z1, Z2, S_star = S_star, Y_star = Y_star)
}

model_matrices_Heckman <- function(formSel, formOut, data, selVar="U", outVar="Yobs") {
  selVars <- all.vars(formSel); outVars <- all.vars(formOut)
  allNeeded <- unique(c(selVars, outVars, selVar, outVar))
  dt <- data[, allNeeded, drop = FALSE]
  trS <- stats::terms(formSel, data=dt); trY <- stats::terms(formOut, data=dt)
  Xsel <- stats::model.matrix(stats::delete.response(trS), data=dt, na.action=stats::na.pass)
  Xout <- stats::model.matrix(stats::delete.response(trY), data=dt, na.action=stats::na.pass)
  U <- dt[[selVar]]; Y <- dt[[outVar]]
  stopifnot(all(U %in% c(0,1)), nrow(Xsel)==nrow(Xout), nrow(Xsel)==length(U))
  list(Xsel=Xsel, Xout=Xout, U=U, Y=Y)
}

# -----------------------------------------------------------------------------#
# 5) Log-likelihoods and fits (ESN/EST)
# -----------------------------------------------------------------------------#
init_ESN <- function(Xsel, Xout, Y) {
  y_obs <- Y[!is.na(Y)]; sdy <- suppressWarnings(stats::sd(y_obs))
  if (!is.finite(sdy) || sdy<=1e-6) sdy <- 1
  c(betaSel=rep(0, ncol(Xsel)),
    betaOut=rep(0, ncol(Xout)),
    logSigmaY=log(sdy), rho_t=atanh(0.1),
    lamS_t=0, lamY_t=0, tau_t=0)
}
init_EST <- function(Xsel, Xout, Y) c(init_ESN(Xsel, Xout, Y), log_nu_minus2 = log(4)) # nu≈6

logLik_HeckmanESN <- function(params, Xsel, Xout, U, Y, nB, nG, L_eta=2, L_alpha=2) {
  tiny <- .Machine$double.xmin
  betaS <- params[1:nB]
  betaY <- params[(nB+1):(nB+nG)]
  logSigY <- params[nB+nG+1]; kappa <- params[nB+nG+2]
  eS_t <- params[nB+nG+3]; eY_t <- params[nB+nG+4]; a_t <- params[nB+nG+5]
  sigmaS <- 1; sigmaY <- exp(logSigY); rho <- tanh(kappa)
  etaS <- L_eta * tanh(eS_t); etaY <- L_eta * tanh(eY_t); alpha <- L_alpha * tanh(a_t)
  ll <- 0
  for (i in seq_along(U)) {
    muS <- sum(Xsel[i,]*betaS); muY <- sum(Xout[i,]*betaY)
    if (U[i]==0L) {
      p0 <- pS_le0_ESE_gauss(muS, sigmaS, muY, sigmaY, rho, etaS, etaY, alpha)
      ll <- ll + log(pmax(p0, tiny))
    } else if (!is.na(Y[i])) {
      joint <- int_bivar_spos_ESE_gauss(Y[i], muS, muY, sigmaS, sigmaY, rho, etaS, etaY, alpha)
      ll <- ll + log(pmax(joint, tiny))
    } else return(1e12)
  }
  -ll
}

logLik_HeckmanEST <- function(params, Xsel, Xout, U, Y, nB, nG, L_eta=2, L_alpha=2) {
  tiny <- .Machine$double.xmin
  betaS <- params[1:nB]
  betaY <- params[(nB+1):(nB+nG)]
  logSigY <- params[nB+nG+1]; kappa <- params[nB+nG+2]
  eS_t <- params[nB+nG+3]; eY_t <- params[nB+nG+4]; a_t <- params[nB+nG+5]
  lognu2 <- params[nB+nG+6]
  sigmaS <- 1; sigmaY <- exp(logSigY); rho <- tanh(kappa)
  etaS <- L_eta * tanh(eS_t); etaY <- L_eta * tanh(eY_t); alpha <- L_alpha * tanh(a_t)
  nu   <- 2 + exp(lognu2)
  ll <- 0
  for (i in seq_along(U)) {
    muS <- sum(Xsel[i,]*betaS); muY <- sum(Xout[i,]*betaY)
    if (U[i]==0L) {
      p0 <- pS_le0_ESE_t(muS, sigmaS, muY, sigmaY, rho, etaS, etaY, alpha, nu)
      ll <- ll + log(pmax(p0, tiny))
    } else if (!is.na(Y[i])) {
      joint <- int_bivar_spos_ESE_t(Y[i], muS, muY, sigmaS, sigmaY, rho, etaS, etaY, alpha, nu)
      ll <- ll + log(pmax(joint, tiny))
    } else return(1e12)
  }
  -ll
}

fit_Heckman_ESN <- function(formSel, formOut, data, start, L_eta=2, L_alpha=2,
                            method="L-BFGS-B", maxit=2000, trace=0) {
  mm <- model_matrices_Heckman(formSel, formOut, data)
  Xs <- mm$Xsel; Xy <- mm$Xout; U <- mm$U; Y <- mm$Y
  nB <- ncol(Xs); nG <- ncol(Xy)
  pnames <- c(paste0("betaSel_", seq_len(nB)),
              paste0("betaOut_", seq_len(nG)),
              "logSigmaY","rho_t","lamS_t","lamY_t","tau_t")
  stopifnot(length(start)==length(pnames)); names(start) <- pnames
  lower <- rep(-Inf, length(pnames)); upper <- rep(Inf, length(pnames))
  names(lower) <- pnames; names(upper) <- pnames
  lower["logSigmaY"] <- log(1e-2); upper["logSigmaY"] <- log(1e2)
  lower["rho_t"]     <- atanh(-0.995); upper["rho_t"] <- atanh(0.995)
  lower[c("lamS_t","lamY_t","tau_t")] <- -3; upper[c("lamS_t","lamY_t","tau_t")] <- 3
  start[!is.finite(start)] <- 0
  start <- pmin(pmax(start, lower), upper)
  nLL <- function(par) {
    val <- logLik_HeckmanESN(par, Xs, Xy, U, Y, nB, nG, L_eta, L_alpha)
    i0 <- nB+nG
    sigmaY <- exp(par[i0+1]); rho <- tanh(par[i0+2])
    pen <- 0
    if (abs(rho)>0.98) pen <- pen + 50*(abs(rho)-0.98)^2 * length(U)
    if (sigmaY>20)     pen <- pen + 1e-3*(sigmaY-20)^2   * length(U)
    if (is.finite(val)) val + pen else 1e12
  }
  fit <- optim(par=start, fn=nLL, method=method,
               lower=lower, upper=upper,
               control=list(maxit=maxit, trace=trace, factr=1e7, pgtol=1e-6),
               hessian=TRUE)
  parf <- fit$par
  loglik_pure <- -logLik_HeckmanESN(parf, Xs, Xy, U, Y, nB, nG, L_eta, L_alpha)
  est <- list(
    betaSel = parf[1:nB],
    betaOut = parf[(nB+1):(nB+nG)],
    sigmaS  = 1,
    sigmaY  = exp(parf["logSigmaY"]),
    rho     = tanh(parf["rho_t"]),
    etaS    = L_eta   * tanh(parf["lamS_t"]),
    etaY    = L_eta   * tanh(parf["lamY_t"]),
    alpha   = L_alpha * tanh(parf["tau_t"])
  )
  list(par=parf, est=est,
       logLik_pen = -fit$value,
       logLik_pure = loglik_pure,
       convergence=fit$convergence, hessian=fit$hessian, message=fit$message)
}

fit_Heckman_EST <- function(formSel, formOut, data, start, L_eta=2, L_alpha=2,
                            method="L-BFGS-B", maxit=2000, trace=0) {
  mm <- model_matrices_Heckman(formSel, formOut, data)
  Xs <- mm$Xsel; Xy <- mm$Xout; U <- mm$U; Y <- mm$Y
  nB <- ncol(Xs); nG <- ncol(Xy)
  pnames <- c(paste0("betaSel_", seq_len(nB)),
              paste0("betaOut_", seq_len(nG)),
              "logSigmaY","rho_t","lamS_t","lamY_t","tau_t","log_nu_minus2")
  stopifnot(length(start)==length(pnames)); names(start) <- pnames
  lower <- rep(-Inf, length(pnames)); upper <- rep(Inf, length(pnames))
  names(lower) <- pnames; names(upper) <- pnames
  lower["logSigmaY"] <- log(1e-2); upper["logSigmaY"] <- log(1e2)
  lower["rho_t"]     <- atanh(-0.995); upper["rho_t"] <- atanh(0.995)
  lower[c("lamS_t","lamY_t","tau_t")] <- -3; upper[c("lamS_t","lamY_t","tau_t")] <- 3
  lower["log_nu_minus2"] <- log(0.05); upper["log_nu_minus2"] <- log(98)
  start[!is.finite(start)] <- 0
  start <- pmin(pmax(start, lower), upper)
  nLL <- function(par) {
    val <- logLik_HeckmanEST(par, Xs, Xy, U, Y, nB, nG, L_eta, L_alpha)
    i0 <- nB+nG
    sigmaY <- exp(par[i0+1]); rho <- tanh(par[i0+2])
    pen <- 0
    if (abs(rho)>0.98) pen <- pen + 50*(abs(rho)-0.98)^2 * length(U)
    if (sigmaY>20)     pen <- pen + 1e-3*(sigmaY-20)^2   * length(U)
    if (is.finite(val)) val + pen else 1e12
  }
  fit <- optim(par=start, fn=nLL, method=method,
               lower=lower, upper=upper,
               control=list(maxit=maxit, trace=trace, factr=1e7, pgtol=1e-6),
               hessian=TRUE)
  parf <- fit$par
  loglik_pure <- -logLik_HeckmanEST(parf, Xs, Xy, U, Y, nB, nG, L_eta, L_alpha)
  est <- list(
    betaSel = parf[1:nB],
    betaOut = parf[(nB+1):(nB+nG)],
    sigmaS  = 1,
    sigmaY  = exp(parf["logSigmaY"]),
    rho     = tanh(parf["rho_t"]),
    etaS    = L_eta   * tanh(parf["lamS_t"]),
    etaY    = L_eta   * tanh(parf["lamY_t"]),
    alpha   = L_alpha * tanh(parf["tau_t"]),
    nu      = 2 + exp(parf["log_nu_minus2"])
  )
  list(par=parf, est=est,
       logLik_pen = -fit$value,
       logLik_pure = loglik_pure,
       convergence=fit$convergence, hessian=fit$hessian, message=fit$message)
}

# -----------------------------------------------------------------------------#
# 6) Scenarios, seeds, metrics, and CSV I/O
# -----------------------------------------------------------------------------#
sigmaS_impl <- function(beta_true, beta_hat, j_fixed = NULL) {
  if (!is.null(j_fixed)) return(beta_true[j_fixed] / beta_hat[j_fixed])
  num <- sum(beta_true * beta_hat); den <- sum(beta_hat^2)
  if (den<=1e-12) return(NA_real_)
  s <- num/den; if (is.finite(s) && s>0) s else NA_real_
}

.make_rep_seeds <- function(n_reps, master_seed) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(as.integer(master_seed))
  as.integer(sample.int(1e9, size = n_reps, replace = FALSE))
}

.sel_intercepts <- c(LOW = -3.0, MID = -5.0, HIGH = -7.0)
.fmt <- function(x) gsub("\\+","p", gsub("-","m", gsub("\\.","", as.character(x))))

.make_one_scenario <- function(family, rho, eta_mag, sign_pair, alpha, nu, sel_tag,
                               n = N_PER_DATASET, reps = REPS_PER_SCEN,
                               betaSel_base = c(-5.0, 0.15, 1.0),
                               betaOut_base = c( 5.0, 1.50, 0.8),
                               sigmaS = 4.0, sigmaY = 3.0, L_eta = 1.5, L_alpha = 1.5) {
  sp <- sign_pair
  etaS <- sp[1]*eta_mag; etaY <- sp[2]*eta_mag
  betaSel <- betaSel_base; betaSel[1] <- .sel_intercepts[[sel_tag]]
  gen <- if (family=="ESN") "ESN" else "EST"
  fit <- if (family=="ESN") "ESN" else "EST"
  nm <- sprintf("%s_r%s_eS%s_eY%s_a%s_sel%s",
                tolower(family),
                .fmt(sprintf("%.2f", rho)),
                .fmt(sprintf("%+.1f", etaS)),
                .fmt(sprintf("%+.1f", etaY)),
                .fmt(sprintf("%.1f", alpha)),
                sel_tag)
  list(
    name = nm, reps = reps, n = n, generator = gen, fit_model = fit,
    betaSel = betaSel, betaOut = betaOut_base,
    sigmaS = sigmaS, sigmaY = sigmaY, rho = rho,
    etaS = etaS, etaY = etaY, alpha = alpha,
    nu = if (family=="EST") nu else NA_real_,
    L_eta = L_eta, L_alpha = L_alpha
  )
}

define_scenarios <- function() {
  rhos <- c(0.0, 0.3, 0.7)
  eta_mags <- c(0.3, 0.6, 1.0)
  sign_pairs <- list(c(-1, +1), c(+1, -1))
  sel_tags <- c("LOW", "MID", "HIGH")
  nus_est <- list(NU4=4, NU6=6, NU10=10, NUINF=Inf)
  alpha <- 0.4
  grid <- list()
  # ESN
  for (rho in rhos)
    for (e in eta_mags)
      for (sp in sign_pairs)
        for (sel in sel_tags) {
          sc <- .make_one_scenario("ESN", rho, e, sp, alpha, nu=NA, sel_tag=sel)
          grid[[sc$name]] <- sc
        }
  # EST
  for (rho in rhos)
    for (e in eta_mags)
      for (sp in sign_pairs)
        for (nu_nm in names(nus_est))
          for (sel in sel_tags) {
            sc <- .make_one_scenario("EST", rho, e, sp, alpha, nu=nus_est[[nu_nm]], sel_tag=sel)
            sc$name <- paste0(sc$name, "_", tolower(nu_nm))
            grid[[sc$name]] <- sc
          }
  grid
}

make_start <- function(scen, df) {
  mm <- model_matrices_Heckman(U ~ X1 + X2, Yobs ~ Z1 + Z2, df)
  st <- if (identical(scen$fit_model, "EST")) init_EST(mm$Xsel, mm$Xout, mm$Y) else
    init_ESN(mm$Xsel, mm$Xout, mm$Y)
  st["lamS_t"] <- atanh((scen$etaS)/scen$L_eta)
  st["lamY_t"] <- atanh((scen$etaY)/scen$L_eta)
  st["tau_t"]  <- atanh((scen$alpha)/scen$L_alpha)
  st["rho_t"]  <- atanh(scen$rho)
  st[is.na(st)] <- 0
  st
}

fit_once <- function(scen, df, start, trace=0) {
  if (identical(scen$fit_model, "EST")) {
    fit <- fit_Heckman_EST(U ~ X1 + X2, Yobs ~ Z1 + Z2, df,
                           start=start, L_eta=scen$L_eta, L_alpha=scen$L_alpha, trace=trace)
  } else {
    fit <- fit_Heckman_ESN(U ~ X1 + X2, Yobs ~ Z1 + Z2, df,
                           start=start, L_eta=scen$L_eta, L_alpha=scen$L_alpha, trace=trace)
  }
  fit
}

compute_metrics_row <- function(scen, gen_info, fit_info) {
  df <- gen_info$data; tr <- gen_info$truth; est <- fit_info$est
  nm_bS <- paste0("betaSel_", seq_along(tr$betaSel))
  nm_bY <- paste0("betaOut_", seq_along(tr$betaOut))
  row <- list(
    scenario    = scen$name,
    generator   = scen$generator,
    fit_family  = scen$fit_model,
    n           = nrow(df),
    L_eta       = scen$L_eta,
    L_alpha     = scen$L_alpha,
    accept_rate = gen_info$accept_rate,
    sel_rate    = gen_info$select_rate,
    PE_theo     = gen_info$PE_theo,
    n_selected  = sum(df$U==1),
    logLik_pure = fit_info$logLik_pure,
    logLik_pen  = fit_info$logLik_pen,
    convergence = fit_info$convergence
  )
  # selection betas (identified)
  for (j in seq_along(tr$betaSel)) {
    true_id <- tr$betaSel[j] / tr$sigmaS
    row[[paste0(nm_bS[j], "_true")]] <- true_id
    row[[paste0(nm_bS[j], "_hat") ]] <- est$betaSel[j]
    row[[paste0(nm_bS[j], "_bias")]] <- est$betaSel[j] - true_id
  }
  # outcome betas
  for (j in seq_along(tr$betaOut)) {
    row[[paste0(nm_bY[j], "_true")]] <- tr$betaOut[j]
    row[[paste0(nm_bY[j], "_hat") ]] <- est$betaOut[j]
    row[[paste0(nm_bY[j], "_bias")]] <- est$betaOut[j] - tr$betaOut[j]
  }
  # implicit selection scale
  row[["sigmaS_true"]]      <- tr$sigmaS
  sig_impl <- sigmaS_impl(tr$betaSel, est$betaSel)
  row[["sigmaS_impl_hat"]]  <- sig_impl
  row[["sigmaS_impl_err"]]  <- if (is.finite(sig_impl)) sig_impl - tr$sigmaS else NA_real_
  # key parameters
  add_par <- function(name, hat, tru) {
    row[[paste0(name, "_true")]] <- tru
    row[[paste0(name, "_hat") ]] <- hat
    row[[paste0(name, "_bias")]] <- hat - tru
  }
  add_par("sigmaY", est$sigmaY, tr$sigmaY)
  add_par("rho",    est$rho,    tr$rho)
  add_par("etaS",   est$etaS,   tr$etaS)
  add_par("etaY",   est$etaY,   tr$etaY)
  add_par("alpha",  est$alpha,  tr$alpha)
  add_par("nu",     if (!is.null(est$nu)) est$nu else NA_real_,
          if (!is.na(tr$nu)) tr$nu else NA_real_)
  row[["error"]] <- ifelse(isTRUE(fit_info$convergence == 0), "ok", "fit_failed")
  as.data.frame(row, check.names = FALSE)
}

write_row_csv_atomic <- function(row_df, out_dir, prefix, rep_id) {
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  tmp  <- file.path(out_dir, sprintf("%s_rep_%05d.tmp.csv", prefix, rep_id))
  path <- file.path(out_dir, sprintf("%s_rep_%05d.csv",     prefix, rep_id))
  utils::write.csv(row_df, tmp, row.names=FALSE)
  file.rename(tmp, path)
  invisible(path)
}

summarize_by_scenario <- function(df_rep) {
  rmse <- function(x) sqrt(mean(x^2, na.rm = TRUE))
  has <- function(nm) nm %in% names(df_rep)
  safe_stat <- function(true_nm, hat_nm, fun) {
    if (has(true_nm) && has(hat_nm)) return(fun(df_rep[[hat_nm]] - df_rep[[true_nm]]))
    NA_real_
  }
  df_rep %>%
    mutate(convergence = suppressWarnings(as.integer(convergence))) %>%
    group_by(scenario, generator, fit_family, n, L_eta, L_alpha) %>%
    summarise(
      reps = n(),
      bias_betaSel_1 = safe_stat("betaSel_1_true","betaSel_1_hat",    mean),
      rmse_betaSel_1 = safe_stat("betaSel_1_true","betaSel_1_hat",    rmse),
      bias_betaSel_2 = safe_stat("betaSel_2_true","betaSel_2_hat",    mean),
      rmse_betaSel_2 = safe_stat("betaSel_2_true","betaSel_2_hat",    rmse),
      bias_betaSel_3 = safe_stat("betaSel_3_true","betaSel_3_hat",    mean),
      rmse_betaSel_3 = safe_stat("betaSel_3_true","betaSel_3_hat",    rmse),
      bias_betaOut_1 = safe_stat("betaOut_1_true","betaOut_1_hat",    mean),
      rmse_betaOut_1 = safe_stat("betaOut_1_true","betaOut_1_hat",    rmse),
      bias_betaOut_2 = safe_stat("betaOut_2_true","betaOut_2_hat",    mean),
      rmse_betaOut_2 = safe_stat("betaOut_2_true","betaOut_2_hat",    rmse),
      bias_betaOut_3 = safe_stat("betaOut_3_true","betaOut_3_hat",    mean),
      rmse_betaOut_3 = safe_stat("betaOut_3_true","betaOut_3_hat",    rmse),
      bias_sigmaY    = safe_stat("sigmaY_true", "sigmaY_hat", mean),
      rmse_sigmaY    = safe_stat("sigmaY_true", "sigmaY_hat", rmse),
      bias_rho       = safe_stat("rho_true",    "rho_hat",    mean),
      rmse_rho       = safe_stat("rho_true",    "rho_hat",    rmse),
      bias_etaS      = safe_stat("etaS_true",   "etaS_hat",   mean),
      rmse_etaS      = safe_stat("etaS_true",   "etaS_hat",   rmse),
      bias_etaY      = safe_stat("etaY_true",   "etaY_hat",   mean),
      rmse_etaY      = safe_stat("etaY_true",   "etaY_hat",   rmse),
      bias_alpha     = safe_stat("alpha_true",  "alpha_hat",  mean),
      rmse_alpha     = safe_stat("alpha_true",  "alpha_hat",  rmse),
      bias_nu        = safe_stat("nu_true",     "nu_hat",     mean),
      rmse_nu        = safe_stat("nu_true",     "nu_hat",     rmse),
      mean_PE_theo   = mean(PE_theo,    na.rm=TRUE),
      mean_acc_emp   = mean(accept_rate,na.rm=TRUE),
      mean_sel_rate  = mean(sel_rate,   na.rm=TRUE),
      conv_rate      = mean(convergence == 0, na.rm=TRUE),
      time_sec_avg   = if ("seconds" %in% names(df_rep)) mean(seconds, na.rm=TRUE) else NA_real_
    ) %>% ungroup()
}

# -----------------------------------------------------------------------------#
# 7) Per-scenario runner (parallel within, sequential across scenarios)
# -----------------------------------------------------------------------------#
run_local_parallel <- function(scen, out_dir, n_cores) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("== running scenario: ", scen$name, " ==")
  message("   using ", n_cores, " workers (psock, load-balanced); blas/lapack threads per worker: 1")
  rep_ids   <- seq_len(scen$reps)
  rep_seeds <- .make_rep_seeds(max(rep_ids), MASTER_SEED)
  run_rep <- function(rep_id) {
    t0 <- proc.time()[3]
    row <- tryCatch({
      set.seed(rep_seeds[rep_id])
      # generate
      if (identical(scen$generator, "ESN")) {
        gen_df <- generate_data_ESE_ESN(scen$n, scen$betaSel, scen$betaOut,
                                        scen$sigmaS, scen$sigmaY, scen$rho,
                                        scen$etaS, scen$etaY, scen$alpha)
        PE_theo <- PE_gauss_canonical(scen$alpha, c(scen$etaS, scen$etaY))
      } else {
        gen_df <- generate_data_ESE_EST(scen$n, scen$betaSel, scen$betaOut,
                                        scen$sigmaS, scen$sigmaY, scen$rho,
                                        scen$etaS, scen$etaY, scen$alpha, scen$nu)
        PE_theo <- PE_t_canonical(scen$alpha, c(scen$etaS, scen$etaY), scen$nu)
      }
      gen_info <- list(
        data = gen_df,
        PE_theo = PE_theo,
        accept_rate = attr(gen_df$U, "acceptance_rate"),
        select_rate = mean(gen_df$U==1),
        truth = list(
          betaSel=scen$betaSel, betaOut=scen$betaOut,
          sigmaS=scen$sigmaS, sigmaY=scen$sigmaY, rho=scen$rho,
          etaS=scen$etaS, etaY=scen$etaY, alpha=scen$alpha, nu=scen$nu
        )
      )
      st  <- make_start(scen, gen_df)
      fit <- fit_once(scen, gen_df, start=st, trace=0)
      rw  <- compute_metrics_row(scen, gen_info, fit)
      rw$rep     <- rep_id
      rw$seconds <- as.numeric(proc.time()[3] - t0)
      rw
    }, error = function(e) {
      data.frame(error="fit_failed", scenario=scen$name, generator=scen$generator,
                 fit_family=scen$fit_model, n=scen$n, L_eta=scen$L_eta, L_alpha=scen$L_alpha,
                 rep=rep_id, seconds=as.numeric(proc.time()[3] - t0), message=as.character(e),
                 check.names=FALSE)
    })
    write_row_csv_atomic(row, out_dir, scen$name, rep_id)
    TRUE
  }
  cl <- parallel::makeCluster(n_cores, type="PSOCK")
  on.exit(parallel::stopCluster(cl), add=TRUE)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(mvtnorm); library(dplyr) })
    Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    TRUE
  })
  parallel::clusterExport(cl, varlist = c(
    "scen","MASTER_SEED",
    "generate_data_ESE_ESN","generate_data_ESE_EST",
    "Omega_mat","lambda_eff_from_eta",".clamp01",".df_mvt",
    "PE_gauss_canonical","pS_le0_ESE_gauss","int_bivar_spos_ESE_gauss",
    "PE_t_canonical","pS_le0_ESE_t","int_bivar_spos_ESE_t",
    "model_matrices_Heckman","logLik_HeckmanESN","logLik_HeckmanEST",
    "init_ESN","init_EST","fit_Heckman_ESN","fit_Heckman_EST",
    "make_start","fit_once","compute_metrics_row",".make_rep_seeds",
    "write_row_csv_atomic","sigmaS_impl"
  ), envir = environment())
  # load-balanced
  parallel::parLapplyLB(cl, X=rep_ids, fun=run_rep)
  # aggregate
  csvs <- list.files(out_dir, pattern = paste0("^", scen$name, "_rep_\\d+\\.csv$"), full.names = TRUE)
  df_rep <- dplyr::bind_rows(lapply(csvs, function(p) suppressMessages(readr::read_csv(p, show_col_types = FALSE))))
  readr::write_csv(df_rep, file.path(out_dir, "replicates_all.csv"))
  df_sum <- summarize_by_scenario(df_rep)
  readr::write_csv(df_sum, file.path(out_dir, "summary.csv"))
  message(".. scenario done: ", normalizePath(out_dir))
}

# -----------------------------------------------------------------------------#
# 8) Post-processing: scan results_root, build mc_report (tables + figures)
# -----------------------------------------------------------------------------#
clean_names <- function(x) {
  x %>% stringr::str_replace_all('^\\ufeff', '') %>% stringr::str_replace_all('^ï\\.\\.', '') %>% stringr::str_trim()
}
read_any_csv <- function(path) {
  out <- try(suppressMessages(readr::read_csv(path, show_col_types = FALSE)), silent = TRUE)
  if (inherits(out, "try-error")) out <- suppressMessages(readr::read_delim(path, delim = ";", show_col_types = FALSE))
  names(out) <- clean_names(names(out)); out
}
read_all_scenarios <- function(root = RESULTS_ROOT) {
  scen_dirs <- fs::dir_ls(root, type="directory", recurse = FALSE)
  if (!length(scen_dirs)) stop("no scenario subfolders under: ", root)
  message(">> scanning scenarios under: ", root)
  dfs <- purrr::map(scen_dirs, function(dir_scen) {
    p_all <- fs::path(dir_scen, "replicates_all.csv")
    if (!fs::file_exists(p_all)) return(NULL)
    df <- suppressMessages(readr::read_csv(p_all, show_col_types = FALSE))
    if (!nrow(df)) return(NULL)
    df$scenario <- fs::path_file(dir_scen)
    df
  })
  keep <- purrr::map_lgl(dfs, ~ !is.null(.x) && nrow(.x) > 0)
  dfs <- dfs[keep]
  if (!length(dfs)) stop("no readable CSVs in: ", root)
  reps <- dplyr::bind_rows(dfs)
  message(">> scenarios found: ", length(dfs), " | total rows: ", nrow(reps))
  reps
}

theme_mc <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top",
                   legend.title = ggplot2::element_text(face = "plain"),
                   strip.text = ggplot2::element_text(face = "plain"))
}
save_all <- function(p, file_base, w=10, h=6) {
  fs::dir_create(fs::path_dir(file_base), recurse = TRUE)
  ggplot2::ggsave(paste0(file_base, ".pdf"), p, width = w, height = h, device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(file_base, ".png"), p, width = w, height = h, dpi = 300)
  if (svg_ok) { svglite::svglite(paste0(file_base, ".svg"), width = w, height = h); print(p); grDevices::dev.off() }
}

make_article_tables <- function(sum_by_scen, outdir = file.path(REPORT_ROOT,"tables")) {
  fs::dir_create(outdir, recurse = TRUE)
  tab_meta <- sum_by_scen %>%
    mutate(family = tolower(fit_family)) %>%
    transmute(scenario, family,
              rho = NA_real_, eta_mag = NA_real_, alpha = NA_real_, nu = NA_real_,
              conv_rate, mean_sel_rate, time_min = time_sec_avg/60)
  readr::write_csv(tab_meta, file.path(outdir, "summary_meta_rates.csv"))
  pct <- function(x) sprintf("%.1f\\%%", 100*x)
  fm1 <- function(x) ifelse(is.na(x),"", sprintf("%.1f", x))
  tab_meta_fmt <- tab_meta %>% mutate(
    conv_rate = pct(conv_rate),
    mean_sel_rate = pct(mean_sel_rate),
    time_min = fm1(time_min)
  )
  has_kable <- requireNamespace("kableExtra", quietly = TRUE)
  if (has_kable) {
    ltx1 <- kable(tab_meta_fmt, "latex", booktabs = TRUE,
                  caption = "scenario meta and rates (conversion/selection/time).",
                  col.names = c("scenario","family","rho","|eta|","alpha","nu","conv.","sel.","time (min)")) %>%
      kableExtra::kable_styling(latex_options = c("striped","hold_position","scale_down"))
  } else {
    ltx1 <- knitr::kable(tab_meta_fmt, format = "latex", booktabs = TRUE,
                         caption = "scenario meta and rates (conversion/selection/time).",
                         col.names = c("scenario","family","rho","|eta|","alpha","nu","conv.","sel.","time (min)"))
  }
  writeLines(ltx1, file.path(outdir, "summary_meta_rates.tex"))
  cols_err <- names(sum_by_scen)[stringr::str_detect(names(sum_by_scen), "^(bias|rmse)_(beta(Sel|Out)_[123])$")]
  if (length(cols_err)) {
    tab_err <- sum_by_scen %>%
      dplyr::select(scenario, fit_family, dplyr::all_of(cols_err)) %>%
      mutate(family = tolower(fit_family)) %>% dplyr::select(-fit_family)
    readr::write_csv(tab_err, file.path(outdir, "summary_errors_main.csv"))
    fm2 <- function(x) ifelse(is.na(x),"", sprintf("%.3f", x))
    tab_err_fmt <- tab_err %>% mutate(across(-c(scenario, family), fm2))
    if (has_kable) {
      ltx2 <- kable(tab_err_fmt, "latex", booktabs = TRUE,
                    caption = "bias/rmse for selected coefficients.",
                    col.names = stringr::str_replace_all(names(tab_err_fmt), "_", "\\\\_")) %>%
        kableExtra::kable_styling(latex_options = c("striped","hold_position","scale_down"))
    } else {
      ltx2 <- knitr::kable(tab_err_fmt, format = "latex", booktabs = TRUE,
                           caption = "bias/rmse for selected coefficients.",
                           col.names = stringr::str_replace_all(names(tab_err_fmt), "_", "\\\\_"))
    }
    writeLines(ltx2, file.path(outdir, "summary_errors_main.tex"))
  }
}

plot_bias_violin <- function(reps, outdir = file.path(REPORT_ROOT,"figs")) {
  params <- c("betaOut_2","betaOut_3","betaSel_2","betaSel_3")
  have <- params[params %in% unique(stringr::str_remove(
    names(reps)[stringr::str_detect(names(reps), "_(true|hat)$")], "_(true|hat)$"))]
  if (!length(have)) return(invisible(NULL))
  reps2 <- if ("rep" %in% names(reps)) reps else dplyr::mutate(reps, rep = dplyr::row_number())
  df_long <- reps2 %>%
    dplyr::select(scenario, rep, fit_family, dplyr::matches(paste0("^(", paste(have, collapse="|"), ")_(true|hat)$"))) %>%
    tidyr::pivot_longer(cols = dplyr::matches(paste0("^(", paste(have, collapse="|"), ")_(true|hat)$")),
                        names_to = c("param","what"), names_pattern = "^(.*)_(true|hat)$",
                        values_to = "val") %>%
    dplyr::mutate(val = suppressWarnings(as.numeric(val))) %>%
    tidyr::drop_na(val) %>%
    tidyr::pivot_wider(names_from = what, values_from = val) %>%
    tidyr::drop_na(true, hat) %>%
    dplyr::mutate(bias = hat - true, family = tolower(fit_family))
  if (!nrow(df_long)) return(invisible(NULL))
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = forcats::fct_inorder(param), y = bias, fill = family)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.85) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95, color = "grey20") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(x = NULL, y = "bias (hat − true)", fill = "family") +
    theme_mc()
  save_all(p, file.path(outdir, "bias_violin"))
}

plot_intercepts <- function(reps, outdir = file.path(REPORT_ROOT,"figs")) {
  need <- c("betaSel_1_true","betaSel_1_hat","betaOut_1_true","betaOut_1_hat")
  if (!all(need %in% names(reps))) return(invisible(NULL))
  reps_num <- reps %>% mutate(across(all_of(need), ~ suppressWarnings(as.numeric(.x))))
  dfI <- reps_num %>%
    transmute(scenario, family = tolower(fit_family),
              bias_sel = betaSel_1_hat - betaSel_1_true,
              bias_out = betaOut_1_hat - betaOut_1_true) %>%
    tidyr::pivot_longer(c(bias_sel, bias_out), names_to = "which", values_to = "bias") %>%
    mutate(which = dplyr::recode(which, "bias_sel" = "selection intercept", "bias_out" = "outcome intercept")) %>%
    tidyr::drop_na(bias)
  if (!nrow(dfI)) return(invisible(NULL))
  p <- ggplot2::ggplot(dfI, ggplot2::aes(x = bias, fill = family)) +
    ggplot2::geom_histogram(bins = 40, alpha = 0.9) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::facet_grid(which ~ ., scales = "free_x") +
    ggplot2::labs(x = "bias", y = "count", fill = "family") +
    theme_mc()
  save_all(p, file.path(outdir, "intercepts_bias_hist"))
}

plot_accept_vs_select <- function(sum_by_scen, outdir = file.path(REPORT_ROOT,"figs")) {
  if (!nrow(sum_by_scen)) return(invisible(NULL))
  df <- sum_by_scen %>%
    transmute(family = tolower(fit_family), mean_acc_emp, mean_sel_rate)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_acc_emp, y = mean_sel_rate, shape = family)) +
    ggplot2::geom_point(alpha = 0.9, size = 2.6) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    ggplot2::labs(x = "acceptance rate (hidden truncation)", y = "selection rate (u=1)", shape = "family") +
    theme_mc()
  save_all(p, file.path(outdir, "accept_vs_select"))
}

plot_time_convergence <- function(sum_by_scen, outdir = file.path(REPORT_ROOT,"figs")) {
  if (!nrow(sum_by_scen)) return(invisible(NULL))
  df <- sum_by_scen %>% mutate(family = tolower(fit_family))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = family, y = time_sec_avg/60, size = conv_rate)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::labs(x = "avg time (min)", y = "avg time (min)", size = "conv. rate") +
    theme_mc()
  save_all(p, file.path(outdir, "time_convergence"))
}

post_process_all <- function(root = RESULTS_ROOT, report_dir = REPORT_ROOT) {
  message(">> reading results from: ", root)
  reps <- read_all_scenarios(root)
  message(">> summarizing by scenario ...")
  sum_by <- summarize_by_scenario(reps)
  fs::dir_create(report_dir, recurse = TRUE)
  fs::dir_create(file.path(report_dir,"tables"), recurse = TRUE)
  fs::dir_create(file.path(report_dir,"figs"),   recurse = TRUE)
  readr::write_csv(sum_by, file.path(report_dir, "tables", "summary_by_scenario.csv"))
  message(">> creating article tables (csv + latex) ...")
  make_article_tables(sum_by, outdir = file.path(report_dir, "tables"))
  message(">> generating figures (ggplot, no titles, english, lowercase) ...")
  plot_bias_violin(reps, outdir = file.path(report_dir, "figs"))
  plot_intercepts(reps,    outdir = file.path(report_dir, "figs"))
  plot_accept_vs_select(sum_by, outdir = file.path(report_dir, "figs"))
  plot_time_convergence(sum_by, outdir = file.path(report_dir, "figs"))
  message(">> done. artefacts in: ", normalizePath(report_dir))
  invisible(list(replicates = reps, summary = sum_by, outdir = normalizePath(report_dir)))
}

# -----------------------------------------------------------------------------#
# 9) Main pipeline — detect cores, run scenarios, post-proc
# -----------------------------------------------------------------------------#
pipeline_main <- function() {
  cores_logical  <- parallel::detectCores(logical = TRUE)
  cores_physical <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
  use_cores      <- min(USE_CORES_MAX, if (is.na(cores_logical)) 1L else cores_logical)
  
  cat("====================================================================\n")
  cat("ESE Monte Carlo pipeline — starting\n")
  print_hw_banner(use_cores)
  cat(sprintf("results root: %s\n", RESULTS_ROOT))
  cat(sprintf("report root : %s\n", REPORT_ROOT))
  cat(sprintf("reps per scenario: %d | n per dataset (full-grid): %d\n", REPS_PER_SCEN, N_PER_DATASET))
  cat("====================================================================\n")
  
  fs::dir_create(RESULTS_ROOT, recurse = TRUE)
  
  scs <- define_scenarios()
  scenario_names <- names(scs)
  
  # --- Sentinel n-sweep branch ---
  if (isTRUE(RUN_SENTINELS_NSWEEP)) {
    keep <- names(scs) %in% c(SENTINEL_ESN_ID, SENTINEL_EST_ID)
    scs  <- scs[keep]
    if (!length(scs)) stop("Sentinel IDs not found. Check SENTINEL_*_ID strings.")
    message(">> running sentinel n-sweep for: ", paste(names(scs), collapse=", "))
    for (nm in names(scs)) {
      for (nval in N_SWEEP) {
        scen <- scs[[nm]]
        scen$reps <- REPS_PER_SCEN
        scen$n    <- as.integer(nval)
        scen$name <- paste0(scen$name, "_n", sprintf("%04d", nval))
        out_dir   <- file.path(RESULTS_ROOT, scen$name)
        message("--------------------------------------------------------------------")
        message("scenario: ", scen$name)
        message("generator/fit: ", scen$generator, "/", scen$fit_model,
                " | rho=", scen$rho, " | etaS=", scen$etaS, " | etaY=", scen$etaY,
                " | alpha=", scen$alpha, " | nu=", ifelse(is.na(scen$nu),"NA",scen$nu),
                " | selection=", names(.sel_intercepts)[which(.sel_intercepts==scen$betaSel[1])],
                " | n=", scen$n)
        run_local_parallel(scen, out_dir, n_cores = use_cores)
      }
    }
    message("--------------------------------------------------------------------")
    post_process_all(root = RESULTS_ROOT, report_dir = REPORT_ROOT)
    cat("====================================================================\n")
    cat("Sentinel n-sweep completed.\n")
    cat("Check per-scenario outputs under: ", normalizePath(RESULTS_ROOT), "\n")
    cat("Combined report (tables+figures) under: ", normalizePath(REPORT_ROOT), "\n")
    cat("====================================================================\n")
    return(invisible(NULL))
  }
  
  # --- Full-grid branch (unused when RUN_SENTINELS_NSWEEP=TRUE) ---
  if (!RUN_FULL_GRID) {
    scenario_names <- scenario_names[grepl(TEST_SUBSET_REGEX, scenario_names)]
    message(">> RUN_FULL_GRID is FALSE — running a small subset of scenarios (", length(scenario_names), ")")
  } else {
    message(">> RUN_FULL_GRID is TRUE — running the full grid (", length(scenario_names), " scenarios)")
  }
  for (nm in scenario_names) {
    scen <- scs[[nm]]
    scen$reps <- REPS_PER_SCEN
    scen$n    <- N_PER_DATASET
    out_dir   <- file.path(RESULTS_ROOT, nm)
    message("--------------------------------------------------------------------")
    message("scenario: ", nm)
    message("generator/fit: ", scen$generator, "/", scen$fit_model,
            " | rho=", scen$rho, " | etaS=", scen$etaS, " | etaY=", scen$etaY,
            " | alpha=", scen$alpha, " | nu=", ifelse(is.na(scen$nu),"NA",scen$nu),
            " | selection=", names(.sel_intercepts)[which(.sel_intercepts==scen$betaSel[1])])
    run_local_parallel(scen, out_dir, n_cores = use_cores)
  }
  message("--------------------------------------------------------------------")
  post_process_all(root = RESULTS_ROOT, report_dir = REPORT_ROOT)
  cat("====================================================================\n")
  cat("Pipeline completed.\n")
  cat("Check per-scenario outputs under: ", normalizePath(RESULTS_ROOT), "\n")
  cat("Combined report (tables+figures) under: ", normalizePath(REPORT_ROOT), "\n")
  cat("====================================================================\n")
}

# -----------------------------------------------------------------------------#
# 10) Auto-execute when opened in RStudio and Run is pressed
# -----------------------------------------------------------------------------#
if (interactive() && length(commandArgs(trailingOnly = TRUE)) == 0) {
  pipeline_main()
}

# If run via Rscript:
if (!interactive()) pipeline_main()

