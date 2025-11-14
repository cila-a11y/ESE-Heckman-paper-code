#!/usr/bin/env Rscript
## =============================================================================
##  Synthetic EST/ESN/HN Heckman (parallel, reproducible, figures)
##  - DGP: extended skew-t (nu = 5) via hidden truncation
##  - Fits: EST (t-core), ESN (Gaussian-core), HN (sampleSelection, opcional)
##  - Outputs: results_synth/rep_level_estimates.csv (+ summary, figs)
## =============================================================================

## ---------------------------- user toggles -----------------------------------
N_REP        <- 64      # number of parallel replications
N_CORES_MAX  <- 64      # cap workers here (<= 64)
BASE_SEED    <- 123     # reproducibility
MAKE_PLOTS   <- TRUE    # make 3 figures from an extra run (seed = BASE_SEED)
DO_HECKMAN_NORMAL <- TRUE

OUT_DIR   <- "results_synth"
FIG_DIR   <- "figures_synth"

## --------------------- threading hygiene (very important) --------------------
Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

## --------------------- dependency helper (auto-install) ----------------------
.need <- function(pkg, quietly_ok = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("[setup] installing '%s' from CRAN ...", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (quietly_ok) suppressPackageStartupMessages(library(pkg, character.only = TRUE)) else library(pkg, character.only = TRUE)
  TRUE
}
.need("mvtnorm")
svg_ok <- TRUE
tryCatch(.need("svglite"), error = function(e) { svg_ok <<- FALSE })
use_future <- TRUE
tryCatch({ .need("future"); .need("future.apply") }, error = function(e) { use_future <<- FALSE })
if (DO_HECKMAN_NORMAL) {
  tryCatch(.need("sampleSelection"), error = function(e) {
    message("[setup] could not load 'sampleSelection'; will skip HN baseline.")
  })
}

## ------------------------ integrator configs ---------------------------------
INT_MAXPTS <- 6e4; INT_ABSEPS <- 1e-6; INT_RELEPS <- 1e-6
ALG_GAUSS  <- mvtnorm::GenzBretz(maxpts = INT_MAXPTS, abseps = INT_ABSEPS, releps = INT_RELEPS)
ALG_T      <- mvtnorm::GenzBretz(maxpts = INT_MAXPTS, abseps = INT_ABSEPS, releps = INT_RELEPS)

## --------------------------- core math helpers -------------------------------
.clamp01 <- function(x, eps=1e-12) pmin(pmax(x, eps), 1 - eps)

Omega_mat <- function(sigmaS, sigmaY, rho) {
  matrix(c( sigmaS^2, rho*sigmaS*sigmaY,
            rho*sigmaS*sigmaY, sigmaY^2 ), 2, 2, byrow=TRUE)
}
lambda_eff_from_eta <- function(Omega, eta) {
  Omega <- (Omega + t(Omega))/2
  R <- chol(Omega)
  drop(backsolve(R, eta)) # R^{-1} eta
}
.df_mvt <- function(nu) { if (!is.finite(nu)) return(100000L); as.integer(max(1L, round(nu))) }

## ---- robust bivariate CDF wrappers (fallback if pmvnorm/pmvt fail) ----------
.safe_biv_pmvnorm <- function(lower, upper, rho, mean = c(0,0)) {
  out <- tryCatch(as.numeric(mvtnorm::pmvnorm(lower=lower, upper=upper,
                                              mean=mean,
                                              corr=matrix(c(1,rho,rho,1),2,2),
                                              algorithm=ALG_GAUSS)),
                  error=function(e) NA_real_)
  if (is.finite(out)) return(out)
  out2 <- tryCatch(as.numeric(mvtnorm::pmvnorm(lower=lower, upper=upper,
                                               mean=mean,
                                               corr=matrix(c(1,rho,rho,1),2,2))),
                   error=function(e) NA_real_)
  if (is.finite(out2)) return(out2)
  if (is.finite(rho) && abs(rho) < 1e-12) {
    a1 <- stats::pnorm(upper[1]) - stats::pnorm(lower[1])
    a2 <- stats::pnorm(upper[2]) - stats::pnorm(lower[2])
    return(max(0,a1)*max(0,a2))
  }
  NA_real_
}
.safe_biv_pmvt <- function(lower, upper, rho, df, mean = c(0,0)) {
  out <- tryCatch(as.numeric(mvtnorm::pmvt(lower=lower, upper=upper,
                                           delta=mean,
                                           corr=matrix(c(1,rho,rho,1),2,2),
                                           df=.df_mvt(df),
                                           algorithm=ALG_T)),
                  error=function(e) NA_real_)
  if (is.finite(out)) return(out)
  # fallback ~ produto das marginais (aprox) se |rho| ~ 0
  if (is.finite(rho) && abs(rho) < 1e-12) {
    a1 <- stats::pt(upper[1], df=df) - stats::pt(lower[1], df=df)
    a2 <- stats::pt(upper[2], df=df) - stats::pt(lower[2], df=df)
    return(max(0,a1)*max(0,a2))
  }
  NA_real_
}

## ------------------------ canonical selection probs --------------------------
PE_gauss_canonical <- function(alpha, eta) .clamp01(stats::pnorm(alpha / sqrt(1 + sum(eta^2))))
PE_t_canonical     <- function(alpha, eta, nu) .clamp01(stats::pt(alpha / sqrt(1 + sum(eta^2)), df = nu))

## ------------------------- observed-data blocks: ESN --------------------------
pS_le0_ESE_gauss <- function(muS, sigmaS, muY, sigmaY, rho, etaS, etaY, alpha) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam <- lambda_eff_from_eta(Omega, c(etaS, etaY)); lamS <- lam[1]; lamY <- lam[2]
  varW <- 1 + (etaS^2 + etaY^2); sdW <- sqrt(max(varW, 1e-12))
  covSW <- - lamS * sigmaS^2 - lamY * rho * sigmaS * sigmaY
  rSW <- covSW / (sigmaS * sdW); if (!is.finite(rSW)) rSW <- 0
  rSW <- max(min(rSW, 0.999999), -0.999999)
  
  hS <- (0 - muS) / sigmaS
  hW <- alpha / sdW
  if (!all(is.finite(c(hS, hW)))) return(NA_real_)
  PE <- PE_gauss_canonical(alpha, c(etaS, etaY))
  if (!is.finite(PE) || PE<=0) return(NA_real_)
  num <- .safe_biv_pmvnorm(lower=c(-Inf,-Inf), upper=c(hS,hW), rho=rSW, mean=c(0,0))
  if (!is.finite(num)) return(NA_real_)
  .clamp01(num / PE)
}
int_bivar_spos_ESE_gauss <- function(y, muS, muY, sigmaS, sigmaY, rho, etaS, etaY, alpha) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam <- lambda_eff_from_eta(Omega, c(etaS, etaY)); lamS <- lam[1]; lamY <- lam[2]
  
  mean_S_y  <- muS + rho * sigmaS/sigmaY * (y - muY)
  varS_cond <- max(sigmaS^2 * (1 - rho^2), 1e-12); sdS <- sqrt(varS_cond)
  
  varW      <- 1 + (etaS^2 + etaY^2)
  covYW     <- - lamS * rho * sigmaS * sigmaY - lamY * sigmaY^2
  varW_cond <- max(varW - (covYW^2)/(sigmaY^2), 1e-12); sdW <- sqrt(varW_cond)
  
  covSW      <- - lamS * sigmaS^2 - lamY * rho * sigmaS * sigmaY
  covSW_cond <- covSW - (rho * sigmaS * sigmaY) * covYW / (sigmaY^2)
  r_cond     <- covSW_cond / (sdS * sdW); if (!is.finite(r_cond)) r_cond <- 0
  r_cond <- max(min(r_cond, 0.999999), -0.999999)
  
  mean_W_y <- - (lamY + lamS * rho * sigmaS/sigmaY) * (y - muY)
  h1 <- (0     - mean_S_y)/sdS
  h2 <- (alpha - mean_W_y)/sdW
  if (!all(is.finite(c(h1,h2)))) return(NA_real_)
  phiY <- stats::dnorm((y - muY)/sigmaY)/sigmaY
  PE   <- PE_gauss_canonical(alpha, c(etaS, etaY))
  if (!is.finite(PE) || PE<=0) return(NA_real_)
  p_sel_cond <- .safe_biv_pmvnorm(lower=c(h1,-Inf), upper=c(Inf,h2), rho=r_cond, mean=c(0,0))
  if (!is.finite(p_sel_cond)) return(NA_real_)
  pmax(phiY * p_sel_cond / PE, .Machine$double.xmin)
}

## ------------------------- observed-data blocks: EST --------------------------
pS_le0_ESE_t <- function(muS, sigmaS, muY, sigmaY, rho, etaS, etaY, alpha, nu) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam <- lambda_eff_from_eta(Omega, c(etaS, etaY)); lamS <- lam[1]; lamY <- lam[2]
  varW <- 1 + (etaS^2 + etaY^2); sdW <- sqrt(max(varW, 1e-12))
  covSW <- - lamS * sigmaS^2 - lamY * rho * sigmaS * sigmaY
  rSW <- covSW / (sigmaS * sdW); if (!is.finite(rSW)) rSW <- 0
  rSW <- max(min(rSW, 0.999999), -0.999999)
  hS <- (0 - muS)/sigmaS
  hW <- alpha/sdW
  if (!all(is.finite(c(hS, hW)))) return(NA_real_)
  PE <- PE_t_canonical(alpha, c(etaS, etaY), nu)
  if (!is.finite(PE) || PE<=0) return(NA_real_)
  num <- .safe_biv_pmvt(lower=c(-Inf,-Inf), upper=c(hS,hW), rho=rSW, df=nu, mean=c(0,0))
  if (!is.finite(num)) return(NA_real_)
  .clamp01(num / PE)
}
## (S,W|Y=y) ~ t_{nu+1} com inflator sqrt((nu+q_y)/(nu+1))
int_bivar_spos_ESE_t <- function(y, muS, muY, sigmaS, sigmaY, rho, etaS, etaY, alpha, nu) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam <- lambda_eff_from_eta(Omega, c(etaS, etaY)); lamS <- lam[1]; lamY <- lam[2]
  
  mean_S_y  <- muS + rho * sigmaS/sigmaY * (y - muY)
  varS_cond <- max(sigmaS^2 * (1 - rho^2), 1e-12)
  
  varW      <- 1 + (etaS^2 + etaY^2)
  covYW     <- - lamS * rho * sigmaS * sigmaY - lamY * sigmaY^2
  covSW     <- - lamS * sigmaS^2 - lamY * rho * sigmaS * sigmaY
  
  varW_cond  <- max(varW - (covYW^2)/(sigmaY^2), 1e-12)
  covSW_cond <- covSW - (rho * sigmaS * sigmaY) * covYW / (sigmaY^2)
  r_cond     <- covSW_cond / (sqrt(varS_cond) * sqrt(varW_cond)); if (!is.finite(r_cond)) r_cond <- 0
  r_cond <- max(min(r_cond, 0.999999), -0.999999)
  
  mean_W_y <- - (lamY + lamS * rho * sigmaS/sigmaY) * (y - muY)
  qy   <- ((y - muY)/sigmaY)^2
  infl <- (nu + qy)/(nu + 1)
  sdS  <- sqrt(varS_cond * infl)
  sdW  <- sqrt(varW_cond  * infl)
  
  h1 <- (0     - mean_S_y)/sdS
  h2 <- (alpha - mean_W_y)/sdW
  if (!all(is.finite(c(h1,h2)))) return(NA_real_)
  
  fY <- stats::dt((y - muY)/sigmaY, df=nu)/sigmaY
  PE <- PE_t_canonical(alpha, c(etaS, etaY), nu)
  if (!is.finite(PE) || PE<=0) return(NA_real_)
  
  p_sel_cond <- .safe_biv_pmvt(lower=c(h1,-Inf), upper=c(Inf,h2), rho=r_cond, df=nu+1, mean=c(0,0))
  if (!is.finite(p_sel_cond)) return(NA_real_)
  pmax(fY * p_sel_cond / PE, .Machine$double.xmin)
}

## ------------------------- model matrices & logliks --------------------------
model_matrices_Heckman <- function(formSel, formOut, data, selVar="U", outVar="Yobs") {
  selVars <- all.vars(formSel); outVars <- all.vars(formOut)
  dt <- data[, unique(c(selVars,outVars,selVar,outVar)), drop=FALSE]
  trS <- stats::terms(formSel, data=dt); trY <- stats::terms(formOut, data=dt)
  Xs <- stats::model.matrix(stats::delete.response(trS), data=dt, na.action=stats::na.pass)
  Xy <- stats::model.matrix(stats::delete.response(trY), data=dt, na.action=stats::na.pass)
  U  <- as.integer(dt[[selVar]])
  Y  <- as.numeric(dt[[outVar]])
  stopifnot(all(U %in% c(0,1)))
  list(Xsel = Xs, Xout = Xy, U = U, Y = Y)
}

logLik_HeckmanESN <- function(params, Xsel, Xout, U, Y, nB, nG, L_eta=1.5, L_alpha=1.5) {
  tiny <- .Machine$double.xmin
  betaS <- params[1:nB]
  betaY <- params[(nB+1):(nB+nG)]
  logSigY <- params[nB+nG+1]; kappa <- params[nB+nG+2]
  eS_t <- params[nB+nG+3]; eY_t <- params[nB+nG+4]; a_t <- params[nB+nG+5]
  sigmaS <- 1
  sigmaY <- exp(logSigY); rho <- tanh(kappa)
  etaS <- L_eta * tanh(eS_t); etaY <- L_eta * tanh(eY_t); alpha <- L_alpha * tanh(a_t)
  
  ll <- 0
  for (i in seq_along(U)) {
    muS <- sum(Xsel[i,] * betaS); muY <- sum(Xout[i,] * betaY)
    if (U[i] == 0L) {
      p0 <- pS_le0_ESE_gauss(muS, sigmaS, muY, sigmaY, rho, etaS, etaY, alpha)
      ll <- ll + log(pmax(p0, tiny))
    } else if (!is.na(Y[i])) {
      joint <- int_bivar_spos_ESE_gauss(Y[i], muS, muY, sigmaS, sigmaY, rho, etaS, etaY, alpha)
      ll <- ll + log(pmax(joint, tiny))
    } else return(1e12)
  }
  -ll
}

logLik_HeckmanEST <- function(params, Xsel, Xout, U, Y, nB, nG, L_eta=1.5, L_alpha=1.5) {
  tiny <- .Machine$double.xmin
  betaS <- params[1:nB]
  betaY <- params[(nB+1):(nB+nG)]
  logSigY <- params[nB+nG+1]; kappa <- params[nB+nG+2]
  eS_t <- params[nB+nG+3]; eY_t <- params[nB+nG+4]; a_t <- params[nB+nG+5]
  lognu2 <- params[nB+nG+6]
  sigmaS <- 1
  sigmaY <- exp(logSigY); rho <- tanh(kappa)
  etaS <- L_eta * tanh(eS_t); etaY <- L_eta * tanh(eY_t); alpha <- L_alpha * tanh(a_t)
  nu   <- 2 + exp(lognu2)
  
  ll <- 0
  for (i in seq_along(U)) {
    muS <- sum(Xsel[i,] * betaS); muY <- sum(Xout[i,] * betaY)
    if (U[i] == 0L) {
      p0 <- pS_le0_ESE_t(muS, sigmaS, muY, sigmaY, rho, etaS, etaY, alpha, nu)
      ll <- ll + log(pmax(p0, tiny))
    } else if (!is.na(Y[i])) {
      joint <- int_bivar_spos_ESE_t(Y[i], muS, muY, sigmaS, sigmaY, rho, etaS, etaY, alpha, nu)
      ll <- ll + log(pmax(joint, tiny))
    } else return(1e12)
  }
  -ll
}

## ------------------------- initializations -----------------------------------
make_start_ESN <- function(formSel, formOut, data) {
  mm <- model_matrices_Heckman(formSel, formOut, data, "U","Yobs")
  Xs <- mm$Xsel; Xy <- mm$Xout; U <- mm$U; Y <- mm$Y
  pS <- ncol(Xs); pY <- ncol(Xy)
  # probit + OLS as robust defaults
  st <- rep(0, pS + pY + 5)
  names(st) <- c(paste0("betaSel_",1:pS), paste0("betaOut_",1:pY),
                 "logSigmaY","rho_t","lamS_t","lamY_t","tau_t")
  prob <- try(stats::glm(U ~ Xs - 1, family=stats::binomial(link="probit")), silent=TRUE)
  if (!inherits(prob,"try-error")) {
    b <- try(as.numeric(stats::coef(prob)), silent=TRUE)
    if (length(b)==pS && all(is.finite(b))) st[1:pS] <- b
  }
  idx_sel <- which(U==1L)
  if (length(idx_sel) > 5) {
    ols <- try(stats::lm(Y ~ Xy - 1, subset=(U==1L)), silent=TRUE)
    if (!inherits(ols,"try-error")) {
      b <- try(as.numeric(stats::coef(ols)), silent=TRUE)
      if (length(b)==pY && all(is.finite(b))) st[pS + 1:pY] <- b
      sdY <- try(as.numeric(summary(ols)$sigma), silent=TRUE)
      if (is.finite(sdY) && sdY>0) st["logSigmaY"] <- log(sdY)
    }
  }
  if (!is.finite(st["logSigmaY"])) st["logSigmaY"] <- 0 # sd=1 fallback
  st["rho_t"] <- atanh(0.1)
  st[c("lamS_t","lamY_t","tau_t")] <- 0
  st
}
make_start_EST <- function(formSel, formOut, data) {
  c(make_start_ESN(formSel, formOut, data), log_nu_minus2 = log(3))  # nu ~ 5
}

## ----------------------------- fit wrappers ----------------------------------
fit_Heckman_ESN <- function(formSel, formOut, data, start,
                            L_eta=1.5, L_alpha=1.5,
                            fix_rho=FALSE, rho_fix_val=0.0, rho_cap=0.95,
                            method="L-BFGS-B", maxit=3000, trace=0) {
  mm <- model_matrices_Heckman(formSel, formOut, data, "U","Yobs")
  Xs <- mm$Xsel; Xy <- mm$Xout; U <- mm$U; Y <- mm$Y
  nB <- ncol(Xs); nG <- ncol(Xy)
  pnames <- c(paste0("betaSel_",1:nB), paste0("betaOut_",1:nG),
              "logSigmaY","rho_t","lamS_t","lamY_t","tau_t")
  if (length(start)!=length(pnames)) stop("Start length mismatch.")
  names(start) <- pnames
  
  lower <- rep(-Inf, length(pnames)); upper <- rep( Inf, length(pnames))
  names(lower) <- pnames; names(upper) <- pnames
  lower["logSigmaY"] <- log(1e-2); upper["logSigmaY"] <- log(1e2)
  lower["rho_t"] <- atanh(-rho_cap); upper["rho_t"] <- atanh(rho_cap)
  lower[c("lamS_t","lamY_t","tau_t")] <- -3; upper[c("lamS_t","lamY_t","tau_t")] <- 3
  
  start[!is.finite(start)] <- 0
  start <- pmin(pmax(start, lower), upper)
  
  nLL <- function(par) {
    par_use <- par
    if (fix_rho) par_use["rho_t"] <- atanh(max(min(rho_fix_val, rho_cap), -rho_cap))
    val <- logLik_HeckmanESN(par_use, Xs, Xy, U, Y, nB, nG, L_eta, L_alpha)
    # mild penalties near numeric corners
    i0 <- nB + nG
    sigmaY <- exp(par_use[i0 + 1]); rho <- tanh(par_use[i0 + 2])
    pen <- 0
    if (abs(rho) > 0.98) pen <- pen + 50 * (abs(rho) - 0.98)^2 * length(U)
    if (sigmaY > 20)     pen <- pen + 1e-3 * (sigmaY - 20)^2 * length(U)
    if (is.finite(val)) val + pen else 1e12
  }
  
  fit <- stats::optim(par = start, fn = nLL, method = method,
                      lower = lower, upper = upper,
                      control = list(maxit = maxit, trace = trace, factr = 1e7, pgtol = 1e-6))
  list(par = fit$par, mm = mm)
}
fit_Heckman_ESN_case2 <- function(formSel, formOut, data, start,
                                  L_eta=1.5, L_alpha=1.5,
                                  rho_fix_val=0.1, rho_cap=0.95,
                                  maxit=3500) {
  out1 <- fit_Heckman_ESN(formSel, formOut, data, start,
                          L_eta, L_alpha, fix_rho=TRUE, rho_fix_val=rho_fix_val, rho_cap=rho_cap,
                          maxit=floor(maxit/2))
  out2 <- fit_Heckman_ESN(formSel, formOut, data, out1$par,
                          L_eta, L_alpha, fix_rho=FALSE, rho_cap=rho_cap,
                          maxit=maxit - floor(maxit/2))
  out2
}

fit_Heckman_EST <- function(formSel, formOut, data, start,
                            L_eta=1.5, L_alpha=1.5,
                            method="L-BFGS-B", maxit=3000, trace=0) {
  mm <- model_matrices_Heckman(formSel, formOut, data, "U","Yobs")
  Xs <- mm$Xsel; Xy <- mm$Xout; U <- mm$U; Y <- mm$Y
  nB <- ncol(Xs); nG <- ncol(Xy)
  pnames <- c(paste0("betaSel_",1:nB), paste0("betaOut_",1:nG),
              "logSigmaY","rho_t","lamS_t","lamY_t","tau_t","log_nu_minus2")
  if (length(start)!=length(pnames)) stop("Start length mismatch.")
  names(start) <- pnames
  
  lower <- rep(-Inf, length(pnames)); upper <- rep( Inf, length(pnames))
  names(lower) <- pnames; names(upper) <- pnames
  lower["logSigmaY"] <- log(1e-2); upper["logSigmaY"] <- log(1e2)
  lower["rho_t"] <- atanh(-0.995); upper["rho_t"] <- atanh(0.995)
  lower[c("lamS_t","lamY_t","tau_t")] <- -3; upper[c("lamS_t","lamY_t","tau_t")] <- 3
  lower["log_nu_minus2"] <- log(0.05); upper["log_nu_minus2"] <- log(98)  # nu in (2.05, ~100)
  
  start[!is.finite(start)] <- 0
  start <- pmin(pmax(start, lower), upper)
  
  nLL <- function(par) {
    val <- logLik_HeckmanEST(par, Xs, Xy, U, Y, nB, nG, L_eta, L_alpha)
    i0 <- nB + nG
    sigmaY <- exp(par[i0 + 1]); rho <- tanh(par[i0 + 2])
    pen <- 0
    if (abs(rho) > 0.98) pen <- pen + 50 * (abs(rho) - 0.98)^2 * length(U)
    if (sigmaY > 20)     pen <- pen + 1e-3 * (sigmaY - 20)^2 * length(U)
    if (is.finite(val)) val + pen else 1e12
  }
  
  fit <- stats::optim(par = start, fn = nLL, method = method,
                      lower = lower, upper = upper,
                      control = list(maxit = maxit, trace = trace, factr = 1e7, pgtol = 1e-6))
  list(par = fit$par, mm = mm)
}

## ---------------------------- DGP (EST) --------------------------------------
## IMPORTANT: geramos (epsS, epsY0) com var=1 e corr=rho;
## depois usamos Y* = gamma + sigmaY_true * epsY0 (evita "dupla" escala).
sim_ESE_t_errors <- function(n, rho, lambda, alpha, nu) {
  Omega0 <- matrix(c(1, rho, rho, 1), 2, 2, byrow=TRUE)  # base-scale
  Sig3   <- rbind(cbind(Omega0, c(0,0)), c(0,0,1))
  out <- matrix(NA_real_, n, 2)
  got <- 0
  while (got < n) {
    m <- max(2*(n-got), 2000)
    Z <- mvtnorm::rmvt(m, sigma = Sig3, df = nu)   # (epsS0, epsY0, T)
    X <- Z[,1:2, drop=FALSE]; Tvar <- Z[,3]
    keep <- Tvar <= alpha + as.numeric(X %*% lambda)
    if (any(keep)) {
      k <- sum(keep); take <- min(k, n-got)
      out[(got+1):(got+take), ] <- X[which(keep)[1:take], , drop=FALSE]
      got <- got + take
    }
  }
  colnames(out) <- c("epsS","epsY0")
  out
}

## --------------------------- one replication ---------------------------------
run_one_rep <- function(seed) {
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
             VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
  set.seed(seed)
  
  ## truth
  beta0 <- -0.7; beta1 <- 0.15; beta2 <- 1.0
  gamma0 <- 5.0; gamma1 <- 1.50; gamma2 <- 0.8
  sigmaY_true <- 3.0; rho_true <- 0.4
  lambda_true <- c(lambdaS = -0.7, lambdaY = 0.7)
  alpha_true  <- 0.5; nu_true <- 5
  n <- 1000
  
  ## covariates
  X1 <- rnorm(n, 0, sqrt(3))
  X2 <- rbinom(n, 1, 0.5)
  Z1 <- rnorm(n, 0, 1)
  Z2 <- rpois(n, 2)
  
  ## errors (EST hidden truncation) base-scale
  E <- sim_ESE_t_errors(n, rho = rho_true, lambda = lambda_true, alpha = alpha_true, nu = nu_true)
  epsS <- E[,1]; epsY0 <- E[,2]
  
  S   <- beta0 + beta1*X1 + beta2*X2 + epsS
  U   <- as.integer(S > 0)
  Yst <- gamma0 + gamma1*Z1 + gamma2*Z2 + sigmaY_true * epsY0
  Yobs<- ifelse(U==1, Yst, NA_real_)
  
  DF <- data.frame(U=U, Yobs=Yobs, X1=X1, X2=X2, Z1=Z1, Z2=Z2)
  sel_formula <- U    ~ X1 + X2
  out_formula <- Yobs ~ Z1 + Z2
  
  ## ESN (two-stage: fix rho then free)
  st_esn <- make_start_ESN(sel_formula, out_formula, DF)
  fit1   <- fit_Heckman_ESN(sel_formula, out_formula, DF, start=st_esn,
                            fix_rho=TRUE, rho_fix_val=0.1, rho_cap=0.95, maxit=2000)
  fit_esn<- fit_Heckman_ESN(sel_formula, out_formula, DF, start=fit1$par,
                            fix_rho=FALSE, rho_cap=0.95, maxit=2000)
  mm  <- fit_esn$mm
  ll_esn <- -logLik_HeckmanESN(fit_esn$par, mm$Xsel, mm$Xout, mm$U, mm$Y,
                               ncol(mm$Xsel), ncol(mm$Xout), 1.5, 1.5)
  pS <- ncol(mm$Xsel); pY <- ncol(mm$Xout); k_esn <- pS + pY + 5
  aic_esn <- -2*ll_esn + 2*k_esn
  bic_esn <- -2*ll_esn + k_esn*log(nrow(DF))
  
  par_esn <- fit_esn$par
  est_esn <- list(
    betaOut = par_esn[(pS+1):(pS+pY)],
    sigmaY  = exp(par_esn["logSigmaY"]),
    rho     = tanh(par_esn["rho_t"]),
    etaS    = 1.5*tanh(par_esn["lamS_t"]),
    etaY    = 1.5*tanh(par_esn["lamY_t"]),
    alpha   = 1.5*tanh(par_esn["tau_t"])
  )
  
  ## EST
  st_est <- make_start_EST(sel_formula, out_formula, DF)
  fit_est<- fit_Heckman_EST(sel_formula, out_formula, DF, start=st_est, maxit=2500)
  ll_est <- -logLik_HeckmanEST(fit_est$par, mm$Xsel, mm$Xout, mm$U, mm$Y,
                               ncol(mm$Xsel), ncol(mm$Xout), 1.5, 1.5)
  k_est  <- pS + pY + 6
  aic_est<- -2*ll_est + 2*k_est
  bic_est<- -2*ll_est + k_est*log(nrow(DF))
  
  par_est <- fit_est$par
  est_est <- list(
    betaOut = par_est[(pS+1):(pS+pY)],
    sigmaY  = exp(par_est["logSigmaY"]),
    rho     = tanh(par_est["rho_t"]),
    etaS    = 1.5*tanh(par_est["lamS_t"]),
    etaY    = 1.5*tanh(par_est["lamY_t"]),
    alpha   = 1.5*tanh(par_est["tau_t"]),
    nu      = 2 + exp(par_est["log_nu_minus2"])
  )
  
  ## HN baseline (optional)
  ll_hn <- NA_real_; aic_hn <- NA_real_; bic_hn <- NA_real_
  rho_hn <- NA_real_; sigmaY_hn <- NA_real_; hn_gamma0 <- NA_real_
  if (DO_HECKMAN_NORMAL && "sampleSelection" %in% .packages(all.available = TRUE)) {
    hN <- try(sampleSelection::selection(selection = sel_formula, outcome = out_formula,
                                         data = DF, method = "ml"),
              silent = TRUE)
    if (!inherits(hN,"try-error")) {
      ll_hn <- as.numeric(logLik(hN))
      k_hn  <- length(coef(hN)) + 1L   # +1 sigma
      aic_hn <- -2*ll_hn + 2*k_hn
      bic_hn <- -2*ll_hn + k_hn*log(nrow(DF))
      nm <- names(coef(hN))
      idx0 <- which(nm %in% c("outcome_(Intercept)","outcome:(Intercept)"))[1]
      if (length(idx0)==1) hn_gamma0 <- as.numeric(coef(hN)[idx0])
      rho_hn <- tryCatch(as.numeric(coef(hN)["rho"]), error=function(e) NA_real_)
      sigmaY_hn <- tryCatch(as.numeric(summary(hN)$estimate["sigma","Estimate"]),
                            error=function(e) tryCatch(as.numeric(hN$sigma), error=function(e2) NA_real_))
    }
  }
  
  data.frame(
    seed = seed,
    n = nrow(DF),
    sel_rate = mean(DF$U),
    
    ## truth
    truth_gamma0 = gamma0, truth_gamma1=gamma1, truth_gamma2=gamma2,
    truth_rho = rho_true, truth_sigmaY = sigmaY_true,
    
    ## ESN
    esn_gamma0 = est_esn$betaOut[1], esn_gamma1 = est_esn$betaOut[2], esn_gamma2 = est_esn$betaOut[3],
    esn_rho = est_esn$rho, esn_sigmaY = est_esn$sigmaY,
    esn_etaS = est_esn$etaS, esn_etaY = est_esn$etaY, esn_alpha = est_esn$alpha,
    logLik_esn = ll_esn, AIC_esn = aic_esn, BIC_esn = bic_esn,
    
    ## HN
    logLik_hn = ll_hn, AIC_hn = aic_hn, BIC_hn = bic_hn,
    hn_rho = rho_hn, hn_sigmaY = sigmaY_hn, hn_gamma0 = hn_gamma0,
    
    ## EST
    est_gamma0 = est_est$betaOut[1], est_gamma1 = est_est$betaOut[2], est_gamma2 = est_est$betaOut[3],
    est_rho = est_est$rho, est_sigmaY = est_est$sigmaY, est_nu = est_est$nu,
    logLik_est = ll_est, AIC_est = aic_est, BIC_est = bic_est,
    
    stringsAsFactors = FALSE
  )
}

## ----------------------- run parallel replications ---------------------------
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

if (use_future) {
  det <- parallel::detectCores(logical = TRUE)
  N_CORES <- min(N_CORES_MAX, ifelse(is.finite(det) && det > 0, det, 1))
  message(sprintf("[parallel] detected cores: %s | using workers: %s", det, N_CORES))
  future::plan(future::multisession, workers = N_CORES)
  seeds <- BASE_SEED + seq_len(N_REP)
  message("[parallel] running replications ...")
  t0 <- Sys.time()
  rep_list <- future.apply::future_lapply(seeds, run_one_rep, future.seed = TRUE)
  t1 <- Sys.time()
  message(sprintf("[parallel] done in %.1f seconds.", as.numeric(difftime(t1, t0, units="secs"))))
} else {
  message("[serial] 'future' not available; running serially.")
  seeds <- BASE_SEED + seq_len(N_REP)
  t0 <- Sys.time()
  rep_list <- lapply(seeds, run_one_rep)
  t1 <- Sys.time()
  message(sprintf("[serial] done in %.1f seconds.", as.numeric(difftime(t1, t0, units="secs"))))
}

rep_tab <- do.call(rbind, rep_list)

## save per-rep results
fn_rep <- file.path(OUT_DIR, "rep_level_estimates.csv")
utils::write.csv(rep_tab, fn_rep, row.names = FALSE)
message(sprintf("[save] per-rep estimates -> %s", tryCatch(normalizePath(fn_rep), error=function(e) fn_rep)))

## summary stats (mean/sd/median)
num_cols <- sapply(rep_tab, is.numeric)
summ <- data.frame(
  metric = c("mean","sd","median"),
  rbind(
    sapply(rep_tab[, num_cols, drop=FALSE], mean,   na.rm=TRUE),
    sapply(rep_tab[, num_cols, drop=FALSE], sd,     na.rm=TRUE),
    sapply(rep_tab[, num_cols, drop=FALSE], median, na.rm=TRUE)
  ),
  check.names = FALSE
)
fn_sum <- file.path(OUT_DIR, "summary_stats.csv")
utils::write.csv(summ, fn_sum, row.names = FALSE)
message(sprintf("[save] summary stats -> %s", tryCatch(normalizePath(fn_sum), error=function(e) fn_sum)))

## ------------------------------ one run for figures --------------------------
if (MAKE_PLOTS) {
  message("[figures] generating figures from a fresh replicate (seed = BASE_SEED) ...")
  set.seed(BASE_SEED)
  ## truth
  beta0 <- -0.7; beta1 <- 0.15; beta2 <- 1.0
  gamma0 <- 5.0; gamma1 <- 1.50; gamma2 <- 0.8
  sigmaY_true <- 3.0; rho_true <- 0.4
  lambda_true <- c(lambdaS = -0.7, lambdaY = 0.7)
  alpha_true  <- 0.5; nu_true <- 5
  n <- 1000
  
  X1 <- rnorm(n, 0, sqrt(3)); X2 <- rbinom(n, 1, 0.5)
  Z1 <- rnorm(n, 0, 1);        Z2 <- rpois(n, 2)
  E  <- sim_ESE_t_errors(n, rho=rho_true, lambda=lambda_true, alpha=alpha_true, nu=nu_true)
  epsS <- E[,1]; epsY0 <- E[,2]
  S   <- beta0 + beta1*X1 + beta2*X2 + epsS
  U   <- as.integer(S > 0)
  Yst <- gamma0 + gamma1*Z1 + gamma2*Z2 + sigmaY_true*epsY0
  Yobs<- ifelse(U==1, Yst, NA_real_)
  DFp <- data.frame(U=U, Yobs=Yobs, X1=X1, X2=X2, Z1=Z1, Z2=Z2)
  
  sel_formula <- U    ~ X1 + X2
  out_formula <- Yobs ~ Z1 + Z2
  
  ## fits
  st_esn <- make_start_ESN(sel_formula, out_formula, DFp)
  fit1   <- fit_Heckman_ESN(sel_formula, out_formula, DFp, start=st_esn,
                            fix_rho=TRUE, rho_fix_val=0.1, rho_cap=0.95, maxit=2000)
  fit_esn<- fit_Heckman_ESN(sel_formula, out_formula, DFp, start=fit1$par,
                            fix_rho=FALSE, rho_cap=0.95, maxit=2000)
  par_esn<- fit_esn$par
  mm <- fit_esn$mm
  est_esn <- list(
    betaSel = par_esn[1:ncol(mm$Xsel)],
    betaOut = par_esn[(ncol(mm$Xsel)+1):(ncol(mm$Xsel)+ncol(mm$Xout))],
    sigmaY  = exp(par_esn["logSigmaY"]),
    rho     = tanh(par_esn["rho_t"]),
    etaS    = 1.5*tanh(par_esn["lamS_t"]),
    etaY    = 1.5*tanh(par_esn["lamY_t"]),
    alpha   = 1.5*tanh(par_esn["tau_t"])
  )
  
  st_est <- make_start_EST(sel_formula, out_formula, DFp)
  fit_est<- fit_Heckman_EST(sel_formula, out_formula, DFp, start=st_est, maxit=2500)
  par_est<- fit_est$par
  est_est <- list(
    betaSel = par_est[1:ncol(mm$Xsel)],
    betaOut = par_est[(ncol(mm$Xsel)+1):(ncol(mm$Xsel)+ncol(mm$Xout))],
    sigmaY  = exp(par_est["logSigmaY"]),
    rho     = tanh(par_est["rho_t"]),
    etaS    = 1.5*tanh(par_est["lamS_t"]),
    etaY    = 1.5*tanh(par_est["lamY_t"]),
    alpha   = 1.5*tanh(par_est["tau_t"]),
    nu      = 2 + exp(par_est["log_nu_minus2"])
  )
  
  sel_idx <- which(DFp$U==1)
  y_obs <- DFp$Yobs[sel_idx]
  grid  <- seq(quantile(y_obs,.01), quantile(y_obs,.99), length.out=200)
  
  muS_esn <- as.numeric(mm$Xsel %*% est_esn$betaSel)
  muY_esn <- as.numeric(mm$Xout %*% est_esn$betaOut)
  p0_esn  <- vapply(seq_len(nrow(DFp)), function(i)
    pS_le0_ESE_gauss(muS_esn[i], 1, muY_esn[i], est_esn$sigmaY, est_esn$rho,
                     est_esn$etaS, est_esn$etaY, est_esn$alpha), 0.0)
  p_sel_esn <- pmax(1 - p0_esn, 1e-12)
  
  dens_esn <- sapply(grid, function(yy) {
    m <- vapply(sel_idx, function(i)
      int_bivar_spos_ESE_gauss(yy, muS_esn[i], muY_esn[i], 1, est_esn$sigmaY, est_esn$rho,
                               est_esn$etaS, est_esn$etaY, est_esn$alpha) / p_sel_esn[i], 0.0)
    mean(m[is.finite(m)])
  })
  
  muS_est <- as.numeric(mm$Xsel %*% est_est$betaSel)
  muY_est <- as.numeric(mm$Xout %*% est_est$betaOut)
  p0_est  <- vapply(seq_len(nrow(DFp)), function(i)
    pS_le0_ESE_t(muS_est[i], 1, muY_est[i], est_est$sigmaY, est_est$rho,
                 est_est$etaS, est_est$etaY, est_est$alpha, est_est$nu), 0.0)
  p_sel_est <- pmax(1 - p0_est, 1e-12)
  
  dens_est <- sapply(grid, function(yy) {
    m <- vapply(sel_idx, function(i)
      int_bivar_spos_ESE_t(yy, muS_est[i], muY_est[i], 1, est_est$sigmaY, est_est$rho,
                           est_est$etaS, est_est$etaY, est_est$alpha, est_est$nu) / p_sel_est[i], 0.0)
    mean(m[is.finite(m)])
  })
  
  dens_hn <- rep(NA_real_, length(grid))
  hn_ok <- FALSE
  if (DO_HECKMAN_NORMAL && "sampleSelection" %in% .packages(all.available = TRUE)) {
    hN <- try(sampleSelection::selection(selection = sel_formula, outcome = out_formula,
                                         data = DFp, method = "ml"), silent = TRUE)
    if (!inherits(hN,"try-error")) {
      hn_ok <- TRUE
      bS <- as.numeric(coef(hN)[paste0("selection_", colnames(mm$Xsel))]); bS[is.na(bS)] <- 0
      bY <- as.numeric(coef(hN)[paste0("outcome_",   colnames(mm$Xout))]); bY[is.na(bY)] <- 0
      rho_hn <- tryCatch(as.numeric(coef(hN)["rho"]), error=function(e) 0)
      sigmaY_hn <- tryCatch(as.numeric(summary(hN)$estimate["sigma","Estimate"]),
                            error=function(e) tryCatch(as.numeric(hN$sigma), error=function(e2) 1))
      muS_hn <- as.numeric(mm$Xsel %*% bS)
      muY_hn <- as.numeric(mm$Xout %*% bY)
      p0_hn  <- vapply(seq_len(nrow(DFp)), function(i)
        pS_le0_ESE_gauss(muS_hn[i], 1, muY_hn[i], sigmaY_hn, rho_hn, 0,0,0), 0.0)
      p_sel_hn <- pmax(1 - p0_hn, 1e-12)
      dens_hn <- sapply(grid, function(yy) {
        m <- vapply(sel_idx, function(i)
          int_bivar_spos_ESE_gauss(yy, muS_hn[i], muY_hn[i], 1, sigmaY_hn, rho_hn, 0,0,0) / p_sel_hn[i], 0.0)
        mean(m[is.finite(m)])
      })
    }
  }
  
  kde <- density(y_obs, n=200, from=min(grid), to=max(grid))
  pdf(file.path(FIG_DIR, "synthetic_density_overlay.pdf"), width=6, height=5)
  plot(kde$x, kde$y, type="l", lwd=2, xlab="y | u=1", ylab="density")
  lines(grid, dens_est, lwd=2, lty=1)  # EST
  lines(grid, dens_esn, lwd=2, lty=2)  # ESN
  if (hn_ok) lines(grid, dens_hn, lwd=2, lty=3)  # HN
  legend("topright", bty="n",
         legend=c("observed (kde)","EST","ESN", if (hn_ok) "HN" else NULL),
         lty=c(1,1,2,3)[seq_len(3 + as.integer(hn_ok))], lwd=2, col="black")
  dev.off()
  if (svg_ok) {
    svglite::svglite(file.path(FIG_DIR, "synthetic_density_overlay.svg"), width=6, height=5)
    plot(kde$x, kde$y, type="l", lwd=2, xlab="y | u=1", ylab="density")
    lines(grid, dens_est, lwd=2, lty=1); lines(grid, dens_esn, lwd=2, lty=2)
    if (hn_ok) lines(grid, dens_hn, lwd=2, lty=3)
    legend("topright", bty="n",
           legend=c("observed (kde)","EST","ESN", if (hn_ok) "HN" else NULL),
           lty=c(1,1,2,3)[seq_len(3 + as.integer(hn_ok))], lwd=2, col="black")
    dev.off()
  }
  
  ## calibration (by deciles) using EST scores
  muS <- muS_est; muY <- muY_est
  p0  <- p0_est; p_sel <- p_sel_est
  ord <- order(p_sel); b <- cut(seq_along(p_sel)[ord], breaks=10, labels=FALSE)
  pred_b <- tapply(p_sel[ord],  b, function(z) mean(z, na.rm=TRUE))
  obs_b  <- tapply(DFp$U[ord],  b, function(z) mean(z, na.rm=TRUE))
  pdf(file.path(FIG_DIR, "synthetic_calibration.pdf"), width=6, height=5)
  plot(pred_b, obs_b, pch=19, xlab="predicted p(u=1) [EST]", ylab="observed frequency")
  abline(0,1,lty=2); lines(lowess(pred_b, obs_b), lwd=2)
  dev.off()
  if (svg_ok) {
    svglite::svglite(file.path(FIG_DIR, "synthetic_calibration.svg"), width=6, height=5)
    plot(pred_b, obs_b, pch=19, xlab="predicted p(u=1) [EST]", ylab="observed frequency")
    abline(0,1,lty=2); lines(lowess(pred_b, obs_b), lwd=2)
    dev.off()
  }
  
  ## ROC (EST vs HN if available)
  roc_points <- function(score, y) {
    ok <- is.finite(score) & is.finite(y)
    score <- score[ok]; y <- y[ok]
    if (!length(score)) return(rbind(c(0,0), c(1,1)))
    thr <- sort(unique(score), decreasing = TRUE)
    TPR <- FPR <- numeric(length(thr)); P <- sum(y==1); N0 <- sum(y==0)
    for (i in seq_along(thr)) {
      pred <- as.integer(score >= thr[i])
      TP <- sum(pred==1 & y==1); FP <- sum(pred==1 & y==0)
      TPR[i] <- if (P>0) TP/P else 0; FPR[i] <- if (N0>0) FP/N0 else 0
    }
    rbind(c(0,0), cbind(FPR,TPR), c(1,1))
  }
  auc_trapz <- function(curve) { x<-curve[,1]; y<-curve[,2]; sum(diff(x)*(head(y,-1)+tail(y,-1))/2) }
  
  roc_est <- roc_points(p_sel_est, DFp$U); auc_est <- auc_trapz(roc_est)
  roc_esn <- roc_points(p_sel_esn, DFp$U); auc_esn <- auc_trapz(roc_esn)
  
  pdf(file.path(FIG_DIR, "synthetic_roc.pdf"), width=6, height=5)
  plot(roc_est[,1], roc_est[,2], type="l", lwd=2, xlab="fpr", ylab="tpr")
  lines(roc_esn[,1], roc_esn[,2], lty=2, lwd=2)
  abline(0,1,lty=2)
  leg <- c(sprintf("EST (auc=%.3f)", auc_est), sprintf("ESN (auc=%.3f)", auc_esn), "random")
  if (hn_ok) {
    p0_hn <- vapply(seq_len(nrow(DFp)), function(i)
      pS_le0_ESE_gauss(muS_hn[i], 1, muY_hn[i], sigmaY_hn, rho_hn, 0,0,0), 0.0)
    roc_hn <- roc_points(1 - p0_hn, DFp$U); auc_hn <- auc_trapz(roc_hn)
    lines(roc_hn[,1], roc_hn[,2], lty=3, lwd=2)
    leg <- c(sprintf("EST (auc=%.3f)", auc_est),
             sprintf("ESN (auc=%.3f)", auc_esn),
             sprintf("HN  (auc=%.3f)", auc_hn), "random")
  }
  legend("bottomright", bty="n", legend=leg, lty=c(1,2,3,2)[seq_along(leg)], lwd=2, col="black")
  dev.off()
  if (svg_ok) {
    svglite::svglite(file.path(FIG_DIR, "synthetic_roc.svg"), width=6, height=5)
    plot(roc_est[,1], roc_est[,2], type="l", lwd=2, xlab="fpr", ylab="tpr")
    lines(roc_esn[,1], roc_esn[,2], lty=2, lwd=2)
    abline(0,1,lty=2)
    if (hn_ok) lines(roc_hn[,1], roc_hn[,2], lty=3, lwd=2)
    legend("bottomright", bty="n", legend=leg, lty=c(1,2,3,2)[seq_along(leg)], lwd=2, col="black")
    dev.off()
  }
  message("[figures] saved to: ", tryCatch(normalizePath(FIG_DIR), error=function(e) FIG_DIR))
}

message("[done] outputs in: ", tryCatch(normalizePath(OUT_DIR), error=function(e) OUT_DIR))


