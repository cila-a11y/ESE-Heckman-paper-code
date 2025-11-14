#!/usr/bin/env Rscript
## ============================================================================
## Mroz87 — ESN vs Heckman-Normal with guaranteed nesting (λ=0, α=0)
## - Apenas Mroz87 (sem PNADC)
## - Wrappers robustos de CDF bivariada (evitam NA)
## - ESN(λ=α=0) ~ HN; LRT (df=3) para ESN livre
## ============================================================================

## ------------------- tiny installer (CRAN) ----------------------------------
.need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("[setup] installing '%s' from CRAN ...", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  TRUE
}
.need("mvtnorm")
.need("sampleSelection")

## ------------------- threading hygiene --------------------------------------
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

## ------------------- integradores (Genz–Bretz) ------------------------------
INT_MAXPTS <- 2e5; INT_ABSEPS <- 1e-8; INT_RELEPS <- 1e-8
ALG_GAUSS  <- mvtnorm::GenzBretz(maxpts = INT_MAXPTS, abseps = INT_ABSEPS, releps = INT_RELEPS)

## ------------------- helpers básicos ----------------------------------------
.clamp01 <- function(x, eps = 1e-12) pmin(pmax(x, eps), 1 - eps)

Omega_mat <- function(sigmaS, sigmaY, rho) {
  matrix(c(
    sigmaS^2,              rho * sigmaS * sigmaY,
    rho * sigmaS * sigmaY,        sigmaY^2
  ), 2, 2, byrow = TRUE)
}

## λ_eff = Ω^{-1/2} η  (via Cholesky Ω = R'R -> R^{-1} η)
lambda_eff_from_eta <- function(Omega, eta) {
  Omega <- (Omega + t(Omega))/2
  R <- chol(Omega)                    # upper triangular
  drop(backsolve(R, eta))             # R^{-1} eta
}

PE_gauss_canonical <- function(alpha, eta) {
  .clamp01(stats::pnorm(alpha / sqrt(1 + sum(eta^2))))
}

## ------------------- wrappers robustos p/ CDF bivariada ----------------------
.safe_biv_pmvnorm <- function(lower, upper, rho, mean = c(0,0)) {
  # 1) tentativa com parâmetros “fortes”
  out <- tryCatch(as.numeric(mvtnorm::pmvnorm(
    lower = lower, upper = upper, mean = mean,
    corr  = matrix(c(1, rho, rho, 1), 2, 2),
    algorithm = ALG_GAUSS
  )), error = function(e) NA_real_)
  if (is.finite(out)) return(out)
  # 2) fallback com defaults do pacote
  out2 <- tryCatch(as.numeric(mvtnorm::pmvnorm(
    lower = lower, upper = upper, mean = mean,
    corr  = matrix(c(1, rho, rho, 1), 2, 2)
  )), error = function(e) NA_real_)
  if (is.finite(out2)) return(out2)
  # 3) fallback aproximado p/ |rho| ~ 0 (produto das marginais)
  if (is.finite(rho) && abs(rho) < 1e-12) {
    a1 <- stats::pnorm(upper[1]) - stats::pnorm(lower[1])
    a2 <- stats::pnorm(upper[2]) - stats::pnorm(lower[2])
    return(max(0, a1) * max(0, a2))
  }
  # 4) último recurso: devolver algo muito pequeno mas >0 (evita NA)
  1e-300
}

## ------------------- blocos observáveis ESN (núcleo gaussiano) --------------
## P(S <= 0 | A_α)
pS_le0_ESE_gauss <- function(muS, sigmaS, muY, sigmaY, rho, etaS, etaY, alpha) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam   <- lambda_eff_from_eta(Omega, c(etaS, etaY)); lamS <- lam[1]; lamY <- lam[2]
  
  varW <- 1 + (etaS^2 + etaY^2); sdW <- sqrt(max(varW, 1e-12))
  covSW <- - lamS * sigmaS^2 - lamY * rho * sigmaS * sigmaY
  rSW <- covSW / (sigmaS * sdW); if (!is.finite(rSW)) rSW <- 0
  rSW <- max(min(rSW, 0.999999), -0.999999)
  
  hS <- (0 - muS) / sigmaS
  hW <- alpha / sdW
  if (!all(is.finite(c(hS, hW)))) return(0.5)  # fallback neutro
  PE <- PE_gauss_canonical(alpha, c(etaS, etaY))
  if (!is.finite(PE) || PE <= 0) PE <- 0.5
  
  num <- .safe_biv_pmvnorm(lower = c(-Inf, -Inf), upper = c(hS, hW), rho = rSW, mean = c(0,0))
  .clamp01(num / PE)
}

## f_Y(y) * P(S>0, W<=α | Y=y) / PE
int_bivar_spos_ESE_gauss <- function(y, muS, muY, sigmaS, sigmaY, rho, etaS, etaY, alpha) {
  Omega <- Omega_mat(sigmaS, sigmaY, rho)
  lam   <- lambda_eff_from_eta(Omega, c(etaS, etaY)); lamS <- lam[1]; lamY <- lam[2]
  
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
  if (!all(is.finite(c(h1, h2)))) return(.Machine$double.xmin)
  
  phiY <- stats::dnorm((y - muY)/sigmaY) / sigmaY
  PE   <- PE_gauss_canonical(alpha, c(etaS, etaY)); if (!is.finite(PE) || PE <= 0) PE <- 0.5
  
  p_sel_cond <- .safe_biv_pmvnorm(lower = c(h1, -Inf), upper = c(Inf, h2), rho = r_cond, mean = c(0,0))
  max(phiY * p_sel_cond / PE, .Machine$double.xmin)
}

## ------------------- matrizes do modelo & logLik ESN -------------------------
model_matrices_Heckman <- function(formSel, formOut, data, selVar="U", outVar="Yobs") {
  selVars <- all.vars(formSel); outVars <- all.vars(formOut)
  dt <- data[, unique(c(selVars, outVars, selVar, outVar)), drop = FALSE]
  trS <- stats::terms(formSel, data = dt)
  trY <- stats::terms(formOut, data = dt)
  Xs <- stats::model.matrix(stats::delete.response(trS), data = dt, na.action = stats::na.pass)
  Xy <- stats::model.matrix(stats::delete.response(trY), data = dt, na.action = stats::na.pass)
  U  <- as.integer(dt[[selVar]])
  Y  <- as.numeric(dt[[outVar]])
  stopifnot(all(U %in% c(0,1)))
  list(Xsel = Xs, Xout = Xy, U = U, Y = Y)
}

## retorna NEGATIVE log-likelihood (p/ minimização)
logLik_ESN_gauss <- function(params, Xs, Xy, U, Y, nB, nG, L_eta=1.5, L_alpha=1.5) {
  tiny <- .Machine$double.xmin
  betaS <- params[1:nB]
  betaY <- params[(nB+1):(nB+nG)]
  logSigY <- params[nB+nG+1]; kappa <- params[nB+nG+2]
  eS_t <- params[nB+nG+3]; eY_t <- params[nB+nG+4]; a_t <- params[nB+nG+5]
  
  sigmaS <- 1; sigmaY <- exp(logSigY); rho <- tanh(kappa)
  etaS <- L_eta * tanh(eS_t); etaY <- L_eta * tanh(eY_t); alpha <- L_alpha * tanh(a_t)
  
  ll <- 0
  for (i in seq_along(U)) {
    muS <- sum(Xs[i,] * betaS); muY <- sum(Xy[i,] * betaY)
    if (U[i] == 0L) {
      p0 <- pS_le0_ESE_gauss(muS, 1, muY, sigmaY, rho, etaS, etaY, alpha)
      if (!is.finite(p0) || p0 <= 0) return(1e12)
      ll <- ll + log(max(p0, tiny))
    } else {
      if (is.na(Y[i])) return(1e12) # seleção com Y faltante: inválido
      joint <- int_bivar_spos_ESE_gauss(Y[i], muS, muY, 1, sigmaY, rho, etaS, etaY, alpha)
      if (!is.finite(joint) || joint <= 0) return(1e12)
      ll <- ll + log(max(joint, tiny))
    }
  }
  -ll
}

fit_ESN_free <- function(Xs, Xy, U, Y, start, L_eta=1.5, L_alpha=1.5, maxit=4000, trace=0) {
  nB <- ncol(Xs); nG <- ncol(Xy)
  pnames <- c(paste0("betaSel_", 1:nB), paste0("betaOut_", 1:nG),
              "logSigmaY", "rho_t", "lamS_t", "lamY_t", "tau_t")
  names(start) <- pnames
  
  lower <- rep(-Inf, length(start)); upper <- rep( Inf, length(start))
  names(lower) <- pnames; names(upper) <- pnames
  lower["logSigmaY"] <- log(1e-2); upper["logSigmaY"] <- log(1e2)
  lower["rho_t"]     <- atanh(-0.995); upper["rho_t"] <- atanh(0.995)
  lower[c("lamS_t","lamY_t","tau_t")] <- -3; upper[c("lamS_t","lamY_t","tau_t")] <-  3
  
  nLL <- function(par) {
    val <- logLik_ESN_gauss(par, Xs, Xy, U, Y, nB, nG, L_eta, L_alpha)
    if (is.finite(val)) val else 1e12
  }
  
  fit <- stats::optim(par = start, fn = nLL, method = "L-BFGS-B",
                      lower = lower, upper = upper,
                      control = list(maxit = maxit, trace = trace, factr = 1e7, pgtol = 1e-6))
  list(par = fit$par,
       ll_pure = -logLik_ESN_gauss(fit$par, Xs, Xy, U, Y, nB, nG, L_eta, L_alpha),
       conv = fit$convergence)
}

decode_ESN_params <- function(par, L_eta=1.5, L_alpha=1.5, nB, nG) {
  list(
    betaSel = par[1:nB],
    betaOut = par[(nB+1):(nB+nG)],
    sigmaY  = exp(par[nB+nG+1]),
    rho     = tanh(par[nB+nG+2]),
    etaS    = L_eta   * tanh(par[nB+nG+3]),
    etaY    = L_eta   * tanh(par[nB+nG+4]),
    alpha   = L_alpha * tanh(par[nB+nG+5])
  )
}

## ------------------- 1) Carregar Mroz87 e construir variáveis ----------------
data("Mroz87", package = "sampleSelection")
DF <- transform(
  Mroz87,
  U    = as.integer(is.finite(wage) & wage > 0),
  Yobs = ifelse(is.finite(wage) & wage > 0, log(wage), NA_real_)
)

sel_form <- U    ~ nwifeinc + educ + exper + I(exper^2) + age + kids5 + kids618
out_form <- Yobs ~ educ + exper + I(exper^2)

## Eliminar apenas NAs nos preditores (Yobs pode ser NA quando U=0)
vars_need <- unique(c(all.vars(sel_form), all.vars(out_form), "U", "Yobs"))
keep <- stats::complete.cases(DF[, setdiff(vars_need, "Yobs"), drop = FALSE])
DF <- DF[keep, , drop = FALSE]

mm  <- model_matrices_Heckman(sel_form, out_form, DF, "U", "Yobs")
Xs <- mm$Xsel; Xy <- mm$Xout; U <- mm$U; Y <- mm$Y
nB <- ncol(Xs); nG <- ncol(Xy)
cat("\n[Mroz87] n =", nrow(DF), " | sel_rate =", round(mean(U), 4), "\n")

## ------------------- 2) Heckman-Normal baseline (sampleSelection) ------------
hN <- sampleSelection::selection(selection = sel_form, outcome = out_form, data = DF, method = "ml")
ll_HN_pkg <- as.numeric(stats::logLik(hN))

## ------------------- 3) ESN restrito no canto HN (λ=0, α=0) -----------------
coef_hN <- stats::coef(hN)
betaSel_start <- as.numeric(coef_hN[paste0("selection_", colnames(Xs))])
betaOut_start <- as.numeric(coef_hN[paste0("outcome_",   colnames(Xy))])
rho_start     <- suppressWarnings(as.numeric(coef_hN["rho"])); if (!is.finite(rho_start)) rho_start <- 0
sigmaY_start  <- tryCatch({
  sm <- summary(hN); as.numeric(sm$estimate["sigma","Estimate"])
}, error = function(e) if (!is.null(hN$sigma)) as.numeric(hN$sigma) else 1)

start_corner <- c(
  betaSel_start, betaOut_start,
  log(max(sigmaY_start, 1e-6)),
  atanh(max(min(rho_start, 0.995), -0.995)),
  0, 0, 0  # lamS_t=0, lamY_t=0, tau_t=0  => etaS=etaY=alpha=0 (canto HN)
)

## Avaliar LL ESN no canto HN com wrappers seguros (deve ≍ HN)
ll_ESN_restr <- -logLik_ESN_gauss(start_corner, Xs, Xy, U, Y, nB, nG, L_eta=1.5, L_alpha=1.5)

cat("\n[CHECK nesting] logLik_HN(pkg) =", round(ll_HN_pkg, 6),
    " | logLik_ESN(λ=α=0) =", round(ll_ESN_restr, 6),
    " | diff =", round(ll_ESN_restr - ll_HN_pkg, 6), "\n")

## ------------------- 4) ESN livre (a partir do canto) + LRT ------------------
fit_free    <- fit_ESN_free(Xs, Xy, U, Y, start = start_corner, L_eta=1.5, L_alpha=1.5, maxit=4000)
ll_ESN_free <- fit_free$ll_pure
LRT <- 2 * (ll_ESN_free - ll_HN_pkg)
p_value <- stats::pchisq(LRT, df = 3, lower.tail = FALSE)

cat("\n[ESN livre] ll =", round(ll_ESN_free, 3),
    " | LRT (df=3) =", round(LRT, 3),
    " | p-value =", signif(p_value, 3), "\n")

## ------------------- 5) guardar resultados ----------------------------------
OUTDIR <- file.path("results_real","Mroz87_ESN_HN")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

sum_tab <- data.frame(
  n = nrow(DF), sel_rate = mean(DF$U),
  logLik_HN = ll_HN_pkg,
  logLik_ESN_restr = ll_ESN_restr,
  logLik_ESN_free = ll_ESN_free,
  LRT_df3 = LRT, p_value = p_value
)
utils::write.csv(sum_tab, file.path(OUTDIR, "summary_nesting_check.csv"), row.names = FALSE)

saveRDS(list(hN = hN,
             start_corner = start_corner,
             fit_ESN_free = fit_free,
             mm = mm,
             session = sessionInfo()),
        file.path(OUTDIR, "fit_objects.rds"))

cat("\n[Saved] ", tryCatch(normalizePath(OUTDIR), error=function(e) OUTDIR), "\n")

