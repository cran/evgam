## Setup and fitting functions
#
# .setup.formula
# --- takes evgam formula and returns something useful for fitting
#
# .setup.family
# --- identifies family, its # parameters, their names and its functions
#
# .setup.data
# --- takes data and evgam argument and return something useful for fitting
#
# .setup.pp.data
# --- takes data and converts to something useful for point process model with quadrature points
#
# .sandwich.C
# --- forms C matrix from curvature adjustment as in Chandler & Bate (2007)
#
# .setup.inner.inits
# --- gives initial basis coefficients, having first fitted constant model
#
# .guess
# --- guesses initial log smoothing parameters
#
# .sandwich
# --- updates fitting data if sandwich correction is to be used
#
# .outer
# --- performs outer iterations, i.e. smoothin parameter estimation
#
# .VpVc
# --- computes variance-covariance matrices for basis coefficients
#
# .edf
# --- computes effective degrees of freedom
#
# .swap
# --- swaps elements in early GAM objects, created for fitting, with final estimates
#
# .finalise
# --- adds other useful things to final objects for user functions
#

############ .setup.formula ##########################

.setup.formulae <- function(formula, npar, npar2, data) {
if (inherits(formula, "formula")) {
formula <- list(formula)
}
resp1 <- resp2 <- as.character(formula[[1]])[2]
if (substr(resp2, 1, 5) == "cbind") {
tfmla <- as.formula(paste(c(as.character(formula[[1]])[2:1], "1"), collapse=" "))
resp1 <- colnames(model.extract(model.frame(tfmla, data[1,]), "response"))
resp2 <- resp1[1]
}
for (i in seq_along(formula)) {
charform <- as.character(formula[[i]])
if (attr(terms(formula[[i]]), "response") == 0) charform <- c(charform[1], resp2, charform[-1])
if (substr(charform[2], 1, 5) == "cbind") charform[2] <- resp2
if (length(charform) > 3) charform[3] <- paste(charform[-(1:2)], collapse = " + ")
formula[[i]] <- as.formula(paste(charform[c(2, 1, 3)], collapse=" "))
}
attr(formula, "response.name") <- resp1
if (length(resp1) > 1) message("Censored models currently just average convential log-likelihood. In future cdfs will be used.")
formula
}

.setup.family <- function(family) {
if (family == "gev") {
  lik.fns <- .gevfns
  npar <- 3
  nms <- c("mu", "lpsi", "xi")
} else {
if (family == "gpd") {
  lik.fns <- .gpdfns
  npar <- 2
  nms <- c("lpsi", "xi")
} else {
if (family == "modgpd") {
stop("'family='modgpd'' will return; in the mean time use `family='gpd''")
# gone, but not forgotten
  lik.fns <- NULL
  npar <- 2
  nms <- c("lmodpsi", "xi")
} else {
if (family == "pp") {
  lik.fns <- .ppfns
  npar <- 3
  nms <- c("mu", "lpsi", "xi")
} else {
if (family == "weibull") {
  lik.fns <- .weibfns
  npar <- 2
  nms <- c("llambda", "lk")
} else {
if (family == "exi") {
  lik.fns <- .exifns
  npar <- 1
  nms <- c("location")
} else { 
if (family == "ald") {
  lik.fns <- .aldfns
  npar <- 2
  nms <- c("mu", "lsigma")
} else {
if (family == "gamma") {
  lik.fns <- NULL#gammafns
  npar <- 2
  nms <- c("ltheta", "lk")
} else {
if (family == "orthoggpd") {
stop("'family='orthoggpd'' may not return return")
# gone, and possibly forgotten
  lik.fns <- NULL#ogpdfns
  npar <- 2
  nms <- c("lnu", "xi")
} else {
if (family == "transxigpd") {
stop("'family='transxigpd'' may not return")
# gone, and possibly forgotten
  lik.fns <- NULL#txigpdfns
  npar <- 2
  nms <- c("lpsi", "xi")
} else {
if (family == "transgev") {
stop("'family='transgev'' may not return")
# gone, and possibly forgotten
  lik.fns <- NULL#transgevfns
  npar <- 6
  nms <- c("mu", "lpsi", "xi", "A", "lB", "C")
} else {
if (family == "exponential") {
  lik.fns <- .expfns
  npar <- 1
  nms <- c("llambda")
} else {
if (family == "gauss") {
  lik.fns <- .gaussfns
  npar <- 2
  nms <- c("mu", "logsigma")
}
}
}
}
}
}
}
}
}
}
}
}
}
out <- list(npar=npar, npar2=npar, lik.fns=lik.fns, nms=nms)
}

############ .setup.data ##########################

.setup.data <- function(data, responsename, formula, family, nms, removeData, 
exiargs, aldargs, pp, knots, maxdata, maxspline, compact, sargs, 
outer, trace) {

## data
for (i in seq_along(responsename)) data <- data[!is.na(data[,responsename[i]]),]

if (nrow(data) > maxdata) {
    id <- sort(sample(nrow(data), maxdata))
    data <- data[id,]
    message("`data' truncated to `maxdata' rows. Re-supply `data' to, e.g., `predict.evgam'")
}

if  (compact) {
data.undup <- as.list(data[,unique(unlist(lapply(formula, function(y) unlist(lapply(mgcv::interpret.gam(y)$smooth.spec, function(x) x$term))))), drop=FALSE])
data.undup <- lapply(data.undup, function(x) as.integer(as.factor(x)))
if (length(data.undup) > 1) for (i in 2:length(data.undup)) data.undup[[1]] <- paste(data.undup[[1]], data.undup[[i]], sep=":")
data.undup <- data.undup[[1]]
gc()
unq.id <- which(!duplicated(data.undup))
data.unq <- data.undup[unq.id]
dup.id <- match(data.undup, data.unq)
}

subsampling <- FALSE

## gams
gams <- list()

if (family == "pp") data <- .setup.pp.data(data, responsename, pp)

if (nrow(data) > maxspline) {
for (i in seq_along(formula)) {
id <- sample(nrow(data), maxspline)
gams[[i]] <- mgcv::gam(formula[[i]], data=data[id,], fit=FALSE, knots=knots)[c("smooth", "pterms", "sp")]
}
} else {
for (i in seq_along(formula)) {
gams[[i]] <- mgcv::gam(formula[[i]], data=data, fit=FALSE, knots=knots)[c("smooth", "pterms", "sp")]
}
}
gc()

## likelihood
lik.data <- list()
lik.data$control <- list()
lik.data$outer <- outer
lik.data$control$outer <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=1e-2, stepmax=3)
lik.data$control$inner <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=1e-4, stepmax=1e2)
lik.data$y <- as.matrix(data[,responsename, drop=FALSE])
lik.data$Mp <- sum(unlist(sapply(gams, function(y) c(1, sapply(y$smooth, function(x) x$null.space.dim)))))
lik.data$const <- .5 * lik.data$Mp * log(2 * pi)
lik.data$nobs <- nrow(lik.data$y)
if (family == "exi") {
if (is.null(exiargs$id)) stop("no `id' in `exi.args'.")
if (is.null(exiargs$nexi)) {
    message("`exiargs$nexi' assumed to be 2.")
    exiargs$nexi <- 2
}
if (is.null(exiargs$link)) {
    message("`exiargs$link' assumed to be `logistic'.")
    exiargs$link <- "logistic"
}
lik.data$exiname <- exiargs$id
lik.data$y <- list(lik.data$y, data[,exiargs$id])
lik.data$nexi <- exiargs$nexi
if (exiargs$link == "cloglog") {
lik.data$exilink <- 2
lik.data$linkfn <- function(x) 1 - exp(-exp(x))
}
if (exiargs$link == "logistic") {
lik.data$exilink <- 1
lik.data$linkfn <- function(x) 1 / (1 + exp(-x))
}
if (exiargs$link == "probit") {
lik.data$exilink <- 0
lik.data$linkfn <- function(x) pnorm(x)
}
}
## These will be re-introduced ##
if (family == "pp") {
lik.data$ppw <- attr(data, "weights") # point process quadrature weights
lik.data$y <- as.matrix(rbind(as.matrix(attr(data, "quad")[,responsename]), lik.data$y))
lik.data$ppq <- rep(as.logical(1:0), c(nrow(attr(data, "quad")), nrow(data))) # identify quadrature points
lik.data$cens <- attr(data, "cens")
lik.data$weights <- attr(data, "cweights")
}
if (family == "ald") {
if (is.null(aldargs$tau)) aldargs$tau <- .5
if (is.null(aldargs$C)) aldargs$C <- .5
lik.data$tau <- aldargs$tau
lik.data$C <- aldargs$C
}
lik.data$sandwich <- !is.null(sargs$id)
if (lik.data$sandwich) lik.data$sandwich.split <- data[,sargs$id]
if (!compact) {
lik.data$X <- lapply(gams, function(x) lapply(x$smooth, function(y) mgcv::PredictMat(y, data)))
lik.data$X <- lapply(lik.data$X, function(x) do.call(cbind, x))
for (i in seq_along(gams)) lik.data$X[[i]] <- cbind(model.matrix(gams[[i]]$pterms, data), lik.data$X[[i]])
if (family == "pp") {
ppX <- lapply(gams, function(x) lapply(x$smooth, function(y) mgcv::PredictMat(y, attr(data, "quad"))))
ppX <- lapply(ppX, function(x) do.call(cbind, x))
for (i in seq_along(gams)) {
    ppX[[i]] <- cbind(model.matrix(gams[[i]]$pterms, attr(data, "quad")), ppX[[i]])
    lik.data$X[[i]] <- rbind(ppX[[i]], lik.data$X[[i]])
    }
}
} else {
lik.data$X <- lapply(gams, function(x) lapply(x$smooth, function(y) mgcv::PredictMat(y, data[unq.id,])))
lik.data$X <- lapply(lik.data$X, function(x) do.call(cbind, x))
for (i in seq_along(gams)) lik.data$X[[i]] <- cbind(model.matrix(gams[[i]]$pterms, data[unq.id,]), lik.data$X[[i]])
}
for (i in seq_along(gams)) if (length(gams[[i]]$smooth) == 0) gams[[i]]$smooth <- list(list(first.para=NULL, last.para=ncol(lik.data$X[[i]])))
for (i in seq_along(gams)) {
class(gams[[i]]) <- "gamlist"
if (removeData) gams[[i]]$y <- NULL
}
if (compact) {
lik.data$dupid <- dup.id - 1
lik.data$duplicate <- 1
} else {
lik.data$dupid <- -1
lik.data$duplicate <- 0
}
if (length(lik.data$X) == 1 & length(nms) > 1) {
for (i in 2:length(nms)) {
lik.data$X[[i]] <- lik.data$X[[1]]
gams[[i]] <- gams[[1]]
}
}
nbk <- sapply(lik.data$X, ncol)
lik.data$nb <- sum(nbk)
lik.data$idpars <- rep(seq_along(lik.data$X), nbk)
lik.data$LAid <- lik.data$idpars > 0
lik.data$subsampling <- subsampling
gotsmooth <- which(sapply(gams, function(x) length(x$sp)) > 0)
lik.data$k <- 1
if (is.null(sargs$id)) {
    lik.data$adjust <- 0
} else {
if (is.null(sargs$method)) sargs$method <- "magnitude"
    if (sargs$method == "curvature") {
        if (trace > 0) message(paste("Sandwich adjustment method: curvature"))
        lik.data$adjust <- 2
    } else {
        if (trace > 0) message(paste("Sandwich adjustment method: magnitude"))
        lik.data$adjust <- 1
    }
}
if (is.null(sargs$force)) sargs$force <- FALSE
lik.data$force <- sargs$force
list(lik.data=lik.data, gotsmooth=gotsmooth, data=data, gams=gams, sandwich=lik.data$adjust > 0)
}

.setup.pp.data <- function(data, responsename, pp) {
data$row <- seq_len(nrow(data))
ds <- split(data, data[,pp$id])
wts <- pp$ny
if (length(wts) == 1) {
  wts <- rep(wts, length(ds))
} else {
  wts <- wts[match(names(ds), names(wts))]
}
nobs2 <- sapply(ds, nrow)
enough <- nobs2 > pp$r
if (any(!enough)) warning(paste(sum(!enough), "unique pp.args$id removed for having fewer than r observations."))
ds <- ds[enough]
wts <- wts[enough]
nid <- sum(enough)
data.quad <- do.call(rbind, lapply(ds, function(x) x[1,]))
if (pp$r != -1) {
    du <- sapply(ds, function(x) x[order(x[,responsename], decreasing=TRUE)[pp$r], responsename])
} else {
    du <- sapply(ds, function(x) x[, responsename])
}
data.quad[,responsename] <- du
ds <- lapply(seq_len(nid), function(i) subset(ds[[i]], ds[[i]][,responsename] >= du[i]))
out <- dfbind(ds)
attr(out, "weights") <- wts
attr(out, "quad") <- data.quad
if (is.null(pp$cens)) {
  attr(out, "cens") <- NULL
} else {
  attr(out, "cens") <- data[out$row, pp$cens]
}
if (is.null(pp$weights)) {
  attr(out, "cweights") <-rep(1, length(out$row))
} else {
  attr(out, "cweights") <- pp$weights[out$row]
}
out
}

############ .sandwich.C ##########################

.sandwich.C <- function(H, J) {
iJ <- pinv(J)
HA <- crossprod(H, crossprod(iJ, H))
sH <- svd(H)
M <- sqrt(sH$d) * t(sH$v)
sHA <- svd(HA)
MA <- sqrt(sHA$d) * t(sHA$v)
solve(M, MA)
}

############ .setup.inner.inits ##########################

.setup.inner.inits <- function(inits, likdata, likfns, npar, family) {

likdata0 <- likdata
likdata0$X <- lapply(seq_along(likdata$X), function(x) matrix(1, nrow=nrow(likdata$X[[x]]), ncol=1))
likdata0$pp$X <- lapply(seq_along(likdata$pp$X), function(x) matrix(1, nrow=nrow(likdata$pp$X[[x]]), ncol=1))
likdata0$S <- diag(0, npar)
likdata0$idpars <- seq_len(npar)

if (is.null(inits)) {
if (npar == 1) {
inits <- 2
}
if (npar == 2) {
inits <- c(log(mean(likdata$y[,1])), .05)
if (family == "transxigpd") inits[2] <- .9
if (family == "ald") {
    inits[1] <- quantile(likdata0$y[,1], likdata0$tau)
    inits[2] <- log(sd(likdata0$y[,1]))
}
}
if (npar == 3) {
inits <- c(sqrt(6) * sd(likdata0$y[,1]) / pi, .05)
inits <- c(mean(likdata0$y[,1]) - .5772 * inits[1], log(inits[1]), inits[2])
}
if (npar == 6) {
inits <- c(sqrt(6) * sd(likdata0$y[,1]) / pi, .05)
inits <- c(mean(likdata0$y[,1]) - .5772 * inits[1], log(inits[1]), inits[2])
inits <- c(inits, 0, 0, 1)
}
likdata0$CH <- diag(length(inits))
likdata0$compmode <- numeric(length(inits))
beta0 <- .newton_step_inner(inits, .nllh.nopen, .search.nopen, likdata=likdata0, likfns=likfns, control=likdata$control$inner)$par
} else {
if (is.list(inits)) {
betamat <- expand.grid(inits)
betanllh <- numeric(nrow(betamat))
for (i in seq_len(nrow(betamat))) {
beta0 <- unlist(betamat[i,])
betanllh[i] <- likfns$nllh(beta0, likdata0)
}
beta0 <- betamat[which.min(betanllh),]
print(beta0)
} else {
beta0 <- inits
}
}
beta0 <- unlist(lapply(seq_len(npar), function(i) c(beta0[i], rep(0, ncol(likdata$X[[i]]) - 1))))
compmode <- 0 * beta0
CH <- diag(compmode + 1)
k <- 1
likdata[c("k", "CH", "compmode")] <- list(k, CH, compmode)
diagH <- diag(.gH.nopen(beta0, likdata=likdata, likfns=likfns)[[2]])
if (likdata$sandwich) {
beta0 <- .newton_step(beta0, .nllh.nopen, .search.nopen, likdata=likdata, likfns=likfns, control=likdata$control$inner)$par
H <- .gH.nopen(beta0, likdata=likdata, likfns=likfns, sandwich=TRUE)
if (family == "pp") {
J0 <- H[[1]]
J <- J0[,!likdata$ppq]
J0 <- rowSums(J0[,likdata$ppq])
J <- split(as.data.frame(t(J)), likdata$sandwich.split)
wts <- sapply(J, nrow)
wts <- wts / sum(wts)
J <- sapply(J, colSums)
J <- J + J0 %o% wts
J <- tcrossprod(J)
} else {
J <- split(as.data.frame(t(H[[1]])), likdata$sandwich.split)
J <- sapply(J, colSums)
J <- tcrossprod(J)
}
H <- H[[2]]
diagH <- diag(H)
cholH <- try(chol(H), silent=TRUE)
if (inherits(cholH, "try-error")) {
    if (!likdata$force) {
        stop("Hessian of unpenalised MLE not positive definite.\n  Supply `force=TRUE' to `sandwich.args' to perturb it to be positive definite.")
    } else {
        message("Hessian perturbed to be positive definite for sandwich adjustment.")
        iH <- pinv(H)
    }
} else {
    iH <- chol2inv(cholH)
}
if (likdata$adjust == 2) {
  cholJ <- try(chol(J), silent=TRUE)
  if (inherits(cholJ, "try-error") & likdata$adjust == 2) {
    HA <- crossprod(backsolve(cholJ, H, transpose=TRUE))
  } else {
    iHA <- tcrossprod(crossprod(iH, J), iH)
    choliHA <- try(chol(iHA), silent=TRUE)
    if (inherits(choliHA, "try-error")) {
      if (!likdata$force) {
        stop("Sandwich variance not positive definite.\n  Supply `force=TRUE' to `sandwich.args' to perturb it to be positive definite.")
      } else {
        message("Sandwich variance perturbed to be positive definite.")
        HA <- pinv(iHA)
      }
    } else {
      HA <- chol2inv(choliHA)
    }
  }
sH <- svd(H)
M <- sqrt(sH$d) * t(sH$v)
sHA <- svd(HA)
MA <- sqrt(sHA$d) * t(sHA$v)
CH <- solve(M, MA)
compmode <- beta0
} else {
k <- 1 / mean(diag(crossprod(iH, J)))
}
}
attr(beta0, "k") <- k
attr(beta0, "CH") <- CH
attr(beta0, "compmode") <- compmode
attr(beta0, "diagH") <- diagH
beta0
}

############ .guess ##########################

.guess <- function(x, d, s) {
okay <- s != 0
val <- d / (d + exp(x) * s)
mean(val[okay]) - .4
}

############ .sandwich ##########################

.sandwich <- function(likdata, beta) {
likdata$k <- attr(beta, "k")
likdata$CH <- attr(beta, "CH")
likdata$compmode <- attr(beta, "compmode")
bigX <- do.call(cbind, likdata$X)
CHX <- bigX %*% likdata$CH
CHX <- lapply(unique(likdata$idpars), function(i) CHX[,likdata$idpars == i])
likdata$CHX <- CHX
likdata
}

############ .outer ##########################

.outer <- function(rho0, beta, likfns, likdata, Sdata, control, correctV, outer, trace) {

attr(rho0, "beta") <- beta

if (outer == "newton") {
  fit.reml <- .newton_step_inner(rho0, .reml0, .search.reml, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
} else {
  if (outer == "fd") {
    fit.reml <- .BFGS(rho0, .reml0, .reml1.fd, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
  } else {
    fit.reml <- .BFGS(rho0, .reml0, .reml1, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
  }
  rho1 <- fit.reml$par
  attr(rho1, "beta") <- fit.reml$beta
  fit.reml$Hessian <- .reml12(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata)[[2]]
}

fit.reml$invHessian <- .solve_evgam(fit.reml$Hessian)

fit.reml$trace <- trace

if (trace == 1) {
    report <- "\n Final max(|grad|))"
    likdata$S <- .makeS(Sdata, exp(fit.reml$par))
    report <- c(report, paste("   Inner:", signif(max(abs(.gH.pen(fit.reml$beta, likdata, likfns)[[1]])), 3)))
    report <- c(report, paste("   Outer:", signif(max(abs(fit.reml$gradient)), 3)))
    report <- c(report, "", "")
    cat(paste(report, collapse="\n"))
}

fit.reml

}

############ .outer.nosmooth ##########################

.outer.nosmooth <- function(beta, likfns, likdata, control, trace) {

fit.inner <- .newton_step(beta, .nllh.nopen, .search.nopen, likdata=likdata, likfns=likfns, control=likdata$control$inner)

list(beta=fit.inner$par)

}

############ .VpVc ##########################

.VpVc <- function(fitreml, likfns, likdata, Sdata, correctV, sandwich, smooths) {
lsp <- fitreml$par
H0 <- .gH.nopen(fitreml$beta, likdata, likfns)[[2]]
if (smooths) {
  sp <- exp(lsp)
  H <- H0 + likdata$S
} else {
  H <- H0
}
cholH <- try(chol(H), silent=TRUE)
if (inherits(cholH, "try-error"))
  message("Final Hessian of negative penalized log-likelihood not numerically positive definite.")
Vc <- Vp <- pinv(H)
if (smooths) {
if (correctV) {
cholVp <- try(chol(Vp), silent=TRUE)
if (inherits(cholVp, "try-error")) {
    cholVp <- attr(.perturb(Vp), "chol")
}
spSl <- lapply(seq_along(sp), function(i) sp[i] * attr(Sdata, "Sl")[[i]])
attr(lsp, "beta") <- fitreml$beta
dbeta <- .dbeta(lsp, spSl, .Hdata(H), fitreml$beta, likdata, likfns, deriv=1)$d1
Vrho <- fitreml$invHessian
Vbetarho <- tcrossprod(dbeta %*% Vrho, dbeta)
# eps <- 1e-4
# R0 <- .grad.R(fitreml$par, Sdata=Sdata, R=0, eps=1, likfns=likfns, likdata=likdata, H0=H0)
# dR <- lapply(seq_along(sp), function(i) grad.R(replace(fitreml$par, i, fitreml$par[i] + eps), Sdata=Sdata, R=R0, eps=eps, likfns=likfns, likdata=likdata, H0=H0))
VR <- matrix(0, nrow=likdata$nb, ncol=likdata$nb)
# for (k in seq_along(sp)) for (l in seq_along(sp)) VR <- VR + crossprod(dR[[k]] * Vrho[k, l], dR[[l]])
# VR <- .5 * (VR + t(VR))
Vc <- .perturb(Vp + Vbetarho + VR)
}
} else {
  Vrho <- 0
}
list(Vp=Vp, Vc=Vc, Vlsp=Vrho, H0=H0, H=H)
}

############ .edf ##########################

.edf <- function(beta, likfns, likdata, VpVc, sandwich) {
diag(crossprod(VpVc$Vp, VpVc$H0))
}

############ .swap ##########################

.swap <- function(fitreml, gams, likdata, VpVc, gotsmooth, edf, smooths) {
Vp <- VpVc$Vp
Vc <- VpVc$Vc
if (smooths) {
  sp <- split(exp(fitreml$par), unlist(sapply(seq_along(gams), function(x) rep(x, length(gams[[x]]$sp)))))
}
for (i in seq_along(gams)) {
idi <- likdata$idpars == i
gams[[i]]$coefficients <- fitreml$beta[idi]
gams[[i]]$Vp <- Vp[idi, idi, drop = FALSE]
gams[[i]]$Vc <- Vc[idi, idi, drop = FALSE]
if (i %in% gotsmooth) gams[[i]]$sp <- sp[[i]]
gams[[i]]$edf <- edf[idi]
}
gams
}

############ .finalise ##########################

.finalise <- function(gams, data, likfns, likdata, Sdata, fitreml, VpVc, family, gotsmooth,
formula, responsenm, removeData, edf) {
nms <- c("location", "logscale", "shape")
if (length(gams) == 2) {
  if (family %in% c("ald", "gauss")) {
    nms <- nms[1:2]
  } else {
    nms <- nms[-1]
}}
if (family == "exponential") nms <- "lograte"
if (family == "weibull") nms[2] <- "logshape"
names(gams) <- nms
smooths <- length(gotsmooth) > 0
Vp <- VpVc$Vp
Vc <- VpVc$Vc
if (smooths) gams$sp <- exp(fitreml$par)
gams$nobs <- likdata$nobs
gams$logLik <- -1e20
fit.lik <- list(convergence=0)
if (fit.lik$convergence == 0) {
gams$logLik <- -.nllh.nopen(fitreml$beta, likdata, likfns)
gams$logLik <- gams$logLik - likdata$const
}
gams$AIC <- gams$BIC <- -2 * gams$logLik
if (fit.lik$convergence != 0) gams$AIC <- gams$BIC <- 1e20
attr(gams, "df") <- sum(edf)
gams$AIC <- gams$AIC + 2 * attr(gams, "df")
gams$BIC <- gams$BIC + attr(gams, "df") * log(gams$nobs)
gams$simulate <- list(mu=fitreml$beta, Sigma=Vp)
gams$family <- family
gams$idpars <- likdata$idpars
names(formula) <- names(gams)[seq_along(formula)]
gams$call <- formula
gams$response.name <- responsenm
gams$gotsmooth <- gotsmooth
if (!removeData) {
    if (family == "pp") {
        gams$data <- attr(data, "quad")
    } else {
        gams$data <- data
    }
}
gams$Vc <- Vc
gams$Vp <- Vp
gams$Vlsp <- VpVc$Vlsp
gams$negREML <- fitreml$objective
gams$coefficients <- fitreml$beta
if (family == "ald") gams$tau <- likdata$tau
if (family == "exi") {
gams$linkfn <- likdata$linkfn
gams$exi.name <- likdata$exiname
}
for (i in seq_along(likdata$X)) {
gams[[i]]$X <- likdata$X[[i]]
if (likdata$dupid[1] != -1) gams[[i]]$X <- gams[[i]]$X[likdata$dupid + 1,]
gams[[i]]$fitted <- as.vector(likdata$X[[i]] %*% gams[[i]]$coefficients)
}
gams$likdata <- likdata
gams$likfns <- likfns
if (smooths) gams$Sdata <- Sdata
gams$formula <- formula
gams$compacted <- likdata$dupid[1] != -1
if (gams$compacted) gams$compactid <- likdata$dupid + 1
smooth.terms <- unique(lapply(lapply(gams[gotsmooth], function(x) x$smooth), function(y) lapply(y, function(z) z$term)))
smooth.terms <- unique(unlist(smooth.terms, recursive=FALSE))
gams$plotdata <- lapply(smooth.terms, function(x) unique(data[,x, drop=FALSE]))
if (family == "weibull") names(gams)[2] <- "logshape"
if (family == "exponential") names(gams)[1] <- "lograte"
gams$ngam <- length(formula)
class(gams) <- "evgam"
return(gams)
}

