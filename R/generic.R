#' Fitting generalised additive extreme-value family models
#'
#' Function \code{evgam} fits generalised additive extreme-value models. It allows 
#' the fitting of various extreme-value models, including the generalised 
#' extreme value and Pareto distributions. It can also perform quantile regression
#' via the asymmetric Laplace dsitribution.
#'
#' @param formula a list of formulae for location, scale and shape parameters, as in \link[mgcv]{gam}
#' @param data a data frame
#' @param family a character string giving the type of family to be fitted; defaults to \code{"gev"}
#' @param correctV logicial: should the variance-covariance matrix include smoothing parameter uncertainty? Defaults to \code{TRUE}
#' @param rho0 a scalar or vector of initial log smoothing parameter values; a scalar will be repeated if there are multiple smoothing terms
#' @param inits a vector or list giving initial values for constant basis coefficients; if a list, a grid is formed using \link[base]{expand.grid}, and the `best' used; defaults to \code{NULL}, so initial values are automatically found
#' @param outer a character string specifying the outer optimiser is full \code{"Newton"}, \code{"BFGS"} or uses finite differences, \code{"FD"}; defaults to \code{"BFGS"}
#' @param control a list of lists of control parameters to pass to inner and outer optimisers; defaults to \code{evgam.control()}
#' @param removeData logical: should \code{data} be removed from \code{evgam} object? Defaults to \code{FALSE}
#' @param trace an integer specifying the amount of information supplied about fitting; defaults to \code{1}
#' @param knots passed to \link[mgcv]{s}; defaults to \code{NULL}
#' @param maxdata an integer specifying the maximum number of \code{data} rows. \code{data} is sampled if its number of rows exceeds \code{maxdata}; defaults to \code{1e20}
#' @param maxspline an integer specifying the maximum number of \code{data} rows used for spline construction; defaults to \code{1e20}
#' @param compact logical: should duplicated \code{data} rows be compacted? Defaults to \code{FALSE}
#' @param ald.args a list of arguments for \code{family="ald"}; see Details
#' @param exi.args a list of arguments for \code{family="exi"}; see Details
#' @param pp.args a list of arguments for \code{family="pp"}; see Details
#' @param sandwich.args a list of arguments for sandwich adjustment; see Details
#' 
#' @details
#' 
#' The following families are currently available: \code{"ald"}, the asymmetric Laplace distribution,
#' primarily intended for quantile regression, as in Yu & Moyeed (2001); \code{"gev"} (default), the
#' generalised extreme valued distribution; \code{"exp"}, the exponential distribution; \code{"gpd"},
#' the generalised Pareto distribution; \code{"gauss"}, the Gaussian distribution; \code{"pp"}, the 
#' point process model for extremes, implemented through \eqn{r}-largest order statistics; \code{"weibull"}, the Weibull distribution; \code{"exi"}, estimation if the
#' extremal index, as in Schlather & Tawn (2003).
#' 
#' Arguments for the asymmetric Laplace distribution are given by \code{ald.args}. A 
#' scalar \code{tau} defines the quantile sought, which has no default. The scalar
#' \code{C} specifies the curvature parameter of Oh et al. (2011).
#'
#' Arguments for extremal index estimation are given by \code{exi.args}. A character
#' string \code{id} specifies the variable in \code{data}over which an \code{nexi}
#' (default 2) running max. has been taken. The \code{link} is specified as a character string,
#' which is one of \code{"logistic"}, \code{"probit"}, \code{"cloglog"}; defaults to \code{"logistic"}.
#' 
#' Arguments for the point process model are given by \code{pp.args}. An integer \code{r}
#' specifies the number of order statistics from which the model will be estimated.
#' If \code{r = -1}, all \code{data} will be used. The character string \code{id} specifies the variable 
#' in \code{data} over which the point process isn't integrated; e.g. if a map 
#' of parameter estimates related to extremes over time is sought, integration isn't 
#' over locations. The scalar \code{nper} number of data per period of interest; scalar or
#' integer vector \code{ny} specifies the number of periods; if \code{length(ny) > 1}
#' then \code{names(ny)} must ne supplied and must match to every unique \code{id}. 
#' logical \code{correctny} specifies whether \code{ny} is 
#' corrected to adjust proportionally for data missingness.
#'
#' Arguments for the sandwich adjustment are given by \code{sandwich.args}. A character
#' string \code{id} can be supplied to the list, which identifies the name of the 
#' variable in \code{data} such that independence will be assumed between its values. The
#' \code{method} for the adjustment is supplied as \code{"magnitude"} (default) or \code{"curvature"};
#' see Chandler & Bate (2007) for their definitions.
#' 
#' @references 
#' 
#' Chandler, R. E., & Bate, S. (2007). Inference for clustered data
#' using the independence loglikelihood. Biometrika, 94(1), 167-183.
#'
#' Oh, H. S., Lee, T. C., & Nychka, D. W. (2011). Fast nonparametric 
#' quantile regression with arbitrary smoothing methods. Journal of 
#' Computational and Graphical Statistics, 20(2), 510-526.
#'
#' Schlather, M., & Tawn, J. A. (2003). A dependence measure for multivariate and 
#' spatial extreme values: Properties and inference. Biometrika, 90(1), 139-156.
#'
#' Wood, S. N., Pya, N., & Safken, B. (2016). Smoothing parameter and model 
#' selection for general smooth models. Journal of the American Statistical 
#' Association, 111(516), 1548-1563.
#'
#' Yu, K., & Moyeed, R. A. (2001). Bayesian quantile regression. 
#' Statistics & Probability Letters, 54(4), 437-447.
#' 
#' @examples
#'
#' library(evgam)
#' 
#' data(COprcp)
#'
#' ## fit generalised Pareto distribution to excesses on 20mm
#'
#' COprcp <- cbind(COprcp, COprcp_meta[COprcp$meta_row,])
#' threshold <- 20
#' COprcp$excess <- COprcp$prcp - threshold
#' COprcp_gpd <- subset(COprcp, excess > 0)
#' fmla_gpd <- list(excess ~ s(lon, lat, k=12) + s(elev, k=5, bs="cr"), ~ 1)
#' m_gpd <- evgam(fmla_gpd, data=COprcp_gpd, family="gpd")
#'
#' \donttest{
#'
#' ## fit generalised extreme value distribution to annual maxima
#'
#' COprcp$year <- format(COprcp$date, "%Y")
#' COprcp_gev <- aggregate(prcp ~ year + meta_row, COprcp, max)
#' COprcp_gev <- cbind(COprcp_gev, COprcp_meta[COprcp_gev$meta_row,])
#' fmla_gev <- list(prcp ~ s(lon, lat, k=30) + s(elev, bs="cr"), ~ s(lon, lat, k=20), ~ 1)
#' m_gev <- evgam(fmla_gev, data=COprcp_gev, family="gev")
#' summary(m_gev)
#' plot(m_gev)
#' predict(m_gev, newdata=COprcp_meta, type="response")
#' 
#' ## fit point process model using r-largest order statistics
#'
#' # we have `ny=30' years' data and use top 45 order statistics
#' pp_args <- list(id="id", ny=30, r=45)
#' m_pp <- evgam(fmla_gev, COprcp, family="pp", pp.args=pp_args)
#'
#' ## estimate 0.98 quantile using asymmetric Laplace distribution
#'
#' fmla_ald <- prcp ~ s(lon, lat, k=15) + s(elev, bs="cr")
#' m_ald <- evgam(fmla_ald, COprcp, family="ald", ald.args=list(tau=.98))
#'
#' }
#'
#' @seealso \link{predict.evgam}
#'
#' @return An object of class \code{evgam}
#' 
#' @export
#' 
evgam <- function(formula, data, family="gev", correctV=TRUE, rho0=0, 
inits=NULL, outer="bfgs", control=NULL, removeData=FALSE, trace=0, 
knots=NULL, maxdata=1e20, maxspline=1e20, compact=FALSE, 
ald.args=list(), exi.args=list(), pp.args=list(), sandwich.args=list()) {

## setup family
family.info <- .setup.family(family)

## setup formulae
formula <- .setup.formulae(formula, family.info$npar, family.info$npar2, data)
response.name <- attr(formula, "response.name")

## setup mgcv objects and data
temp.data <- .setup.data(data, response.name, formula, family, family.info$nms, 
  removeData, exi.args, ald.args, pp.args, knots, maxdata, 
  maxspline, compact, sandwich.args, tolower(outer), trace)
data <- temp.data$data

## initialise inner iteration
beta <- .setup.inner.inits(inits, temp.data$lik.data, family.info$lik.fns, family.info$npar, family)
lik.data <- .sandwich(temp.data$lik.data, beta)
if (trace > 0 & lik.data$adjust > 0) cat(paste("\n Sandwich correct lambda =", signif(lik.data$k, 3), "\n"))

## check whether any smoothing parameters need estimating

smooths <- length(temp.data$gotsmooth) > 0

if (smooths) {

## initialise outer iteration
S.data <- .joinSmooth(lapply(temp.data$gams, function(x) x$smooth))
nsp <- length(attr(S.data, "Sl"))
if (is.null(rho0)) {
    diagSl <- sapply(attr(S.data, "Sl"), diag)
    rho0 <- apply(diagSl, 2, function(y) uniroot(.guess, c(-1e2, 1e2), d=attr(beta, "diagH"), s=y)$root)
} else {
    if (length(rho0) == 1) rho0 <- rep(rho0, nsp)
}

lik.data$S <- .makeS(S.data, exp(rho0))

## perform outer iteration
fit.reml <- .outer(rho0, beta, family.info$lik.fns, lik.data, S.data, control, correctV, lik.data$outer, trace)

sp <- exp(fit.reml$par)
lik.data$S <- .makeS(S.data, sp)

} else {

S.data <- NULL
fit.reml <- .outer.nosmooth(beta, family.info$lik.fns, lik.data, control, trace)

}

## covariance matrices
VpVc <- .VpVc(fit.reml, family.info$lik.fns, lik.data, S.data, correctV=correctV, sandwich=temp.data$sandwich, smooths=smooths)

## effective degrees of freedom
edf <- .edf(fit.reml$beta, family.info$lik.fns, lik.data, VpVc, temp.data$sandwich)

## update mgcv objects
names(temp.data$gams) <- family.info$nms
gams <- .swap(fit.reml, temp.data$gams, lik.data, VpVc, temp.data$gotsmooth, edf, smooths)

## add extra things that make an evgam object
## differ from a list of mgcv objects

gams <- .finalise(gams, data, family.info$lik.fns, lik.data, S.data, fit.reml, VpVc, family, temp.data$gotsmooth, formula, response.name, removeData, edf)

return(gams)
}

#' @rdname evgam
#' @name fevgam
#' @export
NULL
fevgam <- function(...) {
message("`fevgam' will soon be deprecated: please migrate to `evgam'.")
evgam(...)
}

#' Predictions from a fitted \code{evgam} object
#'
#' @param object a fitted \code{evgam} object
#' @param newdata a data frame
#' @param type a character string giving the type of prediction sought; see Details. Defaults to \code{"link"}
#' @param probs a scalar or vector of probabilities for quantiles to be estimated if \code{type == "quantile"}; defaults to 0.5
#' @param se.fit a logical: should estimated standard errors be returned? Defaults to FALSE
#' @param marginal a logical: should uncertainty estimates integrate out smoothing parameter uncertainty? Defaults to TRUE
#' @param ... unused
#'
#' @details
#'
#' There are five options for \code{type}: 1) \code{"link"} distribution parameters 
#' transformed to their model fitting scale; 2) \code{"response"} as 1), but on their 
#' original scale; 3) "lpmatrix" a list of design matrices; 4) "quantile"
#' estimates of distribution quantile(s); and 5) "qqplot" a quantile-quantile
#' plot.
#'
#' @return A data frame or list of predictions, or a plot if \code{type == "qqplot"}
#'
#' @examples
#'
#' \donttest{
#'
#' library(evgam)
#' data(COprcp)
#' COprcp$year <- format(COprcp$date, "%Y")
#' COprcp_gev <- aggregate(prcp ~ year + meta_row, COprcp, max)
#' COprcp_gev <- cbind(COprcp_gev, COprcp_meta[COprcp_gev$meta_row,])
#' fmla_gev <- list(prcp ~ s(lon, lat, k=30) + s(elev, bs="cr"), ~ s(lon, lat, k=20), ~ 1)
#' m_gev <- evgam(fmla_gev, data=COprcp_gev, family="gev")
#' predict(m_gev, COprcp_meta)
#' predict(m_gev, COprcp_meta, type="response")
#' predict(m_gev, COprcp_meta, probs=.99)
#' COprcp_qq1 <- subset(COprcp_gev, name == "BOULDER")
#' predict(m_gev, COprcp_qq1, type="qqplot")
#' COprcp_qq2 <- subset(COprcp_gev, name %in% c("BOULDER", "FT COLLINS"))
#' predict(m_gev, COprcp_qq2, type="qqplot")
#' fitted(m_gev)
#'
#' }
#'
#' @name predict.evgam
#' 
#' @export
#'
predict.evgam <- function(object, newdata=NULL, type="link", probs=NULL, se.fit=FALSE, marginal=TRUE, ...) {
if (!is.null(probs)) type <- "quantile"
if (type == "quantile" & is.null(probs)) stop("non-NULL `probs' required if `type' is `response'")
family <- object$family
ndat <- nrow(object$data)
if (marginal) {
  V.type <- "Vc" 
} else {
  V.type <- "Vp"
}
if (type != "qqplot") {
if (type == "quantile" & !(family %in% c("gev", "gpd", "weibull", "pp")))
  stop(paste("Use of `probs' not implemented for family `", family, "'", sep=""))
gotsmooth <- replace(logical(length(object)), object$gotsmooth, TRUE)
compacted <- object$compacted
if (compacted) compactid <- object$compactid
if (family == "pp") {
    family <- "gev"
    if (is.null(newdata)) {
        message("Predictions for point process model given for quadrature points, not original data frame")
        newdata <- object$data
    }
}
nms <- names(object)[seq_along(object$call)]
if (family == "exi") linkfn <- object$linkfn
conf.pars <- list(object$coefficients, object[[V.type]], object$idpars)
object <- subset(object, sapply(object, class) == "gamlist")
hasformula <- sapply(object, function(x) length(x$coefficients) > 0)
object <- subset(object, hasformula)
hasformula <- hasformula[seq_along(object)]
out <- list()
if (!is.null(newdata)) {
ndat <- nrow(newdata)
if (ndat > 1e5) {
newdata <- split(newdata, ceiling(seq_len(nrow(newdata))/1e5))
} else {
newdata <- list(newdata)
}
}# else ndat <- nrow(object$data)
for (i in seq_along(object)) {
if (is.null(newdata)) {
out[[i]] <- object[[i]]$X
if (compacted) out[[i]] <- out[[i]][compactid, , drop=FALSE]
} else {
if (gotsmooth[[i]]) {
out[[i]] <- lapply(newdata, function(x) matrix(unlist(lapply(object[[i]]$smooth, mgcv::PredictMat, data=x)), nrow(x)))
out[[i]] <- do.call(rbind, out[[i]])
} else {
out[[i]] <- lapply(newdata, function(x) matrix(0, nrow(x), 0))
out[[i]] <- do.call(rbind, out[[i]])
}
pfmla <- as.formula(as.character(object[[i]]$pterms)[-2])
out[[i]] <- cbind(do.call(rbind, lapply(newdata, function(x) model.matrix(pfmla, x))), out[[i]])
}
}
if (type != "lpmatrix") {
if (se.fit & type != "quantile") {
std.err <- as.data.frame(lapply(seq_along(object), function(i) sqrt(rowSums(out[[i]] * (out[[i]] %*% object[[i]][[V.type]])))))
}
  outX <- out
  out <- lapply(seq_along(out), function(i) out[[i]] %*% object[[i]]$coefficients)
  out <- as.data.frame(lapply(out, function(x) x[,1]))
}
if (type %in% c("response", "quantile")) {
unlink <- which(substr(nms, 1, 3) == "log")
for (i in unlink) {
  out[,i] <- exp(out[,i])
  if (se.fit & type == "response") {
    std.err[,i] <- out[,i] * std.err[,i]
  }
}
if (family == "exi") out[,1] <- 2 * linkfn(out[,1]) - 1
nms <- gsub("log", "", nms)
}
names(out) <- nms
if (se.fit & type != "quantile") names(std.err) <- nms
if (type == "quantile") {
if (se.fit) out0 <- out
qnms <- probs
probs <- matrix(probs, nrow=ndat, ncol=length(probs), byrow=TRUE)
pars <- matrix(NA, nrow=ndat, ncol=length(nms), byrow=TRUE)
pars[,hasformula] <- sapply(out, c)
if (family == "gpd") {
out <- apply(probs, 2, function(x) .qgpd(x, 0, pars[,1], pars[,2]))
} else {
if (family == "gev") {
out <- apply(probs, 2, function(x) .qgev(x, pars[,1], pars[,2], pars[,3]))
} else {
if (family == "weibull") {
out <- apply(probs, 2, function(x) .qweibull(x, scale=pars[,1], shape=pars[,2]))
} else {
stop("invalid family")
}}}
if (se.fit) {
ny <- nrow(outX[[1]])
np <- length(outX)
Sigma <- array(NA, dim=c(ny, np, np))
idp <- conf.pars[[3]]
for (i in seq_len(ny)) {
  for (j in seq_len(np)) {
    for (k in j:np) {
      xj <- outX[[j]][i,]
      xk <- outX[[k]][i,]
      V <- conf.pars[[2]][idp == j, idp == k, drop=FALSE]
      Sigma[i, j, k] <- sum(xj * (V %*% xk))
      if (k != j) Sigma[i, k, j] <- Sigma[i, j, k]
    }
  }
}
std.err <- matrix(NA, ny, length(qnms))
for (j in seq_along(qnms)) {
  if (family == "gev") {
  jac <- .dqgev(qnms[j], out0[,1], log(out0[,2]), out0[,3])
  }
  if (family == "ggpd") {
  jac <- .dqgpd(qnms[j], log(out0[,2]), out0[,3])
  }
  if (family == "weibull") {
  jac <- .dqweibull(qnms[j], log(out0[,2]), out0[,3])
  }
  for (i in seq_len(ny)) {
    std.err[i, j] <- sum(jac[i,] * (Sigma[i,,] %*% jac[i,]))
  }
}
std.err <- as.data.frame(sqrt(std.err))
names(std.err) <- paste("q", round(qnms, 3), sep=":")
}
out <- as.data.frame(out)
names(out) <- paste("q", round(qnms, 3), sep=":")
}
if (se.fit) out <- list(fitted = out, se.fit = std.err)
return(out)
} else {
if (is.null(newdata))
  newdata <- object$data
pars <- predict(object, newdata, type="response")
pit <- !all(apply(pars, 2, function(x) all(diff(x) < 1e-12)))
x <- ppoints(nrow(newdata))
y <- newdata[,object$response.name]
if (is.null(y))
  stop("No response in `newdata'")
if (!pit) {
if (!(family %in% c("gev", "gpd", "weibull")))
  stop("Unsuitable `family' for `type == 'qqplot''")
if (family == "gev")
  x <- .qgev(x, pars[,1], pars[,2], pars[,3])
if (family == "gpd")
  x <- .qgpd(x, 0, pars[,1], pars[,2])
if (family == "weibull") 
  x <- .qweibull(x, pars[,1], pars[,2])
} else {
message("Margins converted to unit exponential by probability integral transformation.")
x <- qexp(x)
if (family == "gev")
  y <- .pgev(y, pars[,1], pars[,2], pars[,3])
if (family == "gpd")
  y <- .pgpd(y, 0, pars[,1], pars[,2])
if (family == "weibull") 
  y <- .pweibull(y, pars[,1], pars[,2])
y <- qexp(y)
}
qqplot(x, y)
abline(0, 1)
}
}

#' @rdname predict.evgam
#' 
#' @export
#' 
fitted.evgam <- function(object, ...) {
object <- subset(object, sapply(object, class) == "gamlist")
out <- lapply(object, function(x) x$fitted)
names(out) <- names(object)
return(out)
}

#' Simulations from a fitted \code{evgam} object
#'
#' @param object a fitted \code{evgam} object
#' @param nsim an integer giving the number of simulations
#' @param seed an integer giving the seed for simulations
#' @param newdata a data frame
#' @param type a character string, as in \code{predict.evgam}; defaults to \code{"quantile"}
#' @param probs a scalar or vector of probabilities for quantiles; defaults to NULL
#' @param threshold a scalar, vector or matrix, which is added to each simulation if \code{family == "gpd"}; defaults to 0
#' @param marginal a logical: should simulations integrate out smoothing parameter uncertainty? Defaults to TRUE
#' @param ... arguments to be passed to \code{predict.evgam}
#'
#' @return Simulations of parameters or quantiles
#'
#' @seealso \link{predict.evgam}
#'
#' @examples
#'
#' \donttest{
#'
#' library(evgam)
#' data(COprcp)
#' COprcp$year <- format(COprcp$date, "%Y")
#' COprcp_gev <- aggregate(prcp ~ year + meta_row, COprcp, max)
#' COprcp_gev <- cbind(COprcp_gev, COprcp_meta[COprcp_gev$meta_row,])
#' fmla_gev <- list(prcp ~ s(lon, lat, k=30) + s(elev, bs="cr"), ~ s(lon, lat, k=20), ~ 1)
#' m_gev <- evgam(fmla_gev, data=COprcp_gev, family="gev")
#' simulate(m_gev)
#' simulate(m_gev, probs=c(.95, .99))
#'
#' }
#'
#' @export
#' 
simulate.evgam <- function(object, nsim=1e3, seed=NULL, newdata=NULL, 
type="link", probs=NULL, threshold=0, marginal=TRUE, ...) {
if (!is.null(probs)) type <- "quantile"
if (is.null(newdata)) newdata <- object$data
if (type %in% c("link", "response")) {
if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)                     # initialize the RNG if necessary
if(is.null(seed)) {
  RNGstate <- get(".Random.seed", envir = .GlobalEnv)
} else {
  R.seed <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(seed)
  RNGstate <- structure(seed, kind = as.list(RNGkind()))
  on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
}
family <- object$family
if (marginal) {
  V.type <- "Vc" 
} else {
  V.type <- "Vp"
}
B <- .pivchol_rmvn(nsim, object$coefficients, object[[V.type]])
idpars <- object$idpars
X <- predict.evgam(object, newdata, type="lpmatrix")
nms <- names(X)
B <- lapply(seq_along(X), function(i) B[idpars == i, , drop=FALSE])
X <- lapply(seq_along(X), function(i) X[[i]] %*% B[[i]])
names(X) <- nms
if (type %in% c("response", "quantile")) {
unlink <- which(substr(names(X), 1, 3) == "log")
nms <- gsub("log", "", nms)
for (i in unlink) X[[i]] <- exp(X[[i]])
names(X) <- nms
}
}
if (type == "quantile") {
X <- simulate.evgam(object, nsim=nsim, seed=seed, newdata=newdata, type="response")
out <- list()
for (i in seq_along(probs)) {
if (object$family == "gpd") {
out[[i]] <- threshold + .qgpd(probs[i], 0, X[[1]], X[[2]])
} else {
out[[i]] <- .qgev(probs[i], X[[1]], X[[2]], X[[3]])
}
}
names(out) <- paste("q", probs, sep=":")
if (length(probs) == 1) {
X <- out[[1]]
} else {
if (nrow(newdata) == 1) {
X <- t(sapply(out, c))
} else {
X <- out
}
}
}
return(X)
}

#' Log-likelihood, AIC and BIC from a fitted \code{evgam} object
#'
#' @param object a fitted \code{evgam} object
#' @param ... not used
#' @param k numeric, the penalty per parameter to be used; the default \code{k = 2} is the classical AIC
#'
#' @return A scalar
#'
#' @examples
#'
#' \donttest{
#'
#' library(evgam)
#' data(COprcp)
#' COprcp$year <- format(COprcp$date, "%Y")
#' COprcp_gev <- aggregate(prcp ~ year + meta_row, COprcp, max)
#' COprcp_gev <- cbind(COprcp_gev, COprcp_meta[COprcp_gev$meta_row,])
#' fmla_gev <- list(prcp ~ s(lon, lat, k=30) + s(elev, bs="cr"), ~ s(lon, lat, k=20), ~ 1)
#' m_gev <- evgam(fmla_gev, data=COprcp_gev, family="gev")
#' logLik(m_gev)
#' AIC(m_gev)
#' BIC(m_gev)
#'
#' }
#'
#' @name logLik
#'
#' @export
#' 
logLik.evgam <- function(object, ...) {
if (!missing(...)) warning("extra arguments discarded")
object$logLik
}

#' @rdname logLik
#' @export
#' 
AIC.evgam <- function(object, ..., k = 2) {
if (!missing(...)) warning("extra arguments discarded")
-2 * object$logLik + k * attr(object, "df")
}

#' @rdname logLik
#' @export
#' 
BIC.evgam <- function(object, ...) {
if (!missing(...)) warning("extra arguments discarded")
-2 * object$logLik + log(nobs(object)) * attr(object, "df")
}

#' Plot a fitted \code{evgam} object
#'
#' @param x a fitted \code{evgam} object
#' @param given.vals a list specifying variables values that are fixed.
#' @param add.map logical: should a map outline be added to any two-dimensional smooths using \link[maps]{map}? Defaults to \code{FALSE}
#' @param use.image logical: should \link[graphics]{image} be used to represent two-dimensional smooths, as opposed to \link[graphics]{contour}? Defaults to \code{FALSE}
#' @param map.env a character string identifying the map to superimpose via \link[maps]{map}; defaults to \code{"world"}
#' @param ... unused
#'
#' @return Plots representing all one- or two-dimensional smooths
#'
#' @examples
#'
#' \donttest{
#'
#' library(evgam)
#' data(COprcp)
#' COprcp$year <- format(COprcp$date, "%Y")
#' COprcp_gev <- aggregate(prcp ~ year + meta_row, COprcp, max)
#' COprcp_gev <- cbind(COprcp_gev, COprcp_meta[COprcp_gev$meta_row,])
#' fmla_gev <- list(prcp ~ s(lon, lat, k=30) + s(elev, bs="cr"), ~ s(lon, lat, k=20), ~ 1)
#' m_gev <- evgam(fmla_gev, data=COprcp_gev, family="gev")
#' plot(m_gev)
#'
#' }
#'
#' @export
#' 
plot.evgam <- function(x, given.vals=NULL, use.image=FALSE, add.map=FALSE, map.env="world", ...) {
plotdata <- x$plotdata
nms <- names(x)[x$gotsmooth]
x <- lapply(x$gotsmooth, function(i) x[[i]])
names(x) <- nms
nplot <- sapply(x, function(y) length(y$smooth))
willplot <- nplot > 0
x <- subset(x, willplot)
nplot <- nplot[willplot]
prefixes <- names(nplot)
if (prod(par("mfcol")) < sum(nplot)) par(mfrow=rev(n2mfrow(sum(nplot))))
for (i in seq_along(x)) for (j in seq_along(x[[i]]$smooth)) {
smthij <- x[[i]]$smooth[[j]]
dij <- subset(plotdata, sapply(lapply(plotdata, names), function(x) all(smthij$vn %in% x)))[[1]]
if (ncol(dij) > 2) {
if (is.null(given.vals)) {
given.vals <- lapply(3:ncol(dij), function(j) quantile(dij[,j], .5))
names(given.vals) <- names(dij)[3:ncol(dij)]
} else {
if (inherits(given.vals, "character")) {
given.vals <- as.vector(na.omit(match(names(dij), given.vals)))
}
if (inherits(given.vals, "integer")) {
id.given <- given.vals
given.vals <- lapply(id.given, function(i) quantile(dij[,j], .5))
names(given.vals) <- names(dij)[id.given]
} else {
given.vals <- given.vals
}}
}
.plotSmooth(x[[i]], j, dij, prefixes[i], given.vals, add.map=add.map, use.image=use.image, map.env=map.env)
}
}

#' Summary method for a fitted \code{evgam} object
#'
#' @param object a fitted \code{evgam} object
#' @param ... not used
#'
#' @details
#' 
#' The key part of summary.evgam is p-values for smooths.
#' The tests use code directly taken from \code{mgcv 1.8-14}. This is 
#' to avoid use of \code{mgcv:::...} . Tests implement the method of
#' Wood (2013).
#'
#' @references
#'
#' Wood, S. N., (2013) On p-values for smooth components of an extended
#' generalized additive model, Biometrika 100(1) 221--228
#'
#' @return A \code{summary.evgam} object
#'
#' @examples
#'
#' \donttest{
#'
#' library(evgam)
#' data(COprcp)
#' COprcp$year <- format(COprcp$date, "%Y")
#' COprcp_gev <- aggregate(prcp ~ year + meta_row, COprcp, max)
#' COprcp_gev <- cbind(COprcp_gev, COprcp_meta[COprcp_gev$meta_row,])
#' fmla_gev <- list(prcp ~ s(lon, lat, k=30) + s(elev, bs="cr"), ~ s(lon, lat, k=20), ~ 1)
#' m_gev <- evgam(fmla_gev, data=COprcp_gev, family="gev")
#' summary(m_gev)
#'
#' }
#'
#' @name summary.evgam
#'
#' @export
#' 
summary.evgam <- function(object, ...) {
if (!missing(...)) warning("extra arguments discarded")
out <- list()
out[[1]] <- .parametric.summary.evgam(object)
out[[2]] <- .smooth.summary.evgam(object)
class(out) <- "summary.evgam"
out
}

#' @param x a \code{summary.evgam} object
#'
#' @rdname summary.evgam
#' 
#' @export
#' 
print.summary.evgam <- function(x, ...) {
if (!missing(...)) warning("extra arguments discarded")
cat("\n")
cat("** Parametric terms **")
tab <- lapply(x[[1]], .tidyParametricTable)
cat("\n")
for (i in seq_along(tab)) {
cat("\n")
cat(names(tab)[[i]])
cat("\n")
print(tab[[i]])
}
cat("\n")
cat("** Smooth terms **")
tab <- lapply(x[[2]], .tidySmoothTable)
cat("\n")
for (i in seq_along(tab)) {
cat("\n")
cat(names(tab)[[i]])
cat("\n")
print(tab[[i]])
}
invisible(x)
}

#' Print a fitted \code{evgam} object
#'
#' @param x a fitted \code{evgam} object
#' @param ... not used
#'
#' @return The call of the \code{evgam} object
#'
#' @examples
#' \donttest{
#'
#' library(evgam)
#' data(COprcp)
#' COprcp$year <- format(COprcp$date, "%Y")
#' COprcp_gev <- aggregate(prcp ~ year + meta_row, COprcp, max)
#' COprcp_gev <- cbind(COprcp_gev, COprcp_meta[COprcp_gev$meta_row,])
#' fmla_gev <- list(prcp ~ s(lon, lat, k=30) + s(elev, bs="cr"), ~ s(lon, lat, k=20), ~ 1)
#' m_gev <- evgam(fmla_gev, data=COprcp_gev, family="gev")
#' print(m_gev)
#'
#' }
#'
#' @export
#' 
print.evgam <- function(x, ...) {
if (!missing(...)) warning("extra arguments discarded")
print(x$call)
invisible(x)
}


#' Bind a list a data frames
#'
#' @param x a list of data frames
#'
#' @return A data frame
#'
#' @examples
#' z <- list(data.frame(x=1, y=1), data.frame(x=2, y=2))
#' dfbind(z)
#'
#' @seealso \link[base]{rbind}
#'
#' @export
#' 
dfbind <- function(x) {
nms <- names(x[[1]])
cls <- sapply(x[[1]], class)
x <- lapply(nms, function(i) unlist(lapply(x, function(y) y[,i])))
x <- as.data.frame(x)
dt <- cls == "Date"
if (any(dt)) {
  for (i in which(dt)) x[,i] <- as.Date(x[,i], origin="1970-01-01")
}
names(x) <- nms
x
}

#' Scatter plot, with variable-based point colours
#'
#' @param x a vector of x coordinates
#' @param y a vector of y coordinates
#' @param z a variable for defining colours
#' @param n an integer giving the number of colour levels, supplied to \link[base]{pretty}
#' @param rev logical: should the palette be reversed? Defaults to \code{TRUE}
#' @param cex a scalar for character expansion, supplied to \link[graphics]{plot}
#' @param pch an integer giving the plotting character, supplied to \link[graphics]{plot}
#' @param add logical: add to an existing plot? Defaults to \code{FALSE}
#' @param breaks a vector or breaks for defining color intervals; defaults to \code{NULL}, so \link[base]{pretty} and \code{n} are used on \code{z}
#'
#' @return A plot
#'
#' @examples
#'
#' x <- runif(50)
#' y <- runif(50)
#' colplot(x, y, x * y)
#'
#' @export
#' 
colplot <- function(x, y, z, n=20, rev=TRUE, cex=1, pch=21, add=FALSE, breaks=NULL) {
brks <- pretty(z, n)
pal <- heat.colors(length(brks[-1]))
if (rev) pal <- rev(pal)
col <- pal[as.integer(cut(z, brks))]
if (!add) {
if (pch %in% 21:25) {
  plot(x, y, bg=col, pch=pch, cex=cex)
} else {
  plot(x, y, col=col, pch=pch, cex=cex)
}
} else {
if (pch %in% 21:25) {
  points(x, y, bg=col, pch=pch, cex=cex)
} else {
  points(x, y, col=col, pch=pch, cex=cex)
}
}
}

#' Moore-Penrose pseudo-inverse of a matrix
#'
#' @param x a matrix
#' @param tol a scalar
#' 
#' @details
#' 
#' This function is merely a wrapper for Armadillo's pinv function with its
#' default settings, which, in particular uses the divide-and-conquer
#' method. If \code{tol} isn't provided Armadillo's default for pinv is used.

#' \code{ginv.evgam} mimics \link[MASS]{ginv} using Armadillo's pinv.
#'
#' @return A matrix
#' 
#' @references
#'
#' http://arma.sourceforge.net/docs.html#pinv
#'
#' @seealso \link[MASS]{ginv}
#'
#' @export
#' 
pinv <- function(x, tol=-1) {
armapinv(x, tol)
}

#' @rdname pinv
#'
#' @export
#' 
ginv.evgam <- function(x, tol=sqrt(.Machine$double.eps)) {
armaginv(x, tol)
}

#' More Sequence Generation
#'
#' Generate a sequence of values between a range.
#'
#' @param x a 2-vector
#' @param length an integer
#'
#' @return A vector
#'
#' @seealso \link[base]{seq}, \link[base]{seq_len}, \link[base]{seq_along}
#'
#' @examples
#'
#' seq_between(c(1, 9))
#' seq_between(range(runif(10)), 5)
#'
#' @export
#' 
seq_between <- function(x, length=NULL) {
if (is.null(length)) {
    return(seq(x[1], x[2]))
    } else {
    return(seq(x[1], x[2], length=length))
    }
}
