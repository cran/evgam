## Generalised Pareto negative log-likelihood functions

.gpd.d0 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) gpdd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.gpd.d12 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) gpdd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.gpd.d34 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) gpdd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.gpdfns <- list(d0=.gpd.d0, d120=.gpd.d12, d340=.gpd.d34)

.qgpd <- function(p, loc, scale, shape, zeta=1, theta=1, m=1) {
shape <- sign(shape) * pmax(abs(shape), 1e-6)
out <- 1 - p ^ (1 / (m * theta))
loc + scale * ((out / zeta)^(-shape) - 1) / shape
}

.pgpd <- function(x, loc, scale, shape, tau, NAOK=FALSE, log=FALSE) {
# function to evaluate daily cdf, e.g. Coles (2001, pp.138)
below <- x < loc
x <- pmax(x, loc)
temp <- 1 + shape * (x - loc) / scale
if (!NAOK) temp <- pmax(temp, 0)
out <- temp ^ (-1/shape)
out <- 1 - (1 - tau) * out
if (log) out <- log(out)
if (NAOK) out[below] <- NA
out
}

.dqgpd <- function(p, lscale, shape) {
shape <- sign(shape) * pmax(abs(shape), 1e-6)
.e1 <- 1 - p
.e2 <- .e1^shape
.e4 <- 1/.e2 - 1
.e5 <- exp(lscale)
d1 <- .e4 * .e5/shape
d2 <- -((.e4/shape + log(.e1)/.e2) * .e5/shape)
cbind(d1, d2)
}
