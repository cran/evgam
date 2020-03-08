## Generalised Pareto negative log-likelihood functions

.pp.d0 <- function(pars, likdata) {
ny <- ncol(likdata$y)
prs <- split(pars, likdata$idpars)
X <- lapply(likdata$X, function(x) x[likdata$ppq, , drop=FALSE])
y <- likdata$y[likdata$ppq, , drop=FALSE]
out1 <- lapply(seq_len(ny), function(i) pp1d0(prs, X[[1]], X[[2]], X[[3]], y[, i], likdata$ppw))
out1 <- Reduce("+", out1) / ny
X <- lapply(likdata$X, function(x) x[!likdata$ppq, , drop=FALSE])
y <- likdata$y[!likdata$ppq, , drop=FALSE]
if (is.null(likdata$cens)) {
  out2 <- lapply(seq_len(ny), function(i) pp2d0(prs, X[[1]], X[[2]], X[[3]], y[,i], likdata$weights))
  out2 <- Reduce("+", out2) / ny
} else {
  wts <- likdata$weights[!likdata$cens]
  out21 <- lapply(seq_len(ny), function(i) pp2d0(prs, X[[1]][!likdata$cens, , drop=FALSE], X[[2]][!likdata$cens, , drop=FALSE], X[[3]][!likdata$cens, , drop=FALSE], y[!likdata$cens, i], wts))
  out21 <- Reduce("+", out21) / ny
  wts <- likdata$weights[likdata$cens]
  out22 <- lapply(seq_len(ny), function(i) ppcd0(prs, X[[1]][likdata$cens, , drop=FALSE], X[[2]][likdata$cens, , drop=FALSE], X[[3]][likdata$cens, , drop=FALSE], y[likdata$cens, i], wts))
  out22 <- Reduce("+", out22) / ny
  out2 <- out21 + out22
}
out <- out1 + out2
if (!is.finite(out)) out <- 1e20
out
}

.pp.d12 <- function(pars, likdata) {
ny <- ncol(likdata$y)
prs <- split(pars, likdata$idpars)
X <- lapply(likdata$X, function(x) x[likdata$ppq, , drop=FALSE])
y <- likdata$y[likdata$ppq, , drop=FALSE]
out1 <- lapply(seq_len(ny), function(i) pp1d12(prs, X[[1]], X[[2]], X[[3]], y[, i], likdata$ppw))
out1 <- Reduce("+", out1) / ny
X <- lapply(likdata$X, function(x) x[!likdata$ppq, , drop=FALSE])
y <- likdata$y[!likdata$ppq, , drop=FALSE]
if (is.null(likdata$cens)) {
  out2 <- lapply(seq_len(ny), function(i) pp2d12(prs, X[[1]], X[[2]], X[[3]], y[,i]))
  out2 <- Reduce("+", out2) / ny
} else {
  out21 <- lapply(seq_len(ny), function(i) pp2d12(prs, X[[1]][!likdata$cens, , drop=FALSE], X[[2]][!likdata$cens, , drop=FALSE], X[[3]][!likdata$cens, , drop=FALSE], y[!likdata$cens, i]))
  out21 <- Reduce("+", out21) / ny
  out22 <- lapply(seq_len(ny), function(i) ppcd12(prs, X[[1]][likdata$cens, , drop=FALSE], X[[2]][likdata$cens, , drop=FALSE], X[[3]][likdata$cens, , drop=FALSE], y[likdata$cens, i]))
  out22 <- Reduce("+", out22) / ny
  out2 <- matrix(NA, length(likdata$cens), ncol(out21))
  out2[!likdata$cens,] <- out21
  out2[likdata$cens,] <- out22
}
out2 <- likdata$weights * out2
rbind(out1, out2)
}

.pp.d34 <- function(pars, likdata) {
ny <- ncol(likdata$y)
prs <- split(pars, likdata$idpars)
X <- lapply(likdata$X, function(x) x[likdata$ppq, , drop=FALSE])
y <- likdata$y[likdata$ppq, , drop=FALSE]
out1 <- lapply(seq_len(ny), function(i) pp1d34(prs, X[[1]], X[[2]], X[[3]], y[, i], likdata$ppw))
out1 <- Reduce("+", out1) / ny
X <- lapply(likdata$X, function(x) x[!likdata$ppq, , drop=FALSE])
y <- likdata$y[!likdata$ppq, , drop=FALSE]
if (is.null(likdata$cens)) {
  out2 <- lapply(seq_len(ny), function(i) pp2d34(prs, X[[1]], X[[2]], X[[3]], y[,i]))
  out2 <- Reduce("+", out2) / ny
} else {
  out21 <- lapply(seq_len(ny), function(i) pp2d34(prs, X[[1]][!likdata$cens, , drop=FALSE], X[[2]][!likdata$cens, , drop=FALSE], X[[3]][!likdata$cens, , drop=FALSE], y[!likdata$cens, i]))
  out21 <- Reduce("+", out21) / ny
  out22 <- lapply(seq_len(ny), function(i) ppcd34(prs, X[[1]][likdata$cens, , drop=FALSE], X[[2]][likdata$cens, , drop=FALSE], X[[3]][likdata$cens, , drop=FALSE], y[likdata$cens, i]))
  out22 <- Reduce("+", out22) / ny
  out2 <- matrix(NA, length(likdata$cens), ncol(out21))
  out2[!likdata$cens,] <- out21
  out2[likdata$cens,] <- out22
}
out2 <- likdata$weights * out2
rbind(out1, out2)
}

.ppfns <- list(d0=.pp.d0, d120=.pp.d12, d340=.pp.d34)
