#' @title Long memory parameter estimation
#'
#' @description
#' Let \eqn{\theta_h} be the copula parameter associated to
#' \eqn{(X_t,X_{t+h})} and \eqn{\hat\theta_h} be an estimate of \eqn{\theta_h}
#' based on pseudo observations. The long memory parameter
#' \eqn{d} is estimated by
#'
#' \deqn{\hat d:=\underset{|d|<0.5}{\mathrm{argmin}}\bigg\{\sum_{h=s}^m \bigg|{\hat K_1}(\hat\theta_h-a)-\frac{\Gamma(1-d)}{\Gamma(d)}h^{2d-1}\bigg|^r\bigg\}, \quad r > 0}.
#'
#' @param xt a vector with the observed time series. Missing observations are allowed.
#'
#' @param copula an object of class \sQuote{copula}. Readily available options
#' are \code{frank}, \code{amh}, \code{fgm} and  \code{gauss}. Other
#' copulas can be used but the user must provide the corresponding \code{dCdtheta}.
#' Default is \code{gauss}.
#'
#' @param dCdtheta a two parameter function that returns the limit of the copula
#' derivative, with respect to \eqn{\theta}, as \eqn{\theta} goes to \eqn{a}, where \eqn{a} is
#' such that \eqn{C_a(u,v)=uv}. Readily available
#' options are \code{dCtheta_frank}, \code{dCtheta_amh}, \code{dCtheta_fgm}
#' and  \code{dCtheta_gauss}. Default is  \code{dCtheta_gauss}.
#'
#' @param theta.lower the lower bound for \eqn{\theta}. Default is -1.
#'
#' @param theta.upper the upper bound for \eqn{\theta}. Default is 1.
#'
#' @param optim.method a character string specifying the optimization method.
#' For all available options see \code{\link[stats]{optim}}.
#' Default is \sQuote{Brent}. See \cite{\link[copula]{fitCopula}} for
#' more details.
#'
#' @param method a character string specifying the copula parameter
#' estimator used. This can be one of: \sQuote{mpl}, \sQuote{itau}, \sQuote{irho},
#' \sQuote{itau.mpl} or \sQuote{ml}. See \cite{\link[copula]{fitCopula}} for details.
#' Default is \sQuote{mpl}.
#'
#' @param s integer. The smallest lag \eqn{h} considered in the estimation. Default is 1.
#'
#' @param m integer. The  largest lag \eqn{h} considered in the estimation. Default is 24.
#'
#' @param theta.start starting value for the parameter optimization via \code{\link[stats]{optim}}.
#'
#' @param empirical logical. If \code{TRUE}, the sample estimators for the density
#' and quantile functions are considered. Otherwhise, the gaussian density and
#' quantile functions are used. Default is \code{TRUE}
#'
#' @param r the exponent used in the Minkowski distance used to calculate \eqn{\hat d}.
#' Default is 2, the Euclidean distance.
#'
#' @param a the value of \eqn{\theta} such that \eqn{\lim_{\theta \to a} C_\theta(u,v)=uv},
#' is the product (or independence) copula. Default is 0, which is the common value for
#' the available copulas, namely, \code{frank}, \code{amh}, \code{fgm} and  \code{gauss}.
#'
#' @param d.interval a vector of size 2 giving the lower and upper bound for the
#' long memory parameter \eqn{d}. Default is \code{c(-0.5,0.5)}.
#'
#' @return \eqn{\hat d}, the estimated value of \eqn{d}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' #-------------------------
#' # ARFIMA(0,d,0) process
#' #-------------------------
#' trunc <- 50000
#' n = 1000
#' cks <- arfima.coefs(d = 0.25, trunc = trunc)
#' eps <- rnorm(trunc+n)
#' x <- sapply(1:n, function(t) sum(cks*rev(eps[t:(t+trunc)])))
#'
#' #----------------------
#' # Original time series
#' #-----------------------
#' # For Frank copula, -Inf < theta < Inf. However, "Brent" requires
#' # finite lower and upper bounds so we use c(-100, 100) here
#' d_frank <- d.fit(xt = x, copula = frank, dCdtheta = dCtheta_frank,
#'                  theta.lower = -100, theta.upper = 100)
#' d_amh <- d.fit(xt = x, copula = amh, dCdtheta = dCtheta_amh,
#'                  theta.lower = -1, theta.upper = 1)
#' d_fgm <- d.fit(xt = x, copula = fgm, dCdtheta = dCtheta_fgm,
#'                  theta.lower = -1, theta.upper = 1)
#' d_gauss <- d.fit(xt = x, copula = gauss, dCdtheta = dCtheta_gauss,
#'                  theta.lower = -1, theta.upper = 1)
#'
#' c(FRANK = d_frank, AMH = d_amh, FGM = d_fgm, GAUSS = d_gauss)
#'
#' #----------------------------
#' # Adding some missing values
#' #----------------------------
#' index <- sample(1:n, size = round(0.2*n))
#' xt <- x
#' xt[index] <- NA
#'
#' d_frank_m <- d.fit(xt = xt, copula = frank,
#'                    dCdtheta = dCtheta_frank,
#'                    theta.lower = -100, theta.upper = 100)
#' d_amh_m <- d.fit(xt = xt, copula = amh, dCdtheta = dCtheta_amh,
#'                  theta.lower = -1, theta.upper = 1)
#' d_fgm_m <- d.fit(xt = xt, copula = fgm, dCdtheta = dCtheta_fgm,
#'                  theta.lower = -1, theta.upper = 1)
#' d_gauss_m <- d.fit(xt = xt, copula = gauss,
#'                    dCdtheta = dCtheta_gauss,
#'                    theta.lower = -1, theta.upper = 1)
#'
#' data.frame(
#'   series = c("Complete", "Missing"),
#'   FRANK = c(d_frank, d_frank_m),
#'   AMH = c(d_amh, d_amh_m),
#'   FGM = c(d_fgm, d_fgm_m),
#'   GAUSS = c(d_gauss, d_gauss_m))
#'
#' #-------------------------
#' # ARFIMA(1,d,1) process
#' #-------------------------
#' # For a faster algorithm to generate ARFIMA processes,
#' # see the package "arfima"
#' trunc <- 50000
#' cks <- arfima.coefs(d = 0.35, trunc = trunc, ar = -0.2, ma = 0.4)
#' n = 1000
#' eps <- rnorm(trunc+n)
#' x <- sapply(1:n, function(t) sum(cks*rev(eps[t:(t+trunc)])))
#'
#' #----------------------
#' # Original time series
#' #-----------------------
#' # For Frank copula, -Inf < theta < Inf. However, "Brent" requires
#' # finite lower and upper bounds so we use c(-100, 100) here
#' d_frank <- d.fit(xt = x, copula = frank, dCdtheta = dCtheta_frank,
#'                  theta.lower = -100, theta.upper = 100)
#' d_amh <- d.fit(xt = x, copula = amh, dCdtheta = dCtheta_amh,
#'                  theta.lower = -1, theta.upper = 1)
#' d_fgm <- d.fit(xt = x, copula = fgm, dCdtheta = dCtheta_fgm,
#'                  theta.lower = -1, theta.upper = 1)
#' d_gauss <- d.fit(xt = x, copula = gauss, dCdtheta = dCtheta_gauss,
#'                  theta.lower = -1, theta.upper = 1)
#'
#' c(FRANK = d_frank, AMH = d_amh, FGM = d_fgm, GAUSS = d_gauss)
#'
#' #----------------------------
#' # Adding some missing values
#' #----------------------------
#' n = 1000
#' index <- sample(1:n, size = round(0.2*n))
#' xt <- x
#' xt[index] <- NA
#'
#' d_frank_m <- d.fit(xt = xt, copula = frank,
#'                    dCdtheta = dCtheta_frank,
#'                    theta.lower = -100, theta.upper = 100)
#' d_amh_m <- d.fit(xt = xt, copula = amh, dCdtheta = dCtheta_amh,
#'                  theta.lower = -1, theta.upper = 1)
#' d_fgm_m <- d.fit(xt = xt, copula = fgm, dCdtheta = dCtheta_fgm,
#'                  theta.lower = -1, theta.upper = 1)
#' d_gauss_m <- d.fit(xt = xt, copula = gauss,
#'                    dCdtheta = dCtheta_gauss,
#'                    theta.lower = -1, theta.upper = 1)
#'
#' data.frame(
#'   series = c("Complete", "Missing"),
#'   FRANK = c(d_frank, d_frank_m),
#'   AMH = c(d_amh, d_amh_m),
#'   FGM = c(d_fgm, d_fgm_m),
#'   GAUSS = c(d_gauss, d_gauss_m))
#'}
#'
#' @md
d.fit <- function(xt,
                  copula = gauss, dCdtheta = dCtheta_gauss,
                  theta.lower = -1, theta.upper = 1,
                  optim.method = "Brent", method = "mpl",
                  s = 1, m = 24, theta.start = 0.1,
                  empirical = TRUE, r = 2, a = 0,
                  d.interval = c(-0.5, 0.5)){

  #------------------------
  #  STEP 1: estimate theta
  #------------------------
  thetas <- .theta_hat(xt = xt, s = s, m = m,
             copula = copula,
             lower = theta.lower, upper = theta.upper,
             optim.method = optim.method, method = method,
             start = theta.start, empirical = empirical)

  #------------------------
  #  STEP 2: estimate d
  #------------------------
  d <- .d_hat(xt = xt, thetas = thetas,
             s = s, m = m, r = r, a = a,
             dCdtheta = dCdtheta,
             interval = d.interval)

  return(d)
}
