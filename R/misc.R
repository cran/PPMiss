#' @title Coefficients of an ARFIMA(p,d,q) model
#'
#' @description
#' This function calculates the coefficients \eqn{c_k, k \geq 0} corresponding to
#' \eqn{\theta(z)\phi^{-1}(z)(1-z)^{-d} = \sum_{k = 0}^{\infty}c_k z^k},
#' up to a truncation point
#'
#' @param ar the coefficients of the autoregressive polinomial. Default is NULL
#' @param ma the coefficients of the moving average polinomial. Default is null
#' @param d the long memory parameter. Default is 0.
#' @param trunc the truncation point. Default is 1.
#'
#' @return The coefficients values up to order \sQuote{trunc}.
#'
#' @export
#'
#' @examples
#' cks <- arfima.coefs(d = 0.3, trunc = 5)
#' cks
#'
#' cks <- arfima.coefs(d = 0.1, trunc = 5, ar = 0.5, ma = 0.6)
#' cks
#'
arfima.coefs <- function(ar = NULL, ma = NULL, d = 0, trunc = 1){
  out = .Fortran('coefs',
                 p = length(ar), phi = -c(-1, ar),
                 q = length(ma), theta = c(1, ma),
                 d = as.numeric(d),
                 m = as.integer(trunc),
                 cks = numeric(trunc + 1))

  return(out$cks)
}


#' @title Kernel density estimator
#'
#' @description
#' The probability density function \eqn{F'} is estimated using a kernel density
#' approach. More specifically, first \eqn{y_i = \hat{f}(x_i^\ast)}  is estimated
#' using \eqn{T = 512} (default for the function \code{\link[stats]{density}})
#' equally spaced points \eqn{x_i^\ast}, \eqn{1 \leq i \leq T}, in the interval
#' \eqn{[x_{(1)} - 3b, x_{(n)} + 3b]}, where \eqn{b} is the bandwidth for
#' the Gaussian kernel density estimator, chosen by applying the Silverman's
#' rule of thumb (the default procedure in \code{\link[stats]{density}}).
#' A cubic spline interpolation (the default method for \code{\link[stats]{spline}})
#' is then applied to the pairs \eqn{\{(x_i^\ast, y_i)\}_{i=1}^T} to obtain
#' \eqn{\hat F_n'(x)} for all \eqn{x \in [x_{(1)} - 3b, x_{(n)} + 3b]}.
#'
#' @param x the data from which the estimate is to be computed.
#'
#' @return a function that approximates the probability density function.
#'
#' @examples
#' # creating a time series
#' trunc = 50000
#' cks <- arfima.coefs(d = 0.25, trunc = trunc)
#' eps <- rnorm(trunc+1000)
#' x <- sapply(1:1000, function(t) sum(cks*rev(eps[t:(t+trunc)])))
#'
#' # kernel density function
#' dfun <- kdens(x)
#'
#' # plot
#' curve(dfun, from = min(x), to = max(x))
#'
#' @export
#'
#' @md
kdens <- function(x){
  f <- stats::density(x)
  return(stats::splinefun(f$x, f$y))
}


# F'
.dens <- function(x, fun, empirical, mean = 0, sd = 1){
  # fun = kernel density estimator
  if(empirical) return(fun(x))
  # gaussian marginal.
  return(stats::dnorm(x, mean = mean, sd = sd))
}

# F^{-1}
.Finv <- function(u, data, empirical, mean = 0, sd = 1){
  # empirical quantile
  if(empirical) return(stats::quantile(x = data, prob = u))
  # gaussian quantile
  return(stats::qnorm(p = u, mean = mean, sd = sd))
}


#' @title Constant K1
#'
#' @description Calculates an estimate for the constant \eqn{K_1} given by
#' \deqn{K_1 = \iint_{I^2}\frac{1}{l_0(u)l_n(v)}\lim_{\theta\to a}\frac{\partial C_{\theta}(u,v)}{\partial\theta}\,dudv, }
#' where \eqn{l_m(x):= F_m'\big(F_m^{(-1)}(x)\big)}, \eqn{a} is such that \eqn{C_a(u,v)=uv} (the product copula), and
#' \eqn{\{F_n\}_{n\in\mathbb{N}}} is a sequence of absolutely continuous distribution
#' functions.
#'
#' @param dCdtheta a function providing the limit as \eqn{\theta \to a} of the
#' copula derivative with respect to \eqn{\theta}. For the readily available copulas, namely, \code{frank}, \code{amh},
#' \code{fgm} and  \code{gauss}, \eqn{a=0}.
#'
#' @param fun optionally, a function providing an estimator for the probability density function.
#'
#' @param data the observed time series. Only used to obtain the quantile function
#' when \code{empirical = TRUE}.
#'
#' @param empirical logical. If \code{TRUE}, the sample estimators for the density
#' and quantile functions are considered. Otherwise, the gaussian density and
#' quantile functions are used instead.
#'
#' @param mean the mean of the gaussian distribution.
#' Only used if \code{empirical = FALSE}
#'
#' @param sd the standard deviation of the gaussian distribution.
#' Only used if \code{empirical = FALSE}
#'
#' @details
#' Here \eqn{F'} and \eqn{F^{(-1)}} are replaced by sample estimators for these
#' functions or the gaussian density and quantile functions are used, depending
#' on the context.
#'
#' The function \code{\link{kdens}} is used as sample estimator of \eqn{F'} and
#' \code{\link[stats]{quantile}} is the sample estimator of \eqn{F^{(-1)}}.
#'
#' @return
#' The value of \eqn{K_1}.
#'
#' @examples
#' trunc = 50000
#' cks <- arfima.coefs(d = 0.25, trunc = trunc)
#' eps <- rnorm(trunc+1000)
#' x <- sapply(1:1000, function(t) sum(cks*rev(eps[t:(t+trunc)])))
#'
#' # kernel density function
#' dfun <- kdens(x)
#'
#' # calculating K1 using four copulas and empirical estimates for F' and F^{(-1)}
#' K1_frank_e <- k1fun(dCdtheta = dCtheta_frank, fun = dfun,
#'                  data = x, empirical = TRUE)
#' K1_amh_e <- k1fun(dCdtheta = dCtheta_amh, fun = dfun,
#'                  data = x, empirical = TRUE)
#' K1_fgm_e <- k1fun(dCdtheta = dCtheta_fgm, fun = dfun,
#'                  data = x, empirical = TRUE)
#' K1_gauss_e <- k1fun(dCdtheta = dCtheta_gauss, fun = dfun,
#'                  data = x, empirical = TRUE)
#'
#' # calculating K1 using four copulas and gaussian marginals
#' K1_frank_g <- k1fun(dCdtheta = dCtheta_frank, fun = NULL, data = NULL,
#'                   empirical = FALSE, mean = mean(x), sd = sd(x))
#' K1_amh_g <- k1fun(dCdtheta = dCtheta_amh, fun = NULL, data = NULL,
#'                   empirical = FALSE, mean = mean(x), sd = sd(x))
#' K1_fgm_g <- k1fun(dCdtheta = dCtheta_fgm, fun = NULL, data = NULL,
#'                   empirical = FALSE, mean = mean(x), sd = sd(x))
#' K1_gauss_g <- k1fun(dCdtheta = dCtheta_gauss, fun = NULL, data = NULL,
#'                   empirical = FALSE, mean = mean(x), sd = sd(x))
#'
#' # comparing results
#'  data.frame(MARGINAL = c("Empirical", "Gaussian"),
#'             FRANK = c(K1_frank_e, K1_frank_g),
#'             AMH = c(K1_amh_e,  K1_amh_g),
#'             FGM = c(K1_fgm_e, K1_fgm_g),
#'             GAUSS = c(K1_gauss_e, K1_gauss_g))
#'
#' @export
#'
#' @md
k1fun <- function(dCdtheta, fun, data, empirical, mean = 0, sd = 1){

  f <- function(u,v){
    x <- .Finv(u, data = data, empirical = empirical, mean = mean, sd = sd)
    y <- .Finv(v, data = data, empirical = empirical, mean = mean, sd = sd)
    d1 <- .dens(x,fun = fun, empirical = empirical, mean = mean, sd = sd)
    d2 <- .dens(y,fun = fun, empirical = empirical, mean = mean, sd = sd)
    deno <- d1*d2

    return(dCdtheta(u,v)/deno)
  }
  pracma::integral2(f, xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q
}




#' @title Copula functions and the corresponding derivative limit.
#'
#' @name PPMiss.copulas
#' @order 1
#'
#' @description Implemented copulas and the corresponding derivative limit,
#' as \eqn{\theta \to a}, where \eqn{a} is such that \eqn{C_a(u,v)=uv}.
#' An estimate for \eqn{\theta} is obtained based on the copula function used.
#' and the derivatives are used to obtain an estimate for \eqn{K_1}.
#' The functions \sQuote{frank}, \sQuote{amh}, \sQuote{fgm} and \sQuote{gauss}
#' are shortcuts for  \code{\link[copula:frankCopula]{copula::frankCopula()}},
#' \code{\link[copula:amhCopula]{copula::amhCopula()}}, \code{\link[copula:fgmCopula]{copula::fgmCopula()}} and
#' \code{\link[copula:normalCopula]{copula::normalCopula()}} from package \sQuote{copula},
#'  respectively.
#'
#' @details The constant \eqn{K_1} is given by
#'
#' \deqn{K_1 = \iint_{I^2}\frac{1}{l_0(u)l_n(v)}\lim_{\theta\rightarrow a}\frac{\partial C_{\theta}(u,v)}{\partial\theta}\,dudv, }
#'
#' where \eqn{I=[0,1]}, \eqn{l_m(x):= F_m'\big(F_m^{(-1)}(x)\big)} and
#' \eqn{\{F_n\}_{n\in\mathbb{N}}} is a sequence of absolutely continuous distribution
#' functions
#'
#' @return Archimedean copula objects of class \sQuote{frankCopula}, \sQuote{amhCopula} or a
#' Farlie-Gumbel-Morgenstern copula object of class \sQuote{fgmCopula} or an elliptical
#' copula object of class \sQuote{normalCopula}. For details, see
#' \code{\link[copula]{archmCopula}}, \code{\link[copula]{fgmCopula}} and
#' \code{\link[copula]{ellipCopula}}.
#'
#' The derivative functions return the limit, as \eqn{\theta \to 0}, of the
#' derivative with respect to \eqn{\theta}, corresponding to the copula functions.
#'
#' @md
NULL
#> NULL


#' @rdname PPMiss.copulas
#' @order 2
#' @export
#' @md
frank <- copula::frankCopula()
#' @rdname PPMiss.copulas
#' @order 3
#' @param u a real number between 0 and 1.
#' @param v a real number between 0 and 1.
#' @export
#' @md
dCtheta_frank <- function(u,v) 1/2*u*v*(1-u)*(1-v)

#' @rdname PPMiss.copulas
#' @order 4
#' @export
#' @md
amh <- copula::amhCopula()
#' @rdname PPMiss.copulas
#' @order 5
#' @export
#' @md
dCtheta_amh <- function(u,v) u*v*(1-u)*(1-v)

#' @rdname PPMiss.copulas
#' @order 6
#' @export
#' @md
fgm <- copula::fgmCopula()
#' @rdname PPMiss.copulas
#' @order 7
#' @export
#' @md
dCtheta_fgm <- function(u,v) u*v*(1-u)*(1-v)

#' @rdname PPMiss.copulas
#' @order 8
#' @export
#' @md
gauss <- copula::normalCopula()
#' @rdname PPMiss.copulas
#' @order 9
#' @export
#' @md
dCtheta_gauss <- function(u,v) 1/(2*pi)*exp(-stats::qnorm(u)^2/2 - stats::qnorm(v)^2/2)

