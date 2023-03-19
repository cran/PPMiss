# objective function
.obj <- function(d, rhos, a = 0, k1, s = 1, m = 24, r = 2){
  h = s:m
  if(length(rhos) != length(h)) stop("rhos must have length(s:m)")
  out <- k1*(rhos - a) - (gamma(1-d)/gamma(d))*(s:m)^(2*d-1)
  return(sum(abs(out)^r))
}

.d_hat <-
  function(xt, thetas, empirical = TRUE,
           s = 1, m = 24, r = 2, a = 0,
           dCdtheta = dCtheta_gauss,
           interval = c(-0.5, 0.5)){

    # find the position of missing values
    w <- is.na(xt)
    x <- xt[!w]

    # density
    fun <- NULL
    if(empirical) fun <- kdens(x)

    k1 <- k1fun(dCdtheta = dCdtheta, fun = fun,
                data = x, empirical = empirical,
                mean = mean(x), sd = stats::sd(x))

    dhat <- stats::optimize(.obj,
                     interval = interval, a = a,
                     k1 = k1, rhos = thetas,
                     s = s, m = m, r = r)$minimum

    return(dhat)
  }

