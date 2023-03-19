#-------------------------------
# estimation of theta
#-------------------------------
.theta_hat <-
  function(xt, s = 1, m = 24,
           copula = gauss, lower = -1, upper = 1,
           optim.method = "Brent",
           method = "mpl",
           start = 0.1,
           empirical = TRUE){

    # find the position of missing values
    w <- is.na(xt)

    # pseudo observations based on non-missing values
    u <- xt
    if(empirical)
      u[!w] <- copula::pobs(xt[!w])
    else
      u[!w] <- stats::pnorm(xt,mean = mean(xt[!w]), sd = stats::sd(xt[!w]))

    n <- length(u)
    theta_hat <- NULL
    #---------------------------
    # loop - lags
    for(h in s:m){

      # creating the pair (u_t, u_{t+h})
      samp <- cbind(u[1:(n-h)], u[(1+h):n])

      # find missing values position
      w <- apply(samp, 1, function(x) any(is.na(x)))
      if(sum(w) == nrow(samp))
        stop("all pairs have missing values. Try a different lag")

      samp <- samp[!w,]
      out <- copula::fitCopula(copula = copula,
                       data = samp,
                       method = method,
                       start = start,
                       estimate.var = FALSE,
                       optim.method = optim.method,
                       lower = lower, upper = upper)

      if(("try-error" %in% class(out))){
        warning("Fail to fit the copula.
           Try a different optim.method or a different method.",
           immediate. = TRUE)
        theta_hat <- c(theta_hat, NA)
      }
      else theta_hat <- c(theta_hat, out@estimate)
    }

    if(length(theta_hat) == 1) return(theta_hat)

    #-----------------
    # fixing NA
    #-----------------
    w <- which(!is.na(theta_hat))
    if(length(w) == 0) return(theta_hat)
    # First value is NA
    if(is.na(theta_hat[1])) theta_hat[1] <- theta_hat[w[1]]
    # Last value is NA
    if(is.na(theta_hat[length(theta_hat)]))
      theta_hat[length(theta_hat)] <- theta_hat[w[length(w)]]
    # Other values are NA
    if(sum(is.na(theta_hat)) > 0){
      theta_hat <- zoo::na.approx(theta_hat)
      theta_hat <-  sapply(theta_hat, function(x) max(min(upper, x), lower))
    }

    return(theta_hat)
  }
