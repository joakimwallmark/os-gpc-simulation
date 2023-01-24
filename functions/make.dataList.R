# Same as package but allows for specification of lambda and 
# specification of nbasis when non-penalized splines are used
function (U, key, optList, scrrng = NULL, titlestr = NULL, nbin = nbinDefault(N), 
          NumBasis = NULL, Wnorder = 5, WfdPar = NULL, jitterwrd = TRUE, PcntMarkers = c(5, 
                                                                            25, 50, 75, 95), quadwrd = FALSE, verbose = FALSE) 
{
  N <- nrow(U)
  n <- ncol(U)
  if (min(U) < 1) 
    stop("Zero data values encountered in data.  Are they score values?")
  if (length(key) != n && length(key) > 0) {
    stop("length of key is neither n or 0")
  }
  noption <- rep(0, n)
  for (item in 1:n) {
    noption[item] <- length(optList$optScr[[item]])
  }
  grbg <- rep(FALSE, n)
  for (item in 1:n) {
    grbgind <- U[, item] > noption[item]
    if (sum(grbgind) > 0) {
      noption[item] <- noption[item] + 1
      optScri <- optList$optScr[[item]]
      optScri <- c(optScri, 0)
      optList$optScr[[item]] <- optScri
      grbg[item] <- TRUE
      U[grbgind, item] <- noption[item]
    }
  }
  Wdim <- sum(noption)
  scrvec <- matrix(0, N, 1)
  itmvec <- matrix(0, n, 1)
  for (item in 1:n) {
    for (j in 1:N) {
      scoreij <- optList$optScr[[item]][U[j, item]]
      scrvec[j] <- scrvec[j] + scoreij
      itmvec[item] <- itmvec[item] + scoreij
    }
  }
  scrmin <- min(scrvec)
  scrmax <- max(scrvec)
  scrrng <- c(scrmin, scrmax)
  if (is.null(scrrng)) 
    scrrng <- c(scrmin, scrmax)
  nfine <- 101
  scrfine <- seq(scrrng[1], scrrng[2], len = nfine)
  thetaQnt <- seq(0, 100, len = 2 * nbin + 1)
  if (jitterwrd) {
    scrjit <- scrvec + rnorm(N) * 0.1
    scrjit[scrjit < scrmin] <- scrmin
    scrjit[scrjit > scrmax] <- scrmax
  }
  else {
    scrjit = scan("scrjit.txt", 0)
  }
  scrrnk <- matrix(0, N, 1)
  for (j in 1:N) scrrnk[j] <- sum(scrjit <= scrjit[j])
  percntrnk <- 100 * scrrnk/N
  if (is.null(NumBasis)) {
    NumBasis = NumBasisDefault(N)
  }
  if (is.null(WfdPar)) {
    # if (quadwrd) {
    Wnbasis <- NumBasis
    Wbasis <- fda::create.bspline.basis(c(0, 100), Wnbasis, 
                                        Wnorder)
    WfdPar <- fdPar(Wbasis)
    # }
    # else {
    #   Wnbasis <- NumBasis
    #   Wbasis <- fda::create.bspline.basis(c(0, 100), Wnbasis, 
    #                                       Wnorder)
    #   Wnderiv <- 3
    #   Wpenmat <- fda::eval.penalty(Wbasis, Wnderiv)
    #   WfdPar <- fdPar(Wbasis, Wnderiv, Wlambda, TRUE, 
    #                   Wpenmat)
    # }
  }
  WfdList <- Wbinsmth.init(percntrnk, nbin, WfdPar, grbg, 
                           optList, U)
  dataList <- list(U = U, N = N, n = n, titlestr = titlestr, 
                   optList = optList, WfdList = WfdList, key = key, grbg = grbg, 
                   WfdPar = WfdPar, noption = noption, nbin = nbin, scrrng = scrrng, 
                   scrfine = scrfine, scrvec = scrvec, scrjit = scrjit, 
                   itmvec = itmvec, percntrnk = percntrnk, thetaQnt = thetaQnt, 
                   Wdim = Wdim, PcntMarkers = PcntMarkers, titlestr = titlestr, 
                   grbg = grbg)
  return(dataList)
}


make.dataList <- function (U, key, optList, grbg = rep(0, n), scrrng = NULL, 
          titlestr = NULL, nbin = TestGardener:::nbinDefault(N), NumBasis = NULL, Wnorder = 5,
          WfdPar = NULL, jitterwrd = TRUE, PcntMarkers = c(5, 25, 50, 75, 95), 
          verbose = FALSE, Wlambda = 10000) 
{
  N <- nrow(U)
  n <- ncol(U)
  if (min(U) < 1) 
    stop("Zero data values encountered in data.  Are they score values?")
  if (length(key) != n && length(key) > 0) {
    stop("length of key is neither n or 0")
  }
  if (length(grbg) != n && length(grbg) > 0) {
    stop("length of grbg is not n")
  }
  noption <- matrix(0, n, 1)
  for (item in 1:n) {
    noption[item] <- length(optList$optScr[[item]])
  }
  Wdim <- sum(noption)
  scrvec <- matrix(0, N, 1)
  itmvec <- matrix(0, n, 1)
  for (i in 1:n) {
    for (j in 1:N) {
      scoreij <- optList$optScr[[i]][U[j, i]]
      if (!is.null(scoreij)) {
        if (is.na(scoreij)) 
          print(paste("is.na score:", j))
      }
      else {
        print(paste("score of length 0:", j))
      }
      scrvec[j] <- scrvec[j] + scoreij
      itmvec[i] <- itmvec[i] + scoreij
    }
  }
  scrmin <- min(scrvec)
  scrmax <- max(scrvec)
  if (is.null(scrrng)) {
    scrrng <- c(scrmin, scrmax)
  }
  nfine <- 101
  scrfine <- seq(scrrng[1], scrrng[2], len = nfine)
  denscdf <- seq(0, 1, len = nfine)
  thetaQnt <- seq(0, 100, len = 2 * nbin + 1)
  if (jitterwrd) {
    scrjit <- scrvec + rnorm(N) * 0.1
    scrjit[scrjit < scrmin] <- scrmin
    scrjit[scrjit > scrmax] <- scrmax
  }
  else {
    scrjit = scan("scrjit.txt", 0)
  }
  scrrnk <- matrix(0, N, 1)
  for (j in 1:N) {
    scrrnk[j] <- sum(scrjit <= scrjit[j])
  }
  percntrnk <- 100 * scrrnk/N
  if (is.null(NumBasis)) {
    NumBasis = TestGardener:::NumBasisDefault(N)
  }
  if (is.null(WfdPar)) {
    # if (quadwrd) {
      Wnbasis <- NumBasis
      Wbasis <- fda::create.bspline.basis(c(0, 100), Wnbasis, 
                                          Wnorder)
      WfdPar <- fdPar(Wbasis)
    # }
    # else {
    #   Wnbasis <- NumBasis
    #   Wbasis <- fda::create.bspline.basis(c(0, 100), Wnbasis, 
    #                                       Wnorder)
    #   Wnderiv <- 3
    #   Wpenmat <- fda::eval.penalty(Wbasis, Wnderiv)
    #   WfdPar <- fdPar(Wbasis, Wnderiv, Wlambda, TRUE, 
    #                   Wpenmat)
    # }
  }
  WfdList <- TestGardener:::Wbinsmth.init(percntrnk, nbin, WfdPar, grbg, 
                           optList, U)
  dataList <- list(U = U, optList = optList, WfdList = WfdList, 
                   key = key, grbg = grbg, WfdPar = WfdPar, noption = noption, 
                   nbin = nbin, scrrng = scrrng, scrfine = scrfine, scrvec = scrvec, 
                   scrjit = scrjit, itmvec = itmvec, percntrnk = percntrnk, 
                   thetaQnt = thetaQnt, Wdim = Wdim, PcntMarkers = PcntMarkers, 
                   titlestr = titlestr)
  return(dataList)
}
