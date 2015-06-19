myDist <- function(x) {
  x <- x[, colSums(x) > 0]
  dist(x, "binary")
}

library("Rcpp")
library("plyr")
sourceCpp("bDist.cpp")

# Converts a binary integer vector into a packed raw vector,
# padding out at the end to make the input length a multiple of 8
packRow <- function(row) {
  packBits(as.raw(c(row, rep(0, (8 - length(row)) %% 8 ))))
}

makeRandomData <- function(nr, nc, maxBits, packed = FALSE) {
  t(replicate(nr, {
    y <- integer(nc)
    y[sample(nc, sample(maxBits, 1))] <- 1L
    if (packed) {
      packBits(y)
    } else {
      y
    }
  }))
}

y <- makeRandomData(2000, 400, maxBits = 5, packed = TRUE)

padPackMatrix <- function(x) {
  storage.mode(x) <- "raw"
  if (ncol(x) %% 8 != 0) {
    x <- cbind(x, matrix(0, nrow = nrow(x), ncol = 8 - (ncol(x) %% 8)))
  }
  aaply(x, 1, packBits)
}

binaryDist <- function(x) {
  if (storage.mode(x) != "raw") {
    x <- padPackMatrix(x)
  }
  dst <- bDist(x)
  class(dst) <- "dist"
  attr(dst, "Size") <- nrow(x)
  attr(dst, "Diag") <- attr(dst, "Upper") <- "FALSE"
  attr(dst, "method") <- "binary"
  dst
}

system.time(bd <- binaryDist(y))
