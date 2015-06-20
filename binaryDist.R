library("Rcpp")
library("plyr")
sourceCpp("bDist.cpp")

# Converts a binary integer vector into a packed raw vector,
# padding out at the end to make the input length a multiple of packWidth
packRow <- function(row, packWidth = 64L) {
  packBits(as.raw(c(row, rep(0, (packWidth - length(row)) %% packWidth))))
}

as.PackedMatrix <- function(x, packWidth = 64L) {
  UseMethod("as.PackedMatrix")
}

# Converts a binary integer matrix into a packed raw matrix
# padding out at the end to make the input length a multiple of packWidth
as.PackedMatrix.matrix <- function(x, packWidth = 64L) {
  stopifnot(packWidth %% 8 == 0, class(x) %in% c("matrix", "Matrix"))
  storage.mode(x) <- "raw"
  if (ncol(x) %% packWidth != 0) {
    x <- cbind(x, matrix(0L, nrow = nrow(x), ncol = packWidth - (ncol(x) %% packWidth)))
  }
  out <- packBits(t(x))
  dim(out) <- c(ncol(x) %/% 8, nrow(x))
  class(out) <- "PackedMatrix"
  out
}

# Converts back to an integer matrix
as.matrix.PackedMatrix <- function(x) {
  out <- rawToBits(x)
  dim(out) <- c(nrow(x) * 8L, ncol(x))
  storage.mode(out) <- "integer"
  t(out)
}

# Generates random sparse data for testing the main function
makeRandomData <- function(nObs, nVariables, maxBits, packed = FALSE) {
  x <- replicate(nObs, {
    y <- integer(nVariables)
    y[sample(nVariables, sample(maxBits, 1))] <- 1L
    if (packed) {
      packRow(y, 64L)
    } else {
      y
    }
  })
  if (packed) {
    class(x) <- "PackedMatrix"
    x
  } else {
    t(x)
  }
}

# Reads a binary matrix from file or character vector
# Borrows the first bit of code from read.table
readPackedMatrix <- function(file = NULL, text = NULL, packWidth = 64L) {
  if (missing(file) && !missing(text)) {
    file <- textConnection(text)
    on.exit(close(file))
  }
  if (is.character(file)) {
    file <- file(file, "rt")
    on.exit(close(file))
  }
  if (!inherits(file, "connection")) 
    stop("'file' must be a character string or connection")
  if (!isOpen(file, "rt")) {
    open(file, "rt")
    on.exit(close(file))
  }
  lst <- list()
  i <- 1
  while(length(line <- readLines(file, n = 1)) > 0) {
    lst[[i]] <- packRow(as.integer(strsplit(line, "", fixed = TRUE)[[1]]), packWidth = packWidth)
    i <- i + 1
  }
  out <- do.call("cbind", lst)
  class(out) <- "PackedMatrix"
  out
}

# Wrapper for the C++ code which 
binaryDist <- function(x) {
  if (class(x) != "PackedMatrix") {
    x <- as.PackedMatrix(x)
  }
  dst <- bDist(x)
  attr(dst, "Size") <- ncol(x)
  attr(dst, "Diag") <- attr(dst, "Upper") <- FALSE
  attr(dst, "method") <- "binary"
  attr(dst, "call") <- match.call()
  class(dst) <- "dist"
  dst
}

x <- makeRandomData(2000, 400, maxBits = 5, packed = TRUE)

system.time(bd <- binaryDist(x))
