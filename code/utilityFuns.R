
# same as diag, except return x if it isn't a matrix or has 1 column
myDiag = function(x) {
  if(is.matrix(x) && (ncol(x) != 1)) {
    diag(x)
  } else {
    x
  }
}