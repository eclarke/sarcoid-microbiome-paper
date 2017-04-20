params <- new.env()

#' Get and set global options.
#' @export
opts <- list(
  get = function(key) params$key,
  set = function(key, value) params$key <- value
)
