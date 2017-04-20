params <- new.env()

#' @export
opts <- list(
  get = function(key) params$key,
  set = function(key, value) params$key <- value
)
