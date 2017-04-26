
initialize_defaults <- function() {
  opts <- new.env(parent=emptyenv())
  opts$data_fp <- "data_files"

  opts$figure_fp = "figures"

  opts$barchart.max <- 46
  opts$barchart.col <- 1
  opts$barchart.width <- 1.8*opts$barchart.col
  opts$barchart.height <- (opts$barchart.max*0.5)/opts$barchart.col


  opts$barchart.theme <- theme_bw() + theme(
    axis.line.y=element_line(color="black", size=1, linetype =1),
    axis.text.y=element_blank(),
    axis.text.x=element_text(hjust=0, size=14),
    axis.title=element_blank(),
    axis.ticks.x=element_blank(),

    legend.text=element_text(size=14),
    legend.position="top",
    legend.key.height=unit(1, "lines"),
    legend.key.width=unit(1, "lines"),

    panel.margin.x=unit(0.8, "lines"),
    panel.grid=element_blank(),
    panel.border=element_rect(fill=NA, color="black", size=0.8, linetype=1),

    plot.margin=unit(c(0.4,0,0.4,0), "lines"),

    strip.text.y=element_blank(),
    strip.background=element_blank()
  )
  return(opts)
}

#' @export
opts <- initialize_defaults()
