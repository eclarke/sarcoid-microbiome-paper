
initialize_defaults <- function() {
  opts <- new.env(parent=emptyenv())

  opts$data_fp <- "data_files"
  opts$figure_fp = "figures"

  ## Barchart options

  # Max number of barcharts to show in a plot (will pad to this)
  opts$barchart.max <- 46
  opts$barchart.width <- 1.25
  opts$barchart.height <- 9

  opts$barchart.theme <- theme_bw() + theme(
    axis.line=element_blank(),
    axis.line.x=element_blank(),
    axis.line.y=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(hjust=1, size=6, family="Helvetica"),
    axis.title=element_blank(),
    axis.ticks.x=element_blank(),

    legend.text=element_text(size=14),
    legend.position="top",
    legend.key.height=unit(1, "lines"),
    legend.key.width=unit(1, "lines"),

    panel.margin.x=unit(0.05, "lines"),
    panel.grid=element_blank(),
    panel.border=element_rect(fill=NA, color="black", size=0.2, linetype=1),

    plot.margin=unit(c(0.05,0,0.05,0), "lines"),

    strip.text.y=element_blank(),
    # strip.text.y=element_text(size=5),
    strip.background=element_blank()
  )

  opts$barchart.guide <- guides(
    fill=guide_legend(
      ncol=1, title = NULL,
      label.theme = element_text(size=10, angle=0),
      override.aes = list(color="white", size=0.4)))

  # Makes LN and BAL filled circles, and paraffin and prewash open circles
  opts$ln.labels <-  c("lymph_node"="\u25CF", "paraffin"="\u25CB", "water"="\u2205")
  opts$bal.labels <-  c("BAL"="\u25CF", "prewash"="\u25CB", "water"="\u2205")

  ## PCoA plot options

  # Consistent color/fill pairing for PCoA plots
  opts$pcoa.colors <- c(healthy="#a6611a", sarcoidosis="#018571")
  opts$pcoa.fills <- c(healthy="#dfc27d", sarcoidosis="#80cdc1")
  opts$pcoa.height <- 3.7
  opts$pcoa.width <- 3

  opts$pcoa.theme <- theme_classic(base_size = 10) +
    theme(
      axis.ticks=element_blank(),
      axis.text=element_blank(),
      aspect.ratio=10,
      legend.title=element_blank(),
      legend.background=element_blank(),
      legend.position="bottom",
      plot.title=element_blank(),
      plot.subtitle=element_blank())

  return(opts)
}

#' @export
opts <- initialize_defaults()
