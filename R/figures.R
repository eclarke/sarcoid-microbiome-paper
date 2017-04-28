






#' @export
#' @importFrom magrittr %$%
#' @import ggplot2
GenerateBarchart <- function(
  data, transformations, fills=NULL, labels=NULL, ...)
{

  stopifnot(inherits(data, "SarcoidDataset"))

  plot <- with(data, {

    barchart.data <- agg %>%
      transformations %>%
      MakeBarchartData(...) %>%
      mutate(taxa = fct_rev(taxa))

    sarcoid <- barchart.data %>%
      filter(StudyGroup == "sarcoidosis") %>% droplevels
    healthy <- barchart.data %>%
      filter(StudyGroup == "healthy")
    ext_ctl <- barchart.data %>%
      filter(StudyGroup == "extr_ctrl") %>%
      group_by(SubjectID, SampleType, StudyGroup, taxa) %>%
      summarize(count=sum(count))

    fills <- ifelse(!is.na(fills), fills, quick_palette(barchart.data$taxa))

    p <- PlotSingleBarchart(
      barchart.data,
      x.labels = labels,
      fill.values = fills,
      fill.breaks = rev(levels(barchart.data$taxa))) +
      horiz.barchart.guide +
      horiz.barchart.theme

    legend <- get_ggplot_legend(p)
    sarcoid.plots <- make_standard_barcharts(p, sarcoid, "Sarcoidosis")
    healthy.plots <- make_standard_barcharts(p, healthy, "Healthy")
    control.plots <- make_standard_barcharts(p, ext_ctl, "Extraction controls")

    gridExtra::arrangeGrob(
      grobs=list(sarcoid.plots, healthy.plots, control.plots, legend),
      ncol=4, top=sprintf("Set %s/%s", sampleset, kingdom))
  })
}

#' Saves a barchart using global options.
#' @export
SaveBarchartFigure <- function(plot, filename, ncol=4) {
  ggsave(
    filename, plot=plot, height=barcharts.height, width=ncol*barcharts.width,
    device=cairo_pdf)
}
# cohort <- "A"
# kingdom <- "Bacteria"
# figure.name <- sprintf(
#   "Main/Fig_1A_Set%s_%s_Barcharts.pdf", cohort, kingdom)
#
# barchart.data <- agg %>%
#   filter(SampleType %in% c("lymph_node", "paraffin", "extr_blank")) %>%
#   mutate(SampleType = fct_recode(SampleType, "water"="extr_blank")) %>%
#   mutate(SampleType = relevel(SampleType, "paraffin")) %>%
#   MakeBarchartData() %>%
#   mutate(taxa=fct_rev(taxa)) %>%
#   mutate(taxa=fct_recode(taxa, "Other Bacteria"="Other"))
#
# sarcoid <- barchart.data %>%
#   filter(StudyGroup == "sarcoidosis") %>% droplevels
# healthy <- barchart.data %>%
#   filter(StudyGroup == "healthy")
# ext_ctl <- barchart.data %>%
#   filter(StudyGroup == "extr_ctrl") %>%
#   group_by(SubjectID, SampleType, StudyGroup, taxa) %>%
#   summarize(count=sum(count))
#
# p <- PlotSingleBarchart(
#   barchart.data,
#   x.labels = ln.labels,
#   fill.values = bacterial.colors,
#   fill.breaks = rev(levels(barchart.data$taxa))) +
#   horiz.barchart.guide +
#   horiz.barchart.theme
#
# legend <- get_ggplot_legend(p)
# sarcoid.plots <- make_standard_barcharts(p, sarcoid, "Sarcoidosis")
# healthy.plots <- make_standard_barcharts(p, healthy, "Healthy")
# control.plots <- make_standard_barcharts(p, ext_ctl, "Extraction controls")
# # browser()
# combined.plots <- arrangeGrob(
#   grobs=list(sarcoid.plots, healthy.plots, control.plots, legend),
#   ncol=4, top=sprintf("Set %s/%s", cohort, kingdom))
# ggsave(
#   file.path(params$fig_dir, figure.name),
#   plot = combined.plots,
#   height=barcharts.height, width=4*barcharts.width, device=cairo_pdf)
