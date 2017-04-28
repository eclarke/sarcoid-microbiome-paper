
#' Groups taxa into top members and samples by their larger group (if any)
#' in preparation for plotting summary barcharts
#' @param agg an agglomerated dataframe (\code{\link{eclectic::agglomerate}})
#' @param grouping.column a column that defines larger groups of samples (optional; set to NA to omit)
#' @param taxa.level the taxonomic rank to show in the barcharts
#' @param top.n how many colors/levels to display
#' @export
MakeBarchartData <- function(
  agg, sampleset, datatype, grouping.col="SubjectID", taxa.level="Order",
  top.n = 8)
{

  if (missing(sampleset)) sampleset <- get("sampleset", parent.frame())
  if (missing(datatype)) datatype <- get("datatype", parent.frame())
  agg$GroupingTaxa <- as.character(agg[[taxa.level]])
  agg$GroupingTaxa[is.na(agg$GroupingTaxa)] <- "Other"
  agg$GroupingTaxa <- factor(agg$GroupingTaxa)

  # Finds the top three taxa in each sample group, then votes on the top n
  # taxa that appear the most times
  top.taxa <- agg %>% filter(GroupingTaxa != "Other") %>%
    group_by_(.dots=c(grouping.col, "GroupingTaxa")) %>%
    tally(proportion) %>% top_n(3) %>% group_by(GroupingTaxa) %>%
    summarize(n=n_distinct(SubjectID)) %>% top_n(top.n) %>% droplevels
  top.taxa <- top.taxa$GroupingTaxa

  nixed.taxa <- levels(agg$GroupingTaxa)[!levels(agg$GroupingTaxa) %in% top.taxa]

  if (length(nixed.taxa) == 0) nixed.taxa <- c("placeholder")

  data <- agg %>%
    mutate(taxa = as.character(fct_collapse(GroupingTaxa, "Other" = nixed.taxa))) %>%
    mutate(taxa = plyr::revalue(taxa, c("Other"=NA))) %>%
    group_by_(.dots=c("SampleID", "SampleType", "StudyGroup", grouping.col, "taxa")) %>%
    summarize(count=sum(count, na.rm=TRUE)) %>%
    droplevels() %>% ungroup() %>%
    mutate(taxa = fct_explicit_na(taxa, na_level="Other")) %>%
    mutate(taxa = fct_relevel(taxa, as.character(top.taxa))) %>%
    mutate(
      taxa=lvl_toend(taxa, "Fungi (unclassified)"),
      taxa=lvl_toend(taxa, "Other")
    )

  attr(data, "labels") <- list(
    "A"=opts$ln.labels,
    "B"=opts$ln.labels,
    "C"=opts$bal.labels)[[sampleset]]


  attr(data, "fill") <- list(
    "16S"=bacterial.colors,
    "ITS"=fungal.colors)[[datatype]]

  data
}

#' Plots the results of MakeBarchartData as a set of grouped barcharts.
#' @param .data result from \code{\link{MakeBarchartData}}
#' @param x.col the mapping for the x-axis
#' @param fill.col the mapping for the fill aesthetic
#' @param grouping.col the grouping columns specified in MakeBarchartData, if any
#' @param x.labs if desired, what to revalue the x-axis labels to
#' @param fill.values named vector for fill colors
#' @param fill.breaks ordering of the fills (currently a noop)
#' @param theme a theme object to add
#' @return a \code{\link{ggplot2::ggplot}} object
#' @export
PlotBarchart <- function(
  .data, x.col="SampleType", fill.col="taxa", grouping.col="SubjectID",
  x.labs = attributes(.data)$labels,
  fill.values=attributes(.data)$fill,
  fill.breaks = rev(levels(.data$taxa)),
  theme=opts$barchart.theme,
  guide=opts$barchart.guide)
{

  x.scale <- if (is.null(x.labs)) {
    scale_x_discrete(expand=c(0,0))
  } else {
    scale_x_discrete(labels=x.labs, expand=c(0,0))
  }

  if (is.null(fill.values)) {
    fill.values <- quick_palette(levels(.data$taxa)) %>%
      eclpalettes::color_level("Other", "grey85")
  }

  p <- ggplot(.data, aes_string(x=x.col, y="count", fill=fill.col))

  if (missing(grouping.col)) {
    facet_formula <- ". ~ StudyGroup"
  } else {
    facet_formula <- sprintf(". ~ %s + StudyGroup", grouping.col)
  }

  p + geom_bar(stat="identity", position="fill", color="white", size=0.1, width=0.95) +
    x.scale +
    scale_fill_manual(values=fill.values) +#, breaks=fill.breaks) +
    scale_y_continuous("", expand=c(0,0)) +
    facet_wrap(c(grouping.col), switch="y") +
    coord_flip() +
    theme +
    guide
}

#' @export
SplitBarcharts <- function(
  plot, groups=c("sarcoidosis", "healthy", "extr_ctrl"),
  variables=c("StudyGroup", "SubjectID"), pad_to=45, ncol=1)
{

  plots <- lapply(groups, function(group) {
    data <- plot$data %>% filter(StudyGroup == group)

    # browser()
    little.barcharts <- plyr::dlply(data, .variables=variables, function(dat) {
      p <- plot %+% dat
      p + theme(legend.position="none")
    })

    n <- length(little.barcharts)
    if (n < pad_to) {
      little.barcharts[[n + 1]] <- grid::nullGrob()
    } else if (n > pad_to) {
      stop(sprintf("Too many barcharts (%d), adjust pad_to parameter", n))
    }

    padding <- rep(n+1, pad_to-n)
    layout <- matrix(c(1:n, padding), ncol=ncol)

    gridExtra::arrangeGrob(
      grobs=little.barcharts, layout_matrix = layout)
  })
  names(plots) <- groups
  plots$legend <- get_ggplot_legend(plot)
  plots
}

#' @export
SaveBarcharts <- function(
  barcharts, path, height=opts$barchart.height, width=opts$barchart.width)
{
  sapply(names(barcharts), function(name) {
    filename <- sprintf(path, name)
    if (name == "legend") {
      width=2
      height=2
    }
    ggsave(
      filename, barcharts[[name]], width=width, height=height, device=cairo_pdf)
    })
}

pick_attr <- function(data, attr) {
  attributes(data)[[attr]]
}
