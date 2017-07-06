
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
  barcharts, path, height=opts$barchart.height, width=opts$barchart.width, device=opts$device)
{
  sapply(names(barcharts), function(name) {
    filename <- sprintf(path, name)
    if (name == "legend") {
      width=2
      height=2
    }
    ggsave(
      filename, barcharts[[name]], width=width, height=height, device=cairo_failsafe)
    })
}

pick_attr <- function(data, attr) {
  attributes(data)[[attr]]
}


#' Prepares data for making the summary barcharts showing the dominant taxa
#' by group (rather than by sample).
#' @param agg agglomerated data frame, filtered as desired
#' @param s sample metadata
#' @param top top taxa to select from each group (displayed taxa will be union)
#' @param min.rank min taxonomic rank to display (will group at this or higher)
#' @export
MakeTopBarchartData <- function(agg, s, top=6, min.rank="Species") {
  medians <- agg %>%
    SarcoidMicrobiome:::MakeHeatmapData(
      s, min.samples = 1, min.rank.1 = min.rank) %>%
    # Take care of indeterminate fungi labels
    mutate(
      MinRank1 = forcats::fct_recode(
        MinRank1,
        Indeterminate="Fungi",
        Indeterminate="uncultured fungus (phylum)")) %>%
    # Change all NAs to 0s (otherwise medians are wrong)
    mutate(proportion = ifelse(is.na(proportion), 0, proportion)) %>%
    group_by(SampleType, StudyGroup, MinRank1) %>%
    # Error bars are determined by the quantiles of the binomial distribution
    summarize(
      med=median(proportion),
      max = max(proportion),
      min = min(proportion),
      conf.high=sort(proportion)[qbinom(0.975, length(proportion), 0.5)],
      conf.low=sort(proportion)[qbinom(0.025, length(proportion), 0.5)]) %>%
    mutate(rank = row_number(desc(med))) %>%
    mutate(Group = paste(StudyGroup, SampleType)) %>%
    ungroup %>%
    mutate(
      MinRank1 = fct_reorder(MinRank1, med, na.rm=TRUE, .desc=TRUE)) %>%
    mutate(
      Group = fct_recode(
        Group,
        "Healthy lymph node"="healthy lymph_node",
        "Healthy BAL"="healthy BAL",
        "Healthy paraffin"="healthy paraffin",
        "Healthy prewash"="healthy prewash",
        "Sarcoidosis lymph node"="sarcoidosis lymph_node",
        "Sarcoidosis BAL"="sarcoidosis BAL",
        "Sarcoidosis paraffin"="sarcoidosis paraffin",
        "Sarcoidosis prewash"="sarcoidosis prewash"
      )) %>%
    mutate(Group = fct_relevel(
      Group,
      "Sarcoidosis lymph node", "Sarcoidosis paraffin",
      "Sarcoidosis BAL", "Sarcoidosis prewash",
      "Healthy lymph node", "Healthy paraffin",
      "Healthy BAL", "Healthy prewash"
    ))

  # Isolate to only the top X taxa
  top.taxa <- unique((medians %>% filter(
    rank <= top, MinRank1 != "Fungi"))$MinRank1)
  medians %>% filter(MinRank1 %in% top.taxa)
}

#' Plots data prepared by \link{MakeTopBarchartData}.
#' @export
PlotTopBarcharts <- function(dat) {
  ggplot(dat, aes(x=MinRank1, y=med, fill=Group, color=Group)) +
    geom_crossbar(aes(ymin=conf.low, ymax=conf.high), position=position_dodge(0.7), width=0.6) +
    geom_crossbar(aes(ymin=med, ymax=med), color="white", position=position_dodge(0.7), width=0.4, fatten=1) +
    theme_classic(base_size = 14) +
    guides(color="none", fill=guide_legend(override.aes = list(color=NA))) +
    ylab("Median (95% CI) abundance") +
    theme(
      axis.text.x = element_text(angle=-35, hjust=0, vjust=1),
      axis.title.x = element_blank(),
      legend.position=c(0.9,0.9),
      legend.justification=c(1,1),
      legend.spacing = unit(0, "lines"),
      legend.margin=margin(t=0, b=4, r=4, l=4),
      legend.box.background=element_rect(color="black"),
      legend.title=element_blank(),
      legend.key.size=unit(0.8, "lines"),
      plot.margin=unit(c(.4,3,.4,.4), "lines")
    )
}
