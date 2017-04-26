params <- new.env()

#' Get and set global options.
#' @export
opts <- list(
  get = function(key) params$key,
  set = function(key, value) params$key <- value
)



# General functions -------------------------------------------------------

#' Creates a MinRank vector from an agglomerated data frame (based on
#' eclectic::tax_climber). Use inside dplyr::mutate or directly assign to col.
#' MinRank columns are the most specific taxonomic assignments (up to a
#' threshold)
#' @param agg agglomerated data frame (with taxonomic ranks)
#' @param end.rank the most specific taxonomic rank if available
#' @param rank a character vector of taxonomic ranks present in agg
#' @param ... additional arguments passed to eclectic::tax_climber
min_rank <- function(agg, end.rank, ranks=qiimer::taxonomic_ranks, ...) {
  .agg <- filter(agg, !is.na(otu))
  md <- .agg[,c(ranks, "otu")]
  md <- data.frame(distinct(md))
  rownames(md) <- md$otu
  minrank <- tax_climber(md$otu, md, end=end.rank, ranks=ranks, ...)
  left_join(select(agg, otu), data.frame(otu=md$otu, MinRank=minrank))$MinRank
}

#' Creates a counts matrix from an agglomerated data frame.
#' Optionally, metadata columns can be included (useful for some functions), in
#' which case the output is a data.frame, not a matrix.
#' @param agg agglomerated data frame
#' @param additional.columns any columns to include as columns in the matrix
#' @return default: a numeric matrix with otu as rows and samples as columns.
#' If additional columns specified, a data.frame with otus as columns and samples as rows
counts_matrix <- function(agg, additional.columns) {
  cols <- c("otu", "SampleID", "count")
  if (!missing(additional.columns)) cols <- c(cols, additional.columns)
  .mat <- agg %>% select_(.dots=cols) %>%
    spread(SampleID, count, fill=0) %>%
    as.data.frame %>%
    filter(!is.na(otu))
  if (!missing(additional.columns)) {
    return(.mat)
  } else {
    rownames(.mat) <- .mat$otu
    .mat$otu <- NULL
    .mat <- t(as.matrix(.mat))
    return(.mat)
  }
}

#' Root a tree using a random node as the root.
#' @param tree a phylogenetic tree (class "phylo")
#' @param otus vector of otus to keep
#' @return rooted tree
root_tree <- function(tree, otus) {
  # Root the tree using a random node as the root
  dropped.otus <- tree$tip.label[!(tree$tip.label %in% otus)]
  .tree <- drop.tip(tree, dropped.otus)
  .rtree <- root(.tree, sample(.tree$tip.label, 1), resolve.root = TRUE)
  stopifnot(is.rooted(.rtree))
  return(.rtree)
}

#' Non-random version of \code{\link{forcats::fct_anon}}
#' @inheritParams forcats::fct_anon
#' @export
fct_anon_deterministic <- function(f, prefix="") {
  levels <- paste0(prefix, forcats:::zero_pad(seq_len(nlevels(f))))
  lvls_revalue(f, levels)
}

lvl_toend <- function(f, level) {
  forcats:::check_factor(f)
  lvls <- levels(f)[levels(f) != level]
  lvls <- c(lvls, level)
  factor(f, levels=lvls)
}

get_ggplot_legend <- function(p) {
  gp <- ggplot_gtable(ggplot_build(p))
  legend_idx <- which(sapply(gp$grobs, function(x) x$name) == "guide-box")
  gp$grobs[[legend_idx]]
}
# Color themes ------------------------------------------------------------

scale_fill_custom <- function(..., levels, palette="brewer", other="grey85") {
  palettes <- list(
    brewer=RColorBrewer::brewer.pal(length(levels), "Set3"),
    # Monochrome palettes
    redmono = c("#99000D", "#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1", "#FEE0D2", "#FFF5F0"),
    greenmono = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B", "#C7E9C0", "#E5F5E0", "#F7FCF5"),
    bluemono = c("#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF"),
    grey8mono = c("#000000","#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0"),
    grey6mono = c("#242424", "#494949", "#6D6D6D", "#929292", "#B6B6B6", "#DBDBDB"),

    # Qualitative color schemes by Paul Tol
    tol1qualitative=c("#4477AA"),
    tol2qualitative=c("#4477AA", "#CC6677"),
    tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677"),
    tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677"),
    tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"),
    tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499"),
    tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499"),
    tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499"),
    tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"),
    tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
    tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
    tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499"),
    tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
    tol21rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")[1:length(levels)],
    viridis=viridis::plasma(n=length(levels)),
    gdocs=gdocs_pal()(length(levels))
  )
  pal <- palettes[[palette]]
  names(pal) <- levels
  pal["Other"] <- other
  scale_fill_manual(..., values=pal, na.value=other)
}

alt_heatmap_scale <- function(threshold=0.6, ...) {
  scale_fill_gradientn(colors = c("#FFE591", "#EF2D2D", "#EF2D2D"), values=c(0,threshold,1), ...)
  # scale_fill_gradientn(colors = rev(heat.colors(5)[1:3]), values=c(0,threshold,1), ...)
}

# Heatmap functions -------------------------------------------------------


#' Prepare agglomerated dataframe for heatmap display by aggregating counts
#' by minimum rank, completing missing cases, and clustering rows
#' @param min.samples the minimum number of samples required for a taxa to appear
#' @param s sample metadata dataframe
#' @param min.rank.1 the most specific rank to aggregate by (determines y-axis)
#' Uses tax_climber to select next most specific rank if not available.
#' @param min.rank.2 (optional) higher rank to prefix min.rank.1 (e.g. Phylum)
#' Uses tax_climber to select next most specific rank if not available.
#' @param min.samples a row must appear more than this number of samples to be shown
MakeHeatmapData <- function(agg, s, min.rank.1 = "Genus", min.rank.2 = NA, min.samples=1) {
  # Create the y-axes for the data via the desired taxonomic rank level(s)
  .agg <- agg %>% mutate(
    MinRank1 = .makeMinRankFromAgg(agg, end=min.rank.1, label=TRUE, sep="__"))
  if (!is.na(min.rank.2)) {
    .agg <- .agg %>% mutate(
      MinRank2 = .makeMinRankFromAgg(agg, end=min.rank.2, label=TRUE, sep="__"),
      MinRank1 = paste(MinRank2, MinRank1)
    )
  }

  .agg <- .agg %>%
    # Aggregate counts by minimum available rank
    group_by(SampleID, MinRank1) %>%
    summarize(count = sum(count)) %>%
    # Update proportions
    group_by(SampleID) %>%
    mutate(proportion=count/sum(count))

  # Convert to matrix form
  .mat <- reshape2::dcast(
    .agg, MinRank1 ~ SampleID, value.var="proportion", fill = 0) %>%
    filter(!is.na(MinRank1))
  rownames(.mat) <- .mat$MinRank1
  .mat$MinRank1 <- NULL

  # Cluster and pull out row order
  minrank_order <- dist(.mat, method="euclidean") %>%
    hclust(method="average") %$% order

  .agg %>%
    # Reorder by clustering order
    mutate(MinRank1 = factor(MinRank1, levels=rownames(.mat)[minrank_order])) %>%
    # Complete missing cases: fills missing combinations in with NA
    # instead of just omitting the row entirely (necessary to show blank cells)
    ungroup() %>%
    tidyr::complete(SampleID, nesting(MinRank1)) %>%
    select(SampleID, MinRank1, count, proportion) %>%
    distinct() %>%
    # Filter taxa that appear in fewer than required number of samples
    group_by(MinRank1) %>%
    filter(sum(proportion > 0, na.rm=TRUE) > min.samples) %>%
    # Re-add metadata
    left_join(s, by="SampleID", all.y=TRUE) %>%
    ungroup()
}

#' Plot the results of MakeHeatmapData.
#' @param heatmap.data the dataframe from MakeHeatmapData
#' @param use.reads if TRUE, plot the raw readcounts; if FALSE, use proportions
#' @param threshold at what relative proportion the color scale goes to red
#' @param x.text display labels on the x-axis
PlotHeatmap <- function(
  heatmap.data, use.reads=TRUE, threshold=0.6, x.text=TRUE, x.axis='SampleID') {

  if (use.reads) {
    fill.var <- "count"
    scale.name <- "Reads"
    fill.scale <- saturated_rainbow_cts(threshold = threshold, name=scale.name)
  } else {
    fill.var <- "proportion"
    scale.name <- "Abundance"
    fill.scale <- saturated_rainbow_pct(threshold = threshold, name=scale.name)
  }

  p <- ggplot(heatmap.data, aes_string(x.axis, "MinRank1", fill=fill.var)) +
    geom_tile(color="grey40", size=0.4) +
    facet_grid(. ~ StudyGroup + SampleType, space="free", scales="free") +
    fill.scale +
    theme_grey() +
    theme(
      strip.text.y = element_text(angle=0, vjust=0),
      strip.text.x = element_text(angle=90, vjust=0),
      panel.border = element_blank()) +
    if (x.text) {
      p <- p + theme(
        axis.text.x = element_text(angle=45, hjust=1, vjust=1))
    } else{
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
    }
  p
}


# Barchart functions ------------------------------------------------------


#' Groups taxa into top members and samples by their larger group (if any)
#' in preparation for plotting summary barcharts
#' @param agg an agglomerated dataframe
#' @param grouping.column a column that defines larger groups of samples (optional; set to NA to omit)
#' @param taxa.level.1 the primary taxonomic rank to show in the barcharts
#' @param taxa.level.2 the second, higher rank to show
#' @param top.1 how many colors/levels to display
#' @param top.2 how many to display for what didn't fit into top.1
MakeBarchartData <- function(agg, grouping.col="SubjectID",
                             taxa.level.1="Order", taxa.level.2="Order",
                             top.1 = 8, top.2 = 0) {

  .data <- agg

  .data$GroupingTaxa <- as.character(.data[[taxa.level.1]])
  .data$GroupingTaxa[is.na(.data$GroupingTaxa)] <- "Other"
  .data$GroupingTaxa <- factor(.data$GroupingTaxa)

  .data$GroupingTaxa2 <- as.character(.data[[taxa.level.2]])
  .data$GroupingTaxa2[is.na(.data$GroupingTaxa2)] <- "Other"
  .data$GroupingTaxa2 <- factor(.data$GroupingTaxa2)

  # top.taxa <- .data %>% group_by(GroupingTaxa) %>%
  #   summarize(n=sum(proportion)) %>% top_n(top) %$% GroupingTaxa
  # top.taxa <- .data %>% filter(GroupingTaxa != "Other") %>%
  #   group_by_(.dots=c("GroupingTaxa", grouping.col)) %>%
  #   tally(proportion) %>% top_n(1) %>% select_(.dots=c("GroupingTaxa", "SubjectID")) %>% group_by(GroupingTaxa) %>%
  #   summarize(n=n()) %>% top_n(top.1) %>% droplevels %$% GroupingTaxa
  #
  # Finds the top three taxa in each sample group, then votes on the top n
  # taxa that appear the most times
  top.taxa <- .data %>% filter(GroupingTaxa != "Other") %>%
    group_by_(.dots=c(grouping.col, "GroupingTaxa")) %>%
    tally(proportion) %>% top_n(3) %>% group_by(GroupingTaxa) %>%
    summarize(n=n_distinct(SubjectID)) %>% top_n(top.1) %>% droplevels %$%
    GroupingTaxa
  nixed.taxa <- levels(.data$GroupingTaxa)[!levels(.data$GroupingTaxa) %in% top.taxa]
  # Gather the remaining taxa into their higher levels
  if (top.2 > 0) {

    .data2 <- .data[.data$GroupingTaxa %in% nixed.taxa,]
    top.taxa2 <- (.data2 %>% group_by_(.dots=c(grouping.col, taxa.level.2)) %>%
                    tally(proportion) %>% top_n(3) %>% group_by_(taxa.level.2) %>%
                    summarize(n=n_distinct(SubjectID)) %>% top_n(top.2) %>% droplevels)[[taxa.level.2]]
    names(top.taxa2) <- paste("Other", top.taxa2)
    remaining <- levels(.data[[taxa.level.2]])[!levels(.data[[taxa.level.2]]) %in% top.taxa2]
  }

  if (length(nixed.taxa) == 0) nixed.taxa <- c("placeholder")

  .data3 <- .data %>%
    mutate(taxa = as.character(fct_collapse(GroupingTaxa, "Other" = nixed.taxa)))
  if (top.2 > 0) {
    .data3 <- .data3 %>%
      mutate(GroupingTaxa2 = factor(ifelse(
        !GroupingTaxa2 %in% remaining,
        paste("Other", as.character(GroupingTaxa2)),
        as.character(GroupingTaxa2)))) %>%
      mutate(taxa2 = as.character(fct_collapse(GroupingTaxa2, "Other" = remaining))) %>%
      mutate(taxa = factor(ifelse(taxa == "Other", taxa2, taxa)))
  }
  .data3 <- .data3 %>%
    mutate(taxa = plyr::revalue(taxa, c("Other"=NA))) %>%
    group_by_(.dots=c("SampleID", "SampleType", "StudyGroup", grouping.col, "taxa")) %>%
    summarize(count=sum(count, na.rm=TRUE)) %>%
    # filter(count > 0, !is.na(count)) %>%
    droplevels() %>% ungroup() %>%
    mutate(taxa = fct_explicit_na(taxa, na_level="Other")) %>%
    mutate(taxa = fct_relevel(taxa, as.character(top.taxa)))

  return(.data3)
}

#' Plots the results of MakeBarchartData as a set of grouped barcharts.
#' @param .data data.frame from MakeBarchartData
#' @param x.axis.col the mapping for the x-axis
#' @param fill.col the mapping for the fill aesthetic
#' @param x.labels if desired, what to revalue the x-axis labels to
#' @param grouping.col the grouping columns specified in MakeBarchartData, if any
#' @param reorder.taxa a noop (do it outside this function)
#' @return a ggplot object
PlotBarchart <- function(.data, x.axis.col="SampleType", fill.col="taxa",
                         x.labels, grouping.col, reorder.taxa=TRUE) {
  if (missing(x.labels)) {
    x.scale <- scale_x_discrete(expand=c(0,0))
  } else {
    x.scale <- scale_x_discrete(labels=x.labels, expand=c(0,0))
  }

  p <- ggplot(.data, aes_string(x=x.axis.col, y="count", fill=fill.col))

  if (missing(grouping.col)) {
    facet_formula <- ". ~ StudyGroup"
  } else {
    facet_formula <- sprintf(". ~ %s + StudyGroup", grouping.col)
  }

  p +
    geom_bar(stat="identity", position="fill", color="white", size=0.4, width=1) +
    x.scale +
    scale_y_continuous(
      "", expand=c(0,0), labels=scales::percent) +
    guides(fill=guide_legend(nrow=1, title=NULL)) +
    theme_bw() +
    facet_grid(
      facet_formula,
      scales="free_x", space="free", switch="x") +
    theme(
      axis.text=element_text(size=12),
      legend.text=element_text(size=14),
      strip.text.x=element_text(angle=90, vjust=1, hjust=0.5),
      strip.background = element_blank(),
      # plot.background=element_rect(fill="white"),
      legend.position="top",
      panel.margin.x=unit(0.8, "lines"),
      panel.border=element_rect(color="black", size=0.8))
}

#' Based on PlotBarchart, but facets along ID.col. Intended for use with
#' MakeEqualSizePlots function. Not suggested for use on its own.
#' @param fill.values the palette to use for the scale_fill_manual
#' @param ID.col the column with the desired display indices
PlotSingleBarchart <- function(..., fill.values, fill.breaks=waiver(), ID.col="SubjectID") {
  p <- PlotBarchart(...) +
    scale_fill_manual(values=fill.values, breaks=fill.breaks) +
    facet_wrap(c(ID.col), switch = "y") +
    coord_flip()
}

#' Facets on a layout matrix to make individual plot sizes consistent. Removes
#' the legend; grab and display it separately with get_ggplot_legend.
#' @param p a ggplot object (whose data will be replaced by the chunked data)
#' @param data the dataframe to facet
#' @param .variables a character vector listing the columns to facet by in data
#' @param at_most make space for at most this many facets
#' @param ncol the number of columns in the layout matrix (must be multiple of at_most)
MakeEqualSizePlots <- function(p, data, .variables=c("StudyGroup", "SubjectID"),
                               at_most=45, ncol=1, extra_text="") {

  plots <- plyr::dlply(data, .variables=.variables, function(dat) {
    .p <- p %+% dat
    .p + theme(legend.position="none")
  })

  n <- length(plots)

  if (n < at_most) {
    plots[[n + 1]] <- nullGrob()# textGrob(extra_text)
  } else if (n > at_most) {
    stop(sprintf(
      "More plot elements (%d) than at_most parameter (%d).", n-1, at_most))
  }
  padding <- rep(n+1, at_most-n)
  layout <- matrix(c(1:n, padding), ncol=ncol)
  arrangeGrob(grobs=plots, layout_matrix=layout, top=extra_text)
}

# Unifrac functions -------------------------------------------------------

#' Calculates the generalized UniFrac distances for an agglomerated dataframe.
#' This involves a random rooting of the tree and rarefaction, which are
#' non-deterministic. For reproducibility, recommend calling set.seed before this
#' @param agg agglomerated dataframe
#' @param s sample metadata dataframe
#' @param tree a tree (will root using random node)
#' @param alpha the alpha parameter for GUniFrac
MakeUnifracData <- function(agg, s, tree, rarefy.depth, alpha=0.5) {

  # Root tree
  .rtree <- root_tree(tree, agg$otu)
  # Create sample matrix
  .mat <- counts_matrix(agg)
  # Keep only otus that appear in the tree
  .mat <- .mat[, colnames(.mat) %in% .rtree$tip.label]

  # Rarefy to specified depth
  if (!missing(rarefy.depth)) {
    .mat.rf <- vegan::rrarefy(.mat, sample=rarefy.depth)
  } else {
    .mat.rf <- .mat
  }

  # Calculate generalized UniFrac
  .unifrac <- GUniFrac(.mat.rf, .rtree, alpha=alpha)$unifracs
  # Return the generalized UniFrac matrix at specified alpha level
  .unifrac[,,paste0("d_",alpha)]
}

#' Performs principle coordinates analysis on a GUniFrac distance and merges
#' sample metadata for plotting.
#' @param unifrac.data the output of MakeUnifracData
#' @param s sample metadata dataframe
#' @param grouping.col column in s defining centroid groups
MakePCOAData <- function(unifrac.data, s, grouping.col) {
  .pc <- data.frame(ape::pcoa(unifrac.data)$vectors)[,1:4]
  .pc$SampleID <- rownames(.pc)
  .pc <- left_join(.pc, s)
  if (!missing(grouping.col)) {
    .pc <- group_by_(.pc, .dots=list(grouping.col)) %>%
      mutate(med.x=mean(Axis.1), med.y=mean(Axis.2))
  }
  .pc
}

#' Performs PERMANOVA on the given terms using the unifrac distance matrix.
#' @param unifrac.data matrix/dist of unifrac distances (coerced to dist if not)
#' @param s sample metadata
#' @param terms a string giving the RHS of the adonis formula (eg "counts*StudyGroup")
TestPermanova <- function(unifrac.data, s, terms) {

  .s <- filter(s, SampleID %in% rownames(unifrac.data))
  .s <- .s[match(rownames(unifrac.data), .s$SampleID), ]

  if (class(unifrac.data) != "dist") {
    warning("unifrac.data not dist; coerced to dist")
    unifrac.data <- dist(unifrac.data)
  }

  pn <- adonis(
    formula(sprintf("unifrac.data ~ %s", terms)),
    data = .s)

  broom::tidy(pn$aov.tab)
}

#' Standardized PCoA plot.
#' @param pcoa.data dataframe generated by MakePCOAData
#' @param fill column to set as fill aesthetic
#' @param linetype column to set as linetype aesthetic
#' @param fill.palette color palette to use as fill
#' @param color.palette color palette to use for outlines and line segments
PlotPCOA <- function(pcoa.data, fill = "StudyGroup", linetype="SampleType",
                     fill.palette, color.palette) {
  ggplot(
    pcoa.data,
    aes_string(x="Axis.1", y="Axis.2", fill=fill, color=fill)) +
    geom_segment(aes_string(xend="med.x", yend="med.y", linetype=linetype)) +
    geom_point(shape=21, size=2) +
    stat_ellipse(alpha=0.3) +
    scale_fill_manual(values=pcoa.fills) +
    scale_color_manual(values=pcoa.colors) +
    theme_classic(base_size = 16) +
    theme(
      axis.ticks=element_blank(),
      axis.text=element_blank(),
      aspect.ratio=1,
      # legend.position=c(0,1),
      # legend.justification=c(0,1),
      legend.title=element_blank(),
      legend.background=element_blank(),
      plot.subtitle=element_text(size=10))
}


# GLMM functions ----------------------------------------------------------

#' Create and test a mixed effect model between two sample types and study groups
#' (e.g. BAL and prewash, tissue and paraffin), where the interaction term for the
#' model is the coefficient of interest.
#' @param agg an agglomerated data frame; ensure the StudyGroup and SampleType
#'        levels are ordered correctly (i.e. interesting ones last)
#' @param grouping.var the column defining the sample/control grouping
#' @param tax.level the taxonomic level to test; lower levels will be summed
#' @param sparsity the proportions of zeros below which a taxa will be skipped
TestMatchedControlMixedModel <- function(agg, grouping.var="OriginatingSampleID", tax.level="Genus", sparsity=0.1) {
  # Create the data matrix from aggregated data frame, aggregating by the
  # specified taxonomy level
  .mat <- reshape2::dcast(
    agg,
    formula = formula(sprintf(
      "SampleID + %s + SampleType + StudyGroup + total ~ %s",
      grouping.var, tax.level)),
    value.var = "count",
    fill = 0,
    fun.aggregate = sum)

  # This is the testing function that builds the models
  .fitGLMM <- function(test.col) {
    # We return this default in case of error/no fit
    default <- list(null.model=NA, model=NA)

    counts <- .mat[, colnames(.mat)==test.col]
    if (sum(counts > 0)/length(counts) < sparsity) {
      message(sprintf("%s: too sparse, skipping", test.col))
      return(default)
    }

    # Build positive and negative counts for taxa under test
    test.data <- .mat[, c(1:5)]
    test.data$neg.cts <- test.data[["total"]] - counts
    test.data$pos.cts <- counts

    # The model of interest contains an interaction between StudyGroup and
    # SampleType; the null model does not. We are stricter about forcing
    # convergence for this model than the null model. Errors return NA.
    model1 <- plyr::try_default(
      lme4::glmer(
        formula=sprintf(
          "cbind(pos.cts, neg.cts) ~ StudyGroup * SampleType + (1|%s) + (1|SampleID)",
          grouping.var),
        data=test.data,
        family=binomial,
        control=glmerControl(
          optimizer="bobyqa",
          check.conv.singular = "warning",
          check.conv.grad = "warning",
          check.conv.hess = "warning")
      ), NA)

    # The null model
    model0 <- plyr::try_default(
      lme4::glmer(
        formula=sprintf(
          "cbind(pos.cts, neg.cts) ~ StudyGroup + SampleType + (1|%s) + (1|SampleID)",
          grouping.var),
        data=test.data,
        family=binomial,
        control=glmerControl(
          optimizer="bobyqa")
      ), NA)

    # Return NA if the primary model broke
    if (is.na(model1)) {
      message(sprintf("Nonconvergent: %s", test.col))
      return(default)
    } else {
      # browser()
      return(list(null.model=model0, model=model1))
    }
  }

  # Get a list of things to test:
  test.candidates <- data.frame(taxa=colnames(subset(.mat, select=-c(SampleID:total))))

  # Fit the models
  results <- test.candidates %>%
    mutate(res = purrr::map(taxa, .fitGLMM))

  # this function generates ANOVA p-values, returning NA if either of the
  # models is NA
  .safeAnova <- function(x) {
    if (!is.na(x$null.model) & !is.na(x$model)) {
      p <- anova(x$null.model, x$model)$`Pr(>Chisq)`[2]
      if (length(p) != 1 | is.null(p)) return(NA)
      else p
    } else {
      return(NA)
    }
  }

  # Calculate p-values for all taxa, clean up, and return results
  results %>%
    mutate(anova.p = purrr::map_dbl(res, .safeAnova)) %>%
    mutate(model = purrr::map(res, ~ broom::tidy(.x$model))) %>%
    dplyr::select(-res) %>%
    tidyr::unnest() %>%
    dplyr::select(-x)  # An artifact from tidy(NA)
}

#' Extracts the significant results from TestMatchedControlMixedModel and
#' adjusts for multiple testing (p-values corrected are only those for the
#' interaction terms with the correct magnitude sign)
#' @param glmm.results results from TestMatchedControlMixedModel
#' @param direction a binary function such as > or < determining coefficient relationship to 0
#' @param padj.threshold only return adjusted p-values below this value
ExtractSignificantResults <- function(glmm.results, direction = `>`, padj.threshold=0.1) {
  .results <- glmm.results
  .interaction.terms <- filter(.results, grepl(":", term))

  if (nrow(.interaction.terms) == 0) return()

  .interaction.terms %>% filter(direction(estimate, 0)) %>%
    mutate(padj = p.adjust(anova.p)) %>%
    filter(padj <= padj.threshold)
}


# DESeq2 functions --------------------------------------------------------

MakeDESeqFromAgglomerated <- function(agg, tax.level="Genus") {
  agg <- group_by(agg, SampleID) %>% filter(sum(count) > 0) %>% droplevels()
  # browser()
  mat <- reshape2::dcast(
    agg,
    formula=sprintf("SampleID + StudyGroup + SampleType ~ %s", tax.level),
    value.var="count",
    fill = 0,
    fun.aggregate = sum)
  pd <- mat[, c(1:3)]
  mat <- mat[, -c(1:3)]
  mat <- as.matrix(t(mat))
  colnames(mat) <- pd$SampleID
  DESeqDataSetFromMatrix(mat, pd, ~ StudyGroup)
}

TestDESeq <- function(ds, fitType="local") {
  # Referencing McMurdie (http://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html)
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }

  geoMeans <- apply(counts(ds), 1, gm_mean)
  ds <- estimateSizeFactors(ds, geoMeans=geoMeans)
  DESeq(ds, fitType=fitType)
}


# Macroscopic testing functions -------------------------------------------

#' Test differences between study groups using the generalized linear mixed model
#' @param ranks.to.test a vector giving the taxonomic ranks to collapse and test
#' @param agg an appropriately filtered agglomerated data frame
#' @param md the metadata dataframe
#' @param sparsity omit testing taxa that appear in fewer that this pct of samples
RunGLMM <- function(ranks.to.test, agg, md, sparsity=0.1, grouping.var="OriginatingSampleID") {
  plyr::adply(ranks.to.test, .margins=1, function(tax.level) {
    .md <- md[, 1:which(colnames(md) == tax.level)+1]
    results <- TestMatchedControlMixedModel(
      agg, tax.level=tax.level, sparsity=0.1, grouping.var = grouping.var)
    sig.results <- ExtractSignificantResults(results, padj.threshold = 1) %>%
      select(-c(term, group)) %>%
      mutate(tax.level=tax.level,
             taxa=as.character(taxa))
    sig.results[[tax.level]] <- sig.results$taxa
    left_join(sig.results, .md)
  }, .id=NULL) %>% distinct(taxa, .keep_all=TRUE)
}

#' Test differences between study groups using DESeq2.
#' @param ranks.to.test a vector giving the taxonomic ranks to collapse and test
#' @param agg an appropriately filtered agglomerated data frame
#' @param md the metadata dataframe
RunDESeq2 <- function(ranks.to.test, agg, md, sarcoid.name="sarcoid") {
  plyr::adply(ranks.to.test, .margins=1, function(tax.level) {
    .md <- md[, 1:which(colnames(md) == tax.level)+1]
    ds <- MakeDESeqFromAgglomerated(agg, tax.level=tax.level)
    ds <- TestDESeq(ds)
    ds <- ds %>% results(name=sprintf("StudyGroup%s", sarcoid.name)) %>% data.frame
    ds$taxa <- as.character(rownames(ds))
    ds[[tax.level]] <- rownames(ds)
    ds <- ds %>% filter(stat > 0) %>%
      left_join(.md) %>%
      mutate(test.rank=tax.level)
  }, .id=NULL) %>% distinct(taxa, .keep_all=TRUE)
}


# Testing output functions ----------------------------------------

SummarizeResults <- function(results, p.cutoff=0.1, taxa.level="taxa") {
  out <- filter(results, padj < p.cutoff)
  nlineages <- nrow(out)
  if (nlineages > 0) {
    .out <- out[[taxa.level]]

  } else {
    return(NA)
  }
  return(out)
}

# Case-Control Plotting ---------------------------------------------------

#' Links matched samples and (environmental) controls via the grouping.col
#' by completing any missing cases.
#' @param agg agglomerated dataframe
#' @param s sample metadata
#' @param md otu metadata
#' @return the modified agg with a new column, PlotGroup, for line drawing
LinkMatchedControls <- function(agg, s, md, grouping.col) {
  if (nrow(agg) == 0) {
    warning("Empty input")
    return()
  }

  # Complete missing cases
  agg <- ungroup(agg) %>% select_(grouping.col, "SampleType", "otu", "count") %>%
    droplevels()

  agg <- complete_(agg, list(grouping.col, "otu", "SampleType"), fill=list(count=0))
  # Re-add metadata
  agg <- left_join(agg, s) %>%
    left_join(md)
  agg$PlotGroup <- paste(agg[[grouping.col]], agg$otu)
  return(agg)
}


