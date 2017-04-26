#' Load datasets.
#' @param set the sample set in question (of A-C, DE)
#' @param data the data type (of 16S, ITS, virome, or WGS)
#' @param data_fp path to data files
#' @import dplyr
#' @export
LoadData <- function(
  set=c("A", "B", "C", "DE"),
  data=c("16S", "ITS", "virome", "WGS"),
  data_fp=opts$get("data_fp"))
{
  stopifnot(!is.null(data_fp))

  set <- match.arg(set)
  data <- match.arg(data)

  datasets <- list(
    A=list("16S"=.set_a_16s, "ITS"=.set_a_its),
    B=list("16S"=.set_b_16s, "ITS"=.set_b_its),
    C=list("16S"=.set_c_16s, "ITS"=.set_c_its, "virome"=.set_c_vir),
    DE=list("16S"=.set_de_16s, "ITS"=.set_de_its, "WGS"=.set_de_wgs))

  if (!data %in% names(datasets[[set]]))
    stop("Requested data not available for that set")

  suppressMessages(
    suppressWarnings(
      datasets[[set]][[data]](data_fp)))

}

# Helper functions --------------------------------------------------------
.relpath <- function(path, opts) {
  file.path(opts$data_fp, path)
}

.merge_studygroup <- function(x) {
  x <- as.character(x)
  if (any(x == "cut_ctrl")) {
    non.ctrl <- x[which(x != "cut_ctrl")]
    rep(non.ctrl, length(x))
  } else {
    x
  }
}

.set_ab_combined_its <- function(data_fp) {

  .load_its_data <- function(
    otu_table, mapping_file, tree_file, tax_file, map_processing_fun, root)
  {

    o <- qiimer::read_qiime_otu_table(otu_table)

    cts <- o$counts

    s <- qiimer::read_qiime_mapping_file(mapping_file)
    s <- map_processing_fun(s, root)

    tree <- ape::read.tree(tree_file)
    tree <- ape::root(tree, sample(tree$tip.label, 1), resolve.root=TRUE)

    md2 <- read.table(tax_file, sep="\t")
    md2 <- md2 %>% tidyr::separate(
      V2, into=c("Domain", qiimer::taxonomic_ranks),
      sep=";", extra="drop", fill="right")
    md2 <- select(md2, -Domain)
    md2 <- mutate(md2, Kingdom = ifelse(is.na(Kingdom), "Indeterminate", Kingdom))
    md2 <- eclectic::reorder_taxa(md2)
    rownames(md2) <- md2$V1
    md2 <- dplyr::rename(md2, "otu"=V1)
    md2 <- md2[, c(qiimer::taxonomic_ranks, "otu")]
    md2$MinRank1 <- eclectic::tax_climber(md2$otu, md2, end="Genus", label=TRUE, sep="__")
    md2$MinRank <- eclectic::tax_climber(md2$otu, md2, end="Species", label=TRUE, sep="__")
    md = md2

    list(s=s, md=md, cts=cts, tree=tree)
  }

  # Modifies the cut controls to be part of their paired study groups
  # Adds in building, date info
  .pol.map.fun <- function(s, root) {
    s <- group_by(s, OriginatingSampleID) %>%
      mutate(StudyGroup = .merge_studygroup(StudyGroup)) %>%
      ungroup()
    s2.Healthy <- readxl::read_excel(file.path(root, "extra_metadata.xlsx")) %>%
      filter(!is.na(TubeLabel))
    s2.Sarcoid <- readxl::read_excel(file.path(root, "extra_metadata.xlsx"), sheet = "Sarcoid (Box 2)") %>%
      filter(!is.na(TubeLabel))
    s2 <- rbind(s2.Healthy, s2.Sarcoid) %>%
      mutate(
        Building = factor(ifelse(!is.na(Building), paste("Building", Building), NA))) %>%
      mutate(
        OriginatingSampleID =  make.names(tolower(SampleID))) %>%
      mutate(
        OriginatingSampleID = sub("4819", "4810", OriginatingSampleID)) %>%
      select(OriginatingSampleID, Building, Date) %>% distinct(.keep_all=TRUE) %>%
      right_join(s)

    s <- s2
  }

  # Fixes sample names
  .hup.map.fun <- function(s, root) {
    s %>% mutate(SampleID = gsub("-|_", ".", SampleID)) %>%
      dplyr::rename("StudyGroup"=Study_Group,
                    "SampleType"=sample_type)
  }

  .hup.root <- file.path(data_fp, "B_tissue/ITS")
  .pol.root <- file.path(data_fp, "A_tissue/ITS")
  .hup.map <- file.path(.hup.root, "samples.txt")
  .pol.map <- file.path(.pol.root, "samples.txt")
  pol <- .load_its_data(
    otu_table = file.path(.pol.root, "counts.txt"),
    mapping_file = .pol.map,
    tree_file = file.path(.pol.root, "tree.tre"),
    tax_file = file.path(.pol.root, "taxonomy.txt"),
    map_processing_fun = .pol.map.fun,
    root= .pol.root)

  hup <- .load_its_data(
    otu_table = file.path(.hup.root, "counts.txt"),
    mapping_file = .hup.map,
    tree_file = file.path(.hup.root, "tree.tre"),
    tax_file = file.path(.hup.root, "taxonomy.txt"),
    map_processing_fun = .hup.map.fun,
    root= .hup.root)

  hup.s <- hup$s
  levels(hup.s$StudyGroup) = c("extr_ctrl", "healthy", "sarcoidosis")
  levels(hup.s$SampleType) = c("extr_blank", "kveim", "lymph_node", "pcr_h2o", "saline", "spleen")

  hup.s <- hup.s %>%
    mutate_each("as.character", StudyGroup, SampleType)
  pol.s <- pol$s %>%
    mutate_each("as.character", StudyGroup, SampleType)

  hup.s$Site = "HUP"
  pol.s$Site = "Poland"
  s <- plyr::rbind.fill(hup.s, pol.s)
  # There should not be any merging of rows!
  stopifnot(nrow(s) == nrow(hup.s) + nrow(pol.s))
  stopifnot(n_distinct(s$SampleID) == nrow(s))

  # Drop irrelevant columns and convert to factors where necessary
  s <- select(s, -c(
    BarcodeSequence, LinkerPrimerSequence, ReversePrimer:RevBarcode,
    CrudeDNAConc, Extraction_Method, Description)) %>%
    mutate_each("as.factor", StudyGroup, SampleType, Site)

  # Agglomerate
  agg <- eclectic::agglomerate(s, hup$cts, hup$md)

  # Re-add samples with no reads
  agg <- left_join(s, agg)
  stopifnot(all(s$SampleID %in% agg$SampleID))
  agg$count[is.na(agg$count)] <- 0

  agg <- agg %>%
    group_by(SampleID) %>%
    mutate(
      proportion = count/sum(count),
      total = sum(count))

  list(s=s, cts=hup$cts, md=hup$md, agg=agg, tree=hup$tree)
}

.set_a_16s <- function(data_fp) {

  root <- file.path(data_fp, "A_tissue/16S")

  o <- qiimer::read_qiime_otu_table(file.path(root, "counts.txt"))

  cts <- o$counts

  s <- qiimer::read_qiime_mapping_file(file.path(root, "samples.txt"))
  s <- s %>% group_by(OriginatingSampleID) %>%
    mutate(StudyGroup = .merge_studygroup(StudyGroup)) %>%
    ungroup()

  s.healthy <- readxl::read_excel(file.path(root, "extra_metadata.xlsx")) %>%
    filter(!is.na(TubeLabel))
  s.sarcoid <- readxl::read_excel(
    file.path(root, "extra_metadata.xlsx"), sheet = "Sarcoid (Box 2)") %>%
    filter(!is.na(TubeLabel))
  s <- rbind(s.healthy, s.sarcoid) %>%
    mutate(
      Building = factor(ifelse(!is.na(Building), paste("Building", Building), NA))) %>%
    mutate(
      OriginatingSampleID = make.names(tolower(SampleID))) %>%
    mutate(
      OriginatingSampleID = sub("4819", "4810", OriginatingSampleID)) %>%
    select(OriginatingSampleID, Building, Date) %>% distinct(.keep_all=TRUE) %>%
    right_join(s)

  tree <- ape::read.tree(file.path(root, "tree.tre"))
  tree <- ape::root(tree, sample(tree$tip.label, 1), resolve.root=TRUE)


  md <- as.data.frame(o$metadata)
  md$otu <- rownames(md)
  colnames(md)[1] <- "taxonomy"
  md <- md %>% tidyr::separate(
    taxonomy,
    into=qiimer::taxonomic_ranks,
    sep="; ",
    extra = "drop", fill="right") %>%
    mutate_each(funs(sub("[kpcofgs]__", "", .)), Kingdom:Species) %>%
    mutate_each(funs(ifelse(. == "", NA, .)), Kingdom:Species)
  md <- md %>%
    mutate(
      Species = ifelse(
        !is.na(Genus) & !is.na(Species),
        paste(as.character(Genus), as.character(Species)),
        Species))
  md <- eclectic::reorder_taxa(md)
  rownames(md) <- md$otu
  md$MinRank <- eclectic::tax_climber(
    md$otu, md, end='Species', label = TRUE, sep="__")

  agg <- eclectic::agglomerate(s, cts, md)

  # Re-add samples with no reads
  agg <- left_join(s, agg)
  stopifnot(all(s$SampleID %in% agg$SampleID))
  agg$count[is.na(agg$count)] <- 0

  agg <- agg %>%
    filter(Order != "Streptophyta") %>%
    group_by(SampleID) %>%
    mutate(
      total=sum(count),
      proportion=count/total) %>%
    ungroup() %>%
    droplevels()

  list(agg, s, md, tree, cts)
}

.set_a_its <- function(data_fp) {
  full <- .set_ab_combined_its(data_fp)
  within(full, {
    agg <- filter(agg, Site=="Poland")
    s <- filter(s, Site=="Poland")
  })
}

.set_b_16s <- function(data_fp) {
  root <- file.path(data_fp, "B_tissue/16S")

  o <- qiimer::read_qiime_otu_table(file.path(
    root, "counts.txt"))

  s <- read.delim(file.path(root, "samples.txt")) %>%
    filter(ExtractionType=="PowerSoil", Plate=="Plate2") %>%
    droplevels()

  # There are replicates; pick the ones with higher amplicon conc.
  s <- s %>% group_by(DerivingOriginalSampleID) %>%
    filter(PostPCRConc == max(PostPCRConc)) %>%
    mutate(SampleType = forcats::fct_recode(SampleType, "lymph_node"="tissue"),
           StudyGroup = forcats::fct_recode(StudyGroup, "sarcoidosis"="sarcoid"))

  tree <- ape::read.tree(file.path(root, "tree.tre"))
  tree <- ape::root(tree, sample(tree$tip.label, 1), resolve.root=TRUE)

  md <- as.data.frame(o$metadata)
  md$otu <- rownames(md)
  colnames(md)[1] <- "taxonomy"
  md <- md %>%
    tidyr::separate(
      taxonomy,
      into=qiimer::taxonomic_ranks, sep="; ",
      extra = "drop", fill="right") %>%
    mutate_each(funs(sub("[kpcofgs]__", "", .)), Kingdom:Species) %>%
    mutate_each(funs(ifelse(. == "", NA, .)), Kingdom:Species)
  md <- md %>%
    mutate(
      Species = ifelse(!is.na(Genus) & !is.na(Species),
                       paste(as.character(Genus), as.character(Species)),
                       Species))
  md <- eclectic::reorder_taxa(md)
  rownames(md) <- md$otu
  md$MinRank <- eclectic::tax_climber(
    md$otu, md, end="Species", label = TRUE, sep="__")

  agg <- eclectic::agglomerate(s, o$counts, md)

  # Re-add samples with no reads
  agg <- left_join(s, agg)
  stopifnot(all(s$SampleID %in% agg$SampleID))
  agg$count[is.na(agg$count)] <- 0

  agg <- agg %>%
    filter(Order != "Streptophyta") %>%
    group_by(SampleID) %>%
    mutate(
      total=sum(count),
      proportion=count/total) %>%
    ungroup() %>%
    droplevels()

  list(agg, s, md, tree, cts=o$counts)
}

.set_b_its <- function(data_fp) {
  full <- .set_ab_combined_its(data_fp)
  within(full, {
    agg <- filter(agg, Site=="HUP")
    s <- filter(s, Site=="HUP")
  })
}

.set_c_16s <- function(data_fp) {
  root <- file.path(data_fp, "C_bal/16S")

  # Load count data
  o <- qiimer::read_qiime_otu_table(file.path(root, "counts.txt"))

  # Load sample metadata
  s <- read.delim(file.path(root, "samples.txt")) %>%
    filter(ExtractionType=="PowerSoil", Plate=="Plate4") %>%
    mutate(StudyGroup = forcats::fct_recode(StudyGroup, "sarcoidosis"="sarcoid"))

  # Load phylogenetic tree
  tree <- ape::read.tree(file.path(root, "tree.tre"))
  tree <- ape::root(tree, sample(tree$tip.label, 1), resolve.root=TRUE)

  # Load taxonomic data
  md <- as.data.frame(o$metadata)
  md$otu <- rownames(md)
  colnames(md)[1] <- "taxonomy"
  md <- md %>%
    tidyr::separate(
      taxonomy,
      into=qiimer::taxonomic_ranks, sep="; ",
      extra = "drop", fill="right") %>%
    mutate_each(funs(sub("[kpcofgs]__", "", .)), Kingdom:Species) %>%
    mutate_each(funs(ifelse(. == "", NA, .)), Kingdom:Species)
  md <- md %>%
    mutate(
      Species = ifelse(!is.na(Genus) & !is.na(Species),
                       paste(as.character(Genus), as.character(Species)),
                       Species))
  md <- eclectic::reorder_taxa(md)
  rownames(md) <- md$otu
  md$MinRank <- eclectic::tax_climber(
    md$otu, md, end="Genus", label = TRUE, sep="__")

  # Agglomerate sample data, counts, and taxonomy into melted df
  agg <- eclectic::agglomerate(s, o$counts, md)

  # Re-add samples with no reads
  agg <- left_join(s, agg)
  stopifnot(all(s$SampleID %in% agg$SampleID))
  agg$count[is.na(agg$count)] <- 0

  agg <- agg %>%
    filter(Order != "Streptophyta") %>%
    group_by(SampleID) %>%
    mutate(
      total=sum(count),
      proportion=count/total) %>%
    ungroup()

  list(s, agg, md, cts=o$counts, tree)
}

.set_c_its <- function(data_fp) {
  # browser()
  root <- file.path(data_fp, "C_bal/ITS")

  o <- qiimer::read_qiime_otu_table(file.path(root, "counts.txt"))

  cts <- o$counts
  colnames(cts) <- sub("^X", "", colnames(cts))

  s <- qiimer::read_qiime_mapping_file(file.path(root, "samples.txt"))

  s <- filter(s, Plate == "Plate3") %>%
    droplevels() %>%
    mutate(StudyGroup = forcats::fct_recode(StudyGroup, "sarcoidosis"="sarcoid"))

  md2 <- read.table(file.path(root, "taxonomy.txt"), sep="\t")
  md2 <- md2 %>% tidyr::separate(
    V2, into=c("Domain", qiimer::taxonomic_ranks),
    sep=";", extra="drop", fill="right")
  md2 <- select(md2, -Domain)
  md2 <- mutate(md2, Kingdom = ifelse(is.na(Kingdom), "Indeterminate", Kingdom))
  md2 <- eclectic::reorder_taxa(md2)
  rownames(md2) <- md2$V1
  md2 <- dplyr::rename(md2, "otu"=V1)
  md2 <- md2[, c(qiimer::taxonomic_ranks, "otu")]
  md2$MinRank1 <- eclectic::tax_climber(md2$otu, md2, end="Genus", label=TRUE, sep="__")
  md2$MinRank <- eclectic::tax_climber(md2$otu, md2, end="Species", label=TRUE, sep="__")
  md = md2

  agg <- eclectic::agglomerate(s, cts, md)

  # Re-add samples with no reads
  agg <- left_join(s, agg)
  stopifnot(all(s$SampleID %in% agg$SampleID))
  agg$count[is.na(agg$count)] <- 0

  agg <- agg %>%
    group_by(SampleID) %>%
    mutate(
      total=sum(count),
      proportion=count/total) %>%
    ungroup()

  list(s, cts, agg, md)
}

.set_c_vir <- function(data_fp) {
  root <- file.path(data_fp, "C_bal/virome")
  # o <- qiimer::read_qiime_otu_table(file.path(root, "merged_table.tsv"))
  o <- qiimer::read_qiime_otu_table(file.path(root, "counts.txt"))
  s <- read.delim(file.path(root, "samples.txt"))
  cts <- o$counts
  colnames(cts) <- stringr::str_replace(colnames(cts), "\\.taxa", "")

  baltimore <- read.table(file.path(root, "lineages.txt"), sep="\t", header = TRUE) %>% select(tax_id, baltimore) %>%
    mutate(otu=as.factor(tax_id)) %>%
    select(-tax_id)

  s <- s %>% mutate(
    StudyGroup = ifelse(StudyGroup=="control", as.character(SampleType), as.character(StudyGroup)),
    SampleType = forcats::fct_recode(SampleType, "water"="neg_ctrl", "spike"="spike_ctrl"))

  md <- as.data.frame(o$metadata)
  md$otu <- rownames(md)
  colnames(md)[1] <- "taxonomy"
  md <- md %>%
    tidyr::separate(
      taxonomy,
      into=qiimer::taxonomic_ranks, sep="; ",
      extra = "drop", fill="right") %>%
    mutate_each(funs(sub("[kpcofgs]__", "", .)), Kingdom:Species) %>%
    mutate_each(funs(ifelse(. == "", NA, .)), Kingdom:Species)
  # Fix stupidly named species
  md <- md %>% mutate(
    Species = as.character(Species)) %>%
    mutate(
      Species = ifelse(
        grepl("^[0-9]+$", Species),
        paste(as.character(Genus), Species),
        Species)) %>%
    mutate(Species = as.factor(Species))
  # Fix some papillomaviruses having caudovirales for orders
  md <- md %>% mutate(
    Order = ifelse(
      Family == "Papillomaviridae",
      NA, Order))

  md <- eclectic::reorder_taxa(md)
  rownames(md) <- md$otu
  md$MinRank1 <- eclectic::tax_climber(
    md$otu, md, end="Genus", label = TRUE, sep="__")
  md$MinRank <- eclectic::tax_climber(
    md$otu, md, end="Species", label=FALSE
  )
  md <- left_join(md, baltimore)

  agg <- eclectic::agglomerate(s, cts, md)

  agg <- agg %>%
    filter(Kingdom=="Viruses") %>%
    group_by(SampleID) %>%
    mutate(
      total=sum(count),
      proportion=count/total) %>%
    ungroup() %>%
    droplevels()

  list(s, agg, md, cts)
}

.set_de_16s <- function(data_fp) {
  root <- file.path(data_fp, "DE_kveimspleen/16S")

  # Load count data
  o <- qiimer::read_qiime_otu_table(file.path(
    root, "counts.txt"))

  # Load sample metadata
  s <- read.delim(file.path(
    root, "samples.txt")) %>%
    rename(SampleID=X.SampleID) %>%
    filter(!Sample %in% c("PCR_H2O_1", "PCR_H2O_2"))

  s <- s %>% mutate(SampleType = forcats::fct_recode(
    as.factor(Sample), Spleen="S1", Spleen="S2", Spleen="S3", Kveim="KV", Blank="B", Saline="Sal"))
  s <- s %>% mutate(StudyGroup = forcats::fct_recode(
    SampleType, Sarcoidosis="Spleen", Sarcoidosis="Kveim", Control="Blank", Control="Saline"
  ))

  # Load phylogenetic tree
  tree <- ape::read.tree(file.path(
    root, "tree.tre"))
  tree <- ape::root(tree, sample(tree$tip.label, 1), resolve.root=TRUE)

  # Load taxonomic data
  md <- as.data.frame(o$metadata)
  md$otu <- rownames(md)
  colnames(md)[1] <- "taxonomy"
  md <- md %>%
    tidyr::separate(
      taxonomy,
      into=qiimer::taxonomic_ranks, sep="; ",
      extra = "drop", fill="right") %>%
    mutate_each(funs(sub("[kpcofgs]__", "", .)), Kingdom:Species) %>%
    mutate_each(funs(ifelse(. == "", NA, .)), Kingdom:Species)
  md <- md %>%
    mutate(
      Species = ifelse(!is.na(Genus) & !is.na(Species),
                       paste(as.character(Genus), as.character(Species)),
                       Species))
  md <- eclectic::reorder_taxa(md)
  rownames(md) <- md$otu
  md$MinRank <- eclectic::tax_climber(
    md$otu, md, end="Species", label = TRUE, sep="__")

  # Remove otus that appear only one time
  cts <- o$counts
  cts <- cts[rowSums(cts) > 1, ]


  # Agglomerate sample data, counts, and taxonomy into melted df
  agg <- eclectic::agglomerate(s, cts, md)

  # Re-add samples with no reads
  agg <- left_join(s, agg)
  stopifnot(all(s$SampleID %in% agg$SampleID))
  agg$count[is.na(agg$count)] <- 0

  agg <- agg %>%
    filter(Order != "Streptophyta") %>%
    group_by(SampleID) %>%
    mutate(
      total=sum(count),
      proportion=count/total) %>%
    ungroup()

  list(s, agg, md, cts)
}

.set_de_its <- function(data_fp) {
  with(.set_b_its(data_fp), {
    rename.samples <- list(
      Spleen1="spleen1",
      Spleen2="spleen2",
      Spleen3="spleen3",
      Kveim="Kveim",
      Blank="spleen.ext.blank",
      Saline="saline"
    )
    agg <- agg %>%
      filter(SampleID %in% c("spleen1", "spleen2", "spleen3", "Kveim", "saline", "spleen.ext.blank")) %>%
      ungroup()
    agg <- mutate(
      agg,
      StudyGroup = forcats::fct_recode(StudyGroup, Sarcoidosis="sarcoidosis", Control="extr_ctrl")
    )
    # browser()
    s <- s %>% filter(SampleID %in% agg$SampleID)
    agg$SampleID <- do.call(function(...) forcats::fct_recode(factor(agg$SampleID), ...), rename.samples)
    s$SampleID <- do.call(function(...) forcats::fct_recode(factor(s$SampleID), ...), rename.samples)
    md <- md %>% filter(otu %in% agg$otu)

    list(s, agg, md, cts)
  })
}

.set_de_wgs <- function(data_fp) {
  root <- file.path(data_fp, "DE_kveimspleen/wgs")

  # Load count data
  # o <- qiimer::read_qiime_otu_table(file.path(
  #   root, "sunbeam_output/classify/kraken/all_samples.tsv"))
  o <- qiimer::read_qiime_otu_table(file.path(
    root, "counts.txt"))

  md <- as.data.frame(o$metadata)
  md$otu <- rownames(md)
  colnames(md)[1] <- "taxonomy"
  md <- md %>%
    tidyr::separate(
      taxonomy,
      into=qiimer::taxonomic_ranks, sep="; ",
      extra = "drop", fill="right") %>%
    mutate_each(funs(sub("[kpcofgs]__", "", .)), Kingdom:Species) %>%
    mutate_each(funs(ifelse(. == "", NA, .)), Kingdom:Species)
  md <- md %>%
    mutate(
      Species = ifelse(!is.na(Genus) & !is.na(Species),
                       paste(as.character(Genus), as.character(Species)),
                       Species))
  md <- eclectic::reorder_taxa(md)
  rownames(md) <- md$otu
  md$MinRank <- eclectic::tax_climber(
    md$otu, md, end="Species", label = TRUE, sep="__")

  # Remove otus that appear only one time
  cts <- o$counts
  cts <- cts[rowSums(cts) > 1, ]

  s <- data.frame(SampleID = colnames(cts))
  s <- s %>% mutate(SampleID = gsub("Sarc_|_gDNA.taxa", "", SampleID))
  s <- s %>% mutate(SampleType = forcats::fct_recode(
    as.factor(SampleID), Spleen="S1", Spleen="S2", Spleen="S3", Kveim="KV", Blank="B", Saline="Sal"))
  s <- s %>% mutate(StudyGroup = forcats::fct_recode(
    SampleType, Sarcoidosis="Spleen", Sarcoidosis="Kveim", Control="Blank", Control="Saline"
  ))
  colnames(cts) <- s$SampleID

  # Agglomerate sample data, counts, and taxonomy into melted df
  agg <- eclectic::agglomerate(s, cts, md)

  # Re-add samples with no reads
  agg <- left_join(s, agg)
  stopifnot(all(s$SampleID %in% agg$SampleID))
  agg$count[is.na(agg$count)] <- 0

  agg <- agg %>%
    group_by(SampleID) %>%
    mutate(
      total=sum(count),
      proportion=count/total) %>%
    ungroup()

  list(s, agg, md, cts)
}
