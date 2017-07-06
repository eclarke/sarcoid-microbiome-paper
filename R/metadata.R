#' Standardize a metadata table (s) to be compatible with SRA_data requirements
#' @param ds a SarcoidDataset object created by \code{\link{LoadData}}
#' @export
StandardizeMetadata <- function(ds) {

  stopifnot(inherits(ds, "SarcoidDataset"))
  print(sprintf("%s/%s", ds$sampleset, ds$datatype))
  s <- ds$s
  new.s <- data.frame(library_ID=s$SampleID)

  # title: short description
  describe <- function(study.group, sampleset.type, sample.type) {
    description <- ""
    if (!study.group %in% c("sarcoidosis", "healthy"))
      description <- "experimental control"
    else
      description <- sprintf(
        "%s %s (%s)", study.group, sampleset.type, sub("_", " ", sample.type))
  }
  sampleset.type <- list(
    A="FFPE tissue", B="FFPE tissue", C="BAL", DE="Kveim/spleen")[[ds$sampleset]]
  new.s$title <- paste(
    ds$datatype, "sequencing of",
    mapply(describe, s$StudyGroup, sampleset.type, s$SampleType))

  new.s$sample_name <- make.names(paste(s$SubjectID, s$SampleType, ds$datatype, sep="."), unique=TRUE)

  # library_strategy
  new.s$library_strategy <- list(
    "16S" = "AMPLICON",
    "ITS" = "AMPLICON",
    "virome" = "OTHER",
    "WGS" = "WGS")[[ds$datatype]]

  # library_source
  new.s$library_source <- "METAGENOMIC"
  if (ds$datatype == "virome")
    new.s$library_source <- ifelse(
      s$LibraryType == "cDNA", "VIRAL RNA", "METAGENOMIC")

  # library_selection
  new.s$library_selection <- list(
    "16S" = "PCR",
    "ITS" = "PCR",
    "WGS" = "RANDOM")[[ds$datatype]]
  if (ds$datatype == "virome")
    new.s$library_selection <- ifelse(
      s$LibraryType == "Genomiphi", "MDA", "RANDOM"
    )

  # library_layout
  new.s$library_layout <- "paired"

  # platform
  new.s$platform <- "ILLUMINA"

  # instrument_model
  new.s$instrument_model <- list(
    "16S" = "Illumina MiSeq",
    "ITS" = "Illumina MiSeq",
    "virome" = "Illumina HiSeq 2500",
    "WGS" = "Illumina HiSeq 2500")[[ds$datatype]]

  new.s$design_description <- list(
    "16S" = "PCR amplification of bacterial 16S rRNA, Nextera chemistry",
    "ITS" = "PCR amplification of bacterial 16S rRNA, Nextera chemistry",
    "WGS" = "Random fragmentation of genomic DNA, Nextera chemistry")[[ds$datatype]]
  if (ds$datatype == "virome")
    new.s$design_description <- ifelse(
      s$LibraryType == "Genomiphi",
      "Viral purification, followed by Genomiphi amplification of gDNA and Nextera chemistry",
      "Viral purification, followed by reverse transcription of viral RNA to cDNA and Nextera chemistry"
    )
  new.s$filetype <- "fastq"
  new.s$SampleSet <- ds$sampleset
  new.s$SubjectID <- s$SubjectID
  new.s$SampleID <- s$SampleID

  .nrow <- nrow(new.s)
  new.s <- ds$agg %>% group_by(SampleID) %>% summarize(readcount=sum(count)) %>%
    rename(library_ID=SampleID) %>%
    full_join(new.s) %>%
    select(-readcount, everything())
  stopifnot(.nrow == nrow(new.s))

  new.s
}

MakeSRAMetadata <- function(ds) {
  stopifnot(inherits(ds, "SarcoidDataset"))
  print(sprintf("%s/%s", ds$sampleset, ds$datatype))

  metadata <- StandardizeMetadata(ds)
  metadata$bioproject_accession <- "PRJNA392272"

  stopifnot(!anyDuplicated(metadata$sample_name))
  metadata$filetype <- "fastq"
  metadata$filename <- sprintf("%s_1.fastq.gz", metadata$SampleID)
  metadata$filename2 <- sprintf("%s_2.fastq.gz", metadata$SampleID)

  metadata %>%
    select(
      bioproject_accession, sample_name, library_ID, title, library_strategy,
      library_source, library_selection, library_layout, platform,
      instrument_model, design_description, filetype, filename, filename2)
}

MakeBioSample <- function(ds) {
  stopifnot(inherits(ds, "SarcoidDataset"))
  print(sprintf("%s/%s", ds$sampleset, ds$datatype))

  canonical_materials <- list(
    "paraffin" = "paraffin wax",
    "lymph_node" = "lymphoid tissue",
    "BAL" = "lung",
    "prewash" = "medical instrument",
    "extr_blank" = "water",
    "saline" = "medical instrument",
    "Saline" = "medical instrument",
    "pcr_h2o" = "water",
    "kveim" = "spleen",
    "Kveim" = "spleen",
    "Spleen" = "spleen",
    "Blank" = "water",
    "water" = "water",
    "spleen" = "spleen"
  )
  .s <- ds$s %>% filter(SampleType %in% names(canonical_materials))

  s <- data.frame(sample_name = make.names(paste(.s$SubjectID, .s$SampleType, ds$datatype, sep="."), unique=TRUE))
  stopifnot(!anyDuplicated(s$sample_name))

  s$bioproject_accession <- "PRJNA392272"
  s$organism <- "metagenome"
  s$collection_date <- "not collected"
  s$env_biome <- "not applicable"
  s$env_feature <- "not applicable"
  s$env_material <- sapply(as.character(.s$SampleType), function(x) canonical_materials[[x]])

  s$geo_loc_name <- list(
    "A" = "Poland",
    "B" = "USA",
    "C" = "USA",
    "DE" = "USA"
  )[[ds$sampleset]]
  s$geo_loc_name <- ifelse(s$env_material == "water", "USA", s$geo_loc_name)
  s$host <- ifelse(s$env_material %in% c("lung", "lymphoid_tissue", "spleen"), "Homo sapiens", "not applicable")
  s$lat_lon <- "not collected"
  s$source_material_id <- .s$SubjectID
  s$description <- .s$SampleID
  if ('LibraryType' %in% colnames(.s)) {
    s$description <- paste(s$description, .s$LibraryType, sep = " - ")
  }
  s <- s %>% group_by(source_material_id, env_material) %>%
    mutate(replicate = seq_along(source_material_id))
  return(s)
}
