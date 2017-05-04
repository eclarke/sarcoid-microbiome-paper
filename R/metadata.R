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
  new.s$library_layout <- "Paired end"

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

  .nrow <- nrow(new.s)
  new.s <- ds$agg %>% group_by(SampleID) %>% summarize(readcount=sum(count)) %>%
    rename(library_ID=SampleID) %>%
    full_join(new.s) %>%
    select(-readcount, everything())
  stopifnot(.nrow == nrow(new.s))

  new.s
}
