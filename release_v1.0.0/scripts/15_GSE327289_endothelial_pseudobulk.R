suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

cat("=== SCRIPT 15: GSE327289 ENDOTHELIAL PSEUDOBULK ===\n")

input <- "data/processed/GSE327289/GSE327289_Vascular_Cells.rds"
out <- "results/validation_GSE327289"

obj <- readRDS(input)
meta <- obj@meta.data

stopifnot(all(c(
  "Donor", "Diagnosis", "Batch", "MajorCellType",
  "CellTypeFinalAll", "Age", "Sex"
) %in% colnames(meta)))

keep <- meta$Batch == "AD Progression" &
        meta$Diagnosis %in% c("AD", "CTRL") &
        meta$MajorCellType == "Endothelial"

keep[is.na(keep)] <- FALSE

selected <- rownames(meta)[keep]
m <- meta[selected, , drop=FALSE]

cat("\nSelected endothelial cells:", nrow(m), "\n")

cat("\nDiagnosis distribution:\n")
print(table(m$Diagnosis))

cat("\nEndothelial subtypes:\n")
print(table(m$CellTypeFinalAll))

m$Donor <- as.character(m$Donor)
m$Diagnosis <- as.character(m$Diagnosis)

donor_ids <- sort(unique(m$Donor))

cat("\nNumber of donors:", length(donor_ids), "\n")

# Confirm that each donor has a single diagnosis and consistent covariates.
for (d in donor_ids) {
  x <- m[m$Donor == d, , drop=FALSE]

  stopifnot(length(unique(x$Diagnosis)) == 1)
  stopifnot(length(unique(as.character(x$Age))) == 1)
  stopifnot(length(unique(as.character(x$Sex))) == 1)
}

cat("\nLoading raw RNA counts...\n")

counts <- LayerData(
  obj,
  assay="RNA",
  layer="counts"
)

stopifnot(all(selected %in% colnames(counts)))

counts <- counts[, selected, drop=FALSE]

stopifnot(inherits(counts, "sparseMatrix"))

cat("Count matrix dimensions:\n")
print(dim(counts))

# Sparse aggregation: genes x cells multiplied by cells x donors.
group_index <- match(m$Donor, donor_ids)

aggregation <- sparseMatrix(
  i=seq_len(nrow(m)),
  j=group_index,
  x=1,
  dims=c(nrow(m), length(donor_ids))
)

pb <- counts %*% aggregation

colnames(pb) <- donor_ids
rownames(pb) <- rownames(counts)

# Metadata: exactly one row per donor.
donor_meta <- do.call(
  rbind,
  lapply(donor_ids, function(d) {
    x <- m[m$Donor == d, , drop=FALSE]

    data.frame(
      Donor=d,
      Diagnosis=x$Diagnosis[1],
      Age=as.numeric(as.character(x$Age[1])),
      Sex=as.character(x$Sex[1]),
      n_cells=nrow(x),
      n_capillary=sum(x$CellTypeFinalAll == "Capillary EC"),
      n_venous=sum(x$CellTypeFinalAll == "Venous EC"),
      n_arterial=sum(x$CellTypeFinalAll == "Arterial EC"),
      stringsAsFactors=FALSE
    )
  })
)

rownames(donor_meta) <- donor_meta$Donor

stopifnot(identical(colnames(pb), rownames(donor_meta)))

stopifnot(
  all(donor_meta$n_cells ==
      donor_meta$n_capillary +
      donor_meta$n_venous +
      donor_meta$n_arterial)
)

stopifnot(
  all(Matrix::colSums(pb) > 0)
)

cat("\nDONOR DISTRIBUTION:\n")
print(table(donor_meta$Diagnosis))

cat("\nPSEUDOBULK DIMENSIONS:\n")
print(dim(pb))

cat("\nCELLS PER DONOR:\n")
print(donor_meta[, c("Donor", "Diagnosis", "n_cells")],
      row.names=FALSE)

cat("\nAGE AND SEX SUMMARY:\n")
print(donor_meta[, c("Donor", "Diagnosis", "Age", "Sex")],
      row.names=FALSE)

saveRDS(
  pb,
  file.path(out, "GSE327289_endothelial_pseudobulk_counts.rds")
)

write.csv(
  donor_meta,
  file.path(out, "GSE327289_endothelial_donor_metadata.csv"),
  row.names=FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out, "15_sessionInfo.txt")
)

cat("\nSCRIPT 15 COMPLETE\n")
