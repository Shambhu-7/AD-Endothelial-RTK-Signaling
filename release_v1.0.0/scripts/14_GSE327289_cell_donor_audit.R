suppressPackageStartupMessages({
  library(SeuratObject)
})

cat("=== SCRIPT 14: GSE327289 CELL AND DONOR AUDIT ===\n")

path <- "data/processed/GSE327289/GSE327289_Vascular_Cells.rds"

obj <- readRDS(path)
meta <- obj@meta.data

out <- "results/validation_GSE327289"

required <- c(
  "Donor", "Diagnosis", "MajorCellType",
  "CellTypeFinalAll", "Batch", "Age", "Sex",
  "Sample", "predicted_doublets_manual"
)

stopifnot(all(required %in% names(meta)))

cat("\n=== MAJOR CELL TYPES ===\n")
print(sort(table(meta$MajorCellType), decreasing=TRUE))

cat("\n=== DETAILED CELL TYPES ===\n")
print(sort(table(meta$CellTypeFinalAll), decreasing=TRUE))

cat("\n=== DIAGNOSIS DISTRIBUTION ===\n")
print(table(meta$Diagnosis, useNA="ifany"))

cat("\n=== UNIQUE DONORS BY DIAGNOSIS ===\n")

donors <- unique(
  meta[, c("Donor", "Diagnosis", "Batch", "Age", "Sex")]
)

print(table(donors$Diagnosis, useNA="ifany"))

cat("\n=== ENDOTHELIAL SUBTYPES ===\n")

endo <- meta[
  grepl(
    "endothelial|\\bEC\\b",
    meta$MajorCellType,
    ignore.case=TRUE
  ) |
  grepl(
    "endothelial|\\bEC\\b",
    meta$CellTypeFinalAll,
    ignore.case=TRUE
  ),
  ,
  drop=FALSE
]

cat("Endothelial cells:", nrow(endo), "\n")

print(
  sort(
    table(endo$CellTypeFinalAll),
    decreasing=TRUE
  )
)

cat("\n=== ENDOTHELIAL DONORS BY DIAGNOSIS ===\n")

endo_donors <- unique(
  endo[, c("Donor", "Diagnosis")]
)

print(table(endo_donors$Diagnosis, useNA="ifany"))

cat("\n=== ENDOTHELIAL CELLS PER DONOR ===\n")

donor_counts <- aggregate(
  rep(1, nrow(endo)),
  by=list(
    Donor=endo$Donor,
    Diagnosis=endo$Diagnosis,
    Batch=endo$Batch
  ),
  FUN=sum
)

names(donor_counts)[4] <- "n_cells"

donor_counts <- donor_counts[
  order(donor_counts$Diagnosis, donor_counts$Donor),
]

print(donor_counts, row.names=FALSE)

write.csv(
  donor_counts,
  file.path(out, "GSE327289_endothelial_donor_counts.csv"),
  row.names=FALSE
)

write.csv(
  as.data.frame(
    table(
      Diagnosis=endo$Diagnosis,
      CellType=endo$CellTypeFinalAll
    )
  ),
  file.path(out, "GSE327289_endothelial_subtypes.csv"),
  row.names=FALSE
)

cat("\n=== DOUBLET FLAG DISTRIBUTION ===\n")
print(table(
  meta$predicted_doublets_manual,
  useNA="ifany"
))

cat("\n=== RNA ASSAY LAYERS ===\n")
print(Layers(obj[["RNA"]]))

writeLines(
  capture.output(sessionInfo()),
  file.path(out, "14_sessionInfo.txt")
)

cat("\nSCRIPT 14 COMPLETE\n")
