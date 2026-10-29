suppressPackageStartupMessages(library(Matrix))

base <- "data/processed/GSE252921/"
out <- "results/pseudobulk/"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_double_gzip_rds <- function(path) {
  con <- gzcon(gzfile(path, open = "rb"))
  tryCatch(readRDS(con), finally = close(con))
}

cat("Loading metadata...\n")
meta <- read_double_gzip_rds(
  paste0(base, "GSM8010381_metadata_endo.rds.gz")
)

cat("Loading sparse counts...\n")
counts <- read_double_gzip_rds(
  paste0(base, "GSM8010381_counts_endo.rds.gz")
)

stopifnot(
  inherits(counts, "dgCMatrix"),
  identical(colnames(counts), rownames(meta)),
  !anyNA(meta$donor_ID),
  !anyDuplicated(rownames(counts))
)

# Create one aggregation column per biological donor.
donors <- sort(unique(as.character(meta$donor_ID)))
donor_factor <- factor(meta$donor_ID, levels = donors)

membership <- sparseMatrix(
  i = seq_len(nrow(meta)),
  j = as.integer(donor_factor),
  x = rep(1, nrow(meta)),
  dims = c(nrow(meta), length(donors)),
  dimnames = list(colnames(counts), donors)
)

cat("Aggregating cells by donor...\n")
pb <- counts %*% membership
pb <- as(pb, "dgCMatrix")

# Construct one metadata row per donor.
fields <- c("type", "sex", "age", "batch")

for (field in fields) {
  n_values <- tapply(
    meta[[field]], meta$donor_ID,
    function(z) length(unique(z))
  )
  stopifnot(all(n_values == 1))
}

idx <- match(donors, meta$donor_ID)

donor_meta <- data.frame(
  donor_ID = donors,
  disease = as.character(meta$type[idx]),
  sex = as.character(meta$sex[idx]),
  age = as.numeric(as.character(meta$age[idx])),
  batch = as.character(meta$batch[idx]),
  n_cells = as.integer(table(donor_factor)),
  library_size = as.numeric(colSums(pb)),
  stringsAsFactors = FALSE
)

stopifnot(
  identical(colnames(pb), donor_meta$donor_ID),
  ncol(pb) == 20L,
  sum(donor_meta$n_cells) == ncol(counts),
  isTRUE(all.equal(
    sum(pb),
    sum(counts),
    tolerance = 0
  ))
)

saveRDS(
  pb,
  paste0(out, "GSE252921_endothelial_pseudobulk_counts.rds")
)

write.csv(
  donor_meta,
  paste0(out, "GSE252921_endothelial_donor_metadata.csv"),
  row.names = FALSE
)

cat("\n=== PSEUDOBULK AUDIT ===\n")
cat("Genes:", nrow(pb), "\n")
cat("Donors:", ncol(pb), "\n")
cat("Total cells:", sum(donor_meta$n_cells), "\n")
cat("Total counts conserved: TRUE\n")

cat("\n=== DONORS BY DISEASE ===\n")
print(table(donor_meta$disease))

cat("\n=== DONOR QC SUMMARY ===\n")
print(donor_meta, row.names = FALSE)

cat("\n=== OUTPUT FILES ===\n")
print(list.files(out))

cat("\nPSEUDOBULK CONSTRUCTION: PASS\n")
