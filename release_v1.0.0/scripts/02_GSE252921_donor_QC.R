suppressPackageStartupMessages(library(Matrix))

base <- "results/pseudobulk/"
out <- "results/qc/"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pb <- readRDS(paste0(
  base, "GSE252921_endothelial_pseudobulk_counts.rds"
))

meta <- read.csv(
  paste0(base, "GSE252921_endothelial_donor_metadata.csv"),
  stringsAsFactors = FALSE
)

stopifnot(
  identical(colnames(pb), meta$donor_ID),
  all(meta$library_size == as.numeric(colSums(pb))),
  all(meta$n_cells > 0)
)

meta$counts_per_cell <- meta$library_size / meta$n_cells
meta$log10_library_size <- log10(meta$library_size)
meta$low_cells_200 <- meta$n_cells < 200
meta$low_cells_500 <- meta$n_cells < 500

cat("\n=== DONOR QC ===\n")
print(meta, row.names = FALSE)

cat("\n=== GROUP SUMMARY ===\n")
for (group in c("AD", "CT")) {
  d <- meta[meta$disease == group, ]
  cat("\nGroup:", group, "\n")
  cat("Donors:", nrow(d), "\n")
  cat("Cells:", sum(d$n_cells), "\n")
  cat("Median cells:", median(d$n_cells), "\n")
  cat("Median library size:", median(d$library_size), "\n")
  cat("Median age:", median(d$age), "\n")
  cat("Donors below 200 cells:", sum(d$low_cells_200), "\n")
  cat("Donors below 500 cells:", sum(d$low_cells_500), "\n")
}

cat("\n=== AGE BY DISEASE ===\n")
print(tapply(meta$age, meta$disease, summary))

cat("\n=== LIBRARY SIZE RANGE ===\n")
print(range(meta$library_size))

cat("\n=== GENES DETECTED PER DONOR ===\n")
meta$genes_detected <- as.integer(colSums(pb > 0))
print(meta[, c("donor_ID", "disease",
               "n_cells", "library_size",
               "genes_detected")], row.names = FALSE)

write.csv(
  meta,
  paste0(out, "GSE252921_donor_QC.csv"),
  row.names = FALSE
)

pdf(paste0(out, "GSE252921_donor_QC.pdf"),
    width = 10, height = 7)

par(mfrow = c(2, 2))

boxplot(log10_library_size ~ disease, data = meta,
        main = "Library size by disease",
        ylab = "log10(total counts)")

boxplot(n_cells ~ disease, data = meta,
        main = "Cells per donor",
        ylab = "Number of cells")

boxplot(age ~ disease, data = meta,
        main = "Age by disease",
        ylab = "Age (years)")

plot(meta$n_cells, meta$library_size,
     pch = 19, log = "xy",
     xlab = "Cells per donor",
     ylab = "Library size",
     main = "Library size vs cell number")

dev.off()

cat("\nDONOR QC: COMPLETE\n")
