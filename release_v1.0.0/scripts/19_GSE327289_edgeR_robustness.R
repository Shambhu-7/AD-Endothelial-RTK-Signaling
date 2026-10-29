suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
})

cat("=== SCRIPT 19: GSE327289 edgeR ROBUSTNESS ===\n")

out <- "results/GSE327289_robustness"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

input <- "results/validation_GSE327289"

counts <- readRDS(file.path(
  input, "GSE327289_endothelial_pseudobulk_counts.rds"
))

meta <- read.csv(
  file.path(input, "GSE327289_endothelial_donor_metadata.csv"),
  stringsAsFactors=FALSE
)

meta$Donor <- as.character(meta$Donor)

stopifnot(
  !anyDuplicated(meta$Donor),
  !anyDuplicated(colnames(counts)),
  setequal(meta$Donor, colnames(counts))
)

meta <- meta[
  match(colnames(counts), meta$Donor),
  ,
  drop=FALSE
]

stopifnot(identical(meta$Donor, colnames(counts)))

meta$Diagnosis <- factor(meta$Diagnosis, levels=c("CTRL","AD"))
meta$Sex <- factor(meta$Sex)

meta$p_venous <- meta$n_venous / meta$n_cells
meta$p_arterial <- meta$n_arterial / meta$n_cells

stopifnot(!anyNA(meta[, c(
  "Diagnosis","Sex","Age","p_venous","p_arterial"
)]))

fit_model <- function(label, selected, composition=FALSE) {

  cat("\n====================================\n")
  cat("MODEL:", label, "\n")
  cat("====================================\n")

  m <- meta[selected, , drop=FALSE]
  x <- counts[, m$Donor, drop=FALSE]

  m$Diagnosis <- droplevels(m$Diagnosis)
  m$Sex <- droplevels(m$Sex)

  m$Age_z <- as.numeric(scale(m$Age))

  if (composition) {
    m$Venous_z <- as.numeric(scale(m$p_venous))
    m$Arterial_z <- as.numeric(scale(m$p_arterial))

    design <- model.matrix(
      ~ Diagnosis + Age_z + Sex + Venous_z + Arterial_z,
      data=m
    )
  } else {
    design <- model.matrix(
      ~ Diagnosis + Age_z + Sex,
      data=m
    )
  }

  stopifnot(
    !anyNA(design),
    qr(design)$rank == ncol(design)
  )

  cat("Donors:", nrow(m), "\n")
  print(table(m$Diagnosis))

  cat("Design rank:",
      qr(design)$rank, "/", ncol(design), "\n")

  cat("Condition number:", kappa(design), "\n")

  y <- DGEList(counts=x, group=m$Diagnosis)

  keep <- filterByExpr(y, design=design)

  cat("Genes retained:", sum(keep), "\n")

  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)

  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)

  test <- glmQLFTest(fit, coef="DiagnosisAD")

  res <- topTags(test, n=Inf, sort.by="none")$table
  res$gene <- rownames(res)

  res <- res[, c(
    "gene","logFC","logCPM","F","PValue","FDR"
  )]

  res <- res[order(res$PValue), ]

  write.csv(
    res,
    file.path(out, paste0(label, "_edgeR_all_genes.csv")),
    row.names=FALSE
  )

  cat("FDR < 0.05:",
      sum(res$FDR < 0.05, na.rm=TRUE), "\n")

  invisible(res)
}

primary <- fit_model(
  "primary_refit",
  rep(TRUE, nrow(meta))
)

lowcell <- fit_model(
  "exclude_under100",
  meta$n_cells >= 100
)

composition <- fit_model(
  "composition_adjusted",
  rep(TRUE, nrow(meta)),
  composition=TRUE
)

# Compare common genes across models.
common <- Reduce(
  intersect,
  list(primary$gene, lowcell$gene, composition$gene)
)

comparison <- data.frame(gene=sort(common))

for (label in c(
  "primary_refit",
  "exclude_under100",
  "composition_adjusted"
)) {

  tab <- get(
    switch(
      label,
      primary_refit="primary",
      exclude_under100="lowcell",
      composition_adjusted="composition"
    )
  )

  idx <- match(comparison$gene, tab$gene)

  comparison[[paste0("logFC_", label)]] <- tab$logFC[idx]
  comparison[[paste0("FDR_", label)]] <- tab$FDR[idx]
}

comparison$direction_consistent <- (
  sign(comparison$logFC_primary_refit) ==
  sign(comparison$logFC_exclude_under100) &
  sign(comparison$logFC_primary_refit) ==
  sign(comparison$logFC_composition_adjusted)
)

write.csv(
  comparison,
  file.path(out, "19_gene_robustness_comparison.csv"),
  row.names=FALSE
)

cat("\n=== FOCUSED CANDIDATE ROBUSTNESS ===\n")

candidates <- c(
  "GNG12","CCDC18-AS1","CD109","HLF",
  "DECR1","ACADM","HADHA","HADHB","ACAA2",
  "FYB1"
)

print(
  comparison[comparison$gene %in% candidates, ],
  row.names=FALSE
)

cat("\nCommon genes:", nrow(comparison), "\n")
cat("Direction-consistent genes:",
    sum(comparison$direction_consistent), "\n")

writeLines(
  capture.output(sessionInfo()),
  file.path(out, "19_sessionInfo.txt")
)

cat("\nSCRIPT 19 COMPLETE\n")
