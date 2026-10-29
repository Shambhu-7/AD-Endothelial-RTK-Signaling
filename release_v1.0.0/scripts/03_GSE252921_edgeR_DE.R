suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(Matrix)
})

set.seed(20261028)

input <- "results/pseudobulk/"
out <- "results/manuscript_final/reproducibility/dispersion_check/primary/"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

counts <- readRDS(paste0(
  input, "GSE252921_endothelial_pseudobulk_counts.rds"
))

meta <- read.csv(paste0(
  input, "GSE252921_endothelial_donor_metadata.csv"
), stringsAsFactors=FALSE)

stopifnot(
  identical(colnames(counts), meta$donor_ID),
  all(meta$library_size == colSums(counts)),
  !anyNA(meta[,c("disease","age","sex","batch")])
)

meta$disease <- factor(meta$disease, levels=c("CT","AD"))
meta$sex <- factor(meta$sex)
meta$batch <- factor(meta$batch)
meta$age_centered <- meta$age - mean(meta$age)

cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")
cat("limma version:", as.character(packageVersion("limma")), "\n")

run_de <- function(ids, label) {
  m <- droplevels(meta[ids,,drop=FALSE])
  x <- counts[,ids,drop=FALSE]

  design <- model.matrix(
    ~ disease + age_centered + sex + batch,
    data=m
  )

  cat("\n=== ", label, " ===\n", sep="")
  cat("Donors:", nrow(m), "\n")
  cat("Design rank:", qr(design)$rank, "\n")
  cat("Design columns:", ncol(design), "\n")
  cat("Residual df:", nrow(design)-qr(design)$rank, "\n")
  print(table(m$disease))

  stopifnot(
    qr(design)$rank == ncol(design),
    "diseaseAD" %in% colnames(design),
    nrow(design) > ncol(design)
  )

  y <- DGEList(counts=as.matrix(x))

  keep <- filterByExpr(y, design=design)
  cat("Genes retained:", sum(keep), "/", length(keep), "\n")
  stopifnot(sum(keep) > 1000)

  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y, method="TMM")

  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  test <- glmQLFTest(fit, coef="diseaseAD")

  tab <- topTags(test, n=Inf, sort.by="PValue")$table
  tab$gene <- rownames(tab)
  tab <- tab[,c("gene",setdiff(colnames(tab),"gene"))]

  write.csv(
    tab,
    paste0(out,label,"_edgeR_results.csv"),
    row.names=FALSE
  )

  cat("FDR < 0.05:", sum(tab$FDR < 0.05), "\n")
  cat("FDR < 0.05, positive logFC:",
      sum(tab$FDR < 0.05 & tab$logFC > 0), "\n")
  cat("FDR < 0.05, negative logFC:",
      sum(tab$FDR < 0.05 & tab$logFC < 0), "\n")

  cat("\nTop 10 genes:\n")
  print(head(tab[,c("gene","logFC","PValue","FDR")],10),
        row.names=FALSE)

  write.csv(
    data.frame(
      donor_ID=m$donor_ID,
      disease=m$disease,
      library_size=y$samples$lib.size,
      norm_factor=y$samples$norm.factors,
      effective_library_size=
        y$samples$lib.size*y$samples$norm.factors
    ),
    paste0(out,label,"_normalization.csv"),
    row.names=FALSE
  )

  saveRDS(
    list(design=design, retained_genes=rownames(y)),
    paste0(out,label,"_design_audit.rds")
  )

  tab
}

primary <- run_de(
  seq_len(nrow(meta)),
  "GSE252921_primary"
)

# Prespecified sensitivity: remove the 61-cell control donor.
sensitivity <- run_de(
  which(meta$donor_ID != "NP063/2019"),
  "GSE252921_exclude_61cell"
)

shared <- merge(
  primary[,c("gene","logFC","FDR")],
  sensitivity[,c("gene","logFC","FDR")],
  by="gene",
  suffixes=c("_primary","_sensitivity")
)

sig <- shared$FDR_primary < 0.05

cat("\n=== SENSITIVITY COMPARISON ===\n")
cat("Primary significant genes:",sum(sig),"\n")
cat("Primary genes retaining FDR < 0.05:",
    sum(sig & shared$FDR_sensitivity < 0.05),"\n")
cat("Primary genes with concordant direction:",
    sum(sig & sign(shared$logFC_primary) ==
        sign(shared$logFC_sensitivity)),"\n")

write.csv(
  shared,
  paste0(out,"GSE252921_primary_sensitivity_comparison.csv"),
  row.names=FALSE
)

writeLines(
  capture.output(sessionInfo()),
  paste0(out,"03_R_sessionInfo.txt")
)

cat("\nSCRIPT 03 COMPLETE\n")
