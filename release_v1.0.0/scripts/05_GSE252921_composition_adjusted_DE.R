suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
})

input <- "results/pseudobulk/"
out <- "results/manuscript_final/reproducibility/dispersion_check/composition/"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

counts <- readRDS(paste0(
  input, "GSE252921_endothelial_pseudobulk_counts.rds"
))

meta <- read.csv(paste0(
  input, "GSE252921_endothelial_donor_metadata.csv"
), stringsAsFactors=FALSE)

comp <- read.csv(
  "results/composition/GSE252921_donor_composition.csv",
  stringsAsFactors=FALSE
)

stopifnot(
  identical(colnames(counts), meta$donor_ID),
  !anyDuplicated(comp$donor_ID)
)

comp <- comp[match(meta$donor_ID, comp$donor_ID), ]
stopifnot(identical(meta$donor_ID, comp$donor_ID))

meta$disease <- factor(meta$disease, levels=c("CT","AD"))
meta$sex <- factor(meta$sex)
meta$batch <- factor(meta$batch)
meta$age_centered <- meta$age - mean(meta$age)

meta$cluster6 <- comp$fraction_cluster_6
meta$prepP <- comp$fraction_prep_P

design_base <- model.matrix(
  ~ disease + age_centered + sex + batch,
  data=meta
)

design_adjusted <- model.matrix(
  ~ disease + age_centered + sex + batch +
    cluster6 + prepP,
  data=meta
)

cat("=== DESIGN AUDIT ===\n")
cat("Base rank:", qr(design_base)$rank, "\n")
cat("Adjusted rank:", qr(design_adjusted)$rank, "\n")
cat("Adjusted columns:", ncol(design_adjusted), "\n")
cat("Adjusted residual df:",
    nrow(meta)-qr(design_adjusted)$rank, "\n")

stopifnot(
  qr(design_base)$rank == ncol(design_base),
  qr(design_adjusted)$rank == ncol(design_adjusted),
  nrow(meta) > ncol(design_adjusted)
)

y <- DGEList(counts=as.matrix(counts))
keep <- filterByExpr(y, design=design_base)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")

cat("Common genes tested:", nrow(y), "\n")

run_model <- function(design, label) {
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  test <- glmQLFTest(fit, coef="diseaseAD")
  tab <- topTags(test, n=Inf)$table
  tab$gene <- rownames(tab)

  write.csv(
    tab,
    paste0(out, label, "_results.csv"),
    row.names=FALSE
  )

  cat(label, "FDR < 0.05:",
      sum(tab$FDR < 0.05), "\n")
  tab
}

base <- run_model(design_base, "base")
adjusted <- run_model(design_adjusted, "composition_adjusted")

genes <- c(
  "GNG12", "CCDC18-AS1", "CD109",
  "HLF", "AC097528.1"
)

comparison <- merge(
  base[,c("gene","logFC","PValue","FDR")],
  adjusted[,c("gene","logFC","PValue","FDR")],
  by="gene",
  suffixes=c("_base","_adjusted")
)

comparison <- comparison[
  match(genes, comparison$gene), ,
  drop=FALSE
]

cat("\n=== FIVE CANDIDATE GENES ===\n")
print(comparison, row.names=FALSE)

write.csv(
  comparison,
  paste0(out,"five_candidate_comparison.csv"),
  row.names=FALSE
)

writeLines(
  capture.output(sessionInfo()),
  paste0(out,"05_sessionInfo.txt")
)

cat("\nSCRIPT 05 COMPLETE\n")
