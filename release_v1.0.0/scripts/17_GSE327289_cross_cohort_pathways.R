suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
})

cat("=== SCRIPT 17: CROSS-COHORT PATHWAY REPLICATION ===\n")

out <- "results/cross_cohort_replication"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

discovery <- read.csv(
  "results/manuscript_final/reproducibility/dispersion_check/primary/GSE252921_primary_edgeR_results.csv",
  stringsAsFactors=FALSE
)

validation <- read.csv(
  "results/validation_GSE327289/GSE327289_edgeR_all_genes.csv",
  stringsAsFactors=FALSE
)

make_ranks <- function(tab) {
  stopifnot(all(c("gene","logFC","F") %in% names(tab)))

  tab <- tab[
    is.finite(tab$F) &
    is.finite(tab$logFC) &
    !is.na(tab$gene) &
    nzchar(tab$gene),
  ]

  tab$score <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))

  tab <- tab[order(-abs(tab$score)), ]
  tab <- tab[!duplicated(tab$gene), ]

  ranks <- tab$score
  names(ranks) <- tab$gene

  sort(ranks, decreasing=TRUE)
}

discovery_ranks <- make_ranks(discovery)
validation_ranks <- make_ranks(validation)

common <- intersect(
  names(discovery_ranks),
  names(validation_ranks)
)

cat("\nDISCOVERY GENES:", length(discovery_ranks), "\n")
cat("VALIDATION GENES:", length(validation_ranks), "\n")
cat("SHARED GENES:", length(common), "\n")

stopifnot(length(common) >= 1000)

discovery_shared <- sort(
  discovery_ranks[common], decreasing=TRUE
)

validation_shared <- sort(
  validation_ranks[common], decreasing=TRUE
)

hallmark <- msigdbr(
  species="Homo sapiens",
  collection="H"
)

reactome <- msigdbr(
  species="Homo sapiens",
  collection="C2",
  subcollection="CP:REACTOME"
)

all_sets <- c(
  split(hallmark$gene_symbol, hallmark$gs_name),
  split(reactome$gene_symbol, reactome$gs_name)
)

all_sets <- lapply(all_sets, unique)

focus <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_HYPOXIA",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
  "REACTOME_FATTY_ACID_METABOLISM",
  "REACTOME_CELL_CELL_JUNCTION_ORGANIZATION",
  "REACTOME_TIGHT_JUNCTION_INTERACTIONS",
  "HALLMARK_ANGIOGENESIS"
)

stopifnot(all(focus %in% names(all_sets)))

cat("\nGENE SETS LOADED:", length(all_sets), "\n")

run_fgsea <- function(stats, label) {
  cat("\nRUNNING:", label, "\n")

  set.seed(20261028)

  res <- fgseaMultilevel(
    pathways=all_sets,
    stats=stats,
    minSize=10,
    maxSize=500,
    eps=0
  )

  res <- as.data.frame(res)

  res$leadingEdge <- vapply(
    res$leadingEdge,
    paste,
    collapse=";",
    FUN.VALUE=character(1)
  )

  res <- res[order(res$padj, res$pval), ]

  write.csv(
    res,
    file.path(out, paste0(label, "_all_pathways.csv")),
    row.names=FALSE
  )

  cat("Pathways tested:", nrow(res), "\n")
  cat("FDR < 0.05:", sum(res$padj < 0.05, na.rm=TRUE), "\n")

  res
}

discovery_res <- run_fgsea(
  discovery_shared,
  "GSE252921_shared_genes"
)

validation_res <- run_fgsea(
  validation_shared,
  "GSE327289_shared_genes"
)

validation_full_res <- run_fgsea(
  validation_ranks,
  "GSE327289_full_genes"
)

# Compare pathways tested in both shared-gene analyses.
comparison <- merge(
  discovery_res[, c("pathway","NES","pval","padj","size")],
  validation_res[, c("pathway","NES","pval","padj","size")],
  by="pathway",
  suffixes=c("_discovery","_validation")
)

comparison$same_direction <- (
  sign(comparison$NES_discovery) ==
  sign(comparison$NES_validation)
)

comparison$both_FDR05 <- (
  comparison$padj_discovery < 0.05 &
  comparison$padj_validation < 0.05
)

comparison$replicated_FDR05 <- (
  comparison$same_direction &
  comparison$both_FDR05
)

comparison <- comparison[
  order(comparison$padj_validation),
]

write.csv(
  comparison,
  file.path(out, "cross_cohort_all_pathway_comparison.csv"),
  row.names=FALSE
)

focused <- comparison[
  match(focus, comparison$pathway),
]

write.csv(
  focused,
  file.path(out, "cross_cohort_focused_pathways.csv"),
  row.names=FALSE
)

cat("\n=== FOCUSED PATHWAY COMPARISON ===\n")

print(
  focused[, c(
    "pathway",
    "NES_discovery",
    "padj_discovery",
    "NES_validation",
    "padj_validation",
    "same_direction",
    "replicated_FDR05"
  )],
  row.names=FALSE
)

cat("\n=== CROSS-COHORT SUMMARY ===\n")

cat(
  "Pathways tested in both cohorts:",
  nrow(comparison), "\n"
)

cat(
  "Same-direction pathways:",
  sum(comparison$same_direction, na.rm=TRUE), "\n"
)

cat(
  "Same-direction FDR-significant in both:",
  sum(comparison$replicated_FDR05, na.rm=TRUE), "\n"
)

# Candidate-gene comparison
candidates <- c(
  "GNG12","CCDC18-AS1","CD109","HLF","AC097528.1",
  "DECR1","ACADM","HADHA","ACAA2","HADHB"
)

gene_comparison <- merge(
  discovery[, c("gene","logFC","PValue","FDR")],
  validation[, c("gene","logFC","PValue","FDR")],
  by="gene",
  suffixes=c("_discovery","_validation")
)

gene_comparison <- gene_comparison[
  gene_comparison$gene %in% candidates,
]

gene_comparison$same_direction <- (
  sign(gene_comparison$logFC_discovery) ==
  sign(gene_comparison$logFC_validation)
)

write.csv(
  gene_comparison,
  file.path(out, "cross_cohort_candidate_genes.csv"),
  row.names=FALSE
)

cat("\n=== CANDIDATE GENE COMPARISON ===\n")
print(gene_comparison, row.names=FALSE)

writeLines(
  capture.output(sessionInfo()),
  file.path(out, "17_sessionInfo.txt")
)

cat("\nSCRIPT 17 COMPLETE\n")
