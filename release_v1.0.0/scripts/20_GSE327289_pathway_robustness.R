suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
})

cat("=== SCRIPT 20: GSE327289 PATHWAY ROBUSTNESS ===\n")

input <- "results/GSE327289_robustness"
out <- input

model_files <- c(
  primary="primary_refit_edgeR_all_genes.csv",
  exclude_under100="exclude_under100_edgeR_all_genes.csv",
  composition_adjusted="composition_adjusted_edgeR_all_genes.csv"
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

tables <- lapply(
  model_files,
  function(f) read.csv(
    file.path(input, f),
    stringsAsFactors=FALSE
  )
)

ranks <- lapply(tables, make_ranks)

common <- Reduce(intersect, lapply(ranks, names))

ranks <- lapply(
  ranks,
  function(x) sort(x[common], decreasing=TRUE)
)

cat("Common genes:", length(common), "\n")

stopifnot(length(common) >= 1000)

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
  "HALLMARK_MTORC1_SIGNALING",
  "REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES",
  "HALLMARK_HYPOXIA",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "REACTOME_FATTY_ACID_METABOLISM",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

stopifnot(all(focus %in% names(all_sets)))

cat("Gene sets loaded:", length(all_sets), "\n")

run_model <- function(stats, label) {

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
    file.path(out, paste0(
      "20_", label, "_all_pathways.csv"
    )),
    row.names=FALSE
  )

  cat("Pathways tested:", nrow(res), "\n")
  cat(
    "FDR < 0.05:",
    sum(res$padj < 0.05, na.rm=TRUE),
    "\n"
  )

  res
}

results <- lapply(
  names(ranks),
  function(n) run_model(ranks[[n]], n)
)

names(results) <- names(ranks)

# Restrict comparison to pathways present in all models.
common_pathways <- Reduce(
  intersect,
  lapply(results, function(x) x$pathway)
)

comparison <- data.frame(
  pathway=sort(common_pathways)
)

for (model in names(results)) {

  tab <- results[[model]]

  idx <- match(comparison$pathway, tab$pathway)

  comparison[[paste0("NES_", model)]] <- tab$NES[idx]
  comparison[[paste0("FDR_", model)]] <- tab$padj[idx]
}

nes_columns <- paste0("NES_", names(results))
fdr_columns <- paste0("FDR_", names(results))

comparison$direction_consistent <- apply(
  comparison[, nes_columns, drop=FALSE],
  1,
  function(x) {
    all(is.finite(x)) &&
      all(sign(x) == sign(x[1]))
  }
)

comparison$FDR05_all_models <- apply(
  comparison[, fdr_columns, drop=FALSE],
  1,
  function(x) {
    all(is.finite(x)) && all(x < 0.05)
  }
)

comparison$robust <- (
  comparison$direction_consistent &
  comparison$FDR05_all_models
)

write.csv(
  comparison,
  file.path(out, "20_all_pathway_robustness.csv"),
  row.names=FALSE
)

focused <- comparison[
  match(focus, comparison$pathway),
  ,
  drop=FALSE
]

write.csv(
  focused,
  file.path(out, "20_focused_pathway_robustness.csv"),
  row.names=FALSE
)

cat("\n=== FOCUSED PATHWAY ROBUSTNESS ===\n")

print(
  focused[, c(
    "pathway",
    nes_columns,
    fdr_columns,
    "direction_consistent",
    "robust"
  )],
  row.names=FALSE
)

cat("\n=== OVERALL ROBUSTNESS ===\n")

cat(
  "Pathways compared:",
  nrow(comparison), "\n"
)

cat(
  "Direction-consistent:",
  sum(comparison$direction_consistent),
  "\n"
)

cat(
  "FDR-significant in all three models:",
  sum(comparison$robust),
  "\n"
)

cat("\n=== ROBUST PATHWAYS ===\n")

robust <- comparison[
  comparison$robust,
  c("pathway", nes_columns, fdr_columns),
  drop=FALSE
]

print(
  head(robust, 30),
  row.names=FALSE
)

write.csv(
  robust,
  file.path(out, "20_robust_pathways.csv"),
  row.names=FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out, "20_sessionInfo.txt")
)

cat("\nSCRIPT 20 COMPLETE\n")
