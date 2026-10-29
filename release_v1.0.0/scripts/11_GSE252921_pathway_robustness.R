suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
})

out <- "results/manuscript_final/reproducibility/dispersion_check/pathway_robustness"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

cat("=== SCRIPT 11: PATHWAY ROBUSTNESS ===\n")

primary <- read.csv(
  "results/manuscript_final/reproducibility/dispersion_check/primary/GSE252921_primary_edgeR_results.csv"
)

sensitivity <- read.csv(
  "results/manuscript_final/reproducibility/dispersion_check/primary/GSE252921_exclude_61cell_edgeR_results.csv"
)

adjusted <- read.csv(
  "results/manuscript_final/reproducibility/dispersion_check/composition/composition_adjusted_results.csv"
)

make_ranks <- function(tab) {
  stopifnot(all(c("gene", "logFC", "F") %in% names(tab)))

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

ranks <- list(
  primary=make_ranks(primary),
  exclude_61cell=make_ranks(sensitivity),
  composition_adjusted=make_ranks(adjusted)
)

common <- Reduce(intersect, lapply(ranks, names))

ranks <- lapply(
  ranks,
  function(x) sort(x[common], decreasing=TRUE)
)

cat("Common genes across models:", length(common), "\n")

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

cat("Gene sets loaded:", length(all_sets), "\n")

run_model <- function(stats, model_name) {
  cat("\nRunning:", model_name, "\n")

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
    file.path(out, paste0(model_name, "_all_pathways.csv")),
    row.names=FALSE
  )

  cat("Pathways tested:", nrow(res), "\n")
  cat(
    "Significant at FDR < 0.05:",
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

comparison <- data.frame(pathway=focus)

for (model_name in names(results)) {
  x <- results[[model_name]]

  idx <- match(comparison$pathway, x$pathway)

  comparison[[paste0("NES_", model_name)]] <- x$NES[idx]
  comparison[[paste0("FDR_", model_name)]] <- x$padj[idx]
}

comparison$direction_consistent <-
  apply(
    comparison[
      ,
      paste0("NES_", names(results)),
      drop=FALSE
    ],
    1,
    function(x) all(is.finite(x)) &&
      length(unique(sign(x))) == 1
  )

comparison$significant_all_three <-
  comparison$FDR_primary < 0.05 &
  comparison$FDR_exclude_61cell < 0.05 &
  comparison$FDR_composition_adjusted < 0.05

comparison <- comparison[
  order(comparison$FDR_primary),
]

write.csv(
  comparison,
  file.path(out, "focused_pathway_robustness.csv"),
  row.names=FALSE
)

cat("\n=== FOCUSED PATHWAY ROBUSTNESS ===\n")

print(
  comparison,
  row.names=FALSE,
  digits=4
)

cat("\n=== SUMMARY ===\n")

cat(
  "Direction consistent across all models:",
  sum(comparison$direction_consistent, na.rm=TRUE),
  "/", nrow(comparison), "\n"
)

cat(
  "Significant across all three models:",
  sum(comparison$significant_all_three, na.rm=TRUE),
  "/", nrow(comparison), "\n"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out, "11_sessionInfo.txt")
)

cat("\nSCRIPT 11 COMPLETE\n")
