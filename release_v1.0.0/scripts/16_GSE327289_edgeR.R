suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
})

cat("=== SCRIPT 16: GSE327289 ENDOTHELIAL edgeR ===\n")

out <- "results/validation_GSE327289"

counts <- readRDS(
  file.path(out, "GSE327289_endothelial_pseudobulk_counts.rds")
)

meta <- read.csv(
  file.path(out, "GSE327289_endothelial_donor_metadata.csv"),
  stringsAsFactors=FALSE
)

meta$Donor <- as.character(meta$Donor)
stopifnot(setequal(colnames(counts), meta$Donor))
meta <- meta[match(colnames(counts), meta$Donor), , drop=FALSE]
stopifnot(identical(colnames(counts), meta$Donor))
cat("\nDONOR ALIGNMENT: PASS\n")

meta$Diagnosis <- factor(meta$Diagnosis, levels=c("CTRL","AD"))
meta$Sex <- factor(meta$Sex)
meta$Age_centered <- meta$Age - mean(meta$Age)

stopifnot(!anyNA(meta[,c("Diagnosis","Sex","Age_centered")]))

design <- model.matrix(
  ~ Diagnosis + Age_centered + Sex,
  data=meta
)

cat("\nDESIGN MATRIX:\n")
print(design)

stopifnot(qr(design)$rank == ncol(design))

cat("\nGROUP SIZES:\n")
print(table(meta$Diagnosis))

cat("\nLIBRARY SIZES:\n")
print(summary(Matrix::colSums(counts)))

y <- DGEList(counts=counts, group=meta$Diagnosis)

keep <- filterByExpr(y, design=design)

cat("\nGENES RETAINED:",sum(keep),"\n")

y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

y <- estimateDisp(y, design, robust=TRUE)

fit <- glmQLFit(y, design, robust=TRUE)

test <- glmQLFTest(fit, coef="DiagnosisAD")

res <- topTags(test, n=Inf, sort.by="none")$table

res$gene <- rownames(res)

res <- res[,c(
  "gene","logFC","logCPM","F","PValue","FDR"
)]

res <- res[order(res$PValue),]

write.csv(
  res,
  file.path(out,"GSE327289_edgeR_all_genes.csv"),
  row.names=FALSE
)

cat("\nFDR < 0.05 GENES:\n")
print(sum(res$FDR < 0.05))

cat("\nTOP 20 GENES:\n")
print(head(res,20),row.names=FALSE)

cat("\nDISCOVERY CANDIDATE GENES:\n")

candidates <- c(
  "GNG12","CCDC18-AS1","CD109","HLF","AC097528.1",
  "DECR1","ACADM","HADHA","ACAA2","HADHB"
)

print(
  res[res$gene %in% candidates,],
  row.names=FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out,"16_sessionInfo.txt")
)

cat("\nSCRIPT 16 COMPLETE\n")
