f <- "data/processed/GSE252921/GSM8010381_metadata_endo.rds.gz"
con <- gzcon(gzfile(f, open="rb"))
x <- tryCatch(readRDS(con), finally=close(con))

out <- "results/composition/"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

cat("=== DISEASE BY PREPARATION ===\n")
print(table(x$type, x$prep))

cat("\n=== DISEASE BY CLUSTER ===\n")
print(table(x$type, x$seurat_clusters))

cat("\n=== DONOR-LEVEL COMPOSITION ===\n")

donors <- sort(unique(as.character(x$donor_ID)))

d <- do.call(rbind, lapply(donors, function(id) {
  z <- x[x$donor_ID == id, , drop=FALSE]
  data.frame(
    donor_ID=id,
    disease=as.character(z$type[1]),
    n_cells=nrow(z),
    fraction_cluster_0=mean(as.character(z$seurat_clusters)=="0"),
    fraction_cluster_6=mean(as.character(z$seurat_clusters)=="6"),
    fraction_prep_P=mean(as.character(z$prep)=="P"),
    fraction_prep_V=mean(as.character(z$prep)=="V")
  )
}))

print(d, row.names=FALSE)

write.csv(
  d,
  paste0(out,"GSE252921_donor_composition.csv"),
  row.names=FALSE
)

cat("\n=== GROUP MEDIAN COMPOSITION ===\n")
for (v in c("fraction_cluster_6","fraction_prep_P")) {
  cat(v,"\n")
  print(tapply(d[[v]],d$disease,summary))
}

cat("\nCOMPOSITION AUDIT COMPLETE\n")
