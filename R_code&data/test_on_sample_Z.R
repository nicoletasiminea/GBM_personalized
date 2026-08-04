library(edgeR)
library(limma)

options(timeout = 300) 

# load and format counts
counts <- read.csv("./data/htseq.csv", row.names = 1)
counts <- as.matrix(counts)

colnames(counts) <- gsub("\\.", "-", colnames(counts))
colnames(counts) <- gsub("^X", "", colnames(counts))
rownames(counts) <- sub("\\..*", "", rownames(counts))


run_diff_expr_analysis <- function(
    counts,
    sample_map_file,
    output_dir,
    outlier = NULL,
    n_baseline = 5
) {
  
  
  # load sample map
  sample_map <- read.csv(sample_map_file)
  
  targets <- data.frame(
    sample_id = trimws(sample_map[, 1]),
    group     = c(rep("baseline", n_baseline),
                  rep("case", nrow(sample_map) - n_baseline)),
    stringsAsFactors = FALSE
  )
  rownames(targets) <- targets$sample_id
  
  # subset counts
  counts_sub <- counts[, targets$sample_id]
  stopifnot(all(colnames(counts_sub) == targets$sample_id))
  
  if (!is.null(outlier)) {
    keep <- !targets$sample_id %in% outlier
    targets   <- targets[keep, ]
    counts_sub <- counts_sub[, keep]
  }
  
  # edgeR pipeline
  y <- DGEList(counts = counts_sub, samples = targets)
  keep_genes <- filterByExpr(y, group = targets$group)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ group, data = targets)
  rownames(design) <- targets$sample_id
  
  y <- estimateDisp(y, design, robust = TRUE)
  
  logCPM <- cpm(y, log = TRUE)
  
  baseline_samples <- rownames(targets)[targets$group == "baseline"]
  baseline_mean <- rowMeans(logCPM[, baseline_samples])
  baseline_sd <- apply(logCPM[, baseline_samples], 1, sd)
  
  fit <- glmFit(y, design)
  tr  <- glmTreat(fit, lfc = 2)
  
  fn <- decideTests(tr, lfc = 2, p.value = 0.05)
  summary(fn)
  
  tab <- topTags(tr, n = Inf, sort = "none")$table
  
  upGenes   <- tab[fn == 1, ]
  downGenes <- tab[fn == -1, ]
  
  #save cohort-level results
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.table(upGenes,
              file = file.path(output_dir, "all_up.txt"),
              sep = "\t", quote = FALSE)
  
  write.table(downGenes,
              file = file.path(output_dir, "all_down.txt"),
              sep = "\t", quote = FALSE)
  
  
  
  # output folders
  
  for (p in c("ensembl/up", "ensembl/down")) {
    dir.create(file.path(output_dir, p), recursive = TRUE, showWarnings = FALSE)
  }
  
  # per-sample classification
  case_samples <- rownames(targets)[targets$group == "case"]
  
  for (sample in case_samples) {
    
    #sample_expr <- logCPM[, sample]
    
    
    sample_z_up <- (logCPM[rownames(upGenes), sample] -
                      baseline_mean[rownames(upGenes)]) /
      baseline_sd[rownames(upGenes)]
    
    up_in_sample <- names(sample_z_up)[sample_z_up > 2]
  
    
    write.table(
      up_in_sample,
      file = file.path(output_dir, "ensembl/up", paste0(sample, "_up_ensembl.txt")),
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    
    sample_z_down <- (logCPM[rownames(downGenes), sample] -
                        baseline_mean[rownames(downGenes)]) /
      baseline_sd[rownames(downGenes)]
    
    down_in_sample <- names(sample_z_down)[sample_z_down < -2]
    
    write.table(
      down_in_sample,
      file = file.path(output_dir, "ensembl/down", paste0(sample, "_down_ensembl.txt")),
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }
  
  invisible(list(up = upGenes, down = downGenes))
}


#primary cohort
run_diff_expr_analysis(
  counts = counts,
  sample_map_file = "./data/cases_id_143.v1.csv",
  output_dir = "./results_on_sample_var1/primary",
  outlier = "63a30223-e2e9-45bc-bc42-c00f9202b493"
)


# recurrent cohort
run_diff_expr_analysis(
  counts = counts,
  sample_map_file = "./data/cases_id_13.v1.csv",
  output_dir = "./results_on_sample_var1/recurrent",
  outlier = "63a30223-e2e9-45bc-bc42-c00f9202b493"
)
