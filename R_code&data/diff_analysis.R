library(edgeR)
library(limma)
library(biomaRt)

options(timeout = 300) 

# load and format counts
counts <- read.csv("./data/htseq.csv", row.names = 1)
counts <- as.matrix(counts)

colnames(counts) <- gsub("\\.", "-", colnames(counts))
colnames(counts) <- gsub("^X", "", colnames(counts))
rownames(counts) <- sub("\\..*", "", rownames(counts))

# function to convert to hugo
convert_ensembl_to_hugo <- function(ensembl_ids) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  res <- getBM(
    attributes = "hgnc_symbol",
    filters    = "ensembl_gene_id",
    values     = unique(ensembl_ids),
    mart       = mart
  )
  res[res$hgnc_symbol != "", , drop = FALSE]
}


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
  
  par(family = "Times New Roman", font.axis = 2, font.lab =2, cex=1.3, cex.lab = 1.3, cex.axis = 1.3)
  par(mar = c(5, 5, 4, 2))
  
  plotBCV(y)
  
  plotMDS(
    logCPM,
    col = as.numeric(factor(targets$group)),
    pch = 19,
    text.font=2
  )

  
  legend(
    "bottomleft",
    legend = levels(factor(targets$group)),
    col = seq_along(levels(factor(targets$group))),
    pch = 19,
    cex = 0.9,
    text.font = 2
  )
  
  fit <- glmFit(y, design)
  tr  <- glmTreat(fit, lfc = 2)
  
  plotMD(tr)
  abline(h = c(-2, 2), col = "blue")
  
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
  
  
  hugo_up <- convert_ensembl_to_hugo(rownames(upGenes))
  
  write.table(
    hugo_up,
    file = file.path(output_dir,"all_up_hugo.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  hugo_down <- convert_ensembl_to_hugo(rownames(downGenes))
  
  write.table(
    hugo_down,
    file = file.path(output_dir,"all_down_hugo.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # cohort-level thresholds
  
  up_thresholds   <- upGenes$logCPM
  down_thresholds <- downGenes$logCPM
  names(up_thresholds)   <- rownames(upGenes)
  names(down_thresholds) <- rownames(downGenes)
  
  # output folders
  
  for (p in c("ensembl/up", "ensembl/down", "hugo/up", "hugo/down")) {
    dir.create(file.path(output_dir, p), recursive = TRUE, showWarnings = FALSE)
  }
  
  # per-sample classification
  
  case_samples <- rownames(targets)[targets$group == "case"]
  
  for (sample in case_samples) {
    
    sample_expr <- logCPM[, sample]
    
    up_expr <- sample_expr[names(up_thresholds)]
    up_in_sample <- names(up_expr)[up_expr > up_thresholds[names(up_expr)]]
    
    write.table(
      up_in_sample,
      file = file.path(output_dir, "ensembl/up", paste0(sample, "_up_ensembl.txt")),
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    
    down_expr <- sample_expr[names(down_thresholds)]
    down_in_sample <- names(down_expr)[down_expr < down_thresholds[names(down_expr)]]
    
    
    write.table(
      down_in_sample,
      file = file.path(output_dir, "ensembl/down", paste0(sample, "_down_ensembl.txt")),
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    
    if (length(up_in_sample) > 0) {
      write.table(
        convert_ensembl_to_hugo(up_in_sample),
        file = file.path(output_dir, "hugo/up", paste0(sample, "_up_hugo.txt")),
        quote = FALSE, row.names = FALSE, col.names = FALSE
      )
    }
    
    if (length(down_in_sample) > 0) {
      write.table(
        convert_ensembl_to_hugo(down_in_sample),
        file = file.path(output_dir, "hugo/down", paste0(sample, "_down_hugo.txt")),
        quote = FALSE, row.names = FALSE, col.names = FALSE
      )
    }
  }
  
  invisible(list(up = upGenes, down = downGenes))
}


#primary cohort
run_diff_expr_analysis(
  counts = counts,
  sample_map_file = "./data/cases_id_143.v1.csv",
  output_dir = "./results/primary",
  outlier = "63a30223-e2e9-45bc-bc42-c00f9202b493"
)


# recurrent cohort
run_diff_expr_analysis(
  counts = counts,
  sample_map_file = "./data/cases_id_13.v1.csv",
  output_dir = "./results/recurrent",
  outlier = "63a30223-e2e9-45bc-bc42-c00f9202b493"
)
