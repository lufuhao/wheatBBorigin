library(ape)

# -------------------------------
# 参数
# -------------------------------
input_file <- "Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.tree"
output_normal <- "Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.multi2di_root.tree"
error_log <- "tree_processing_errors.log"

tree_strings <- readLines(input_file)
total_trees <- length(tree_strings)
cat("Found", total_trees, "trees to process\n")

cat("", file = output_normal)
cat("", file = error_log)

success_count <- 0
failure_count <- 0


# ------------------------------------------------
# ID 前缀 → 类别映射
# ------------------------------------------------
classify_taxa <- function(tips) {
    list(
        Os   = tips[grep("^LOC", tips)],
        Bd   = tips[grep("^Bradi", tips)],
        StSt = tips[grep("^Es", tips)],
        SS   = tips[grep("^GWHPBFXR", tips)],
        BB   = tips[grep("^TraesCS", tips)]
    )
}


# ------------------------------------------------
# 逐棵树处理
# ------------------------------------------------
for (i in seq_along(tree_strings)) {

    cat("Processing tree", i, "/", total_trees, "\n")

    tryCatch({

        tree_raw <- read.tree(text = tree_strings[i])

        taxa <- classify_taxa(tree_raw$tip.label)

        # 必须包含 Os（用于定根）
        if (length(taxa$Os) == 0) {
            stop("Tree missing Os (LOC*) for rooting")
        }

        # --- ① 非二叉树 → 转二叉 ---
        if (!is.binary(tree_raw)) {
            set.seed(123)
            tree_bin <- multi2di(tree_raw)
        } else {
            tree_bin <- tree_raw
        }

        # --- ② 定根（Os 作为外群，取第一个 LOC*） ---
        tree_root <- root(
            tree_bin,
            outgroup = taxa$Os[1],
            resolve.root = TRUE
        )

        # --- ③ 输出处理后的树 ---
        cat(write.tree(tree_root), "\n",
            file = output_normal, append = TRUE)

        success_count <- success_count + 1

    }, error = function(e) {

        err <- paste("Tree", i, "ERROR:", conditionMessage(e))
        cat(err, "\n", file = error_log, append = TRUE)
        cat(err, "\n")
        failure_count <- failure_count + 1
    })
}


# ------------------------------------------------
# 总结
# ------------------------------------------------
cat("\n================ DONE ================\n")
cat("Total trees:                 ", total_trees, "\n")
cat("Successfully processed:      ", success_count, "\n")
cat("Failed:                      ", failure_count, "\n")
cat("Output saved to:             ", output_normal, "\n")
cat("Error log saved to:          ", error_log, "\n")
