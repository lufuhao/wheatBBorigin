library(ape)

# 设置参数
input_file <- "Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.multi2di_root.tree"
output_success <- "Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.multi2di_root.chronos.tree"
output_failure <- "time_calibrated_failure.nwk"
error_log <- "calibration_errors.log"

# 时间校准参数
calibration_params <- list(
  age.min = c(58, 33, 5.7),
  age.max = c(62.0, 35.5, 10)
)

# 读取树
tree_strings <- readLines(input_file)
total_trees <- length(tree_strings)

cat("", file = output_success)
cat("", file = output_failure)
cat("时间校准错误日志\n", file = error_log)
cat("=================\n\n", file = error_log, append = TRUE)

success_count <- 0
failure_count <- 0


# ------------------------------------------------
# ID 前缀 → 类别映射
# ------------------------------------------------
classify_taxa <- function(tips) {
  list(
    Os   = tips[grep("^LOC", tips)],
    Bd   = tips[grep("^Bradi", tips)],
    SS   = tips[grep("^GWHPBFXR", tips)],
    StSt = tips[grep("^Es", tips)]
  )
}


# ------------------------------------------------
# 逐棵树进行时间校准
# ------------------------------------------------
for (i in seq_along(tree_strings)) {

  current_tree_string <- tree_strings[i]

  tryCatch({

    # 1. 读取树
    tree <- read.tree(text = current_tree_string)

    # 2. 检查是否二叉 + 定根
    if (!is.binary(tree)) {
      stop("树未二叉化")
    }
    if (!is.rooted(tree)) {
      stop("树未定根")
    }

    taxa <- classify_taxa(tree$tip.label)

    # 必须存在校准所需物种
    if (length(taxa$Os) == 0 ||
        length(taxa$Bd) == 0 ||
        length(taxa$SS) == 0 ||
        length(taxa$StSt) == 0) {
      stop("缺少校准所需的 Os / Bd / SS / StSt")
    }

    # 3. 确定校准节点（逻辑不变）
    mrca_nodes <- c(
      getMRCA(tree, c(taxa$Os[1],   taxa$Bd[1])),
      getMRCA(tree, c(taxa$Bd[1],   taxa$SS[1])),
      getMRCA(tree, c(taxa$SS[1],   taxa$StSt[1]))
    )

    if (any(is.na(mrca_nodes))) {
      stop("无法找到有效的 MRCA 校准节点")
    }

    # 4. 创建校准表
    calibration <- makeChronosCalib(
      phy = tree,
      node = mrca_nodes,
      age.min = calibration_params$age.min,
      age.max = calibration_params$age.max
    )

    # 5. 执行时间校准
    calibrated_tree <- chronos(
      phy = tree,
      calibration = calibration,
      model = "correlated",
      quiet = TRUE
    )

    # 6. 校准结果检查
    if (any(calibrated_tree$edge.length < 0, na.rm = TRUE)) {
      stop("出现负分支长度")
    }
    if (any(is.na(calibrated_tree$edge.length))) {
      stop("分支长度包含 NA")
    }

    # 7. 成功输出
    cat(write.tree(calibrated_tree), "\n",
        file = output_success, append = TRUE)

    success_count <- success_count + 1

  }, error = function(e) {

    err_msg <- paste("第", i, "棵树校准失败：", conditionMessage(e))
    cat(err_msg, "\n", file = error_log, append = TRUE)

    cat(paste0("# 第", i, "棵树失败：", conditionMessage(e), "\n"),
        file = output_failure, append = TRUE)
    cat(current_tree_string, "\n",
        file = output_failure, append = TRUE)

    failure_count <<- failure_count + 1
  })
}


# ------------------------------------------------
# 总结
# ------------------------------------------------
cat("时间校准处理完成！\n")
cat("总树数量：", total_trees, "\n")
cat("成功校准的树数量：", success_count, "\n")
cat("校准失败的树数量：", failure_count, "\n")
cat("成功输出文件：", output_success, "\n")
cat("失败输出文件：", output_failure, "\n")
cat("错误日志：", error_log, "\n")

if (success_count + failure_count != total_trees) {
  warning("计数不一致：成功 + 失败 != 总数")
}

