library(ape)

input_file  <- "Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.tree"
output_file <- "Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.tree"

tree_strings <- readLines(input_file)

cat("", file = output_file)

get_taxa <- function(tips) {
  list(
    StSt = tips[grep("^Es", tips)],
    SS   = tips[grep("^GWHPBFXR", tips)],
    BB   = tips[grep("^TraesCS", tips)]
  )
}

# 检查 StSt + SS + BB 是否完全聚成一个 clade
is_stst_ss_bb_clade <- function(tree) {

  taxa <- get_taxa(tree$tip.label)

  # 三类必须都存在
  if (length(taxa$StSt) == 0 ||
      length(taxa$SS) == 0 ||
      length(taxa$BB) == 0) {
    return(FALSE)
  }

  st <- taxa$StSt[1]
  ss <- taxa$SS[1]
  bb <- taxa$BB[1]

  m1 <- getMRCA(tree, c(st, ss))
  m2 <- getMRCA(tree, c(st, bb))
  m3 <- getMRCA(tree, c(ss, bb))

  # 三者 MRCA 相同，即完全聚成一个 clade
  return(m1 == m2 && m2 == m3)
}

kept <- 0
removed <- 0

for (i in seq_along(tree_strings)) {

  tree <- read.tree(text = tree_strings[i])

  # 排除 StSt + SS + BB 完全聚成一个 clade 的树
  if (!is_stst_ss_bb_clade(tree)) {

    cat(write.tree(tree), "\n",
        file = output_file, append = TRUE)

    kept <- kept + 1

  } else {
    removed <- removed + 1
  }
}

cat("Total:", length(tree_strings), "\n")
cat("Kept (not StSt-SS-BB clade):", kept, "\n")
cat("Removed (StSt-SS-BB clade):", removed, "\n")
