#!/bin/bash
#SBATCH -o %x-%j.mafft批量比对.log
#SBATCH -e %x-%j.mafft批量比对.err
#SBATCH --job-name=mafft批量比对
#SBATCH -p fat
#SBATCH --nodelist=LFH
#SBATCH --mail-user=15225754599@163.com
#SBATCH --mail-type=END,FAIL 
#SBATCH --mem=20G
#SBATCH -c 8
#SBATCH -t 24:00:00

# ----------------------------
# 固定输入输出目录
# ----------------------------
INPUT_DIR="./Single_Copy_Orthologue_Sequences_mRNA"
OUTPUT_DIR="./Single_Copy_Orthologue_Sequences_mRNA_MAFFT"
mkdir -p "$OUTPUT_DIR"

# ----------------------------
# 计数器初始化
# ----------------------------
TOTAL=0
SUCCESS=0
FAIL=0
FAIL_LIST=()

START_TIME=$(date +"%Y-%m-%d %H:%M:%S")

# ----------------------------
# 检查 MAFFT
# ----------------------------
if ! command -v mafft &> /dev/null; then
    echo "错误: MAFFT 未安装或不在 PATH 中"
    exit 1
fi

# ----------------------------
# 检查输入目录
# ----------------------------
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 输入目录不存在: $INPUT_DIR"
    exit 1
fi

# ----------------------------
# 获取 FA 文件列表
# ----------------------------
FA_FILES=("$INPUT_DIR"/*.fa)

if [ ${#FA_FILES[@]} -eq 0 ]; then
    echo "错误: 在 $INPUT_DIR 中未找到 .fa 文件"
    exit 1
fi

TOTAL=${#FA_FILES[@]}
COUNT=0

echo "======================================"
echo "MAFFT 批量比对开始"
echo "开始时间: $START_TIME"
echo "输入文件总数: $TOTAL"
echo "使用线程数: $SLURM_CPUS_PER_TASK"
echo "======================================"

# ----------------------------
# 遍历文件执行 MAFFT
# ----------------------------
for file in "${FA_FILES[@]}"; do
    COUNT=$((COUNT+1))
    FILENAME=$(basename "$file")
    OUTPUT_FILE="$OUTPUT_DIR/${FILENAME%.fa}.aligned.fa"

    echo "[$COUNT/$TOTAL] 正在处理: $FILENAME"

    mafft --auto --thread "$SLURM_CPUS_PER_TASK" "$file" > "$OUTPUT_FILE"

    if [ $? -eq 0 ] && [ -s "$OUTPUT_FILE" ]; then
        SUCCESS=$((SUCCESS+1))
        echo "  ✔ 成功: $OUTPUT_FILE"
    else
        FAIL=$((FAIL+1))
        FAIL_LIST+=("$FILENAME")
        echo "  ✘ 失败: $FILENAME"
        rm -f "$OUTPUT_FILE"
    fi
done

END_TIME=$(date +"%Y-%m-%d %H:%M:%S")

# ----------------------------
# 汇总统计
# ----------------------------
echo
echo "======================================"
echo "MAFFT 批量比对完成"
echo "开始时间: $START_TIME"
echo "结束时间: $END_TIME"
echo "--------------------------------------"
echo "输入文件总数 : $TOTAL"
echo "成功比对数   : $SUCCESS"
echo "失败文件数   : $FAIL"
echo "输出目录     : $OUTPUT_DIR"
echo "--------------------------------------"

if [ $FAIL -gt 0 ]; then
    echo "失败文件列表:"
    for f in "${FAIL_LIST[@]}"; do
        echo "  - $f"
    done
else
    echo "所有文件均成功完成比对"
fi

echo "======================================"

