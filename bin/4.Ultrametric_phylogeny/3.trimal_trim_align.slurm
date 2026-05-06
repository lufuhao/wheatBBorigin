#!/bin/bash
#SBATCH -o %x-%j.trimal批量剪切.log
#SBATCH -e %x-%j.trimal批量剪切.err
#SBATCH --job-name=trimal批量剪切
#SBATCH -p fat
#SBATCH --nodelist=LFH
#SBATCH --mail-user=15225754599@163.com
#SBATCH --mail-type=END,FAIL,START
#SBATCH --mem=10000
#SBATCH -c 6
#SBATCH -t 6:00:00

# ===============================
# 默认值
# ===============================
INPUT_DIR=""
OUTPUT_DIR=""
FORCE_OVERWRITE=0
TRIMAL_PARAMS="-gt 0.4 -cons 70"

# ===============================
# 显示用法信息
# ===============================
usage() {
    echo "用法: sbatch $0 -i <输入目录> [-o 输出目录] [-f] [-p 'trimal参数']"
    echo "选项:"
    echo "  -i  指定输入目录 (必需)"
    echo "  -o  指定输出目录 (默认为: $OUTPUT_DIR)"
    echo "  -f  强制覆盖已存在的输出文件"
    echo "  -p  自定义 trimal 参数"
    echo "  -h  显示此帮助信息"
    exit 1
}

# ===============================
# 解析参数
# ===============================
while getopts ":i:o:fp:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        f) FORCE_OVERWRITE=1 ;;
        p) TRIMAL_PARAMS="$OPTARG" ;;
        h) usage ;;
        \?) echo "无效选项: -$OPTARG" >&2; usage ;;
        :) echo "选项 -$OPTARG 需要参数" >&2; usage ;;
    esac
done
shift $((OPTIND - 1))

# ===============================
# 参数检查
# ===============================
if [ -z "$INPUT_DIR" ]; then
    echo "错误: 必须使用 -i 指定输入目录"
    usage
fi

#生成输出目录名称
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${INPUT_DIR%/}_trimal"
fi

if ! command -v trimal &> /dev/null; then
    echo "错误: trimal 未找到! 请确保 trimal 已安装并在 PATH 中"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 输入目录不存在: $INPUT_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR" || {
    echo "错误: 无法创建输出目录 '$OUTPUT_DIR'"
    exit 1
}

echo "=== Trimal 批处理开始 ==="
echo "SLURM 任务: $SLURM_JOB_ID"
echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "Trimal 参数: $TRIMAL_PARAMS"
echo "覆盖模式: $FORCE_OVERWRITE"
echo "==================================="

# ===============================
# 查找文件
# ===============================
shopt -s nullglob
fa_files=("$INPUT_DIR"/*.{fa,fasta,aln,afa})

if [ ${#fa_files[@]} -eq 0 ]; then
    fa_files=("$INPUT_DIR"/*.{faa,aligned.fa,aligned.afa})
    if [ ${#fa_files[@]} -eq 0 ]; then
        echo "错误: 在 $INPUT_DIR 中未找到 FASTA 文件"
        exit 1
    fi
fi

total=${#fa_files[@]}
count=0
success_count=0

# ===============================
# 批处理
# ===============================
for file in "${fa_files[@]}"
do
    count=$((count+1))
    filename=$(basename "$file")
    base_name="${filename%.*}"
    output_file="$OUTPUT_DIR/${base_name}.trimmed.fa"

    if [ -f "$output_file" ] && [ $FORCE_OVERWRITE -eq 0 ]; then
        echo "跳过 ($count/$total): $filename (输出已存在)"
        continue
    fi

    echo "处理 ($count/$total): $filename"

    trimal -in "$file" -out "$output_file" $TRIMAL_PARAMS 2>&1
    if [ $? -ne 0 ]; then
        echo "错误: $filename 处理失败"
        rm -f "$output_file"
    else
        echo "成功: $output_file"
        success_count=$((success_count+1))
    fi
done

echo "==================================="
echo "处理完成: 成功 $success_count / $total 个文件"
echo "输出目录: $OUTPUT_DIR"
echo "SLURM 任务结束: $SLURM_JOB_ID"
echo "==================================="
