#!/bin/bash
#SBATCH -o %x-%j.muscle批量比对.log
#SBATCH -e %x-%j.muscle批量比对.err
#SBATCH --job-name=muscle批量比对
#SBATCH -p fat
#SBATCH --nodelist=LFH
#SBATCH --mail-user=15225754599@163.com
#SBATCH --mail-type=END,FAIL 
#SBATCH --mem=20000  
#SBATCH -c 6  
#SBATCH -t 04:00:00

# 输入目录
INPUT_DIR="./Single_Copy_Orthologue_Sequences"
# 输出目录
OUTPUT_DIR="./Single_Copy_Orthologue_Sequences_aligned"
mkdir -p "$OUTPUT_DIR"

# 检查输入目录是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 输入目录不存在: $INPUT_DIR"
    exit 1
fi

# 检查是否有.fa文件
shopt -s nullglob  # 若没有匹配文件，让$files为空
files="$INPUT_DIR"/*.fa
if [ -z "$files" ]; then
    echo "错误: 在 $INPUT_DIR 中未找到.fa文件"
    exit 1
fi

# 统计文件总数
total=$(ls -1 "$INPUT_DIR"/*.fa | wc -l)
count=0

# 遍历所有FA文件
for file in "$INPUT_DIR"/*.fa; do
    # 跳过非文件（如子目录）
    [ -f "$file" ] || continue
    
    count=$((count+1))
    filename=$(basename "$file")
    echo "[$(date)] 正在处理文件 ($count/$total): $filename"
    
    # 输出文件路径
    output_file="${OUTPUT_DIR}/${filename%.fa}.aligned.afa"
    
    # 运行Muscle比对（假设muscle在环境变量中）
    muscle -align "$file" -output "$output_file"
    
    # 检查执行结果
    if [ $? -ne 0 ]; then
        echo "[$(date)] 错误: 处理 $filename 失败"
    else
        echo "[$(date)] 完成: $output_file"
    fi
done

echo "[$(date)] 所有文件处理完毕!"
