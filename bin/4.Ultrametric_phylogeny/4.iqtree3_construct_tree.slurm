#!/bin/bash
#SBATCH -o %x-%j.iqtree3批量建树.log
#SBATCH -e %x-%j.iqtree3批量建树.err
#SBATCH --job-name=iqtree3批量建树
#SBATCH -p fat
#SBATCH --nodelist=LFH
#SBATCH --mail-user=15225754599@163.com
#SBATCH --mail-type=END,FAIL,START
#SBATCH --mem=10G
#SBATCH -c 12
#SBATCH -t 24:00:00


# 固定默认参数（无需命令行输入，直接使用这些值）
INPUT_DIR="./Single_Copy_Orthologue_Sequences_aligned_trimal"
OUTPUT_DIR="./Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3"
THREADS=12
FORCE_OVERWRITE=1  # 0=不覆盖已有结果，1=强制覆盖（如需覆盖，手动改为1即可）


# 检查输入目录是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 输入目录不存在: $INPUT_DIR"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR" || {
    echo "错误: 无法创建输出目录 '$OUTPUT_DIR'"
    exit 1
}

echo "输入目录: $INPUT_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "线程数: $THREADS"
echo ""

# 查找所有修剪后的FASTA文件
shopt -s nullglob
trimmed_files=("$INPUT_DIR"/*.trimmed.fa)
if [ ${#trimmed_files[@]} -eq 0 ]; then
    echo "错误: 在 $INPUT_DIR 中未找到修剪后的文件 (*.trimmed.fa)"
    exit 1
fi

total=${#trimmed_files[@]}
count=0
success_count=0

# 遍历所有文件
for file in "${trimmed_files[@]}"
do
    count=$((count+1))
    filename=$(basename "$file")
    
    base_name="${filename%.*}"
    result_dir="$OUTPUT_DIR/${base_name}.iqtree3"
    
    # 检查结果目录是否已存在
    if [ -d "$result_dir" ] && [ $FORCE_OVERWRITE -eq 0 ]; then
        echo "跳过文件 ($count/$total): $filename (结果目录已存在)"
        continue
    fi
    
    # 强制覆盖时删除旧结果
    if [ -d "$result_dir" ] && [ $FORCE_OVERWRITE -eq 1 ]; then
        rm -rf "$result_dir"
    fi
    
    # 创建结果目录
    mkdir -p "$result_dir" || {
        echo "警告: 无法为 $filename 创建结果目录"
        continue
    }
    
    echo "正在处理文件 ($count/$total): $filename"
    
    # 提取外类群（LOC开头序列）
    outgroup_seq=$(grep -m1 -i '^>LOC' "$file" | cut -d'>' -f2 | cut -d'|' -f1 | cut -d' ' -f1)
    
    if [ -z "$outgroup_seq" ]; then
        echo "错误: 在 $filename 中未找到以 'LOC' 开头的序列"
        touch "$result_dir/NO_OUTGROUP"
        continue
    else
        echo "      使用外类群序列: $outgroup_seq"
    fi
    
    # 复制输入文件到结果目录并运行IQ-TREE
    cp "$file" "$result_dir/"
    (
        cd "$result_dir" || { echo "无法进入结果目录"; return 1; }
        echo "开始运行IQ-TREE: $(date)" >> iqtree.log
        # 强制指定蛋白质序列类型，避免误判
        iqtree3 -s "$filename" --seqtype AA -m MFP -B 1000 --bnni -T "$THREADS" -o "$outgroup_seq" >> iqtree.log 2>&1
        
        if [ $? -ne 0 ]; then
            echo "      错误: 建树失败"
            exit 1
        fi
        
        echo "      成功生成系统发育树"
    )
    
    # 统计成功数量
    if [ $? -eq 0 ]; then
        success_count=$((success_count+1))
    else
        touch "$result_dir/FAILED"
    fi
done

echo ""
echo "处理完成! 成功构建 $success_count/$total 个系统发育树"
echo "结果目录: $OUTPUT_DIR"
