#!/bin/bash
#SBATCH -o %x-%j.gmap.log
#SBATCH -e %x-%j.gmap.err
#SBATCH --job-name=flat_gmap
#SBATCH -p fat
#SBATCH --nodelist=LFH
#SBATCH --mail-user=15225754599@163.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=20G
#SBATCH -c 12
#SBATCH -t 6:00:00

# -----------------------------
# 固定参数（可直接改）
# -----------------------------
THREADS=12
MIN_IDENTITY=0.8

# -----------------------------
# 接收终端传入的 4 个参数
# -----------------------------
if [ $# -ne 4 ]; then
    echo "用法：sbatch $0 <输入fasta> <GMAP根目录> <输出根目录> <GMAP数据库名>"
    exit 1
fi

# 从命令行赋值（完全外部指定）
INPUT_FASTA="$1"
GMAP_BASE_DIR="$2"
OUTPUT_DIR="$3"  # 直接作为最终输出目录，不再嵌套
DB_NAME="$4"

# 自动拼接数据库路径
GMAP_DB_DIR="${GMAP_BASE_DIR}/${DB_NAME}"

# 检查输入是否存在
if [ ! -f "$INPUT_FASTA" ]; then echo "错误：输入FASTA不存在 $INPUT_FASTA"; exit 1; fi
if [ ! -d "$GMAP_DB_DIR" ]; then echo "错误：GMAP数据库不存在 $GMAP_DB_DIR"; exit 1; fi

mkdir -p "$OUTPUT_DIR"
echo "========================================"
echo "输入FASTA  : $INPUT_FASTA"
echo "GMAP数据库 : $DB_NAME"
echo "输出目录   : $OUTPUT_DIR"
echo "========================================"

# -----------------------------
# Python 批量拆分每条序列运行 GMAP
# -----------------------------
python3 << END_PYTHON
import os
from pathlib import Path
import subprocess

input_fasta = Path("$INPUT_FASTA")
gmap_db_dir = Path("$GMAP_DB_DIR")
output_dir = Path("$OUTPUT_DIR")
threads = $THREADS
min_id = $MIN_IDENTITY

# 读取FASTA
seqs = {}
with open(input_fasta, 'r') as f:
    seq_id = None
    seq_lines = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq_id:
                seqs[seq_id] = "".join(seq_lines)
            seq_id = line[1:].split()[0]
            seq_lines = []
        else:
            seq_lines.append(line)
    if seq_id:
        seqs[seq_id] = "".join(seq_lines)

# 每条序列单独运行，直接输出到 OUTPUT_DIR
for seq_id, seq in seqs.items():
    # 临时FASTA只在内存/当前目录生成，不写入结果目录
    tmp_fa = Path(f"{seq_id}.tmp.fa")
    gff_out = output_dir / f"{seq_id}.gff3"

    with open(tmp_fa, 'w') as f:
        f.write(f">{seq_id}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

    cmd = [
        "gmap",
        "-d", gmap_db_dir.name,
        "-D", str(gmap_db_dir.parent),
        "-f", "gff3_gene",
        "-t", str(threads),
        "--min-identity", str(min_id),
        str(tmp_fa)
    ]

    print(f"Running: {seq_id} -> {gff_out}")
    with open(gff_out, 'w') as out:
        subprocess.run(cmd, stdout=out, check=True)

    # 删除临时FASTA，不留在结果目录
    tmp_fa.unlink()

END_PYTHON

echo -e "\n✅ 全部完成！结果直接位于：$OUTPUT_DIR"

