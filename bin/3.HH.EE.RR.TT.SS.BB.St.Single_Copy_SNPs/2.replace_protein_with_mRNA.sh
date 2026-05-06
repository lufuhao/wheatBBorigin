#!/bin/bash

# ============================================
# Usage:
#   bash replace_batch.sh <protein_dir> <mrna_dir> <output_dir>
# ============================================

protein_dir=$1
mrna_dir=$2
output_dir=$3

if [ $# -ne 3 ]; then
    echo "Usage: bash replace_batch.sh <protein_dir> <mrna_dir> <output_dir>"
    exit 1
fi

mkdir -p "$output_dir"

main_log="./replace.log"
failed_log="./failed_files.log"

# 清空旧日志
> "$main_log"
> "$failed_log"

python3 - "$protein_dir" "$mrna_dir" "$output_dir" << 'EOF'
import sys, os
from Bio import SeqIO
from datetime import datetime

protein_dir = sys.argv[1]
mrna_dir    = sys.argv[2]
output_dir  = sys.argv[3]

main_log   ="replace.log"
failed_log ="failed_files.log"

# ===============================
# 日志函数
# ===============================
def log(msg):
    with open(main_log, "a") as f:
        f.write(f"{datetime.now()}  {msg}\n")
    print(msg)

def log_failed_group(filename, missing_list):
    with open(failed_log, "a") as f:
        f.write(f"[FILE] {filename}\n")
        for item in missing_list:
            f.write(f"    missing: {item}\n")
        f.write("\n")
    print(f"[FAILED] {filename} 共有 {len(missing_list)} 条序列未找到")

# ===============================
# SECCE 特殊：去掉 .CDS.1
# ===============================
def get_secore(pid):
    if pid.endswith(".CDS.1"):
        return pid[:-6]
    return pid

# ===============================
# 处理的 7 类物种（新增 Ammut）
# ===============================
species_prefix_files = {
    "Ammut":    "Aegilops_mutica.mRNA.TT.primary.fa",
    "GWHPBFXR": "Aegilops_speltoides.mRNA.SS.primary.fa",
    "Es":       "Elymus_sibiricus.mRNA.St.primary.fa",
    "HORVU":    "HvulgareMorex_702_V3.transcript.primary.fa",
    "SECCE":    "Secale_cereale.cdna.RR.primary.fa",
    "GWHPABKY": "Thinopyrum_elongatum.mRNA.EE.primary.fa",
    "TraesCS":  "Triticum_aestivum.mRNA.BB.primary.fa"
}

log("开始加载 mRNA 文件到内存……")

# ===============================
# 加载 mRNA 到字典
# ===============================
mrna_dict = {}

for prefix, filename in species_prefix_files.items():
    path = os.path.join(mrna_dir, filename)
    log(f"加载 {prefix} => {path}")
    mrna_dict[prefix] = {r.id: r for r in SeqIO.parse(path, "fasta")}

log("所有 mRNA 文件加载完成。\n")

# ===============================
# 处理 protein 文件
# ===============================
protein_files = sorted(
    f for f in os.listdir(protein_dir)
    if f.endswith(".fa") or f.endswith(".fasta")
)

log(f"检测到 {len(protein_files)} 个 protein 文件，开始处理。\n")

def find_prefix(pid):
    for p in species_prefix_files:
        if pid.startswith(p):
            return p
    return None

for pfile in protein_files:
    log(f"==== 开始处理：{pfile} ====")

    protein_path = os.path.join(protein_dir, pfile)
    output_path  = os.path.join(
        output_dir,
        pfile.replace(".fa", "_mRNA.fa").replace(".fasta", "_mRNA.fasta")
    )

    out_records = []
    missing = []

    for record in SeqIO.parse(protein_path, "fasta"):
        pid = record.id
        prefix = find_prefix(pid)

        if not prefix:
            log(f"  [FAIL] 未识别前缀：{pid}")
            missing.append(pid)
            continue

        # ---------- SECCE ----------
        if prefix == "SECCE":
            se_core = get_secore(pid)
            matched = None
            for mid, mrec in mrna_dict["SECCE"].items():
                if se_core == mid or se_core in mid:
                    matched = mrec
                    break

            if matched:
                out_records.append(matched)
                log(f"  [OK-SECCE] {pid} => {matched.id}")
            else:
                log(f"  [FAIL-SECCE] 未找到 mRNA：{pid}")
                missing.append(pid)
            continue

        # ---------- 其他 6 类（含 Ammut）：严格 ID 匹配 ----------
        if pid in mrna_dict[prefix]:
            out_records.append(mrna_dict[prefix][pid])
            log(f"  [OK] {pid}")
        else:
            log(f"  [FAIL] 未找到 mRNA：{pid}")
            missing.append(pid)

    # ---------- 输出 ----------
    if out_records:
        SeqIO.write(out_records, output_path, "fasta")
        log(f"[OUTPUT] 已输出 {len(out_records)} 条 mRNA 到 {output_path}")

    if missing:
        log_failed_group(pfile, missing)
        log(f"[SKIP] {pfile} 存在未匹配 ID\n")

log("全部完成！")
EOF

