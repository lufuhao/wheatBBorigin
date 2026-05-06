#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import glob
import argparse


META_RE = re.compile(r"-c\s+(TraesCS\S+)")
CONTIG_RE = re.compile(r"##contig=<ID=([^,>]+)")


def extract_traescs_id(meta_line):
    """
    从 ##msa2vcf.meta 行中提取 TraesCS ID
    """
    m = META_RE.search(meta_line)
    return m.group(1) if m else None


def extract_from_chrom_header(line):
    """
    从 #CHROM header 行中兜底获取 TraesCS
    """
    for col in line.strip().split("\t"):
        if col.startswith("TraesCS"):
            return col
    return None


def process_vcf(input_file, output_file):
    traescs_id = None
    out_lines = []

    with open(input_file) as f:
        for line in f:

            # ---------- meta 行：优先解析 -c ----------
            if line.startswith("##msa2vcf.meta"):
                traescs_id = extract_traescs_id(line)
                out_lines.append(line)
                continue

            # ---------- contig 行：后面需要替换 ----------
            if line.startswith("##contig"):
                out_lines.append(line)
                continue

            # ---------- #CHROM 行：兜底 ----------
            if line.startswith("#CHROM"):
                if traescs_id is None:
                    traescs_id = extract_from_chrom_header(line)
                if traescs_id is None:
                    raise ValueError(
                        f"[ERROR] {input_file} 中未找到 TraesCS ID"
                    )
                out_lines.append(line)
                continue

            # ---------- 真实 SNP 数据行 ----------
            if line.startswith("#"):
                out_lines.append(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if parts[0] == "chrUn":
                parts[0] = traescs_id
            out_lines.append("\t".join(parts) + "\n")

    # ---------- 第二遍：修正 contig ----------
    fixed_lines = []
    for line in out_lines:
        if line.startswith("##contig=<ID=chrUn") and traescs_id:
            line = CONTIG_RE.sub(
                f"##contig=<ID={traescs_id}", line
            )
        fixed_lines.append(line)

    with open(output_file, "w") as out:
        out.writelines(fixed_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Replace chrUn with TraesCS ID in standard VCF files"
    )
    parser.add_argument(
        "-i", "--input_dir", required=True, help="Directory containing VCF files"
    )
    parser.add_argument(
        "-o", "--output_dir", default="updated_vcfs", help="Output directory"
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for vcf in glob.glob(os.path.join(args.input_dir, "*.vcf")):
        out_vcf = os.path.join(
            args.output_dir,
            os.path.basename(vcf).replace(".vcf", ".updated.vcf")
        )
        print(f"[INFO] Processing {vcf} → {out_vcf}")
        process_vcf(vcf, out_vcf)

    print("[DONE] All VCF files processed.")


if __name__ == "__main__":
    main()

