#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re

# ================= 参数解析 =================
if len(sys.argv) not in (4, 5):
    print("Usage: python3 extract_best_gmap_mrna.py <gff3_dir> <output.gff3> <log.txt> [coverage_cutoff]")
    sys.exit(1)

gff_dir = sys.argv[1]
out_gff = sys.argv[2]
log_file = sys.argv[3]

# coverage 阈值（默认 80）
coverage_cutoff = float(sys.argv[4]) if len(sys.argv) == 5 else 80.0

# ================= 输出文件名 =================
prefix = os.path.splitext(out_gff)[0]
cutoff_tag = str(int(coverage_cutoff))

ge_cutoff_gff = prefix + f".mRNA.coverage.{cutoff_tag}.gff3"
ge_cutoff_id_list = prefix + f".mRNA.coverage.{cutoff_tag}.ID.list"

# ================= 正则 =================
coverage_re = re.compile(r'coverage=([\d\.]+)')
id_re = re.compile(r'ID=([^;]+)')
name_re = re.compile(r'Name=([^;]+)')

results = []
failed_logs = []

# ================= 逐文件处理 =================
for fname in sorted(os.listdir(gff_dir)):
    if not fname.endswith(".gff3"):
        continue

    fpath = os.path.join(gff_dir, fname)

    mrna_blocks = {}
    mrna_coverage = {}

    with open(fpath) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                continue

            if fields[2] != "mRNA":
                continue

            attr = fields[8]

            m_id = id_re.search(attr)
            m_name = name_re.search(attr)
            m_cov = coverage_re.search(attr)

            if not m_id:
                continue

            mrna_id = m_id.group(1)
            mrna_name = m_name.group(1) if m_name else mrna_id
            coverage = float(m_cov.group(1)) if m_cov else 0.0

            mrna_blocks[mrna_id] = {
                "line": line,
                "name": mrna_name
            }
            mrna_coverage[mrna_id] = coverage

    mrna_count = len(mrna_blocks)

    if mrna_count == 0:
        failed_logs.append({
            "file": fname,
            "count": 0,
            "best": "NA",
            "all": "NA"
        })
        continue

    all_covs = sorted(mrna_coverage.values(), reverse=True)
    all_cov_str = ",".join(f"{c:.1f}" for c in all_covs)

    best_mrna = max(mrna_coverage, key=mrna_coverage.get)
    best_cov = mrna_coverage[best_mrna]

    results.append({
        "best_cov": best_cov,
        "mrna_id": best_mrna,
        "mrna_name": mrna_blocks[best_mrna]["name"],
        "mrna_line": mrna_blocks[best_mrna]["line"],
        "log": {
            "file": fname,
            "count": mrna_count,
            "best": best_cov,
            "all": all_cov_str
        }
    })

# ================= 排序 =================
results.sort(key=lambda x: x["best_cov"], reverse=True)

# ================= 输出 best mRNA GFF3 =================
with open(out_gff, "w") as out:
    out.write("##gff-version 3\n")
    for r in results:
        out.write(r["mrna_line"] + "\n")
        out.write("###\n")

# ================= 输出 coverage ≥ cutoff 的 GFF3 =================
ge_cutoff_count = 0
with open(ge_cutoff_gff, "w") as outc:
    outc.write("##gff-version 3\n")
    for r in results:
        if r["best_cov"] >= coverage_cutoff:
            outc.write(r["mrna_line"] + "\n")
            ge_cutoff_count += 1

# ================= 输出 coverage ≥ cutoff 的 ID list（Name） =================
with open(ge_cutoff_id_list, "w") as idout:
    for r in results:
        if r["best_cov"] >= coverage_cutoff:
            idout.write(r["mrna_name"] + "\n")

# ================= 输出 log =================
with open(log_file, "w") as log:
    log.write("File\tTotal_mRNA\tBest_coverage\tAll_coverages\n")
    for r in results:
        l = r["log"]
        log.write(f"{l['file']}\t{l['count']}\t{l['best']:.1f}\t{l['all']}\n")
    for l in failed_logs:
        log.write(f"{l['file']}\t0\tNA\tNA\n")

# ================= 终端统计 =================
print("Finished.")
print(f"GFF3 written to: {out_gff}")
print(f"coverage ≥ {coverage_cutoff} GFF3 written to: {ge_cutoff_gff}")
print(f"coverage ≥ {coverage_cutoff} ID list written to: {ge_cutoff_id_list}")
print(f"Total best mRNA lines: {len(results)}")
print(f"coverage ≥ {coverage_cutoff} mRNA lines: {ge_cutoff_count}")
print(f"Log written to: {log_file}")

