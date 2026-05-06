#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
提取所有VCF文件的内容行，生成表格

"""

import sys
import os

indir  = sys.argv[1]
outfile = sys.argv[2]

out_header = ["#CHROM", "BB", "SS", "St", "EE", "RR", "HH", "TT"]

total_valid_lines = 0
total_output_lines = 0
header_written = False

with open(outfile, "w") as fout:

    for fname in sorted(os.listdir(indir)):
        fpath = os.path.join(indir, fname)

        if not os.path.isfile(fpath):
            continue

        with open(fpath) as fin:
            col = {}

            for line in fin:
                line = line.rstrip("\n")
                if not line:
                    continue

                # 处理表头
                if line.startswith("#CHROM"):
                    fields = line.split("\t")
                    for i, name in enumerate(fields):
                        if name == "#CHROM":
                            col["CHROM"] = i
                        elif name.startswith("TraesCS"):
                            col["BB"] = i
                        elif name.startswith("GWHPBFXR"):
                            col["SS"] = i
                        elif name.startswith("Es"):
                            col["St"] = i
                        elif name.startswith("GWHPABKY"):
                            col["EE"] = i
                        elif name.startswith("SECCE"):
                            col["RR"] = i
                        elif name.startswith("HORVU"):
                            col["HH"] = i
                        elif name.startswith("Ammut"):
                            col["TT"] = i

                    # 检查是否齐全
                    missing = [k for k in ["CHROM","BB","SS","St","EE","RR","HH","TT"] if k not in col]
                    if missing:
                        sys.stderr.write(
                            f"[WARN] 文件 {fname} 缺少列: {','.join(missing)}，已跳过该文件\n"
                        )
                        col = {}
                    else:
                        if not header_written:
                            fout.write("\t".join(out_header) + "\n")
                            header_written = True
                    continue

                # 非表头行
                if not col:
                    continue

                total_valid_lines += 1
                fields = line.split("\t")

                try:
                    out = [
                        fields[col["CHROM"]],
                        fields[col["BB"]],
                        fields[col["SS"]],
                        fields[col["St"]],
                        fields[col["EE"]],
                        fields[col["RR"]],
                        fields[col["HH"]],
                        fields[col["TT"]],
                    ]
                except IndexError:
                    continue

                fout.write("\t".join(out) + "\n")
                total_output_lines += 1

# 统计信息输出到 stderr
sys.stderr.write(f"Total valid lines read (excluding headers): {total_valid_lines}\n")
sys.stderr.write(f"Total valid lines written: {total_output_lines}\n")

