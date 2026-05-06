#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import re

def get_variant_type(ref, alts):
    # ALT 是单点，算作 INDEL
    if alts == ".":
        return "INDEL"
    else:
        alt_list = alts.split(",")
    if all(len(a) == 1 and len(ref) == 1 for a in alt_list):
        return "SNP"
    elif all(len(a) == len(ref) and len(ref) > 1 for a in alt_list):
        return "MNP"
    else:
        return "INDEL"

def has_special_char(ref, alts):
    """判断REF或ALT是否含非ACGT字符"""
    combined = ref
    if alts != ".":
        combined += alts
    return not re.fullmatch(r"[ACGTagct,]*", combined)

if len(sys.argv) != 3:
    print("用法: python3 classify_vcf.py <input_dir> <output_dir>")
    sys.exit(1)

input_dir = sys.argv[1]
output_dir = sys.argv[2]

# 创建结果文件夹（不再创建 SPECIAL 文件夹）
for t in ["SNP", "MNP", "INDEL"]:
    os.makedirs(os.path.join(output_dir, t), exist_ok=True)

log_path = os.path.join(output_dir, "classification_log.txt")
special_file_path = os.path.join(output_dir, "special_characters.vcf")

with open(log_path, "w") as log_file, open(special_file_path, "w") as special_file:
    log_file.write("文件名\t总行数\tSNP数\tMNP数\tINDEL数\t特殊字符行数\t是否相等\n")

    for filename in os.listdir(input_dir):
        if not filename.endswith(".vcf"):
            continue
        vcf_in = os.path.join(input_dir, filename)

        snp_path = os.path.join(output_dir, "SNP", filename.replace(".vcf", ".snp.vcf"))
        mnp_path = os.path.join(output_dir, "MNP", filename.replace(".vcf", ".mnp.vcf"))
        indel_path = os.path.join(output_dir, "INDEL", filename.replace(".vcf", ".INDEL.vcf"))

        total_data_lines = 0
        counts = {"SNP":0, "MNP":0, "INDEL":0, "SPECIAL":0}

        with open(vcf_in) as fin, \
             open(snp_path, "w") as snp_out, \
             open(mnp_path, "w") as mnp_out, \
             open(indel_path, "w") as indel_out:

            header_written = False
            for line in fin:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    snp_out.write(line)
                    mnp_out.write(line)
                    indel_out.write(line)
                    header_written = True
                    continue
                if not header_written or line.strip() == "":
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue

                total_data_lines += 1
                ref = parts[3]
                alts = parts[4]

                # 检查特殊字符
                if has_special_char(ref, alts) and alts != ".":
                    special_file.write(f"{filename}\t行号:{total_data_lines}\t{line}")
                    counts["SPECIAL"] += 1
                    # 跳过分类，不写入 SNP/MNP/INDEL 文件
                    continue

                # 分类
                vtype = get_variant_type(ref, alts)
                if vtype == "SNP":
                    snp_out.write(line)
                    counts["SNP"] += 1
                elif vtype == "MNP":
                    mnp_out.write(line)
                    counts["MNP"] += 1
                else:
                    indel_out.write(line)
                    counts["INDEL"] += 1

        sum_counts = counts["SNP"] + counts["MNP"] + counts["INDEL"]
        check = "Equal" if sum_counts + counts["SPECIAL"] == total_data_lines else "Unequal"

        # 写入日志
        log_file.write(f"{filename}\t{total_data_lines}\t{counts['SNP']}\t{counts['MNP']}\t{counts['INDEL']}\t{counts['SPECIAL']}\t{check}\n")

print(f"✅ 分类完成，统计表格保存在: {log_path}")
print(f"✅ 含特殊字符的行已保存至: {special_file_path}")

