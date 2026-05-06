#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import subprocess
import os
import argparse
from multiprocessing import Pool
import glob
import sys


# ------------------------------------------
# 自动查找以 TraesCS 开头的参考序列
# ------------------------------------------
def find_reference_sequence(file_path):
    try:
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith(">TraesCS"):
                    return line[1:].strip()
        return None
    except:
        return None


# ------------------------------------------
# 处理单个文件
# ------------------------------------------
def process_file(file_path, output_dir, jar_path, log_file, err_file):
    base = os.path.splitext(os.path.basename(file_path))[0]
    vcf_output = os.path.join(output_dir, f"{base}.vcf")

    ref_id = find_reference_sequence(file_path)

    command = ["java", "-jar", jar_path, "msa2vcf", file_path, "-o", vcf_output]
    if ref_id:
        command.extend(["-c", ref_id])

    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # 写入总日志
        with open(log_file, "a") as lf:
            lf.write(f"\n\n===== {file_path} =====\n")
            lf.write(result.stdout)
            lf.write(f"[LOG] {file_path} → {vcf_output}  (参考序列: {ref_id if ref_id else '无'})\n")

        # 写入总错误日志
        with open(err_file, "a") as ef:
            ef.write(f"\n\n===== {file_path} =====\n")
            ef.write(result.stderr)

        if result.returncode == 0:
            msg = f"[OK] {file_path} → {vcf_output}  (参考序列: {ref_id if ref_id else '无'})"
            print(msg)
            return True
        else:
            msg = f"[FAILED] {file_path}，查看 {err_file}"
            print(msg)
            return False

    except Exception as e:
        with open(err_file, "a") as ef:
            ef.write(f"\n\n===== {file_path} (Exception) =====\n")
            ef.write(str(e) + "\n")
        print(f"[ERROR] {file_path}: {e}")
        return False


# ------------------------------------------
# 获取 fasta 文件
# ------------------------------------------
def get_fasta_files(input_dir):
    exts = (".fa", ".fasta", ".fas", ".msa", ".aln")
    files = [
        f for f in glob.glob(os.path.join(input_dir, "*"))
        if os.path.isfile(f) and f.lower().endswith(exts)
    ]
    print(f"在目录 {input_dir} 中找到 {len(files)} 个序列文件")
    
    return files


# ------------------------------------------
# 主程序
# ------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="批量msa2vcf，使用TraesCS自动参考序列")
    parser.add_argument("input_dir", type=str, help="输入目录")
    parser.add_argument("output_dir", type=str, help="VCF输出目录")
    parser.add_argument("--threads", type=int, default=4, help="并行数（默认4）")
    parser.add_argument("--jar", type=str, required=True, help="jvarkit.jar路径")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    threads = args.threads
    jar_path = args.jar

    if not os.path.isdir(input_dir):
        print(f"错误：输入目录不存在 {input_dir}")
        return
    if not os.path.isfile(jar_path):
        print(f"错误：jvarkit.jar 不存在 {jar_path}")
        return

    os.makedirs(output_dir, exist_ok=True)

    # 日志文件放到当前路径
    log_file = "msa2vcf_all.log"
    err_file = "msa2vcf_all.err"

    open(log_file, "w").close()
    open(err_file, "w").close()

    files = get_fasta_files(input_dir)
    if not files:
        print("没有找到要处理的序列文件")
        return

    tasks = [(f, output_dir, jar_path, log_file, err_file) for f in files]

    with Pool(threads) as pool:
        results = pool.starmap(process_file, tasks)

    success = sum(1 for r in results if r)
    fail = len(files) - success

    print("\n===== 处理完成 =====")
    print(f"总文件: {len(files)}")
    print(f"成功:   {success}")
    print(f"失败:   {fail}")
    print(f"总日志: {log_file}")
    print(f"错误日志: {err_file}")


if __name__ == "__main__":
    main()

