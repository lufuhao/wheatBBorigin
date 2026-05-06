#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import sys
import argparse
import csv



def extract_gt(gt_str):
	"""Extract GT geno and missing data"""
	if gt_str == './.':
		return "NAN"  # missing data
	return gt_str.split(':', 1)[0].split('/', 1)[0]  # first allele only



def process_vcf(file_path, outpfx):
	### debug mode to print more info
	debug=0
	### initalization
	group_counter = 0
	ss_score = 0
	st_score = 0
	total_lines = 0     # line count for current VCF
	num_filter_out = 0  # lines removed（NOT including BB-specific SNP）
	total_ss_score = 0
	total_st_score = 0
	max_ref_allele_length=7
	### global count
	bb_specific_count = 0        # BB-specific SNP for current VCF
	total_bb_specific = 0        # Global BB-specific SNP
	total_processed_lines = 0    # Global lines
	num_filter_out_sum = 0       # Global lines filtered out (BB-SNP exclusive)
	###
	ss_col = None
	st_col = None
	bb_col = None
	bb_col_name = None
	ss_matches_count = 0
	st_matches_count = 0
	### open file for read or write
	in_file_vcf=open(file_path, 'r')
	out_file=open(output_file, 'w')
	out_file_passed=open(outpfx+".pass.vcf", 'w')          ### passed lines
	out_file_filterout=open(outpfx+".problem.vcf", 'w')    ### filter out lines
	out_file_long=open(outpfx+".long.vcf", 'w')            ### Long Indel lines, useful for find specific genes
	out_file_sum=open(outpfx+".sum", 'w', newline='')      ###
	out_file_sum_csv = csv.writer(out_file_sum, delimiter='\t')
	### parse input
	out_file_sum_csv.writerow(["#CHROM", "POS", "REF", "ALT", "st", "ss", "bb"])
	out_file.write("#ref\ttotal\tSS_var\tSt_var\tBB_var\tfilter_out\n")
	for line in in_file_vcf:
		line = line.strip('\n')
		if line.startswith('##'):
			continue
		parts = line.split('\t')
		### clean spaces
		parts = [p.strip() for p in parts]
		if len(parts) < 12:  ### QC: ncol>=12
			sys.stderr.write(f"Warning: ncol<12: {line}\n")
			continue
		if line.startswith('#CHROM'):
			### write last stat
			if total_lines > 0:
				group_counter += 1
				out_file.write(f"#{group_counter}\t{total_lines}\t{ss_score}\t{st_score}\t{bb_specific_count}\t{num_filter_out}\n")
				### Global count
				total_ss_score += ss_score
				total_st_score += st_score
				total_processed_lines += total_lines
				num_filter_out_sum += num_filter_out
				### reset count for each VCF
				ss_score = 0
				st_score = 0
				total_lines = 0
				num_filter_out = 0
				bb_specific_count = 0
			### identify SS/BB/St column
			for i, col_name in enumerate(parts):
				if "SS" in col_name:
					ss_col = i
					if debug == 1:
						print(f"Info: SS col detected: index {i} col {col_name}")
				elif "St" in col_name:
					st_col = i
					if debug == 1:
						print(f"Info: St col detected: index {i} col {col_name}")
				elif "BB" in col_name:
					bb_col = i
					bb_col_name = col_name
					if debug == 1:
						print(f"Info: BB col detected: index {i} col {col_name}")
			### check if enougn columns
			if ss_col is None or st_col is None or bb_col is None:
				sys.stderr.write(f"Error: missing SS/BB/St column: {file_path}\n")
				sys.stderr.write("    Hint: check VCF for header 'GWH'(SS), 'Es'(St), 'TraesCS'(BB)\n")
				return
			### write header line
			out_file.write(line + "\n")
			out_file_filterout.write(line+"\n")
			out_file_passed.write(line+"\n")
			out_file_long.write(line+"\n")
			continue
		### data processing
		total_lines += 1
		### set ref genome
		if bb_col_name:
			parts[0] = bb_col_name
		### Filter
		ref = parts[3]
		alt = parts[4]
		### Multi-allele in ALT
		if ',' in alt:
			num_filter_out += 1
			out_file_filterout.write(line+"\n")
			continue
		### Long Indel
		if len(ref) > max_ref_allele_length or len(alt) > max_ref_allele_length:
			num_filter_out += 1
			out_file_long.write(line+"\n")
			continue
		### double check col index
		if ss_col >= len(parts) or st_col >= len(parts) or bb_col >= len(parts):
			sys.stderr.write(f"Warnings: out-range col index, ignored: {line}\n")
			out_file_filterout.write(line+"\n")
			num_filter_out += 1
			continue
		### extract geno
		ss_gt = extract_gt(parts[ss_col])
		st_gt = extract_gt(parts[st_col])
		bb_gt = extract_gt(parts[bb_col])
		### missing geno
		if "NAN" in [ss_gt, st_gt, bb_gt]:
			num_filter_out += 1
			out_file_filterout.write(line+"\n")
			continue
		### Grouping
		if ss_gt == bb_gt and st_gt == bb_gt:  ### BB allele is BOTH St or SS, ideally won't happen
			num_filter_out += 1
			out_file_filterout.write(line+"\n")
			continue
		elif ss_gt != bb_gt and st_gt != bb_gt:  ### BB allele is NOT St or SS
			bb_specific_count += 1
			total_bb_specific += 1
			out_file_filterout.write(line+"\n")
			continue
		elif ss_gt == bb_gt:
			ss_score += 1
			ss_matches_count += 1
		elif st_gt == bb_gt:
			st_score += 1
			st_matches_count += 1
		### write passed
		out_file_passed.write(line+"\n")
		### write sum
		out_file_sum_csv.writerow([
			parts[0], parts[1], parts[3], parts[4], st_gt, ss_gt, bb_gt
		])
	### write the last ref
	if total_lines > 0:
		group_counter += 1
		out_file.write(f"#{group_counter}\t{total_lines}\t{ss_score}\t{st_score}\t{bb_specific_count}\t{num_filter_out}\n")
		total_ss_score += ss_score
		total_st_score += st_score
		total_processed_lines += total_lines
		num_filter_out_sum += num_filter_out
	### write sum
	out_file.write(f"#{'='*50}\n")
	out_file.write(f"#sum\t{total_processed_lines}\t{total_ss_score}\t{total_st_score}\t{total_bb_specific}\t{num_filter_out_sum}\n")
	out_file.write(f"#{'='*50}\n")
	out_file_sum_csv.writerow(["", "", "", "", st_matches_count, ss_matches_count, ""])
	### close filehandlers
	close (out_file)
	close (out_file_passed)
	close (out_file_filterout)
	close (out_file_long)
	close (out_file_sum_csv)
	close (out_file_sum)
	### final sum
	print(f"\n{'='*50}\n")
	print(f"Total lines  : {total_processed_lines}")
	print(f"  Num ref    : {group_counter}\n")
	print(f"  num SS     : {total_ss_score}")
	print(f"  num St     : {total_st_score}")
	print(f"  num BB     : {total_bb_specific}")
	print(f"  num ignored: {total_filtered_lines}\n")
	print(f"{'='*50}")



def main():
	parser = argparse.ArgumentParser(description='Stat St/BB/SS SNPs')
	parser.add_argument('vcf', type=str, help='VCF input')
	parser.add_argument('out', type=str, help='Out prefix')
	args = parser.parse_args()
	process_vcf(args.invcf, args.outpfx)

if __name__ == '__main__':
	main()
