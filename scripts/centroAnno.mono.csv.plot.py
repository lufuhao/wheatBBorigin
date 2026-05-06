#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import matplotlib.pyplot as plt
import os
from collections import Counter
import random
# matplotlib.rcParams['font.family'] = 'Times New Roman'
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase
# from matplotlib import colormaps
# cmap = colormaps.get_cmap("viridis")
cmap = plt.colormaps["viridis"]
from mpl_toolkits.axes_grid1 import make_axes_locatable


fig_fontsize = 18
plt.rcParams['font.size'] = fig_fontsize



def find_mode(lst):
    if not lst:
        return None, 0
    counter = Counter(lst)
    mode, count = counter.most_common(1)[0]
    return mode, count



def draw_mono_fig(bed_file, out_fig):
    data = []
    with open(bed_file) as f:
        lines = f.readlines()
        for line in lines:
            info = line.strip('\n').split('\t')
            try:
                this_data = [info[0], int(info[1]), int( info[2] ), int( info[3] )]
                data.append( this_data )
            except:
                print( f"Warning Line: {line}" )
#    print (data)
    hor_lengths = [row[3] for row in data]
#    print (hor_lengths)
    min_len, max_len = min(hor_lengths), max(hor_lengths)
    norm = mcolors.Normalize(vmin=min_len, vmax=max_len)
    if min_len == max_len:
        norm = mcolors.Normalize(vmin=min_len, vmax=max_len + 1)
    all_seqs = sorted(set(row[0] for row in data))
    seq_y_offset = {seq: i * 20 for i, seq in enumerate(all_seqs)}
    fig, ax = plt.subplots(figsize=(24, 9))
    for row in data:
        seq, start, end, hor_len = row
        ypos = seq_y_offset[seq]
        color = cmap(norm(hor_len))
        ax.broken_barh([(start, end - start)], (ypos, 8), facecolors=color)
    for seq in all_seqs:
        ypos = seq_y_offset[seq]
        ax.text(0, ypos + 4, seq, va='center', ha='right', fontsize=12, fontweight='bold')
    ax.set_ylim(-10, max(seq_y_offset.values()) + 20)
    ax.set_xlim(0, max(row[2] for row in data) + 300)
    ax.set_xlabel("Genomic Coordinate")
    ax.set_yticks([])
    ax.set_title("Tandem Repeat Unit Visualization")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cb = ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cb.set_label("Tandem Repeat Unit Length (bp)")
    plt.tight_layout()
    plt.savefig(out_fig, bbox_inches='tight')


def main(incsv, repeat_region_minlen, outsvg):
#    repeat_region_minlen = 100
    out_put_mono_bed = f'{outsvg}.bed'
    w_file = open(out_put_mono_bed, 'w')
    with open (incsv, "r") as f:
        region_nm = ""
        region_st = 0
        region_ed = 0
        rep_len_list = []
        linenum=0
        for line in f:
#            print (line)
            linenum+=1
            if 'name' in line:
                continue
            info = line.strip('\n').split(',')
            if region_nm == "":
                region_nm=info[0]
#                print (linenum)
#                print (line)
                region_st=int(info[2])
                region_ed=int(info[3])
                rep_len_list.append( int( info[-1] ) )
                continue
            elif region_nm != info[0]:
                rep_len = find_mode( rep_len_list )[0]
                w_file.write(f'{region_nm}\t{region_st}\t{region_ed}\t{rep_len}\n')
                rep_len_list = []
                region_nm = info[0]
                region_st=int(info[2])
                region_ed=int(info[3])
                rep_len_list.append( int( info[-1] ) )
                continue
            st = int(info[2])
            ed = int(info[3])
            # iden = float(info[-2])
            if abs(region_ed - st) < 100:
                region_ed = ed
                rep_len_list.append( int( info[-1] ) )
            else:
                rep_len = find_mode( rep_len_list )[0]
                if abs(region_ed - region_st) >= repeat_region_minlen:
                    w_file.write(f'{region_nm}\t{region_st}\t{region_ed}\t{rep_len}\n')
                rep_len_list = []
                region_st = st
                region_ed = ed
                rep_len_list.append( int( info[-1] ) )
        rep_len = find_mode( rep_len_list )[0]
        if abs(region_ed - region_st) >= repeat_region_minlen:
            w_file.write(f'{region_nm}\t{region_st}\t{region_ed}\t{rep_len}\n')
    # sort_bed_file(out_put_mono_bed, out_put_mono_sorted_bed)
    ########################## get repeat regions ##########################
    ########################## Visualization ##########################
    draw_mono_fig(out_put_mono_bed, outsvg)
    ########################## Visualization ##########################



if __name__ == "__main__":
    if len( sys.argv ) != 4:
        print(f"Usage: python {sys.argv[0]} mono.csv min_length mono.svg")
        print("    A script to plot centroAnno's mono csv")
    else:
        main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
