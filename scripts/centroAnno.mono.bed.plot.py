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



if __name__ == "__main__":
    if len( sys.argv ) != 3:
        print(f"Usage: python {sys.argv[0]} mono.csv mono.svg")
        print("    A script to plot centroAnno's mono bed")
    else:
        draw_mono_fig(sys.argv[1], sys.argv[2])
