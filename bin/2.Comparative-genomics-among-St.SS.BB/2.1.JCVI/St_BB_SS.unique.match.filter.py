#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os

### read JCVI syn: 2-3 columns, each ### line will count for block number
def readSyn(input_file, idx):
    id1_dict = {}
    id2_dict = {}
    id1_rm={}
    id2_rm={}
    blockno=0
    linenum1=0
    linenum2=0
    with open(input_file) as f:
        for line in f:
            linenum1+=1
            if not line:
                continue
            if line.startswith("###"):
                blockno+=1
                continue
            linenum2+=1
            line = line.strip()
            parts = line.split("\t")
            if len(parts) < 2:
                print (f"Error: invalid line({linenum1}): {line}", file=sys.stderr)
                continue
            id1, id2 = parts[0], parts[1]
            if id1 in id1_dict:
                id1_rm[id1]=1
                id2_rm[id2]=1
                continue
            if id2 in id2_dict:
                id1_rm[id1]=1
                id2_rm[id2]=1
                continue
            id1_dict[id1]={}
            id1_dict[id1]['no']=blockno
            id1_dict[id1]['syn']=id2
            id2_dict[id2]={}
            id2_dict[id2]['no']=blockno
            id2_dict[id2]['syn']=id1
    numt1=len(id1_dict.keys())
    numt2=len(id2_dict.keys())
    numr1=len(id1_rm.keys())
    numr2=len(id2_rm.keys())
    print (f"Info: file {input_file}")
    print (f"    Total lines: {linenum1}")
    print (f"    Valid lines: {linenum2}")
    print (f"    ID1 dict keys: {numt1}")
    print (f"    ID2 dict keys: {numt2}")
    print (f"    ID1 remove keys: {numr1}")
    print (f"    ID2 remove keys: {numr2}")
    for id3 in id1_rm:
        if id3 in id1_dict:
            idtmp=id1_dict[id3]['syn']
            if idtmp in id2_dict:
                del id2_dict[idtmp]
        if id3 in id1_dict:
            del id1_dict[id3]
    for id4 in id2_rm:
        if id4 in id2_dict:
            idtmp=id2_dict[id4]['syn']
            if idtmp in id1_dict:
                del id1_dict[idtmp]
        if id4 in id2_dict:
            del id2_dict[id4]
    numk1=len(id1_dict.keys())
    numk2=len(id2_dict.keys())
    print (f"    ID1 final keys: {numk1}")
    print (f"    ID2 final keys: {numk2}")
    if numk1 != numk2:
        print (f"Error: {numk1} !=  {numk2}", file=sys.stderr)
        for id5 in sorted(id1_dict.keys()):
            id6=id1_dict[id5]['syn']
            if id6 not in id2_dict:
                print(f"Debug1: check {id6}", file=sys.stderr)
        for id7 in sorted(id2_dict.keys()):
            id8=id2_dict[id7]['syn']
            if id8 not in id1_dict:
                print(f"Debug2: check {id7}", file=sys.stderr)
        sys.exit(100)
    if idx==1:
        return id1_dict
    elif idx==2:
        return id2_dict
    else:
        print ("Error: invalid index for dict return", file=sys.stderr)
        return {}



def writeSyn(bb2st, bb2ss, ss2st, outfile):
    outnum=0
    with open(outfile, "w") as out:
        out.write("St\tBB\tSS\tBB2ST_BLK_NO\tBB2SS_BLK_NO\tSS2ST_BLK_NO\n")
        for bbid in sorted(bb2st.keys()):
            stid=bb2st[bbid]['syn']
            if bbid in bb2ss:
                ssid=bb2ss[bbid]['syn']
            else:
                continue
            if ssid in ss2st and ss2st[ssid]['syn'] == stid:
                out.write(f"{stid}\t{bbid}\t{ssid}\t{bb2st[bbid]['no']}\t{bb2ss[bbid]['no']}\t{ss2st[ssid]['no']}\n")
                outnum+=1
    print (f"Out file: {outfile}")
    print (f"    Total lines: {outnum}")
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("""
Usage: python full_synteny_pipeline_v5.py <SS.BB.anchors><BB.St.anchors>  <SS.St.anchors> <output_file>

1. Input files

- JCVI文件: SS->BB
- JCVI文件: BB->St
- JCVI文件: SS->St 参考文件

2. Process

- 用唯一 BB->St 映射和去重后的 SS->BB 文件
- 仅保留在去重后的 SS->St 参考文件中的配对

""")
        sys.exit(1)
    ss_bb_file = sys.argv[1]          # SS->BB 文件
    bb_st_file = sys.argv[2]          # BB->St 文件
    ss_st_ref_file = sys.argv[3]      # SS->St 参考文件
    final_out = sys.argv[4]           # 输出最终筛选结果
    bb2ss=readSyn(ss_bb_file, 2)
    bb2st=readSyn(bb_st_file, 1)
    ss2st=readSyn(ss_st_ref_file, 1)
    writeSyn(bb2st, bb2ss, ss2st, final_out)
    sys.exit(0)
