# Comparative genomics among St.SS.BB

---

## Overview

This pipeline performs gene preprocessing, ortholog detection, synteny filtering, and unique 1:1:1 ortholog identification across **BB**, **SS**, and **StSt(Elymus sibiricus)** subgenomes using JCVI tools.

---

## Dependencies
| Software | Version | Description |
|------|------|------|
| Python 3 | v3.12.3 | Run JCVI modules and scripts |
| SeqKit | v2.10.1 | Sequence filtering |
| JCVI | v1.5.8 | Ortholog identification and synteny analysis |

## Script Description

|Script | Language | Function |
|------|------|------|
| `St_BB_SS.unique.match.filter.py` | Python | Identify unique 1:1:1 orthologs among BB, SS, and St genomes |


## Pipeline Steps

### Variables
```
rundir=/home/mengdi_wang/Documents/JCVI
```

### 1.File Preprocessing
```bash
# (1).GFF_2_BED
python3 -m jcvi.formats.gff bed --type=mRNA iwgsc_refseqv2.1_annotation.BB.gff3 --key=ID -o BB.bed
python3 -m jcvi.formats.gff bed --type=mRNA SS.GWHBFXR00000000.1.gff --key=ID -o SS.bed
python3 -m jcvi.formats.gff bed --type=mRNA E.sibiricus_St.gff --key=ID -o St.bed

# (2)Deduplicate BED files
python3 -m jcvi.formats.bed uniq BB.bed
python3 -m jcvi.formats.bed uniq SS.bed
python3 -m jcvi.formats.bed uniq St.bed

# (3)Extract protein sequences based on unique BED gene IDs
seqkit grep -f <(cut -f 4 BB.uniq.bed ) iwgsc_refseqv2.1_annotation_HC_pep.BB.fa | seqkit seq -i > BB.uniq.pep
seqkit grep -f <(cut -f 4 SS.uniq.bed ) SS.GWHBFXR00000000.1.Protein.fa | seqkit seq -i > SS.uniq.pep
seqkit grep -f <(cut -f 4 St.uniq.bed ) E.sibiricus_St.pep.fa | seqkit seq -i > St.uniq.pep
```

### 2.JCVI Synteny & Ortholog Analysis
```bash
mkdir -p Run_JCVI && cd Run_JCVI
cp $rundir/BB.uniq.bed ./BB.bed
cp $rundir/BB.uniq.pep ./BB.pep
cp $rundir/SS.uniq.bed ./SS.bed
cp $rundir/SS.uniq.pep ./SS.pep
cp $rundir/St.uniq.bed ./St.bed
cp $rundir/St.uniq.pep ./St.pep

# Run jcvi
python3 -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names SS BB
python3 -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names BB St
python3 -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names SS St

python3 -m jcvi.compara.synteny screen --minspan=30 --simple SS.BB.anchors SS.BB.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple BB.St.anchors BB.St.anchors.new
python3 -m jcvi.compara.synteny screen --minspan=30 --simple SS.St.anchors SS.St.anchors.new

# Identify Unique 1:1:1 Orthologs
python3 St_BB_SS.unique.match.filter.py SS.BB.anchors.new BB.St.anchors.new SS.St.anchors.new BB-SS-St_unique.match.filter_1-1-1.txt
