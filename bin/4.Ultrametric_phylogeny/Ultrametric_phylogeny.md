# Ultrametric Tree Inference Pipeline

---
## Overview

This pipeline infers ultrametric trees from single-copy orthologous genes identified by OrthoFinder across **BB**, **SS**, **StSt(Elymus sibiricus)** subgenome, *Oryza sativa* (Os), and *Brachypodium distachyon* (Bd).
---

## Dependencies

| Software | Version | Description |
|------|------|------|
| Perl | v5.38.2 | Run scripts |
| SLURM | - | Job scheduling and submission |
| OrthoFinder | v3.1.0 | Orthologous group analysis |
| MUSCLE | v5.3 | Sequence alignment |
| trimAl | v1.5.0 | Alignment trimming |
| IQ-TREE | v3.0.1 | Phylogenetic inference |
| R | 4.5.0 | Tree processing |
| APE (R package) | v5.8.1 | Chronos ultrametric conversion |
| DensiTree | v3.1.0 | Tree visualization |

## Script Description

The following 7 scripts are used in this pipeline:
|Script | Language | Function |
|------|------|------|
| `1.extract_primary.pl` | Perl | Extract first isoform per gene from FASTA |
| `2.muscle_align.slurm` | Bash (SLURM) | Batch multiple sequence alignment using MUSCLE |
| `3.trimal_trim_align.slurm` | Bash (SLURM) | Batch trimming of aligned sequences using trimAl |
| `4.iqtree3_construct_tree.slurm` | Bash (SLURM) | Maximum likelihood tree construction using IQ-TREE 3 |
| `5.check_clade_StSt_SS_BB.R` | R | Check and filter trees based on StSt, SS, and BB clades |
| `6.multi2di_root.R` | R | multi2di and outgroup rooting |
| `7.chronos_time_calibrated.R` | R | Generate ultrametric trees with APE::chronos |


## Directory Structure

```bash
/home/mengdi_wang/Documents/Densitree/Orthofinder_Os.Bd.SS.BB.St_ES_Densitree_5/
├── Ultrametric_phylogeny.md
├── 1.extract_primary.pl
├── 2.muscle_align.slurm
├── 3.trimal_trim_align.slurm
├── 4.iqtree3_construct_tree.slurm
├── 5.check_clade_StSt_SS_BB.R
├── 6.multi2di_root.R
└── 7.chronos_time_calibrated.R
```


## Pipeline Steps

### Variables

```
rundir=/home/mengdi_wang/Documents/Densitree/Orthofinder_Os.Bd.SS.BB.St_ES_Densitree_5
```

### 1 Extract Primary pep

```
mkcd Os.Bd.SS.BB.St_ES.pep.primary
perl 1.extract_primary.pl --input Genome.pep.fa -s first --output Genome.pep.primary.fa、
# output directory: Os.Bd.SS.BB.St_ES.pep
# output files:
#  Elymus_sibiricus.pep.St.primary.fa
#  Triticum_aestivum.pep.BB.primary.fa
#  Aegilops_speltoides.pep.SS.primary.fa
#  Oryza_sativa.pep.Os.primary.fa
#  Brachypodium_distachyon.pep.Bd.primary.fa
```

### 2.OrthoFinder Analysis

```
mkcd Orthofinder_Os.Bd.SS.BB.St_ES
orthofinder -f $rundir/Os.Bd.SS.BB.St_ES.pep.primary/Os.Bd.SS.BB.St_ES.pep -S diamond -M msa -A mafft -T fasttree -t 8 -a 8
```

### 3.Single-Copy Orthologue Processing

```
mkcd Single_Copy_2_Densitree
cp -r $rundir/Orthofinder_Os.Bd.SS.BB.St_ES/OrthoFinder/Results/Single_Copy_Orthologue_Sequences ./

# MUSCLE Alignment
sbatch 2.muscle_align.slurm
#output directory: Single_Copy_Orthologue_Sequences_aligned/

# trimAl Trimming
sbatch 3.trimal_trim_align.slurm -i Single_Copy_Orthologue_Sequences_aligned/ -o Single_Copy_Orthologue_Sequences_aligned_trimal/
#output directory: Single_Copy_Orthologue_Sequences_aligned_trimal/

# IQ-TREE3 Phylogenetic tree
sbatch 4.iqtree3_construct_tree.slurm
#output directory: Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3

mkdir -p Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile/
find Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3/ -type f -name "*.treefile" -exec cp {} Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile/ \;

for f in ./Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile/*; do
cat "$f"
done | sed '/^$/d' > Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile.txt

sed -E 's/\)[0-9]+(\.[0-9]+)?/)/g' Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile.txt | sed '/^$/d' > Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.tree

# Check and filter trees
Rscript 5.check_clade_StSt_SS_BB.R
#output files: Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.tree

# multi2di and outgroup rooting
Rscript 6.multi2di_root.R
#output files: Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.multi2di_root.tree

# ultrametric trees
Rscript 7.chronos_time_calibrated.R
#output files: Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.multi2di_root.chronos.tree

# Replace Sequence IDs
sed -E 's/Bradi[[:alnum:]_.]+/Bd/g; s/GWHPBFXR[[:alnum:]_.]+/SS/g; s/Es[[:alnum:]_.]+/StSt/g; s/TraesCS[[:alnum:]_.]+/BB/g; s/LOC[[:alnum:]_.]+/Os/g' Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.multi2di_root.chronos.tree > Single_Copy_Orthologue_Sequences_aligned_trimal_iqtree3_treefile_NO_bootstrap.not-StSt-SS-BB-clade.multi2di_root.chronos.replace_ID.txt
