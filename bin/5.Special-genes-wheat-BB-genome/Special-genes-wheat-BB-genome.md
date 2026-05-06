# Identification Pipeline of St and BB Subgenome-Specific Genes

## Overview

This pipeline identifies BB- and St-subgenome-specific genes in wheat (Triticum aestivum). After initial screening with OrthoFinder, candidate genes are aligned to multiple reference genomes (AA, DD, SS, St, EE, RR, TT) using GMAP. BB and St-specific gene sets are finally filtered using a coverage threshold of ≥95%

## Dependencies

| Software | Version | Description |
|------|------|------|
| Perl | v5.38.2 | Run scripts |
| SLURM | - | Job scheduling and submission |
| OrthoFinder | v3.1.0 | Orthologous group analysis |
| SeqKit | v2.10.1 | Sequence filtering and extraction |
| GMAP/GSNAP | v20250731 | Genomic alignment |
| Python 3 | v3.12.3 | Run scripts |

## Script Description

The following 4 scripts are used in this pipeline:
|Script | Language | Function |
|------|------|------|
| `1.extract_primary.pl` | Perl | Extract first isoform per gene from FASTA |
| `2.run_gmap.slurm` | Bash (SLURM) | Run GMAP alignment for small genomes (<4Gb) |
| `3.run_gmapl.slurm` | Bash (SLURM) | Run GMAP alignment for large genomes (>4Gb) |
| `4.select_mRNA_by_coverage.py` | Python | Filter GMAP output based on coverage threshold |

## Directory Structure

```bash
/home/mengdi_wang/Documents/Orthofinder_BB-St_Special_gene/
├── run_code.md
├── 1.extract_primary.pl
├── 2.run_gmap.slurm
├── 3.run_gmapl.slurm
├── 4.select_mRNA_by_coverage.py
```


## Pipeline Steps

### Variables

```
rundir=/home/mengdi_wang/Documents/Orthofinder_BB-St_Special_gene
```

### 1.Extract primary transcripts

```
mkcd mRNA.pep.primary
perl 1.extract_primary.pl --input Genome.pep.fa -s first --output Genome.pep.primary.fa
#  output directory: Tr_aestivum_AA_BB_DD.pep.primary
#  output files: Triticum_aestivum.pep.AA.primary.fa Triticum_aestivum.pep.BB.primary.fa Triticum_aestivum.pep.DD.primary.fa

perl 1.extract_primary.pl --input Genome.RNA.fa -s first --output Genome.mRNA.primary.fa
#  output directory: Tr_aestivum_BB.mRNA.primary
#  output files: Triticum_aestivum.mRNA.BB.primary.fa
```

### 2.OrthoFinder Analysis

```
mkcd Orthofinder_6AA_6BB_6DD
orthofinder -f $rundir/mRNA.pep.primary/Tr_aestivum_AA_BB_DD.pep.primary -S diamond -M msa -A mafft -T fasttree -t 8 -a 8
```

### 3.Extract BB-specific orthogroup sequences

```
mkcd BB_uniq_gene
cp $rundir/Orthofinder_6AA_6BB_6DD/Results/Orthogroups/Orthogroups.txt ./
awk '!/TraesCS[1-7][AD]/ && /TraesCS[1-7]B/' Orthogroups.txt > BB.uniq.orthogroups.txt
grep -o "TraesCS[^ ]*" BB.uniq.orthogroups.txt > BB.uniq.orthogroups.first_gene_ids.2452.ID
seqkit grep -n -f BB.uniq.orthogroups.first_gene_ids.2452.ID $rundir/mRNA.pep.primary/Tr_aestivum_BB.mRNA.primary.fa/Triticum_aestivum.mRNA.BB.primary.fa > BB.uniq.orthogroups.first_gene_ids.2452.fa
seqkit grep -s -r -v -p "[^ATCGatcg]" BB.uniq.orthogroups.first_gene_ids.2452.fa > BB.uniq.orthogroups.first_gene_ids.mRNA.clean_ATCG.2443.fa
# Input: 2,452 sequences | Output: 2,443 sequences (9 low-quality sequences removed)
```

### 4.Build GMAP Databases

```
mkcd GMAP_Database_Build
gmap_build -d AA_GMAP_Database -D GMAP_Database/ AA_GCF_003073215.2_Tu2.1_genomic.fna
gmap_build -d DD_GMAP_Database -D GMAP_Database/ DD_GCF_002575655.3_Aet_v6.0_genomic.fna
gmap_build -d SS_GMAP_Database -D GMAP_Database/ SS.GWHBFXR00000000.1.genome.fasta
gmap_build -d St_GMAP_Database -D GMAP_Database/ E.sibiricus.St.fa
gmap_build -d EE_GMAP_Database -D GMAP_Database/ GWHABKY00000000.EE.genome.fasta
gmap_build -d RR_GMAP_Database -D GMAP_Database/ Secale.cereale_GCA_902687465.1_Rye_Lo7_2018_v1p1p1_genomic.fa
gmap_build -d TT_GMAP_Database -D GMAP_Database/ lpAmbMuti1_1.curated_primary.TT.fa
```

### 5.GMAP AA DD genomes

```
mkcd BB_uniq_gene_GMAP_AA_DD

# Note: Use gmap for genomes < 4Gb; Use gmapl for genomes > 4Gb
cp $rundir/BB_uniq_gene/BB.uniq.orthogroups.first_gene_ids.mRNA.clean_ATCG.2443.fa ./

# Align against AA genome (use gmapl, as genome size > 4Gb)
sbatch 3.run_gmapl.slurm BB.uniq.orthogroups.first_gene_ids.mRNA.clean_ATCG.2443.fa $rundir/GMAP_Database_Build/GMAP_Database/ BB.uniq.AA_GMAP_Database_gff3 AA_GMAP_Database
# Output directory: BB.uniq.AA_GMAP_Database_gff3

# Align against DD genome (use gmap, as genome size < 4Gb)
sbatch 2.run_gmap.slurm BB.uniq.orthogroups.first_gene_ids.mRNA.clean_ATCG.2443.fa $rundir/GMAP_Database_Build/GMAP_Database/ BB.uniq.DD_GMAP_Database_gff3 DD_GMAP_Database
# Output directory: BB.uniq.DD_GMAP_Database_gff3

# Filter genes with coverage ≥ 95%
python3 4.select_mRNA_by_coverage.py BB.uniq.AA_GMAP_Database_gff3 BB.uniq.AA_GMAP_Database.gff3 BB.uniq.AA_GMAP_Database.log 95
# Output files:
#  BB.uniq.AA_GMAP_Database.gff3
#  BB.uniq.AA_GMAP_Database.log
#  BB.uniq.AA_GMAP_Database.mRNA.coverage.95.gff3
#  BB.uniq.AA_GMAP_Database.mRNA.coverage.95.ID.list
python3 4.select_mRNA_by_coverage.py BB.uniq.DD_GMAP_Database_gff3 BB.uniq.DD_GMAP_Database.gff3 BB.uniq.DD_GMAP_Database.log 95
# Output files:
#  BB.uniq.DD_GMAP_Database.gff3
#  BB.uniq.DD_GMAP_Database.log
#  BB.uniq.DD_GMAP_Database.mRNA.coverage.95.gff3
#  BB.uniq.DD_GMAP_Database.mRNA.coverage.95.ID.list

grep "^>" BB.uniq.orthogroups.first_gene_ids.mRNA.clean_ATCG.2443.fa | sed 's/^>//' > BB.uniq.orthogroups.first_gene_ids.mRNA.clean_ATCG.2443.ID
cat BB.uniq.DD_GMAP_Database.mRNA.coverage.95.ID.list BB.uniq.AA_GMAP_Database.mRNA.coverage.95.ID.list | sort | uniq > AA-DD.GMAP.merged.txt
grep -v -f AA-DD.GMAP.merged.txt BB.uniq.orthogroups.first_gene_ids.mRNA.clean_ATCG.2443.ID | grep -v '^$' > BB_uniq.remove-AA-DD.676.ID
seqkit grep -n -f BB_uniq.remove-AA-DD.676.ID $rundir/mRNA.pep.primary/Tr_aestivum_BB.mRNA.primary.fa/Triticum_aestivum.mRNA.BB.primary.fa > BB_uniq.remove-AA-DD.676.fa
```

### 6.GMAP SS St genomes

```
mkcd BB_uniq_gene_GMAP_SS_St
cp $rundir/BB_uniq_gene_GMAP_AA_DD/BB_uniq.remove-AA-DD.676.fa ./

# Align against SS genome (use gmap, as genome size < 4Gb)
sbatch 2.run_gmap.slurm BB_uniq.remove-AA-DD.676.fa $rundir/GMAP_Database_Build/GMAP_Database/ BB.uniq.SS_GMAP_Database_gff3 SS_GMAP_Database
# Output directory: BB.uniq.SS_GMAP_Database_gff3

# Align against St genome (use gmap, as genome size < 4Gb)
sbatch 2.run_gmap.slurm BB_uniq.remove-AA-DD.676.fa $rundir/GMAP_Database_Build/GMAP_Database/ BB.uniq.St_GMAP_Database_gff3 St_GMAP_Database
# Output directory: BB.uniq.St_GMAP_Database_gff3

# Filter genes with coverage ≥ 95%
python3 4.select_mRNA_by_coverage.py BB.uniq.SS_GMAP_Database_gff3 BB.uniq.SS_GMAP_Database.gff3 BB.uniq.SS_GMAP_Database.log 95
# Output files:
#  BB.uniq.SS_GMAP_Database.gff3
#  BB.uniq.SS_GMAP_Database.log
#  BB.uniq.SS_GMAP_Database.mRNA.coverage.95.gff3
#  BB.uniq.SS_GMAP_Database.mRNA.coverage.95.ID.list
python3 4.select_mRNA_by_coverage.py BB.uniq.St_GMAP_Database_gff3 BB.uniq.St_GMAP_Database.gff3 BB.uniq.St_GMAP_Database.log 95
# Output files:
#  BB.uniq.St_GMAP_Database.gff3
#  BB.uniq.St_GMAP_Database.log
#  BB.uniq.St_GMAP_Database.mRNA.coverage.95.gff3
#  BB.uniq.St_GMAP_Database.mRNA.coverage.95.ID.list
grep -v -f BB.uniq.SS_GMAP_Database.mRNA.coverage.95.ID.list BB.uniq.St_GMAP_Database.mRNA.coverage.95.ID.list | grep -v '^$' > BB_uniq.remove-AA-DD-SS.52.ID
seqkit grep -n -f BB_uniq.remove-AA-DD-SS.52.ID $rundir/mRNA.pep.primary/Tr_aestivum_BB.mRNA.primary.fa/Triticum_aestivum.mRNA.BB.primary.fa > BB_uniq.remove-AA-DD-SS.52.fa
```

### 7.GMAP EE RR TT genomes

```
mkcd BB_uniq_gene_GMAP_EE_RR_TT
cp $rundir/BB_uniq_gene_GMAP_SS_St/BB_uniq.remove-AA-DD-SS.52.fa ./

# Align against EE genome (use gmapl, as genome size > 4Gb)
sbatch 3.run_gmapl.slurm BB_uniq.remove-AA-DD-SS.52.fa $rundir/GMAP_Database_Build/GMAP_Database/ BB.uniq.EE_GMAP_Database_gff3 EE_GMAP_Database
# Output directory: BB.uniq.EE_GMAP_Database_gff3

# Align against RR genome (use gmapl, as genome size > 4Gb)
sbatch 3.run_gmapl.slurm BB_uniq.remove-AA-DD-SS.52.fa $rundir/GMAP_Database_Build/GMAP_Database/ BB.uniq.RR_GMAP_Database_gff3 RR_GMAP_Database
# Output directory: BB.uniq.RR_GMAP_Database_gff3

# Align against TT genome (use gmapl, as genome size > 4Gb)
sbatch 3.run_gmapl.slurm BB_uniq.remove-AA-DD-SS.52.fa $rundir/GMAP_Database_Build/GMAP_Database/ BB.uniq.TT_GMAP_Database_gff3 TT_GMAP_Database
# Output directory: BB.uniq.TT_GMAP_Database_gff3

# Filter genes with coverage ≥ 95%
python3 4.select_mRNA_by_coverage.py BB.uniq.EE_GMAP_Database_gff3 BB.uniq.EE_GMAP_Database.gff3 BB.uniq.EE_GMAP_Database.log 95
# Output files:
#  BB.uniq.EE_GMAP_Database.gff3
#  BB.uniq.EE_GMAP_Database.log
#  BB.uniq.EE_GMAP_Database.mRNA.coverage.95.gff3
#  BB.uniq.EE_GMAP_Database.mRNA.coverage.95.ID.list
python3 4.select_mRNA_by_coverage.py BB.uniq.RR_GMAP_Database_gff3 BB.uniq.RR_GMAP_Database.gff3 BB.uniq.RR_GMAP_Database.log 95
# Output files:
#  BB.uniq.RR_GMAP_Database.gff3
#  BB.uniq.RR_GMAP_Database.log
#  BB.uniq.RR_GMAP_Database.mRNA.coverage.95.gff3
#  BB.uniq.RR_GMAP_Database.mRNA.coverage.95.ID.list
python3 4.select_mRNA_by_coverage.py BB.uniq.TT_GMAP_Database_gff3 BB.uniq.TT_GMAP_Database.gff3 BB.uniq.TT_GMAP_Database.log 95
# Output files:
#  BB.uniq.TT_GMAP_Database.gff3
#  BB.uniq.TT_GMAP_Database.log
#  BB.uniq.TT_GMAP_Database.mRNA.coverage.95.gff3
#  BB.uniq.TT_GMAP_Database.mRNA.coverage.95.ID.list
```

### 8.Summary output

```
mkcd BB_uniq_gene_Summary
cp $rundir/BB_uniq_gene_GMAP_SS_St/BB_uniq.remove-AA-DD-SS.52.ID ./
cp $rundir/BB_uniq_gene_GMAP_EE_RR_TT/BB.uniq.EE_GMAP_Database.mRNA.coverage.95.ID.list ./
cp $rundir/BB_uniq_gene_GMAP_EE_RR_TT/BB.uniq.RR_GMAP_Database.mRNA.coverage.95.ID.list ./
cp $rundir/BB_uniq_gene_GMAP_EE_RR_TT/BB.uniq.TT_GMAP_Database.mRNA.coverage.95.ID.list ./

St="BB_uniq.remove-AA-DD-SS.52.ID"; EE="BB.uniq.EE_GMAP_Database.mRNA.coverage.95.ID.list"; RR="BB.uniq.RR_GMAP_Database.mRNA.coverage.95.ID.list"; TT="BB.uniq.TT_GMAP_Database.mRNA.coverage.95.ID.list"; \
{ echo -e "Gene_ID\tSS\tBB\tSt\tEE\tRR\tTT"; \
grep '^TraesCS' "$St" | sort -u | while read gene; do \
  printf "%s\t0\t1\t1" "$gene"; \
  for g in "$EE" "$RR" "$TT"; do \
    grep -qxF "$gene" "$g" 2>/dev/null && printf "\t1" || printf "\t0"; \
  done; \
  echo; \
done; } > merged.SS_BB_St_EE_RR_TT.tsv
```
