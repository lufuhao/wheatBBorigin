# Single-Copy Ortholog SNP Identification Pipeline

---

## Overview

This pipeline identifies SNPs from single-copy orthologous genes across **H. vulgare** (HH), **StSt(Elymus sibiricus)** , **RR**, **EE**, **TT**, **SS**, and **BB** using OrthoFinder, MAFFT, and jvarkit.

---

## Dependencies

| Software | Version | Description |
|------|------|------|
| Perl | v5.38.2 | Run scripts |
| Bash | - | Run scripts |
| SLURM | - | Job scheduling and submission |
| OrthoFinder | v3.1.0 | Orthologous group analysis |
| MAFFT | v7.526 | Sequence alignment |
| jvarkit |v2024.08.25 | msa2vcf for variant calling |
| Python 3 | v3.12.3 | Run scripts |

## Script Description

The following 7 scripts are used in this pipeline:

|Script | Language | Function |
|------|------|------|
| `1.extract_primary.pl` | Perl | Extract first isoform per gene from FASTA |
| `2.replace_protein_with_mRNA.sh` | Bash | Replace protein sequences with corresponding mRNA sequences |
| `3.mafft_align.slurm` | Bash (SLURM) | Batch mRNA sequence alignment using MAFFT |
| `4.MAFFT_to_VCF-c.py` | Python | Convert alignments to VCF using jvarkit/msa2vcf (BB as reference) |
| `5.replace_chrUn_TraseCS.py` | Python | Replace BB sequence IDs in VCF files |
| `6.Non-ACGT-character.VCF-classify.py` | Python | Classify variants: SNP, MNP, INDEL, and special characters |
| `7.extract_VCF-SNP_to_table.py` | Python | Extract SNPs to summary table |


## Directory Structure

```bash
/home/mengdi_wang/Documents/Orthofinder_HH.TT.EE.RR.BB.SS.StSt_ES_7/
├── 1.extract_primary.pl
├── 2.replace_protein_with_mRNA.sh
├── 3.mafft_align.slurm
├── 4.MAFFT_to_VCF-c.py
├── 5.replace_chrUn_TraseCS.py
├── 6.Non-ACGT-character.VCF-classify.py
├── 7.extract_VCF-SNP_to_table.py
└── HH.EE.RR.TT.SS.BB.St.Single_Copy_SNPs.md
```

## Pipeline Steps

### Variables

```
rundir=/home/mengdi_wang/Documents/Orthofinder_HH.TT.EE.RR.BB.SS.StSt_ES_7
```

### 1.Extract Primary Transcripts (Protein & mRNA)

```
mkcd mRNA.pep.primary
# Extract primary protein sequences
perl 1.extract_primary.pl --input Genome.pep.fa -s first --output Genome.pep.primary.fa
#  output directory: HH.TT.EE.RR.BB.SS.St.pep.primary
#  output files: Aegilops_mutica.pep.TT.primary.fa Elymus_sibiricus.pep.St.primary.fa Secale_cereale.pep.RR.primary.fa Triticum_aestivum.pep.BB.primary.fa Aegilops_speltoides.pep.SS.primary.fa Hordeum_vulgare.pep.HH.primary.fa Thinopyrum_elongatum.pep.EE.primary.fa

# Extract primary mRNA sequences
perl 1.extract_primary.pl --input Genome.RNA.fa -s first --output Genome.mRNA.primary.fa
#  output directory: HH.TT.EE.RR.BB.SS.St.mRNA.primary
#  output files: Aegilops_mutica.mRNA.TT.primary.fa Elymus_sibiricus.mRNA.St.primary.fa Secale_cereale.cdna.RR.primary.fa Triticum_aestivum.mRNA.BB.primary.fa Aegilops_speltoides.mRNA.SS.primary.fa HvulgareMorex_702_V3.transcript.primary.fa Thinopyrum_elongatum.mRNA.EE.primary.fa
```

### 2.OrthoFinder Analysis

```
mkcd Orthofinder_HH.TT.EE.RR.BB.SS.St
orthofinder -f $rundir/mRNA.pep.primary/HH.TT.EE.RR.BB.SS.St.pep.primary -S diamond -M msa -A mafft -T fasttree -t 8 -a 8
```

### 3.SNP Calling from Single-Copy Orthologs

```
mkcd Single_Copy_Orthologue_Sequences_SNP

# Replace protein sequences with mRNA sequences
bash 2.replace_protein_with_mRNA.sh $rundir/Orthofinder_HH.TT.EE.RR.BB.SS.St/Results/Single_Copy_Orthologue_Sequences/ $rundir/mRNA.pep.primary/HH.TT.EE.RR.BB.SS.St.mRNA.primary/ Single_Copy_Orthologue_Sequences_mRNA/
# output directory: Single_Copy_Orthologue_Sequences_mRNA/

# MAFFT alignment
sbatch 3.mafft_align.slurm
# output directory: Single_Copy_Orthologue_Sequences_mRNA_MAFFT

# Convert alignments to VCF using jvarkit/msa2vcf (BB as reference)
python3 4.MAFFT_to_VCF-c.py Single_Copy_Orthologue_Sequences_mRNA_MAFFT Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF --jar ~/Programs/jvarkit/jvarkit.jar
# output directory: Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF/

# Replace BB sequence IDs in VCF
python3 5.replace_chrUn_TraseCS.py -i Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF/ -o Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF_repalce_BB-ID/
# output directory: Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF_repalce_BB-ID

# Classify variants: SNP, MNP, INDEL, special characters
python3 6.Non-ACGT-character.VCF-classify.py Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF_repalce_BB-ID Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF_repalce_BB-ID_Classify
# Output subdirectories:
# INDEL/
# MNP/
# SNP/
# special_characters.vcf
# classification_log.txt

Extract SNPs to summary table
python3 7.extract_VCF-SNP_to_table.py SNP/ Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF_repalce_BB-ID_Classify_SNP.tsv
# output files: Single_Copy_Orthologue_Sequences_mRNA_MAFFT_VCF_repalce_BB-ID_Classify_SNP.tsv
```
