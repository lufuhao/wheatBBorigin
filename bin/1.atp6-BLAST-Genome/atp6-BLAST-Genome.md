# atp6 BLAST genomes

---
## Overview
The atp6 gene was used as the query to search and determine the copy number of its homologs in the following genomes: **Aegilops speltoides (SS)**,**Aegilops tauschii (DD)**, **Elymus nutans (StStYYHH)**, **Elymus sibiricus (StStHH)**, **Hordeum marinum (HH)**, **Pseudoroegneria libanotica (diploid StSt)**, **Secale cereale (RR)**, **Thinopyrum elongatum (EE)**,**Th. intermedium (JrJrJvJvStSt)**, **Thinopyrum  obtusiflorum**, **Triticum aestivum (AABBDD)**, **Triticum urartu (AA)**, using NCBI BLAST (v2.17.0).
---

## Dependencies
| Software | Version | Description |
|------|------|------|
| BLAST+ | v2.17.0 | sequence search |

```
# Build BLAST Database
makeblastdb -in genome.fa -dbtype nucl -out genome
# Note: Replace genome.fa with the corresponding reference genome files

# run BLAST
blastn -query atp6.nucl.fa -out atp6_BLAST_genome -db genome -outfmt '6 qseqid sseqid qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore'
```
