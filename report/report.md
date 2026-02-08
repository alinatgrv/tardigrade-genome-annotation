# Project 4. Tardigrades: from genestealers to space marines

## Structural and functional analysis of chromatin-associated proteins in *Ramazzottius varieornatus*

**Student:** Alina Tagirova  

## 1. Data sources

This project is based on genomic and proteomic data for the tardigrade *Ramazzottius varieornatus* (strain YOKOZUNA-1).

The assembled genome sequence of *R. varieornatus* was obtained from the NCBI database in FASTA format and used as a reference. Structural gene annotation was not performed de novo in this project; instead, precomputed gene predictions generated using the AUGUSTUS gene prediction software were used. These predictions included genomic coordinates of genes and transcripts (GFF format) as well as the corresponding predicted protein sequences (FASTA format).

In addition to genomic data, a set of short peptide sequences obtained by tandem mass spectrometry from the chromatin fraction of *R. varieornatus* cells was provided. These peptides represent fragments of proteins experimentally detected in association with DNA and were used to identify candidate chromatin-associated proteins.

The following commands were used to obtain and unpack the genome assembly:
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/949/185/GCA_001949185.1_Rvar_4.0/GCA_001949185.1_Rvar_4.0_genomic.fna.gz
gunzip -k GCA_001949185.1_Rvar_4.0_genomic.fna.gz
```

The files were renamed for clarity and placed into the project data directory:
```bash
mv augustus.whole.aa  Rvar_augustus_proteins.faa
mv augustus.whole.gff Rvar_augustus.gff
```

## 2. Structural annotation and protein prediction

Structural gene annotation of the *Ramazzottius varieornatus* genome was not performed de novo in this project. Instead, precomputed gene predictions generated using the AUGUSTUS gene prediction software were used. The predicted protein sequences were provided in FASTA format and served as the basis for all subsequent analyses.

The total number of predicted protein-coding genes was estimated by counting the number of entries in the protein FASTA file using the following command:

```bash
grep -c '^>' data/raw/Rvar_augustus_proteins.faa
```

This resulted in a total of **16,435 predicted proteins**, which is consistent with the expected gene count for a multicellular eukaryotic organism.

The predicted proteome contains more than sixteen thousand proteins, which is far beyond what can be experimentally validated on an individual basis. Therefore, it is necessary to apply additional filtering strategies to narrow down the list of potential candidates. In this project, candidate selection is achieved by integrating multiple independent lines of evidence, including experimental proteomics data (chromatin-associated peptides), predicted subcellular localization, functional domain annotation, and sequence similarity to known proteins. This stepwise reduction allows prioritization of a small subset of biologically relevant candidates for further experimental validation.

## 3. Physical localization (chromatin fraction) and mapping peptides to proteins

Considering that DNA is a major target of radiation damage, we hypothesized that *Ramazzottius varieornatus* may have proteins physically associated with chromatin that contribute to DNA protection and/or repair. To explore this, genomic predictions were combined with experimental proteomics data.

Chromatin-associated proteins were enriched experimentally by isolating the chromatin fraction, followed by tandem mass spectrometry. This produced a set of short peptide sequences. Each peptide corresponds to a fragment of a protein that was bound to DNA or co-purified during chromatin extraction. The goal of this step was to identify which predicted proteins from the *R. varieornatus* proteome contain these peptides.

### 3.1 Input peptide set

The peptide sequences were provided in FASTA format:

```bash
head -n 10 data/raw/peptides.fa
grep -c '^>' data/raw/peptides.fa
```
A total of 43 peptides were present in the input file.

### 3.2 Local alignment-based search (DIAMOND)

To match short peptides against the predicted proteome, a local alignment-based search was performed using DIAMOND. First, a local database was created from the predicted protein FASTA:
```bash
diamond makedb \
  --in data/raw/Rvar_augustus_proteins.faa \
  --db data/processed/Rvar_proteome
  ```
Then, the peptide sequences were searched against the local proteome database in tabular output format:

```bash
diamond blastp \
  -d data/processed/Rvar_proteome.dmnd \
  -q data/raw/peptides.fa \
  -f 6 qseqid sseqid evalue pident qcovhsp \
  --very-sensitive \
  -o results/peptide_to_protein.tsv
```
### 3.3 Results: peptide-to-protein mapping

The DIAMOND search produced 4 alignments (4 HSPs), corresponding to 4 peptide queries that matched the predicted proteome:
```bash
wc -l results/peptide_to_protein.tsv
head results/peptide_to_protein.tsv
```
To obtain the set of proteins supported by peptide evidence, unique subject protein IDs were extracted:

```bash
cut -f2 results/peptide_to_protein.tsv | sort | uniq | wc -l
cut -f2 results/peptide_to_protein.tsv | sort | uniq
```

In total, **2 predicted proteins** were supported by peptide matches:

- `g4106.t1`
- `g12510.t1`

All reported matches showed **100% sequence identity** and **100% query coverage (qcovhsp = 100)**, with E-values in the range of approximately *1e-6 to 1e-9*, indicating confident and unambiguous mapping of the experimental peptides to the predicted proteins.

These two proteins were therefore selected as **candidate chromatin-associated proteins** and carried forward for downstream subcellular localization prediction and functional annotation.

## 4. Localization prediction

To further prioritize chromatin-associated candidates, subcellular localization was predicted from protein sequences. Two complementary tools were used: WoLF PSORT (general localization prediction for animal proteins) and TargetP 2.0 (prediction of N-terminal targeting peptides and secretory signal peptides). Predictions were performed using the candidate protein FASTA file (`results/chromatin_candidates.faa`).

### 4.1 WoLF PSORT (Animal)

WoLF PSORT was run in **Animal** mode. The results were interpreted primarily to identify **non-secretory proteins** (i.e., proteins not predicted to follow the ER/Golgi/extracellular secretory pathway), which are more plausible candidates for nuclear localization.

- `g4106.t1`: strong secretory-pathway signal (E.R. 14.5; E.R._golg 9.5; extracellular 7; Golgi 3.5), consistent with a protein entering the ER/Golgi pathway.
- `g12510.t1`: predicted mainly as plasma membrane (plas 29; cyto 3) with no evidence for ER/Golgi/extracellular localization.

**Non-secretory candidates (WoLF PSORT):**
- `g12510.t1`

### 4.2 TargetP 2.0 
![TargetP prediction results for chromatin-associated protein candidates](figures/TargetP.jpg)

TargetP 2.0 was run in **Non-plant** mode to evaluate N-terminal presequences: secretory signal peptide (SP) and mitochondrial targeting peptide (mTP). Proteins classified as **Other** (no SP/mTP) were considered compatible with potential nuclear/cytosolic localization.

- `g12510.t1`: Prediction = Other (Likelihood: Other 0.9997; SP 0.0001; mTP 0.0002)
- `g4106.t1`: Prediction = Other (Likelihood: Other 0.7297; SP 0.2669; mTP 0.0034)

### 4.3 Integrated interpretation

Combining both predictors, `g4106.t1` is likely associated with the secretory pathway (WoLF PSORT: ER/Golgi/extracellular), despite being classified as “Other” by TargetP with a moderate SP probability. Therefore, it was excluded from nuclear localization candidates.

`g12510.t1` showed no evidence of secretory or mitochondrial targeting in TargetP and was predicted as non-secretory by WoLF PSORT. It was retained as the primary candidate for downstream functional annotation.

## 5. Pfam domain prediction
Pfam domain prediction with HMMER (hmmscan against the Pfam database) did not yield any significant domain hits for the initially selected candidates (`g4106.t1` and `g12510.t1`). This result is consistent with the possibility that these proteins are highly divergent or lineage-specific.


### Method update: revision of candidate selection strategy

The initial peptide-to-protein mapping performed with DIAMOND produced only a very small number of peptide matches (4 alignments) and resulted in just 2 candidate proteins. Because the peptide queries are short, this approach was considered overly stringent and potentially insensitive.

Therefore, the peptide-to-protein mapping step was repeated using BLAST+ in the `blastp-short` mode, which is specifically designed for short protein queries. A local BLAST database was created from the predicted *R. varieornatus* proteome, and peptide sequences were aligned against it. To reduce false positives, only matches covering the full peptide length (alignment length = query length) and having at least 90% identity were retained. This resulted in a refined candidate list of **20 proteins** supported by peptide evidence.

All subsequent analyses (subcellular localization prediction and functional annotation) were continued using this BLAST-derived candidate set.

```bash
makeblastdb \
  -in data/raw/Rvar_augustus_proteins.faa \
  -dbtype prot \
  -parse_seqids \
  -out /tmp/loc_db
```
Peptide sequences were then aligned against the local proteome database using parameters optimized for short queries:
```bash
blastp -task blastp-short \
  -db /tmp/loc_db \
  -query data/raw/peptides.fa \
  -evalue 1000 \
  -word_size 2 \
  -outfmt "6 qseqid sseqid pident length qlen evalue bitscore" \
  -out results/prot_pept_ext.tsv
```

To reduce false-positive matches, only alignments covering the full peptide length and having high sequence identity were retained. Specifically, matches with full query coverage (alignment length equal to query length) and at least 90% identity were selected:
```bash
awk '$3>=90 && $4==$5 {print $2}' results/prot_pept_ext.tsv | sort -u > results/chromatin_candidates_blast90_ids.txt
wc -l results/chromatin_candidates_blast90_ids.txt
```
The resulting FASTA file containing 20 candidate proteins was used as the input for subsequent localization prediction and functional annotation analyses.

## 6. Subcellular localization analysis of BLAST-derived candidates

### 6.1 TargetP 2.0 filtering

TargetP 2.0 was run in **Non-plant** mode for all 20 BLAST-derived candidate proteins to identify N-terminal targeting signals, including secretory signal peptides (SP) and mitochondrial targeting peptides (mTP).

The predictions yielded the following distribution:
- **11 proteins** classified as *Other* (non-secretory),
- **9 proteins** predicted to contain a signal peptide,
- **0 proteins** predicted to contain a mitochondrial targeting peptide.

**Table 1. TargetP 2.0 predictions for BLAST-derived candidate proteins**

| Protein ID | Prediction | Other (prob.) | Signal peptide (prob.) | mTP (prob.) |
|-----------|------------|---------------|------------------------|-------------|
| g10513.t1 | Other | 1.0000 | 0.0000 | 0.0000 |
| g10514.t1 | Other | 0.9995 | 0.0003 | 0.0001 |
| g12510.t1 | Other | 0.9997 | 0.0001 | 0.0002 |
| g14472.t1 | Other | 1.0000 | 0.0000 | 0.0000 |
| g15484.t1 | Other | 1.0000 | 0.0000 | 0.0000 |
| g3428.t1  | Other | 0.9999 | 0.0000 | 0.0001 |
| g4106.t1  | Other | 0.7297 | 0.2669 | 0.0034 |
| g5237.t1  | Other | 0.9995 | 0.0003 | 0.0001 |
| g5443.t1  | Other | 0.9529 | 0.0438 | 0.0034 |
| g5510.t1  | Other | 0.9991 | 0.0000 | 0.0009 |
| g12562.t1 | Signal peptide | 0.0001 | 0.9999 | 0.0000 |
| g1285.t1  | Signal peptide | 0.0030 | 0.9968 | 0.0002 |
| g13530.t1 | Signal peptide | 0.1160 | 0.8838 | 0.0002 |
| g15153.t1 | Signal peptide | 0.0000 | 1.0000 | 0.0000 |
| g3679.t1  | Signal peptide | 0.0018 | 0.9980 | 0.0002 |
| g5502.t1  | Signal peptide | 0.0011 | 0.9988 | 0.0000 |
| g5503.t1  | Signal peptide | 0.0012 | 0.9987 | 0.0001 |
| g5616.t1  | Signal peptide | 0.0001 | 0.9999 | 0.0000 |
| g5641.t1  | Signal peptide | 0.0001 | 0.9999 | 0.0000 |
| g702.t1   | Signal peptide | 0.0003 | 0.9997 | 0.0000 |

Proteins predicted to contain a signal peptide were excluded from further analysis. Proteins classified as *Other* were retained for downstream localization refinement using WoLF PSORT.


The remaining **11 non-secretory proteins** were retained for more detailed localization analysis.

---

### 6.2 WoLF PSORT analysis of non-secretory candidates

The non-secretory candidate proteins were further analyzed using **WoLF PSORT** in *Animal* mode to refine subcellular localization predictions.

WoLF PSORT assigns localization scores across multiple cellular compartments. Candidates were evaluated based on the highest-scoring predicted localization, with particular emphasis on nuclear compatibility.

Three proteins exhibited strong and unambiguous nuclear localization signals:
- **g14472.t1** (nucleus score = 28),
- **g10513.t1** (nucleus score = 20),
- **g10514.t1** (nucleus score = 19).

One additional protein (**g15484.t1**) showed mixed nuclear and cytosolic localization scores, suggesting possible partial nuclear association but lower specificity.

Proteins predominantly predicted to localize to the endoplasmic reticulum, Golgi apparatus, extracellular space, plasma membrane, mitochondria, or other non-nuclear compartments were excluded from further consideration.

---

## 6. Pfam domain prediction (HMMER)

To assess the functional properties of the final candidate proteins, Pfam domain prediction was performed using the HMMER web server (hmmscan against the Pfam database). Profile hidden Markov models (HMMs) allow detection of conserved protein domains even in the absence of close homologs identifiable by BLAST.

The following candidate proteins were analyzed:
- g14472.t1  
- g10513.t1  
- g10514.t1  
- g15484.t1  

### Results

For three proteins (**g14472.t1**, **g10513.t1**, **g10514.t1**), no significant Pfam domain matches were detected. This suggests that these proteins may be highly divergent, intrinsically disordered, or lineage-specific, which is consistent with a potential regulatory or structural role in chromatin association.

In contrast, **g15484.t1** showed multiple significant Pfam hits corresponding to protein families involved in vesicular trafficking and membrane-associated complexes, including components of the exocyst and GARP complexes (e.g. Sec5, EXOC6/Sec15, VPS-related domains). These domains are not directly associated with chromatin or DNA repair functions, suggesting that g15484.t1 is less likely to represent a chromatin-specific protein and may instead reflect contamination from membrane-associated proteins during chromatin fraction isolation.

![HMMER hmmscan Pfam results for chromatin-associated protein candidates](figures/HMMER.jpg)

## 7. BLAST search against UniProtKB / Swiss-Prot

To further characterize the selected chromatin-associated protein candidates, BLASTp searches were performed against the **UniProtKB/Swiss-Prot** database using the UniProt web interface. For each protein, the best BLAST hit was recorded, including accession number, E-value, percentage of sequence identity, query coverage, and functional annotation, in accordance with the project requirements.

### BLAST results summary

- **g14472.t1** showed a perfect match to the reviewed Swiss-Prot entry **P0DOW4**, corresponding to the **Damage suppressor protein (Dsup)** from *Ramazzottius varieornatus*.  
- **g10513.t1** and **g10514.t1** matched uncharacterized proteins from *R. varieornatus* with high sequence identity and coverage, but without functional annotation.  
- **g15484.t1** matched a well-characterized vacuolar protein sorting-associated protein (VPS51), a component of the GARP complex involved in vesicular trafficking, and was therefore excluded from further consideration as a chromatin-associated candidate.

---

## 8. Integration of evidence and candidate prioritization

All available evidence peptide support from chromatin fraction proteomics, predicted subcellular localization (TargetP and WoLF PSORT), Pfam domain prediction, and BLAST annotation  was integrated to prioritize proteins potentially involved in chromatin association and radiotolerance.

### Integrated summary of candidate proteins

| Protein ID | Peptide support | TargetP prediction | WoLF PSORT localization | Pfam domains | Best BLAST hit (UniProt) | Annotation | Final status |
|-----------|----------------|--------------------|-------------------------|--------------|--------------------------|------------|--------------|
| **g14472.t1** | Yes | Other | Nuclear | None | P0DOW4 | Damage suppressor protein (Dsup) | Final candidate |
| **g10513.t1** | Yes | Other | Nuclear / cytosolic | None | A0A1D1VM12 | Uncharacterized protein | Candidate |
| **g10514.t1** | Yes | Other | Nuclear / cytosolic | None | A0A1D1VK14 | Uncharacterized protein | Candidate |
| g15484.t1 | Yes | Other | Nuclear / cytosolic | Multiple (GARP / VPS) | A0A1D1W5L8 | Vacuolar protein sorting-associated protein 51 | Excluded |

### Сonclusion

The protein **g14472.t1** was identified as the strongest candidate due to its perfect match to Dsup, a unique chromatin-associated protein experimentally shown to protect DNA from radiation- and ROS-induced damage. Its detection in the chromatin fraction, nuclear localization, lack of conserved Pfam domains, and well-established role in radiotolerance strongly support its functional relevance.

Proteins **g10513.t1** and **g10514.t1** remain promising candidates as potentially novel chromatin-associated proteins. Both lack known functional domains, show nuclear or mixed nuclear/cytosolic localization, and are currently annotated as uncharacterized despite strong peptide and sequence evidence.

In contrast, **g15484.t1** was excluded from the final candidate list due to the presence of conserved domains and a clear functional annotation related to vesicular trafficking, which is inconsistent with a direct role in chromatin protection.
