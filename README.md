![HRP logo](https://github.com/AndolfoG/HRP-2.0/blob/main/LOGO.png)

# Introduction
Welcome to to the full-length **H**omology-based **R**-gene **P**rediction version 2.0 (**HRPv2**) project.
**HRPv2** is a command-line workflow for the genome-wide prediction, classification and filtering of plant NB-LRR resistance genes.

## Table of contents

- [Essential software](#essential-software)
- [Installation checks](#installation-checks)
- [Input files](#input-files)
- [HRPv2 workflow](#hrpv222-workflow)
- [Domain classification](#domain-classification)
- [Redundancy and gene-fusion filtering](#redundancy-and-gene-fusion-filtering)
- [Running HRPv2](#running-hrpv222)
- [Reusing existing outputs](#reusing-existing-outputs)
- [Output files](#output-files)
- [InterProScan configuration](#interproscan-configuration)
- [GenBlastG compatibility wrapper](#genblastg-compatibility-wrapper)
- [Notes and limitations](#notes-and-limitations)
- [Authorship](#authorship)

<a id="essential-software"></a>
## Essential software

Make sure the following programs are correctly installed and available in the
active environment.

### Python

- Python 3.10 or newer.
- HRPv2 uses only the Python standard library.
- No additional `pip` packages are required.

### InterProScan

- InterProScan 5.78-109.0, or a compatible InterProScan 5 release.
- The Pfam, SUPERFAMILY and Coils applications and their data must be installed.
- The executable `interproscan.sh` must be available through `PATH`.

HRPv2 runs InterProScan with:

```text
Pfam,SUPERFAMILY,Coils
```

### Java

- OpenJDK 11 is required by InterProScan 5.78-109.0.

### Species-specific NB-LRR proteins rescue

- HMMER, MEME Suite and MAFFT are required when the species-specific NB-ARC
  rescue is enabled.
- They may be installed in a separate Conda environment selected with
  `--rescue-env`.
- Use `--skip-nb-rescue` only to disable the rescue or when reusing an already
  integrated `with_NB-ARC_rescue.tsv` file.

### GenBlastG and legacy BLAST

The following files must be available:

- `genblastG`
- `formatdb`
- `blastall`
- `alignscore.txt`
- `run_genblastG`

`formatdb` and `blastall` are legacy BLAST programs required by the GenBlastG
release used in this workflow.

### Perl

Perl should be installed in the working Conda environment for compatibility
with the external bioinformatics software used by the workflow.

<a id="installation-checks"></a>
## Installation checks

Activate the environment containing GenBlastG and InterProScan:

```bash
conda activate genblastg_env
```

Check the required software:

```bash
python3 --version
java -version
interproscan.sh --version
which run_genblastG
which genblastG
which formatdb
which blastall
ls "$CONDA_PREFIX/bin/alignscore.txt"
```

For the configuration used during development, the expected InterProScan output
includes:

```text
InterProScan version 5.78-109.0
InterProScan 64-Bit build (requires Java 11)
```

The required InterProScan data directories can be checked with:

```bash
ls /path/to/interproscan-5.78-109.0/data/{pfam,superfamily}
```

Coils must also be available among the installed InterProScan applications.

<a id="input-files"></a>
## Input files

HRPv2 accepts two mandatory inputs and one optional input. Previously
computed InterProScan, HMMER/MAST and GenBlastG results can also be supplied to
reclassify an existing run without repeating the expensive external searches.

### Mandatory inputs

1. **Protein sequences in FASTA format**

   This may be the proteome of the analysed species or a selected protein-query
   dataset.

2. **Genome sequence in FASTA format**

   This is the target genome used by GenBlastG for gene-model prediction.

### Optional input

3. **Genome annotation in GFF3 format**

   This file hould be supplied only when it corresponds to the input proteome and target genome.
   If the proteome is an external or specialised query dataset that does not match the target annotation, omit the GFF3 file. The merging operation will be skipped.

### Data preparation
Download the genome sequence, the protein sequences (encoded by gene set) and genome annotation of interest, for example related to the tomato genome assembly SL2.4
			
	wget -O S_lycopersicum_2.4_genome.fasta https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/Heinz1706/assembly/build_2.40/S_lycopersicum_chromosomes.2.40.fa.gz 
	wget -O ITAG2.3_proteins.fasta https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/Heinz1706/annotation/ITAG2.3_release/ITAG2.3_proteins.fasta 
	wget -O ITAG2.3_gene_models.gff3 https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/Heinz1706/annotation/ITAG2.3_release/ITAG2.3_gene_models.gff3


<a id="hrpv222-workflow"></a>
## HRPv2 workflow

The terminal reports six sections at the start of their corresponding steps:

```text
[1/6] Preflight validation and preparation of input files
  
[2/6] Full-length NB-LRR annotation using a protein motif/domain-based search (PDS)

[3/6] Prediction of NB-LRR gene models

[4/6] Annotation and classification of gene models

[5/6] Filtering of redundant gene models

[6/6] Full NB-LRR resistance gene repertoire
```

The main operations are:

1. Checking input files.
2. Annotate the proteins using a protein motif/domain-based search.
3. Select full-length CNL, TNL, RNL and NL proteins.
4. Use the selected proteins as GenBlastG queries against the target genome.
5. Annotate and classify every protein predicted by GenBlastG.
6. Exclude non-NB-LRR models, proteins shorter than 50 amino acids and gene models
   longer than 20 kb.
7. Resolve overlapping predictions only when they occur on the same chromosome
   and strand.
8. Optionally recover unambiguous, supported annotated partial NB-LRR genes from
   regions not occupied by retained GenBlastG loci.

Opposite-strand predictions are retained independently. The filtering decision
uses NB-LRR class, domain-architecture completeness, LRR evidence, GenBlastG
alignment score, GenBlastG rank, protein length and a deterministic identifier
tie-break.

<a id="domain-classification"></a>
## Domain classification

HRPv2 classifies complete NB-LRR proteins as CNL, TNL, RNL or NL and
retains supported partial architectures.

Functional-description filtering is applied globally to `partial NB-LRR` genes.
Descriptions typical of membrane receptors (RLP: receptor-like proteins and RLK: receptor-like kinase) act as a negative veto.
This prevents generic plasma membrane receptors from being retained as partial NB-LRRs merely because another database reports
an LRR signature. The triggering evidence is recorded in the
`exclusion_evidence` column of the classification TSV.

<a id="redundancy-and-gene-fusion-filtering"></a>
## Redundancy and gene-fusion filtering

Ordinary overlapping models are ranked deterministically, with a complete
CNL/TNL/RNL/NL model preferred to a less complete overlapping model when the
other filtering evidence is compatible.

HRPv2 also prevents a long GenBlastG prediction from incorrectly merging
two adjacent NB-LRR loci:

- `REPEATED_CORE` requires an ordered `NB-LRR` architecture.
- `MULTI_NB_BRIDGE` identifies a multi-NB bridge that does not satisfy the
  strict repeated-core pattern.
- A bridge can be replaced by two models only when both models are on the same
  chromosome and strand as the bridge, are genomically distinct and strictly
  non-overlapping, and each contains a complete NB-LRR core classified as CNL,
  TNL, RNL or NL.
- A complete-core plus partial-model pair cannot justify splitting a bridge.

Every overlap decision and its reason is written to the step-5 overlap report.

<a id="running-hrpv222"></a>
## Running HRPv2

Make the script executable:

```bash
chmod +x HRPv2.py
```

### Interactive execution

```bash
./HRPv2.py
```

or:

```bash
python3 HRPv2.py
```

The program requests:

```text
Proteome FASTA file:
Genome FASTA file:
Genome annotation GFF3 file (optional; press Enter to skip):
```

### Command-line execution with GFF3 merging

```bash
python3 HRPv2.py \
  --proteome species_proteins.fasta \
  --genome species_genome.fasta \
  --annotation-gff species_annotation.gff3 \
  --threads 8
```

Display all options with:

```bash
python3 HRPv2.py --help
```

Display the version with:

```bash
python3 HRPv2.py --version
```

Use `--workdir` to place all six numbered output directories under a separate
root directory. Without this option, they are created in the current directory.

<a id="reusing-existing-outputs"></a>
## Reusing existing outputs

Classification, filtering and merging can be repeated without rerunning
InterProScan, GenBlastG, HMMER or MEME/MAST. Supply the integrated
`with_NB-ARC_rescue.tsv` files generated in steps 2 and 4, together with the
paired GenBlastG GFF and protein files:

```bash
python3 HRPv2.py \
  --proteome species_proteins.fasta \
  --genome species_genome.fasta \
  --annotation-gff species_annotation.gff3 \
  --workdir HRPv2_reclassification \
  --threads 8 \
  --interpro-tsv 02_full-length_NB-LRR_annotation/proteome_NB-LRR.with_NB-ARC_rescue.tsv \
  --genblastg-gff 03_gene_model_prediction/genblastg_output.gff \
  --genblastg-proteins 03_gene_model_prediction/genblastg_output.pro \
  --predicted-interpro-tsv 04_NB-LRR_model_annotation/predicted_gene_models.with_NB-ARC_rescue.tsv \
  --skip-nb-rescue
```

<a id="output-files"></a>
## Output files

The final catalogue is written to:

```text
06_merged_NB-LRR_genes/
```

Its three principal files share the same prefix:

```text
final_NB-LRR-prediction.tsv
final_NB-LRR-prediction.fasta
final_NB-LRR-prediction.gff3
```

- `final_NB-LRR-prediction.tsv` contains the classification, gene identifier,
  protein identifier, source, domain architecture and domain evidence.
- `final_NB-LRR-prediction.fasta` contains the proteins corresponding to the
  final catalogue.
- `final_NB-LRR-prediction.gff3` contains the filtered GenBlastG predictions and,
  when available, the recovered partial genes.

NB-LRR genes retain the original gene `ID` reported in the input GFF3. The same
identifier is recorded in the `gene_id` column of the final TSV.

All workflow results are organised in six numbered directories:

```text
01_preflight_validation_inputs/

02_full-length_NB-LRR_annotation/

03_gene_model_prediction/

04_NB-LRR_model_annotation/

05_filtered_NB-LRR_models/

06_merged_NB-LRR_genes/
```

Diagnostic reports document models excluded as non-NB-LRR, threshold failures,
overlap decisions, unmatched annotation identifiers and ambiguous matches.

<a id="interproscan-configuration"></a>
## InterProScan configuration

HRPv2 executes InterProScan twice:

1. On the input protein dataset to select full-length NB-LRR queries.
2. On all proteins predicted by GenBlastG before filtering.

The equivalent command structure is:

```bash
interproscan.sh \
  -i proteins.fasta \
  -appl Pfam,SUPERFAMILY,Coils \
  -f TSV,GFF3 \
  -b output_prefix \
  -cpu 4 \
  -dp
```

The `-dp` option disables the pre-calculated match lookup service and forces
local match calculation. It is appropriate for a fully local installation.

<a id="genblastg-compatibility-wrapper"></a>
## GenBlastG compatibility wrapper

The `run_genblastG` wrapper used here is a local compatibility script, not an
official GenBlastG executable. It creates temporary symbolic links in the
working directory so that GenBlastG can locate:

```text
formatdb
blastall
alignscore.txt
```

The wrapper removes only the links that it created after GenBlastG terminates.
The executable and supporting files are expected in:

```text
$CONDA_PREFIX/bin/
```

<a id="notes-and-limitations"></a>
## Notes
- HRPv2 is designed for Linux 64-bit environments compatible with the selected
  InterProScan and GenBlastG releases.

<a id="authorship"></a>
## Authorship

**Andolfo Giuseppe, Ph.D.**  
University of Naples "Federico II", Italy
Plant Genetics and Biotechnology Unit   

Please acknowledge the authors and cite the HRP publication when using this
workflow in scientific research: Andolfo *et al.* (2022), *The Plant Journal*.
