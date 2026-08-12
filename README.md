![HRP logo](https://github.com/AndolfoG/HRP/blob/main/LOGO.png)
HRP-2.0
Welcome to The full-length Homology-based R-gene Prediction: HRP-2.0.
HRP-2.0 is a command-line workflow for the genome-wide prediction,
classification and filtering of plant NB-LRR resistance genes. It identifies
full-length NB-LRR proteins, uses them as GenBlastG queries, annotates every
predicted model and removes false-positive or redundant predictions. When a
compatible genome annotation is provided, annotated partial NB-LRR genes missed
by GenBlastG are recovered and merged into the final catalogue.
Table of contents
Essential software
Installation checks
Input files
HRP-2.0 workflow
Running HRP-2.0
Output files
InterProScan configuration
GenBlastG compatibility wrapper
Notes and limitations
Authorship
Essential software
Make sure the following programs are correctly installed and available in the
active environment.
Python
Python 3.10 or newer.
HRP-2.0 uses only the Python standard library.
No additional `pip` packages are required.
InterProScan
InterProScan 5.78-109.0, or a compatible InterProScan 5 release.
The Pfam, SUPERFAMILY and Coils applications and their data must be installed.
The executable `interproscan.sh` must be available through `PATH`.
HRP-2.0 runs InterProScan with:
```text
Pfam,SUPERFAMILY,Coils
```
Java
OpenJDK 11 is required by InterProScan 5.78-109.0.
GenBlastG and legacy BLAST
The following files must be available:
`genblastG`
`formatdb`
`blastall`
`alignscore.txt`
`run_genblastG`
`formatdb` and `blastall` are legacy BLAST programs required by the GenBlastG
release used in this workflow.
Perl
Perl should be installed in the working Conda environment for compatibility
with the external bioinformatics software used by the workflow.
Installation checks
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
Input files
HRP-2.0 accepts two mandatory inputs and one optional input.
Mandatory inputs
Protein sequences in FASTA format
This may be the proteome of the analysed species or a selected protein-query
dataset. Terminal stop symbols (`*`) are removed from a cleaned working copy.
Genome sequence in FASTA format
This is the target genome used by GenBlastG for gene-model prediction.
Optional input
Genome annotation in GFF3 format
This file is used only to recover annotated partial NB-LRR genes that do not
overlap filtered GenBlastG loci. It should be supplied only when it corresponds
to the input proteome and target genome.
`ID=` is the primary identifier used to associate proteins with GFF3 features.
`protein_id=` and `transcript_id=` are accepted as exact-match fallbacks, and
`Parent=` is used when available to reconstruct feature relationships.
If the proteome is an external or specialised query dataset that does not match
the target annotation, omit the GFF3 file. The merging operation will be skipped.
HRP-2.0 workflow
The terminal reports six sections at the start of their corresponding steps:
```text
[1/6] Prediction of NB-LRR resistance genes based on full-length sequence homology
[2/6] Full-length NB-LRR annotation using a protein motif/domain-based search (PDS)
[3/6] Prediction of NB-LRR gene models
[4/6] Annotation and classification of NB-LRR gene models
[5/6] Filtering of redundant NB-LRR gene models
[6/6] Merging of partial and full-length NB-LRR gene-models
```
The main operations are:
Clean terminal stop symbols from the protein input.
Annotate the proteins using Pfam, SUPERFAMILY and Coils.
Select full-length CNL, TNL, RNL and NL proteins containing NB-ARC and LRR.
Use the selected proteins as GenBlastG queries against the target genome.
Annotate and classify every protein predicted by GenBlastG.
Exclude non-NLR models, proteins shorter than 50 amino acids and gene models
longer than 20 kb.
Resolve overlapping predictions only when they occur on the same chromosome
and strand.
Optionally recover unambiguous annotated partial genes containing NB-ARC from
regions not occupied by retained GenBlastG loci.
Opposite-strand predictions are retained independently. The filtering decision
uses NLR class, domain-architecture completeness, LRR evidence, GenBlastG
alignment score, GenBlastG rank, protein length and a deterministic identifier
tie-break.
Running HRP-2.0
Make the script executable:
```bash
chmod +x HRP-2.0.py
```
Interactive execution
```bash
./HRP-2.0.py
```
or:
```bash
python3 HRP-2.0.py
```
The program requests:
```text
Proteome FASTA file:
Genome FASTA file:
Genome annotation GFF3 file (optional; press Enter to skip):
```
Command-line execution with GFF3 merging
```bash
python3 HRP-2.0.py \
  --proteome species_proteins.fasta \
  --genome species_genome.fasta \
  --annotation-gff species_annotation.gff3 \
  --threads 8
```
Command-line execution without GFF3 merging
```bash
python3 HRP-2.0.py \
  --proteome query_proteins.fasta \
  --genome target_genome.fasta \
  --threads 8
```
Display all options with:
```bash
python3 HRP-2.0.py --help
```
Display the version with:
```bash
python3 HRP-2.0.py --version
```
Output files
The final catalogue is written to:
```text
merged_NB-LRR_models/
```
Its three principal files share the same prefix:
```text
Final_NB-LRR-prediction.tsv
Final_NB-LRR-prediction.fasta
Final_NB-LRR-prediction.gff3
```
`Final_NB-LRR-prediction.tsv` contains the classification, gene identifier,
protein identifier, source, domain architecture and domain evidence.
`Final_NB-LRR-prediction.fasta` contains the proteins corresponding to the
final catalogue.
`Final_NB-LRR-prediction.gff3` contains the filtered GenBlastG predictions and,
when available, the recovered partial genes.
Partial genes retain the original gene `ID` reported in the input GFF3. The same
identifier is recorded in the `gene_id` column of the final TSV.
When no annotation GFF3 is supplied, HRP-2.0 still generates the same three final
files, but they contain only the filtered GenBlastG models.
Intermediate prediction results are stored in:
```text
predicted_NB-LRRs/
filtered_NB-LRR_models/
```
Diagnostic reports document models excluded as non-NLR, threshold failures,
overlap decisions, unmatched annotation identifiers and ambiguous matches.
InterProScan configuration
HRP-2.0 executes InterProScan twice:
On the input protein dataset to select full-length NB-LRR queries.
On all proteins predicted by GenBlastG before filtering.
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
GenBlastG compatibility wrapper
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
Notes and limitations
The protein FASTA, genome FASTA and optional GFF3 must use compatible sequence
and feature identifiers when merging is requested.
Identifier matching is exact. Ambiguous or missing associations are reported
and are never merged automatically.
Partial models are recovered only when NB-ARC evidence is present.
Partial models containing only CC, TIR, RPW8 or LRR evidence are not recovered.
GenBlastG redundancy is assessed only for models on the same chromosome and
strand.
HRP-2.0 is designed for Linux 64-bit environments compatible with the selected
InterProScan and GenBlastG releases.
Authorship
Developed by:
Prof. Andolfo Giuseppe  
University of Naples "Federico II"  
Plant Genetics and Biotechnology Unit  
Piazza Carlo di Borbone, 1 - 80055 Portici (Naples), Italy
Please acknowledge the authors and cite the HRP publication when using this
workflow in scientific research:
Andolfo et al. (2022), The Plant Journal.
