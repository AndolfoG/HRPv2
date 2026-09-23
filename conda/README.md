# HRPv2 Conda installation

This directory contains the Conda recipe for HRPv2 2.0.0. The package installs
HRPv2 and its command-line dependencies in a dedicated Linux environment.
InterProScan is not distributed with the package and must be installed
separately.

For workflow usage, input requirements and output descriptions, see the
[main README](../README.md). For a fully manual installation, see the
[Installation Guidelines](../Installation_Guidelines.md).

## Contents

- [Installed components](#installed-components)
- [External requirement: InterProScan](#external-requirement-interproscan)
- [Building the package](#building-the-package)
- [Installing the local build](#installing-the-local-build)
- [Validating the environment](#validating-the-environment)
- [Running HRPv2](#running-hrpv2)
- [Small datasets and NB-ARC rescue](#small-datasets-and-nb-arc-rescue)
- [Removing the environment](#removing-the-environment)
- [Recipe maintenance](#recipe-maintenance)

## Installed components

The Conda package installs the following components in the same environment:

| Component | Version constraint |
|---|---|
| HRPv2 | 2.0.0 |
| Python | `>=3.12,<3.13` |
| OpenJDK | `>=11,<12` |
| GenBlastG | `1.38` |
| BLAST legacy | `2.2.26` |
| HMMER | `>=3.4,<3.5` |
| MEME Suite | `>=5.5.9,<5.6` |
| MAFFT | `>=7.526,<7.527` |

It also installs:

```text
$CONDA_PREFIX/bin/HRPv2
$CONDA_PREFIX/bin/run_genblastG
$CONDA_PREFIX/share/hrpv2/alignscore.txt
```

No separate rescue environment is required. HMMER, MEME/MAST and MAFFT are
executed directly from the active HRPv2 environment.

## External requirement: InterProScan

Install InterProScan separately from its official standalone Linux
distribution. HRPv2 2.0.0 was tested with InterProScan 5.78-109.0 and OpenJDK
11. The required member-database applications are:

```text
Pfam,SUPERFAMILY,Coils,Gene3D,SMART,PANTHER,CDD,FunFam,
PRINTS,ProSiteProfiles,Phobius
```

Phobius is optional. If InterProScan explicitly reports that Phobius is
deactivated, HRPv2 repeats the analysis without it.

Record and validate the absolute path to the installation:

```bash
INTERPROSCAN_BIN=/absolute/path/to/interproscan.sh
test -x "$INTERPROSCAN_BIN"
"$INTERPROSCAN_BIN" --version
```

Supply this path through `--interproscan-bin` whenever existing InterProScan
TSV results are not reused.

## Building the package

### Prerequisites

Install Conda and `conda-build` before building the recipe:

```bash
conda install -n base -c conda-forge conda-build
```

### Release source

The recipe downloads the following GitHub release asset:

```text
Tag:      v2.0.0
Asset:    HRPv2-2.0.0.tar.gz
SHA-256:  eeca6adad28f95f69bad801d5a19cd5e7814f47cab6fb4761d3e998ad7b30856
```

The release must be publicly accessible before the remote-source build can
succeed. The tag, filename, archive contents and checksum must remain
unchanged unless the `source` section of `meta.yaml` is updated accordingly.

Clone the repository and build from its root directory:

```bash
git clone https://github.com/AndolfoG/HRPv2.git
cd HRPv2
conda build conda
```

HRPv2 is packaged as `noarch: generic` because its installed scripts and data
files are platform-independent. Practical execution is currently supported on
Linux because GenBlastG and BLAST legacy are Linux runtime dependencies. During
the build, `run_test.sh` verifies the Python and dependency versions, the HRPv2
command-line interface, required executables and the installed `alignscore.txt`
support file.

## Installing the local build

After a successful build, create a dedicated environment from the local Conda
channel:

```bash
conda create -n hrpv2 \
  -c local \
  -c conda-forge \
  -c bioconda \
  hrpv2=2.0.0
```

Activate it:

```bash
conda activate hrpv2
```

The `-c local` channel is required only for a locally built package. Once
HRPv2 is published in a public Conda channel, use the installation command
provided by that channel instead.

## Validating the environment

Run the following checks after installation:

```bash
HRPv2 --version
HRPv2 --help
python --version
java -version

command -v HRPv2
command -v run_genblastG
command -v genblastG
command -v blastall
command -v formatdb
command -v hmmbuild
command -v hmmsearch
command -v meme
command -v meme2meme
command -v mast
command -v mafft

test -r "$CONDA_PREFIX/share/hrpv2/alignscore.txt"
```

## Running HRPv2

Run an analysis with the command installed by the package:

```bash
HRPv2 \
  --proteome species_proteins.fasta \
  --genome species_genome.fasta \
  --annotation-gff species_annotation.gff3 \
  --interproscan-bin "$INTERPROSCAN_BIN" \
  --workdir HRPv2_results \
  --threads 8
```

The GFF3 annotation is optional, but it should be supplied when it corresponds
to the input proteome and genome. Without it, annotation-aware recovery and
merging are skipped.

## Small datasets and NB-ARC rescue

The species-specific NB-ARC rescue requires at least 10 valid,
non-redundant, complete PF00931 training regions. This threshold concerns the
number of eligible PF00931 regions, not the total number of proteins in the
input FASTA.

For a small dataset that does not satisfy this requirement, run:

```bash
HRPv2 \
  --proteome test_proteins.fasta \
  --genome test_genome.fasta \
  --annotation-gff test_annotation.gff3 \
  --interproscan-bin "$INTERPROSCAN_BIN" \
  --skip-nb-rescue \
  --workdir HRPv2_test \
  --threads 4
```

Do not use `--skip-nb-rescue` when the purpose of the test is to validate the
HMMER plus MEME/MAST rescue itself; use a representative dataset containing at
least 10 eligible PF00931 training regions instead.

## Removing the environment

Deactivate and remove the test environment without affecting other Conda
environments:

```bash
conda deactivate
conda env remove -n hrpv2
```

The locally built package cache is separate from the environment. Remove old
Conda build artifacts only when they are no longer needed:

```bash
conda build purge
```

## Recipe maintenance

Keep these files together in the repository's `conda/` directory:

```text
conda/
├── README.md
├── build.sh
├── meta.yaml
└── run_test.sh
```

For every new HRPv2 release:

1. update the version in `meta.yaml`;
2. create the matching GitHub tag and source archive;
3. calculate the archive SHA-256 checksum;
4. update the checksum in `meta.yaml` and this document;
5. increment the Conda build number when the recipe changes without changing
   the HRPv2 version;
6. rebuild and repeat the package and end-to-end tests.
