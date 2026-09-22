# HRPv2 installation guide for Ubuntu/Linux

This document provides a reproducible, step-by-step installation procedure for
running the complete HRPv2 pipeline on a 64-bit Ubuntu/Linux system.
It is intended as a companion to `README.md`.

The commands were checked against the executable requirements of
`HRPv2.py`.

<a id="contents"></a>
## Contents

- [1. Software used by HRPv2](#software-used)
- [2. Hardware and operating-system recommendations](#requirements)
- [3. Install Ubuntu system packages](#ubuntu-packages)
- [4. Install Miniforge/Conda](#conda)
- [5. Create the HRPv2 Conda environment](#environment)
- [6. Configure GenBlastG](#genblastg)
- [7. Install InterProScan](#interproscan)
- [8. Install HRPv2](#hrpv2)
- [9. Validate the complete installation](#validation)
- [10. Run a smoke test](#smoke-test)
- [11. Running HRPv2](#running-hrpv2)
- [12. Troubleshooting](#troubleshooting)
- [13. Reproducibility and environment export](#reproducibility)
- [14. Official resources](#resources)
- [15. Conda installation](#conda-installation)

<a id="software-used"></a>
## 1. Software used by HRPv2

| Component | Role in HRPv2 | Installation used here |
|---|---|---|
| Python ≥3.10 | Pipeline, classification, filtering and merge | Conda |
| Perl | Compatibility with MEME Suite and legacy tools | Conda |
| GenBlastG 1.38 | Homology-based gene-model prediction | Bioconda |
| NCBI BLAST legacy | Provides `blastall` and `formatdb` for GenBlastG | Bioconda |
| HMMER 3 | Species-specific NB-ARC HMM construction and search | Bioconda |
| MEME Suite | NB-ARC motif discovery; supplies `meme`, `meme2meme`, `mast` | Bioconda |
| MAFFT | Alignment of NB-ARC seed sequences | Bioconda |
| InterProScan 5 | Protein motif/domain annotation | Standalone EBI distribution |
| OpenJDK 11 | Java runtime for the pinned InterProScan release | Conda |

HRPv2 itself uses only the Python standard library. No `pip` packages are
required.

The current script invokes the following InterProScan applications:

```text
Pfam,SUPERFAMILY,Coils,Gene3D,SMART,PANTHER,CDD,FunFam,
PRINTS,ProSiteProfiles,Phobius
```

`Phobius` is optional because its local components require a separate licence.
If InterProScan explicitly reports that Phobius is deactivated, HRPv2 retries
the same analysis without Phobius. The remaining applications must be
available; other InterProScan errors are not ignored.

<a id="requirements"></a>
## 2. Hardware and operating-system recommendations

- Ubuntu 20.04, 22.04 or 24.04, 64-bit (`x86_64`).
- At least 16 GB RAM; 32–64 GB is preferable for a plant proteome.
- At least 100 GB free disk space for InterProScan, temporary files and HRPv2
  results.
- Multiple CPU cores; the examples below use 24 threads.
- Internet access during installation.

The RAM, disk and thread values are operational recommendations rather than
hard-coded HRPv2 limits. Requirements increase with proteome and genome size.

Check the architecture:

```bash
uname -m
```

The expected result for the installation below is:

```text
x86_64
```

<a id="ubuntu-packages"></a>
## 3. Install Ubuntu system packages

```bash
sudo apt update
sudo apt install -y \
  bash build-essential ca-certificates coreutils curl git gzip \
  libx11-6 libxext6 libxrender1 perl tar wget
```

Create a dedicated installation directory in the current location:

```bash
mkdir -p "$PWD/hrpv2-install"
HRPV2_ROOT="$PWD/hrpv2-install"
```

Keep this terminal open while installing, or redefine `HRPV2_ROOT` in every
new terminal. You may replace this path with another absolute path.

<a id="conda"></a>
## 4. Install Miniforge/Conda

Skip this section if a working `conda` installation is already available.

```bash
cd "$HRPV2_ROOT"
wget -O Miniforge3-Linux-x86_64.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

Accept the licence, choose the proposed installation path or another user-owned
path, and allow the installer to initialise the shell. Then close and reopen
the terminal, or load Conda manually using the path selected during setup:

```bash
source /path/to/miniforge3/etc/profile.d/conda.sh
```

Configure channels in the order recommended for Bioconda:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

With Miniforge, `conda-forge` may already be configured. Verify the result:

```bash
conda config --show channels
```

<a id="environment"></a>
## 5. Create the HRPv2 Conda environment

Create one environment containing the pipeline language runtimes, GenBlastG
and the complete NB-ARC rescue toolchain:

```bash
conda create -n hrpv2_env -y \
  python=3.11 \
  perl \
  openjdk=11 \
  genblastg=1.38 \
  blast-legacy \
  hmmer=3.4 \
  meme \
  mafft
```

Activate it:

```bash
conda activate hrpv2_env
```

The HRPv2 option used later must therefore be:

```text
--rescue-env hrpv2_env
```

If the solver cannot produce the environment in one operation, create it in
two stages:

```bash
conda create -n hrpv2_env -y python=3.11 perl openjdk=11
conda activate hrpv2_env
conda install -y -c conda-forge -c bioconda \
  genblastg=1.38 blast-legacy hmmer=3.4 meme mafft
```

Do not replace `blast-legacy` with modern BLAST+: GenBlastG 1.38 expects the
legacy executable names `blastall` and `formatdb`.

<a id="genblastg"></a>
## 6. Configure GenBlastG

### 6.1 Verify installed files

```bash
conda activate hrpv2_env
command -v genblastG
command -v blastall
command -v formatdb
find "$CONDA_PREFIX" -name alignscore.txt -print
```

All four components must exist. Copy `alignscore.txt` into the environment's
`bin` directory only if the package installed it elsewhere:

```bash
ALIGNSCORE_PATH=$(find "$CONDA_PREFIX" -name alignscore.txt -print -quit)
test -n "$ALIGNSCORE_PATH"
if [ "$ALIGNSCORE_PATH" != "$CONDA_PREFIX/bin/alignscore.txt" ]; then
  cp "$ALIGNSCORE_PATH" "$CONDA_PREFIX/bin/alignscore.txt"
fi
```

### 6.2 Install the compatibility wrapper expected by HRPv2

HRPv2 calls `run_genblastG` by default. The wrapper below exposes the legacy
BLAST programs and `alignscore.txt` in GenBlastG's working directory, then
removes only the symbolic links it created.

Create the file:

```bash
nano "$CONDA_PREFIX/bin/run_genblastG"
```

Paste the following content:

```bash
#!/usr/bin/env bash
set -euo pipefail

created_links=()

cleanup() {
  local item
  for item in "${created_links[@]}"; do
    if [ -L "$item" ]; then
      unlink "$item"
    fi
  done
}
trap cleanup EXIT

for item in blastall formatdb alignscore.txt; do
  if [ -e "$item" ] || [ -L "$item" ]; then
    continue
  fi

  if [ "$item" = "alignscore.txt" ]; then
    source_path="${CONDA_PREFIX}/bin/alignscore.txt"
  else
    source_path=$(command -v "$item")
  fi

  if [ ! -e "$source_path" ]; then
    echo "Required GenBlastG support file not found: $source_path" >&2
    exit 1
  fi

  ln -s "$source_path" "$item"
  created_links+=("$item")
done

genblastG "$@"
```

Make it executable:

```bash
chmod 755 "$CONDA_PREFIX/bin/run_genblastG"
```

The wrapper must be on `PATH` whenever the `hrpv2_env` environment is active:

```bash
command -v run_genblastG
```

<a id="interproscan"></a>
## 7. Install InterProScan

InterProScan should be installed from the self-contained EBI Linux archive,
not inside the Conda package set. The example pins InterProScan 5.78-109.0,
the release used during development of this HRPv2 version.

### 7.1 Download and verify the archive

```bash
cd "$HRPV2_ROOT"
mkdir -p interproscan
cd interproscan

wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz.md5

md5sum -c interproscan-5.78-109.0-64-bit.tar.gz.md5
```

Do not continue unless the checksum reports `OK`.

### 7.2 Extract InterProScan

```bash
tar -pxzf interproscan-5.78-109.0-64-bit.tar.gz
cd interproscan-5.78-109.0
```

The standalone archive already contains the supported member-database binaries
and signature data, including PANTHER data. Licensed applications may remain
deactivated.

### 7.3 Initialise and test InterProScan

Ensure that the Conda environment—and therefore Java 11—is active:

```bash
conda activate hrpv2_env
java -version
python3 setup.py -f interproscan.properties
./interproscan.sh -i test_all_appl.fasta -f tsv -dp
```

`setup.py` indexes the HMM models and normally needs to run only once. The EBI
archive may already contain indexed models, but running the official setup step
is a useful validation.

### 7.4 Expose InterProScan in the HRPv2 environment

Create a symbolic link in the active Conda environment:

```bash
ln -s "$PWD/interproscan.sh" "$CONDA_PREFIX/bin/interproscan.sh"
```

If the destination already exists, inspect it rather than overwriting it:

```bash
readlink -f "$CONDA_PREFIX/bin/interproscan.sh"
```

Verify the version and available applications:

```bash
interproscan.sh --version
interproscan.sh 2>&1 | less
```

Confirm that the following names occur under `Available analyses`:

```text
Pfam SUPERFAMILY Coils Gene3D SMART PANTHER CDD PRINTS ProSiteProfiles
```

`FunFam` is produced through the Gene3D/CATH analysis in compatible releases.
Phobius may appear under `Deactivated analyses`; HRPv2 handles this specific
case automatically.

### 7.5 Important local-execution behaviour

HRPv2 passes `-dp` to InterProScan. This disables the EBI pre-calculated match
lookup and performs the searches locally. Consequently, the installation must
contain the relevant member-database data and local binaries.

<a id="hrpv2"></a>
## 8. Install HRPv2

Copy the release into a project directory:

```bash
mkdir -p "$HRPV2_ROOT/pipeline"
cd "$HRPV2_ROOT/pipeline"
cp /path/to/HRPv2.py .
chmod 755 HRPv2.py
```

HRPv2 is a self-contained Python script. Do not split or extract sections from
the file.

Check the command-line interface:

```bash
python3 HRPv2.py --version
python3 HRPv2.py --help
```

<a id="validation"></a>
## 9. Validate the complete installation

Activate the environment and run every check before starting a long analysis:

```bash
conda activate hrpv2_env

python3 --version
perl -v
java -version

interproscan.sh --version
command -v run_genblastG
command -v genblastG
command -v blastall
command -v formatdb
test -s "$CONDA_PREFIX/bin/alignscore.txt"

hmmbuild -h | head
hmmsearch -h | head
meme -version
meme2meme -version
mast -version
mafft --version
```

Validate the rescue environment exactly as HRPv2 will access it:

```bash
conda run -n hrpv2_env sh -c '
  command -v hmmsearch &&
  command -v hmmbuild &&
  command -v meme &&
  command -v meme2meme &&
  command -v mast &&
  command -v mafft
'
```

Check the HRPv2 command without executing external analyses:

```bash
python3 HRPv2.py \
  --proteome /path/to/proteins.fasta \
  --genome /path/to/genome.fasta \
  --annotation-gff /path/to/annotation.gff3 \
  --threads 24 \
  --rescue-env hrpv2_env \
  --dry-run
```

The dry run checks argument handling and displays planned operations. It does
not prove that a biological analysis will complete, so also perform the tests
below.

<a id="smoke-test"></a>
## 10. Run a smoke test

### 10.1 InterProScan test

```bash
mkdir -p "$HRPV2_ROOT/tests/interproscan"
cd "$HRPV2_ROOT/tests/interproscan"
interproscan.sh \
  -i /path/to/interproscan/test_all_appl.fasta \
  -appl Pfam,SUPERFAMILY,Coils,Gene3D,SMART,PANTHER,CDD,FunFam,PRINTS,ProSiteProfiles \
  -f TSV,GFF3 \
  -b ips_test \
  -cpu 4 \
  -dp
```

If the installed InterProScan release does not accept `FunFam` as an explicit
application name, confirm that Gene3D is active and use HRPv2's actual command
to diagnose compatibility. Do not silently remove required applications from
the pipeline script.

### 10.2 GenBlastG test

Use a small protein query and its matching genome:

```bash
mkdir -p "$HRPV2_ROOT/tests/genblastg"
cd "$HRPV2_ROOT/tests/genblastg"

run_genblastG \
  -q query_protein.fasta \
  -t small_genome.fasta \
  -gff -pro \
  -o genblastg_test
```

Confirm that `.gff` and `.pro` outputs were created:

```bash
find . -maxdepth 1 -type f \
  \( -name 'genblastg_test*.gff' -o -name 'genblastg_test*.pro' \) \
  -print
```

<a id="running-hrpv2"></a>
## 11. Running HRPv2

For complete execution commands, output-reuse procedures and output
descriptions, see the corresponding sections of `README.md`.

<a id="troubleshooting"></a>
## 12. Troubleshooting

### `required executable not found in PATH`

Activate the environment and inspect `PATH`:

```bash
conda activate hrpv2_env
command -v python3 interproscan.sh run_genblastG
```

### `rescue environment ... is missing or incomplete`

At least one of these programs is unavailable inside the environment:

```text
hmmsearch hmmbuild meme meme2meme mast mafft
```

Run the exact validation command in section 9 and reinstall the missing
package.

### Phobius is deactivated

This is expected for many standalone InterProScan installations because
Phobius is separately licensed. HRPv2 detects this specific diagnostic and
repeats InterProScan without Phobius. No manual change is required.

### Another InterProScan application is missing

Do not treat this as a Phobius warning. Confirm that the complete standalone
InterProScan archive was extracted and that HRPv2 is resolving the intended
`interproscan.sh`:

```bash
readlink -f "$(command -v interproscan.sh)"
```

### GenBlastG cannot find `blastall`, `formatdb` or `alignscore.txt`

Confirm the files and wrapper:

```bash
command -v blastall formatdb genblastG run_genblastG
ls -l "$CONDA_PREFIX/bin/alignscore.txt"
```

Run GenBlastG through `run_genblastG`, not directly through `genblastG`.

### GenBlastG fails with absolute input paths

HRPv2 stages local, slash-free symbolic links before running GenBlastG because
legacy GenBlastG constructs temporary BLAST filenames from its arguments. Run
GenBlastG through HRPv2 or reproduce the same local-basename behaviour in a
manual test.

### Java errors or unsupported class version

Check that Java comes from `hrpv2_env`:

```bash
conda activate hrpv2_env
command -v java
java -version
```

This guide pins OpenJDK 11 for InterProScan 5.78-109.0. If a different
InterProScan release is selected, follow the Java requirement documented for
that exact release.

### InterProScan is killed or the disk fills

Reduce `--threads`, verify RAM and free space, and inspect temporary storage:

```bash
free -h
df -h
```

<a id="reproducibility"></a>
## 13. Reproducibility and environment export

Record exact Conda packages:

```bash
conda activate hrpv2_env
conda list --explicit > hrpv2_env-linux-64.lock.txt
conda env export --no-builds > hrpv2_env.yml
```

Record external versions:

```bash
{
  python3 --version
  perl -v | head -2
  java -version
  interproscan.sh --version
  hmmbuild -h | head -2
  meme -version
  mast -version
  mafft --version
} > hrpv2_software_versions.txt 2>&1
```

Preserve together with each analysis:

- the HRPv2 script used;
- `HRPv2.log`;
- the HRPv2 manifest JSON;
- `hrpv2_env-linux-64.lock.txt`;
- `hrpv2_software_versions.txt`;
- the InterProScan release number.

<a id="resources"></a>
## 14. Official resources

- [InterProScan download instructions](https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html)
- [InterProScan command-line documentation](https://interproscan-docs.readthedocs.io/en/v5/HowToRun.html)
- [HMMER documentation](https://hmmer.org/documentation.html)
- [HMMER Bioconda recipe](https://bioconda.github.io/recipes/hmmer/README.html)
- [MEME Suite documentation](https://meme-suite.org/meme/doc/overview.html)
- [MEME Suite Bioconda recipe](https://bioconda.github.io/recipes/meme/README.html)
- [MAFFT website](https://mafft.cbrc.jp/alignment/software/)
- [GenBlastG publication and original availability information](https://doi.org/10.1093/bioinformatics/btr342)
- [GenBlastG Bioconda package](https://anaconda.org/bioconda/genblastg)
- [Bioconda usage documentation](https://bioconda.github.io/)
- [Miniforge releases](https://github.com/conda-forge/miniforge/releases)

---

Installation procedure prepared for HRPv2 on Ubuntu/Linux.

## 15. Conda installation

For the recommended Conda-based installation, see the [Conda installation guide](conda/README.md).

