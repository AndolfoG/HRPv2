# HRPv2 Conda recipe
This directory contains the Conda recipe for HRPv2.

## Repository layout
Keep `meta.yaml`, `build.sh`, and `run_test.sh` in the repository's `conda/`
directory. The HRPv2 program, wrapper, licence files, and documentation remain
in the repository root.

## Release asset
Create the GitHub release tag `v2.0.0` and upload the supplied, unmodified file
`HRPv2-2.0.0.tar.gz` as a release asset. Its SHA-256 checksum is:

```text
eeca6adad28f95f69bad801d5a19cd5e7814f47cab6fb4761d3e998ad7b30856
```

The filename, tag, archive contents, and checksum must not be changed unless
the `source` section of `meta.yaml` is updated accordingly.

## Building
Clone the HRPv2 repository and build the Conda package from the repository root:

```bash
git clone https://github.com/AndolfoG/HRPv2.git
cd HRPv2
conda build conda
```

## Local installation
After the build completes successfully, install HRPv2 in a dedicated environment:

```bash
conda create -n hrpv2 \
  -c local \
  -c conda-forge \
  -c bioconda \
  hrpv2=2.0.0
```

Activate the environment with:

```bash
conda activate hrpv2
```

## NOTES
InterProScan is not included in the Conda package. Install it separately and
pass the absolute path to `interproscan.sh` through `--interproscan-bin`.
