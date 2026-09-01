#!/usr/bin/env python3
"""# Introduction

Welcome to to the full-length **H**omology-based **R**-gene **P**rediction version 2.0 (**HRPv2**) project. **HRPv2** is a command-line workflow for the genome-wide prediction, classification and filtering of plant NB-LRR resistance genes.
"""
from __future__ import annotations


import argparse
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


__version__ = "HRPv2"


@dataclass(frozen=True)
class Hit:
    analysis: str
    accession: str
    description: str
    start: int
    end: int
    interpro_accession: str
    interpro_description: str
    kind: str


@dataclass(frozen=True)
class Evidence:
    kind: str
    start: int
    end: int
    sources: tuple[str, ...]
    accessions: tuple[str, ...]


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current: str | None = None
    with path.open() as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"Empty FASTA header at line {number}")
                current = fields[0]
                if current in lengths:
                    raise ValueError(f"Duplicate FASTA identifier: {current}")
                lengths[current] = 0
            elif current is None:
                raise ValueError(f"Sequence before first FASTA header at line {number}")
            else:
                lengths[current] += len(line.replace("*", ""))
    if not lengths:
        raise ValueError("No protein sequences found in FASTA")
    return lengths


def hit_kind(analysis: str, accession: str, description: str,
             interpro_accession: str, interpro_description: str) -> str | None:
    analysis_l = analysis.lower()
    accession_u = accession.upper()
    interpro_u = interpro_accession.upper()
    text = " ".join((description, interpro_accession, interpro_description)).lower()
    if (analysis_l == "coils" or accession_u == "PF18052"
            or "coiled-coil" in text or "coiled coil" in text):
        return "CC"
    if accession_u in {
        "PF01582", "PF13676", "SM00255", "SMART00255", "SSF52200",
        "G3DSA:3.40.50.10140", "CL23801",
    } or interpro_u in {"IPR000157", "IPR035897"}:
        return "TIR"
    if (accession_u == "PF25019"
            or (accession_u.split(":", 1)[0] == "PTHR15140"
                and interpro_u == "IPR032675")
            or (analysis_l == "superfamily"
                and accession_u in {"SSF52047", "SSF52058"})):
        return "LRR"
    if analysis_l == "superfamily" and accession_u == "SSF52540":
        return "NB_SUPPORT"
    if (accession_u in {"PF00931", "PF23559", "PF12061"}
            or interpro_u == "IPR002182" or "nb-arc" in text):
        return "NB"
    if "toll/interleukin" in text or "tir domain" in text:
        return "TIR"
    if (accession_u in {
            "PF05659", "PFAM05659", "PS51153", "CL44420",
        }
            or accession_u.startswith("PTHR36766")
            or interpro_u == "IPR008808"):
        return "RPW8"
    if accession_u.split(":", 1)[0] in PANTHER_PARTIAL_NLR_ACCESSIONS:
        return "LRR"
    if "leucine-rich repeat" in text or "leucine rich repeat" in text or " lrr" in f" {text}":
        return "LRR"
    return None


INCOMPATIBLE_DOMAIN_APPLICATIONS = frozenset({
    "pfam", "superfamily", "gene3d", "smart", "cdd",
    "prositeprofiles",
})

FAMILY_SUPPORT_ANALYSES = frozenset({"panther", "funfam", "prints"})

DIRECT_PANTHER_LRR_SIGNATURES = frozenset({("PTHR15140", "IPR032675")})
PANTHER_PARTIAL_NLR_ACCESSIONS = frozenset({"PTHR11017", "PTHR23155"})
DIRECT_PARTIAL_L_ACCESSIONS = frozenset({"PF25019"})

CORE_MAX_OVERLAP_AA = 50
PRIMARY_NB_MIN_LRR_AA = 60
PRIMARY_NB_ACCESSIONS = frozenset({"PF00931", "PF23559", "PF12061", "HRP_NB_ARC"})
STRUCTURAL_NB_MIN_AA = 90
STRUCTURAL_NB_RESCUE_MIN_AA = 50
STRUCTURAL_NB_RESCUE_MIN_PROTEIN_AA = 300
STRUCTURAL_NB_MIN_LRR_AA = 80
STRUCTURAL_NB_RESCUE_MIN_LRR_SOURCES = 2

PANTHER_LRR_FAMILY_ACCESSIONS = PANTHER_PARTIAL_NLR_ACCESSIONS
PANTHER_LRR_INTERPRO_ACCESSION = "IPR044974"
PANTHER_LRR_MIN_CTERMINAL_AA = 60

PARTIAL_L_NLR_SUPPORT_TERMS = (
    "disease resistance",
    "plant resistance",
    "resistance protein",
    "nbs-lrr",
    "nbs lrr",
    "nb-lrr",
    "nb lrr",
    "nlr protein",
    "nlr-like",
    "nlr like",
)

PARTIAL_L_NON_NLR_TERMS = (
    "f-box",
    "f box",
    "receptor-like",
    "receptor like",
    "receptor-like kinase",
    "receptor like kinase",
    "protein kinase",
    "serine/threonine-protein kinase",
    "serine threonine protein kinase",
    "cysteine-containing",
    "cysteine containing",
    "ve resistance",
    "extracellular",
)


def normalized_function_description(*values: str) -> str:
    """Normalize IPS prose so stable keyword combinations survive DB updates."""
    return " ".join(re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()
                    for value in values if value).strip()


def partial_negative_description_reason(description: str,
                                        interpro_description: str) -> str | None:
    """Return the functional-description rule incompatible with partial NB-LRRs."""
    text = normalized_function_description(description, interpro_description)
    tokens = set(text.split())
    if "f" in tokens and "box" in tokens:
        return "f_box"
    if "extracellular" in tokens:
        return "extracellular"
    if "transmembrane" in tokens:
        return "transmembrane"
    if "signal peptide" in text:
        return "signal_peptide"
    if "non cytoplasmic" in text:
        return "non_cytoplasmic"
    if "ve resistance" in text:
        return "ve_resistance"
    if "receptor" in tokens and "kinase" in tokens:
        return "receptor_kinase"
    if "receptor" in tokens and "like" in tokens:
        return "receptor_like"
    if "receptor" in tokens and "protein" in tokens:
        return "receptor_protein"
    if {"serine", "threonine", "kinase"}.issubset(tokens):
        return "serine_threonine_kinase"
    return None


def is_incompatible_partial_l_hit(analysis: str, description: str,
                                  interpro_description: str) -> bool:
    """Flag non-NLR domains or membrane features in an LRR-only protein."""
    analysis_l = analysis.lower()
    text = f"{description} {interpro_description}".lower()
    if any(term in text for term in (
            "disease resistance", "nbs-lrr", "nb-lrr", "nlr protein",
            "apoptotic protease-activating", "apaf", "winged helix",
            "p-loop containing nucleoside triphosphate hydrolase",
            "p-loop ntpase")):
        return False
    if analysis_l == "phobius":
        return any(term in text for term in (
            "transmembrane", "signal peptide", "non cytoplasmic",
            "non-cytoplasmic",
        ))
    return analysis_l in INCOMPATIBLE_DOMAIN_APPLICATIONS


def partial_l_nlr_support_hits(hits: list[Hit]) -> list[Hit]:
    """Return structural LRR matches explicitly described as resistance/NLR-like."""
    supported: list[Hit] = []
    for hit in hits:
        if hit.kind != "LRR":
            continue
        text = f"{hit.description} {hit.interpro_description}".lower()
        if hit.analysis.lower() in FAMILY_SUPPORT_ANALYSES:
            signature = (hit.accession.upper().split(":", 1)[0],
                         hit.interpro_accession.upper())
            if (signature not in DIRECT_PANTHER_LRR_SIGNATURES
                    and signature[0] not in PANTHER_PARTIAL_NLR_ACCESSIONS):
                continue
            if any(term in text for term in PARTIAL_L_NON_NLR_TERMS):
                continue
            supported.append(hit)
            continue
        if any(term in text for term in PARTIAL_L_NON_NLR_TERMS):
            continue
        if hit.accession.upper() in DIRECT_PARTIAL_L_ACCESSIONS:
            supported.append(hit)
            continue
        if any(term in text for term in PARTIAL_L_NLR_SUPPORT_TERMS):
            supported.append(hit)
    return supported


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def consolidate(hits: list[Hit]) -> list[Evidence]:
    """Merge overlapping calls of the same biological domain across analyses."""
    coordinate_hits = [
        hit for hit in hits
        if (hit.analysis.lower() not in FAMILY_SUPPORT_ANALYSES
            or (hit.accession.upper().split(":", 1)[0],
                hit.interpro_accession.upper()) in DIRECT_PANTHER_LRR_SIGNATURES)
    ]
    direct = [hit for hit in coordinate_hits if hit.kind != "NB_SUPPORT"]
    nb_support = [hit for hit in coordinate_hits if hit.kind == "NB_SUPPORT"]

    lrr_hits = [hit for hit in direct if hit.kind == "LRR"]
    promoted_support = [
        hit for hit in nb_support
        if any(lrr.end > hit.end
               and overlap(hit.start, hit.end, lrr.start, lrr.end)
               <= CORE_MAX_OVERLAP_AA
               for lrr in lrr_hits)
    ]
    for hit in promoted_support:
        direct.append(Hit(
            hit.analysis, hit.accession, hit.description, hit.start, hit.end,
            hit.interpro_accession, hit.interpro_description, "NB",
        ))
    nb_support = [hit for hit in nb_support if hit not in promoted_support]
    evidence: list[Evidence] = []
    for kind in ("TIR", "RPW8", "CC", "NB", "LRR"):
        pending = sorted((hit for hit in direct if hit.kind == kind), key=lambda hit: (hit.start, hit.end))
        groups: list[list[Hit]] = []
        for hit in pending:
            matching = [group for group in groups
                        if any(overlap(hit.start, hit.end, member.start, member.end) for member in group)]
            if not matching:
                groups.append([hit])
            else:
                target = matching[0]
                target.append(hit)
                for extra in matching[1:]:
                    target.extend(extra)
                    groups.remove(extra)
        for group in groups:
            if kind == "NB":
                group.extend(support for support in nb_support
                             if any(overlap(support.start, support.end, member.start, member.end)
                                    for member in group))
            evidence.append(Evidence(
                kind, min(hit.start for hit in group), max(hit.end for hit in group),
                tuple(sorted({hit.analysis for hit in group})),
                tuple(sorted({hit.accession for hit in group})),
            ))
    return sorted(evidence, key=lambda item: (item.start, item.end, item.kind))


def parse_interpro(path: Path, proteins: dict[str, int]) -> dict[str, list[Hit]]:
    hits: dict[str, list[Hit]] = defaultdict(list)
    seen: dict[str, set[tuple[object, ...]]] = defaultdict(set)
    with path.open() as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\n").split("\t")
            if len(columns) < 11:
                raise ValueError(f"InterProScan TSV line {number}: expected at least 11 columns")
            protein_id, analysis = columns[0], columns[3]
            if protein_id not in proteins:
                raise ValueError(f"InterProScan protein absent from FASTA: {protein_id}")
            try:
                reported_length = int(columns[2])
                start, end = int(columns[6]), int(columns[7])
            except ValueError as exc:
                raise ValueError(f"InterProScan TSV line {number}: invalid numeric field") from exc
            if reported_length != proteins[protein_id]:
                raise ValueError(
                    f"Length mismatch for {protein_id}: FASTA={proteins[protein_id]}, TSV={reported_length}")
            ipr_accession = columns[11] if len(columns) > 11 and columns[11] != "-" else ""
            ipr_description = columns[12] if len(columns) > 12 and columns[12] != "-" else ""
            kind = hit_kind(analysis, columns[4], columns[5], ipr_accession, ipr_description)
            if kind is None:
                if not is_incompatible_partial_l_hit(
                        analysis, columns[5], ipr_description):
                    continue
                kind = "INCOMPATIBLE"
            key = (kind, analysis, columns[4], start, end)
            if key in seen[protein_id]:
                continue
            seen[protein_id].add(key)
            hits[protein_id].append(Hit(analysis, columns[4], columns[5], start, end,
                                        ipr_accession, ipr_description, kind))
    return hits


def valid_nb_lrr_core(nb: Evidence, lrr: Evidence, protein_length: int) -> bool:
    """Return True only for an ordered, sufficiently supported NB-LRR core."""
    if (lrr.end <= nb.end
            or overlap(nb.start, nb.end, lrr.start, lrr.end)
            > CORE_MAX_OVERLAP_AA):
        return False
    nb_length = nb.end - nb.start + 1
    lrr_length = lrr.end - lrr.start + 1
    accessions = {value.upper() for value in nb.accessions}

    if PRIMARY_NB_ACCESSIONS & accessions:
        return lrr_length >= PRIMARY_NB_MIN_LRR_AA

    if "SSF52540" in accessions and lrr_length >= STRUCTURAL_NB_MIN_LRR_AA:
        if nb_length >= STRUCTURAL_NB_MIN_AA:
            return True
        return (
            nb_length >= STRUCTURAL_NB_RESCUE_MIN_AA
            and protein_length >= STRUCTURAL_NB_RESCUE_MIN_PROTEIN_AA
            and len(lrr.sources) >= STRUCTURAL_NB_RESCUE_MIN_LRR_SOURCES
        )
    return False


def contextual_panther_lrr(evidence: list[Evidence], hits: list[Hit],
                           protein_length: int) -> Evidence | None:
    """Validate PTHR11017/PTHR23155 plus IPR044974 as contextual LRR evidence."""
    if any(hit.kind == "INCOMPATIBLE" for hit in hits):
        return None
    nb_hits = [item for item in evidence if item.kind == "NB"]
    lrr_hits = [item for item in evidence if item.kind == "LRR"]
    if any(valid_nb_lrr_core(nb, lrr, protein_length)
           for nb in nb_hits for lrr in lrr_hits):
        return None
    family_hits = [
        hit for hit in hits
        if hit.analysis.lower() == "panther"
        and hit.accession.upper().split(":", 1)[0] in PANTHER_LRR_FAMILY_ACCESSIONS
        and hit.interpro_accession.upper() == PANTHER_LRR_INTERPRO_ACCESSION
    ]
    for nb in sorted(nb_hits, key=lambda item: (item.start, item.end)):
        accessions = {value.upper() for value in nb.accessions}
        if not (PRIMARY_NB_ACCESSIONS & accessions):
            continue
        for family in sorted(family_hits, key=lambda hit: (-hit.end, hit.start)):
            cterminal_end = min(protein_length, family.end)
            if protein_length - nb.end < PANTHER_LRR_MIN_CTERMINAL_AA:
                continue
            if cterminal_end - nb.end < PANTHER_LRR_MIN_CTERMINAL_AA:
                continue
            return Evidence(
                "LRR", nb.end + 1, cterminal_end,
                ("PANTHER_CONTEXT",),
                (family.accession.upper().split(":", 1)[0],
                 PANTHER_LRR_INTERPRO_ACCESSION),
            )
    return None


def partial_panther_nlr_lrr(evidence: list[Evidence], hits: list[Hit]) -> Evidence | None:
    """Represent PTHR11017/PTHR23155 alone as partial NB-LRR-associated LRR evidence."""
    if any(item.kind in {"NB", "LRR"} for item in evidence):
        return None
    family_hits = [
        hit for hit in hits
        if hit.analysis.lower() == "panther"
        and hit.accession.upper().split(":", 1)[0] in PANTHER_PARTIAL_NLR_ACCESSIONS
    ]
    if not family_hits:
        return None
    family = min(family_hits, key=lambda hit: (hit.start, -hit.end, hit.accession))
    return Evidence(
        "LRR", family.start, family.end,
        ("PANTHER_PARTIAL_NB-LRR",),
        (family.accession.upper().split(":", 1)[0],),
    )


def classify(evidence: list[Evidence], protein_length: int | None = None) -> tuple[str, str, int]:
    if protein_length is None:
        protein_length = max((item.end for item in evidence), default=0)
    nb_hits = [item for item in evidence if item.kind == "NB"]
    if not nb_hits:
        kinds = {item.kind for item in evidence}
        has_lrr = "LRR" in kinds
        if "TIR" in kinds and has_lrr:
            return "partial_TL", "TIR-LRR without NB-ARC", 1
        if "RPW8" in kinds and has_lrr:
            return "partial_RL", "RPW8-LRR without NB-ARC", 1
        if "CC" in kinds and has_lrr:
            return "partial_CL", "CC-LRR without NB-ARC", 1
        if has_lrr:
            return "partial_L", "LRR without NB-ARC", 1
        if "TIR" in kinds:
            return "partial_T", "TIR without NB-ARC or LRR", 1
        if "RPW8" in kinds:
            return "partial_R", "RPW8 without NB-ARC or LRR", 1
        if any(item.kind == "CC" and "PF18052" in {
                accession.upper() for accession in item.accessions
        } for item in evidence):
            return "partial_C", "PF18052 CC without NB-ARC or LRR", 1
        return "non_NLR", "no NB-LRR-associated domain", 0

    lrr_hits = [item for item in evidence if item.kind == "LRR"]
    valid_pairs = [
        (nb, lrr) for nb in nb_hits for lrr in lrr_hits
        if valid_nb_lrr_core(nb, lrr, protein_length)
    ]
    if valid_pairs:
        nb, _core_lrr = min(
            valid_pairs,
            key=lambda pair: (pair[0].start, pair[1].start, -pair[1].end),
        )
    else:
        nb = min(nb_hits, key=lambda hit: (hit.start, -hit.end))
    upstream = [item for item in evidence
                if item.kind in {"TIR", "RPW8", "CC"} and item.start < nb.start]
    upstream_kinds = {hit.kind for hit in upstream}
    nterm = ("TIR" if "TIR" in upstream_kinds else
             "RPW8" if "RPW8" in upstream_kinds else
             "CC" if "CC" in upstream_kinds else None)

    if valid_pairs:
        if nterm == "TIR":
            return "TNL", "TIR-NB-ARC-LRR", 4
        if nterm == "RPW8":
            return "RNL", "RPW8-NB-ARC-LRR", 4
        if nterm == "CC":
            return "CNL", "CC-NB-ARC-LRR", 4
        return "NL", "NB-ARC-LRR; N-terminal domain unresolved", 3
    if nterm:
        display_nterm = "C" if nterm == "CC" else nterm
        architecture_nterm = "CC" if nterm == "CC" else nterm
        display_nterm = {"TIR": "T", "RPW8": "R"}.get(display_nterm, display_nterm)
        return (f"partial_{display_nterm}N",
                f"{architecture_nterm}-NB-ARC without downstream LRR", 2)
    return "partial_N", "NB-ARC without downstream LRR", 2


def _classification_main() -> int:
    parser = argparse.ArgumentParser(
        description="Classify CNL, TNL, RNL, NL and partial NB-LRR proteins from InterProScan.",
        epilog=AUTHORSHIP,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version="HRPv2")
    parser.add_argument("--proteins", type=Path,
                        help="Protein FASTA submitted to InterProScan")
    parser.add_argument("--interpro-tsv", type=Path,
                        help="InterProScan TSV from the HRP NB-LRR application set")
    parser.add_argument("--output", type=Path, help="Classification TSV")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress interactive banners and summaries")
    parser.add_argument("--summary-only", action="store_true",
                        help="Print only the final classification counts")
    parser.add_argument("--no-authorship", action="store_true",
                        help="Suppress only the concluding authorship block")
    args = parser.parse_args()

    if not args.quiet:
        print("-" * 72, file=sys.stderr)
        print("NB-LRR protein classification", file=sys.stderr)
        print("-" * 72, file=sys.stderr)

    if args.proteins is None or args.interpro_tsv is None:
        print("Required input files:", file=sys.stderr)
        print("  1. Protein sequences in FASTA format", file=sys.stderr)
        print("  2. InterProScan results in TSV format", file=sys.stderr)

    if args.proteins is None:
        value = input("Protein FASTA file: ").strip().strip("'\"")
        if not value:
            raise ValueError("the protein FASTA path cannot be empty")
        args.proteins = Path(value).expanduser()
    if args.interpro_tsv is None:
        value = input("InterProScan TSV file: ").strip().strip("'\"")
        if not value:
            raise ValueError("the InterProScan TSV path cannot be empty")
        args.interpro_tsv = Path(value).expanduser()

    if args.output is None:
        args.output = Path("NB-LRR_protein_classification.tsv")
    metadata = Path(str(args.output) + ".metadata.txt")
    if not args.quiet:
        print("\nOutput files that will be generated:", file=sys.stderr)
        print(f"  Classification table: {args.output}", file=sys.stderr)
        print(f"  Run metadata:         {metadata}", file=sys.stderr)
        print(file=sys.stderr)

    if not args.quiet:
        print("-" * 72, file=sys.stderr)
        print("[1/3] Loading and validating protein inputs...", file=sys.stderr)
    proteins = fasta_lengths(args.proteins)
    by_protein = parse_interpro(args.interpro_tsv, proteins)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = defaultdict(int)

    if not args.quiet:
        print("[2/3] Integrating domain evidence and classifying NB-LRR architectures...",
              file=sys.stderr)
    rows: list[list[object]] = []
    for protein_id, length in proteins.items():
        hits = sorted(by_protein.get(protein_id, []),
                      key=lambda hit: (hit.start, hit.end, hit.kind))
        evidence = consolidate(hits)
        panther_lrr = contextual_panther_lrr(evidence, hits, length)
        if panther_lrr is not None:
            evidence = sorted([*evidence, panther_lrr],
                              key=lambda item: (item.start, item.end, item.kind))
        panther_partial_lrr = partial_panther_nlr_lrr(evidence, hits)
        if panther_partial_lrr is not None:
            evidence = sorted([*evidence, panther_partial_lrr],
                              key=lambda item: (item.start, item.end, item.kind))
        label, architecture, rank = classify(evidence, length)
        incompatible = [hit for hit in hits if hit.kind == "INCOMPATIBLE"]
        family_support = [
            hit for hit in hits
            if hit.analysis.lower() in FAMILY_SUPPORT_ANALYSES
            and hit.kind in {"TIR", "RPW8", "CC", "NB", "LRR"}
        ]
        exclusion_evidence = "."
        partial_l_support: list[Hit] = []
        if label in {"partial_L", "partial_CL"}:
            if label == "partial_L":
                partial_l_support = partial_l_nlr_support_hits(hits)
            negative_descriptions = [
                (hit, reason)
                for hit in hits
                if (reason := partial_negative_description_reason(
                    hit.description, hit.interpro_description)) is not None
            ]
            if incompatible:
                label = "non_NLR"
                rank = 0
                architecture = "Partial LRR model with non-NB-LRR domain or membrane evidence"
                exclusion_evidence = ";".join(sorted({
                    f"{hit.analysis}:{hit.accession}:{hit.description}"
                    for hit in incompatible
                }))
            elif negative_descriptions:
                label = "non_NLR"
                rank = 0
                architecture = "Partial LRR model with non-NB-LRR functional description"
                exclusion_evidence = ";".join(sorted({
                    f"NEGATIVE_DESCRIPTION:{reason}:{hit.analysis}:"
                    f"{hit.accession}:{hit.description}"
                    for hit, reason in negative_descriptions
                }))
            elif label == "partial_L" and not partial_l_support:
                label = "non_NLR"
                rank = 0
                architecture = "LRR without resistance/NB-LRR-like IPS support"
                exclusion_evidence = "no_resistance_or_NB-LRR_like_LRR_signature"
        counts[label] += 1
        kinds = {item.kind for item in evidence}
        rows.append([
            protein_id, length, label, rank, architecture,
            int("CC" in kinds), int("TIR" in kinds), int("RPW8" in kinds), int("NB" in kinds),
            sum(item.kind == "LRR" for item in evidence),
            "-".join(item.kind for item in evidence) or ".",
            ";".join(f"{item.kind}:{item.start}-{item.end}:"
                     f"{','.join(item.sources)}:{','.join(item.accessions)}"
                     for item in evidence) or ".",
            ";".join(sorted({
                f"{hit.kind}:{hit.analysis}:{hit.accession}:{hit.start}-{hit.end}"
                for hit in family_support
            })) or ".",
            ";".join(sorted({
                f"{hit.analysis}:{hit.accession}:{hit.start}-{hit.end}:"
                f"{hit.description.replace(chr(9), ' ')}"
                for hit in partial_l_support
            })) or ".",
            exclusion_evidence,
        ])

    if not args.quiet:
        print("[3/3] Writing classification table and metadata...", file=sys.stderr)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "length_aa", "nlr_class", "class_rank", "architecture",
                         "has_CC", "has_TIR", "has_RPW8", "has_NB_ARC", "LRR_hits",
                         "domain_order", "domain_details", "family_support_details",
                         "partial_L_support_details", "exclusion_evidence"])
        writer.writerows(rows)

    if not args.quiet or args.summary_only:
        print(f"Proteins classified: {len(proteins)}", file=sys.stderr)
        for label in sorted(counts):
            print(f"  {label}: {counts[label]}", file=sys.stderr)
    metadata.write_text(
        f"Program: HRPv2\nVersion: {__version__}\n"
        f"Protein input: {args.proteins}\nInterProScan TSV: {args.interpro_tsv}\n"
        f"Classification output: {args.output}\n\nAuthor:\n{AUTHORSHIP}\n")
    if not args.quiet and not args.no_authorship:
        print("-" * 72, file=sys.stderr)
        print(f"Developed by:\n{AUTHORSHIP}", file=sys.stderr)
        print("-" * 72, file=sys.stderr)
    return 0


"""Filter overlapping GenBlastG predictions using NLR domain classification. The program combines three inputs: 1. the complete GenBlastG GFF; 2. the GenBlastG protein FASTA; 3. the domain-based classification produced by HRPv2.py."""


import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Model:
    model_id: str
    seqid: str
    start: int
    end: int
    strand: str
    genblastg_score: float
    gff_lines: list[str]
    nlr_class: str = ""
    class_rank: int = 0
    architecture: str = ""
    architecture_score: int = 0
    length_aa: int = 0
    lrr_hits: int = 0
    genblastg_rank: int = 0
    query_id: str = ""
    has_nb_arc: bool = False
    domain_order: str = ""

    def is_fusion_like(self) -> bool:
        """Return whether a protein architecture can represent an NB-LRR gene-fusion bridge."""
        return fusion_signature(self.domain_order) is not None

    def priority(self) -> tuple[int, int, int, float, int, int, str]:
        """Return a deterministic biological-first priority tuple."""
        return (
            self.class_rank,
            self.architecture_score,
            self.lrr_hits,
            self.genblastg_score,
            -self.genblastg_rank,
            self.length_aa,
            self.model_id,
        )


ARCHITECTURE_SCORE = {
    "CNL": 4,
    "TNL": 4,
    "RNL": 4,
    "NL": 3,
    "partial_CN": 2,
    "partial_TN": 2,
    "partial_RN": 2,
    "partial_N": 2,
    "partial_CL": 2,
    "partial_TL": 1,
    "partial_RL": 1,
    "partial_T": 1,
    "partial_R": 1,
    "partial_C": 1,
    "partial_L": 1,
    "non_NLR": 0,
}

CLASS_ALIASES = {
    "partial_CCN": "partial_CN",
    "partial_TIRN": "partial_TN",
    "partial_RPW8N": "partial_RN",
    "partial_NB": "partial_N",
    "partial_LRR": "partial_L",
    "partial_TIR": "partial_T",
    "partial_RPW8": "partial_R",
}

MIN_PROTEIN_LENGTH_AA = 50
MAX_GENE_MODEL_SPAN_BP = 20_000
MAX_COMPLETE_NLR_SPAN_BP = 25_000
COMPLETE_NLR_CLASSES = frozenset({"CNL", "TNL", "RNL"})
CORE_NLR_CLASSES = frozenset({"CNL", "TNL", "RNL", "NL"})
FUSION_TERMINAL_MODULES = frozenset({"TIR", "CC", "RPW8"})


def domain_tokens(domain_order: str) -> list[str]:
    """Return normalized domain tokens from a classifier architecture."""
    return [item.strip().upper() for item in domain_order.split("-") if item.strip()]


def fusion_signature(domain_order: str) -> str | None:
    """Identify the REPEATED_NB-LRR_CORE, MULTI_NB_BRIDGE or N/C-TERMINAL_BRIDGE fusion signature."""
    domains = domain_tokens(domain_order)
    nb_positions = [index for index, value in enumerate(domains) if value == "NB"]
    lrr_positions = [index for index, value in enumerate(domains) if value == "LRR"]
    if not nb_positions or not lrr_positions:
        return None
    first_nb = nb_positions[0]
    downstream_lrr = next((index for index in lrr_positions if index > first_nb), None)
    if downstream_lrr is None:
        return None
    for nb1 in nb_positions:
        for first_lrr in (position for position in lrr_positions
                          if position > nb1):
            for second_nb in (position for position in nb_positions
                              if position > first_lrr):
                if any(position > second_nb for position in lrr_positions):
                    return "REPEATED_NB-LRR_CORE"
    leading_modules = [value for value in domains[:first_nb]
                       if value in FUSION_TERMINAL_MODULES]
    trailing_modules = [value for value in domains[downstream_lrr + 1:]
                        if value in FUSION_TERMINAL_MODULES]
    if any(value in trailing_modules for value in leading_modules):
        return "N/C-TERMINAL_BRIDGE"
    if len(nb_positions) >= 2:
        return "MULTI_NB_BRIDGE"
    return None


def nc_terminal_bridge_domain(domain_order: str) -> str | None:
    """Return the repeated CC/TIR/RPW8 domain defining an N/C bridge."""
    domains = domain_tokens(domain_order)
    nb_positions = [index for index, value in enumerate(domains) if value == "NB"]
    lrr_positions = [index for index, value in enumerate(domains) if value == "LRR"]
    for nb in nb_positions:
        lrr = next((position for position in lrr_positions if position > nb), None)
        if lrr is None:
            continue
        leading = [value for value in domains[:nb]
                   if value in FUSION_TERMINAL_MODULES]
        trailing = set(value for value in domains[lrr + 1:]
                       if value in FUSION_TERMINAL_MODULES)
        for value in reversed(leading):
            if value in trailing:
                return value
    return None


def parse_attributes(text: str) -> dict[str, str]:
    return dict(item.split("=", 1) for item in text.split(";") if "=" in item)


def parse_genblastg_rank(model_id: str) -> int:
    """Extract the primary GenBlastG rank from IDs containing '-R<number>-'."""
    match = re.search(r"(?:^|-)R(\d+)(?:-|$)", model_id)
    if not match:
        raise ValueError(f"cannot extract GenBlastG rank from model ID: {model_id}")
    return int(match.group(1))


def parse_query_id(model_id: str) -> str:
    """Recover the GenBlastG query identifier preceding the rank field."""
    match = re.match(r"^(.*)-R\d+(?:-|$)", model_id)
    if not match or not match.group(1):
        raise ValueError(f"cannot extract query ID from GenBlastG model ID: {model_id}")
    return match.group(1)


def parse_gff(path: Path) -> tuple[list[str], dict[str, Model]]:
    """Read transcript records and attach every child feature to its parent."""
    preamble: list[str] = []
    models: dict[str, Model] = {}
    child_lines: list[tuple[str, str]] = []

    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("##gff-version") and line not in preamble:
                    preamble.append(line)
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"invalid GFF line {line_number}: expected 9 columns")
            attrs = parse_attributes(fields[8])
            if fields[2] == "transcript":
                model_id = attrs.get("ID")
                if not model_id:
                    raise ValueError(f"transcript without ID at GFF line {line_number}")
                if model_id in models:
                    raise ValueError(f"duplicate transcript ID in GFF: {model_id}")
                if fields[6] not in {"+", "-"}:
                    raise ValueError(f"unsupported strand for {model_id}: {fields[6]}")
                try:
                    score = float(fields[5])
                except ValueError as exc:
                    raise ValueError(f"invalid GenBlastG score for {model_id}") from exc
                models[model_id] = Model(
                    model_id=model_id,
                    seqid=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    genblastg_score=score,
                    gff_lines=[line],
                    genblastg_rank=parse_genblastg_rank(model_id),
                    query_id=parse_query_id(model_id),
                )
            else:
                parent = attrs.get("Parent")
                if not parent:
                    raise ValueError(f"child feature without Parent at GFF line {line_number}")
                child_lines.append((parent, line))

    for parent, line in child_lines:
        if parent not in models:
            raise ValueError(f"GFF child refers to unknown parent: {parent}")
        models[parent].gff_lines.append(line)
    if not models:
        raise ValueError("the GFF contains no transcript records")
    return preamble or ["##gff-version 3"], models


def add_classification(path: Path, models: dict[str, Model]) -> None:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "length_aa", "nlr_class", "class_rank",
                    "architecture", "LRR_hits", "domain_order"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"classification TSV is missing columns: {', '.join(missing)}")
        seen: set[str] = set()
        for row in reader:
            model_id = row["protein_id"]
            if model_id not in models:
                raise ValueError(f"classification ID absent from GFF: {model_id}")
            if model_id in seen:
                raise ValueError(f"duplicate classification ID: {model_id}")
            seen.add(model_id)
            model = models[model_id]
            model.nlr_class = CLASS_ALIASES.get(row["nlr_class"], row["nlr_class"])
            model.class_rank = int(row["class_rank"])
            model.architecture = row["architecture"]
            model.architecture_score = ARCHITECTURE_SCORE.get(model.nlr_class, model.class_rank)
            model.length_aa = int(row["length_aa"])
            model.lrr_hits = int(row["LRR_hits"])
            model.has_nb_arc = row.get("has_NB_ARC", "") == "1" or "NB" in model.architecture
            model.domain_order = row["domain_order"]
    missing = set(models) - seen
    if missing:
        raise ValueError(f"{len(missing)} GFF models lack classification; example: {sorted(missing)[0]}")


def overlaps(a: Model, b: Model) -> bool:
    return (a.seqid == b.seqid and a.strand == b.strand
            and a.start <= b.end and b.start <= a.end)


def genomic_overlap_bp(a: Model, b: Model) -> int:
    """Return inclusive same-strand genomic overlap, or zero."""
    if a.seqid != b.seqid or a.strand != b.strand:
        return 0
    return max(0, min(a.end, b.end) - max(a.start, b.start) + 1)


def is_complete_nlr(model: Model) -> bool:
    """Return whether a model has a canonical complete NLR classification."""
    return model.nlr_class in COMPLETE_NLR_CLASSES and model.has_nb_arc and model.lrr_hits > 0


def is_core_nlr(model: Model) -> bool:
    """Return whether a model contains a classified NB--LRR core."""
    return (model.nlr_class in CORE_NLR_CLASSES and model.has_nb_arc
            and model.lrr_hits > 0)


def supports_terminal_module(model: Model, module: str) -> bool:
    """Return whether a model supports the specified NLR terminal module."""
    domains = set(domain_tokens(model.domain_order))
    class_support = {
        "TIR": {"TNL", "partial_TN", "partial_TL", "partial_T"},
        "CC": {"CNL", "partial_CN", "partial_CL", "partial_C"},
        "RPW8": {"RNL", "partial_RN", "partial_RL", "partial_R"},
    }
    return module in domains or model.nlr_class in class_support[module]


def distinct_ordered_loci(a: Model, b: Model) -> bool:
    """Require two strictly non-overlapping loci on one sequence and strand."""
    if a.seqid != b.seqid or a.strand != b.strand:
        return False
    left, right = sorted((a, b), key=lambda item: (item.start, item.end))
    return left.end < right.start


def protein_ordered_models(a: Model, b: Model) -> tuple[Model, Model]:
    """Return models in 5'-to-3' protein order on their shared strand."""
    genomic = sorted((a, b), key=lambda item: (item.start, item.end))
    return (genomic[0], genomic[1]) if a.strand == "+" else (genomic[1], genomic[0])


def compatible_fusion_pair(a: Model, b: Model, signature: str,
                           terminal_module: str | None = None) -> bool:
    """Validate two alternatives against the architecture of their bridge."""
    if not distinct_ordered_loci(a, b):
        return False
    if signature in {"REPEATED_NB-LRR_CORE", "MULTI_NB_BRIDGE"}:
        return is_core_nlr(a) and is_core_nlr(b)
    if signature != "N/C-TERMINAL_BRIDGE" or terminal_module is None:
        return False
    required_class = {"CC": "CNL", "TIR": "TNL", "RPW8": "RNL"}[terminal_module]
    core_model, terminal_model = protein_ordered_models(a, b)
    return (core_model.nlr_class == required_class
            and is_core_nlr(core_model)
            and supports_terminal_module(terminal_model, terminal_module))


def preserves_bridge_full_length(bridge: Model, a: Model, b: Model,
                                 signature: str,
                                 terminal_module: str | None = None) -> bool:
    """Validate whether two GenBlastG models biologically resolve an NB-LRR gene fusion."""
    core_count = sum(is_core_nlr(model) for model in (a, b))
    if core_count >= 2:
        return True
    if core_count == 0:
        return False
    if signature in {"REPEATED_NB-LRR_CORE", "MULTI_NB_BRIDGE"}:
        return False
    return compatible_fusion_pair(a, b, signature, terminal_module)


def fusion_query_ids(path: Path | None) -> dict[str, str]:
    """Map fusion-like query IDs to their repeated architectural signature."""
    if path is None:
        return {}
    fusion: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "domain_order"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"query classification is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            signature = fusion_signature(row["domain_order"])
            if signature is not None:
                fusion[row["protein_id"]] = signature
    return fusion


def fusion_bridge_models(models: dict[str, Model], fusion_queries: dict[str, str],
                         bridge_ids: set[str] | None = None
                         ) -> tuple[set[str], list[tuple[Model, Model, Model, int]]]:
    """Replace a fusion-like spanning model with two distinct NB-LRR loci."""
    excluded: set[str] = set()
    decisions: list[tuple[Model, Model, Model, int]] = []
    baseline_winners = ([models[item] for item in bridge_ids]
                        if bridge_ids is not None else [])
    for bridge in sorted(models.values(), key=lambda item: item.priority(), reverse=True):
        if bridge_ids is not None and bridge.model_id not in bridge_ids:
            continue
        if not bridge.is_fusion_like() or bridge.model_id in excluded:
            continue
        signature = (fusion_signature(bridge.domain_order)
                     or "REPEATED_NB-LRR_CORE")
        terminal_module = (nc_terminal_bridge_domain(bridge.domain_order)
                           if signature == "N/C-TERMINAL_BRIDGE" else None)
        candidates = [model for model in models.values()
                      if model.model_id != bridge.model_id
                      and not model.is_fusion_like()
                      and model.nlr_class != "non_NLR"
                      and overlaps(bridge, model)
                      and not any(winner.model_id != bridge.model_id
                                  and overlaps(model, winner)
                                  for winner in baseline_winners)]
        non_overlapping_pairs = [
            (left, right) for index, left in enumerate(candidates)
            for right in candidates[index + 1:]
            if distinct_ordered_loci(left, right)
        ]
        pairs = [pair for pair in non_overlapping_pairs
                 if compatible_fusion_pair(pair[0], pair[1], signature,
                                           terminal_module)
                 and preserves_bridge_full_length(
                     bridge, pair[0], pair[1], signature, terminal_module)]
        if not pairs:
            continue
        left, right = max(
            pairs,
            key=lambda pair: (
                sum(is_complete_nlr(model) for model in pair),
                sum(is_core_nlr(model) for model in pair),
                pair[0].class_rank + pair[1].class_rank,
                pair[0].architecture_score + pair[1].architecture_score,
                pair[0].lrr_hits + pair[1].lrr_hits,
                pair[0].genblastg_score + pair[1].genblastg_score,
                -(pair[0].genblastg_rank + pair[1].genblastg_rank),
                pair[0].length_aa + pair[1].length_aa,
            ))
        if left.start > right.start:
            left, right = right, left
        excluded.add(bridge.model_id)
        overlap_bp = ((min(bridge.end, left.end) - max(bridge.start, left.start) + 1)
                      + (min(bridge.end, right.end) - max(bridge.start, right.start) + 1))
        decisions.append((bridge, left, right, overlap_bp))
    return excluded, decisions


def select_models(models: dict[str, Model], protected: set[str] | None = None,
                  allowed_protected_overlaps: set[frozenset[str]] | None = None
                  ) -> tuple[set[str], list[tuple[Model, Model, int]]]:
    """Keep priority winners, optionally seeding validated replacement models."""
    protected = protected or set()
    allowed_protected_overlaps = allowed_protected_overlaps or set()
    unknown = protected - set(models)
    if unknown:
        raise ValueError(f"protected model is absent: {sorted(unknown)[0]}")
    protected_models = sorted((models[item] for item in protected),
                              key=lambda item: item.priority(), reverse=True)
    for index, model in enumerate(protected_models):
        for other in protected_models[index + 1:]:
            if not overlaps(model, other):
                continue
            pair_key = frozenset({model.model_id, other.model_id})
            if pair_key not in allowed_protected_overlaps:
                raise ValueError("unvalidated overlap between protected fusion replacements")
    kept: list[Model] = list(protected_models)
    removed: list[tuple[Model, Model, int]] = []
    for candidate in sorted((model for key, model in models.items()
                             if key not in protected),
                            key=lambda item: item.priority(), reverse=True):
        blockers = [winner for winner in kept if overlaps(candidate, winner)]
        if not blockers:
            kept.append(candidate)
            continue
        winner = max(blockers, key=lambda item: item.priority())
        overlap_bp = min(candidate.end, winner.end) - max(candidate.start, winner.start) + 1
        removed.append((candidate, winner, overlap_bp))
    return {item.model_id for item in kept}, removed


def read_fasta(path: Path) -> tuple[list[str], dict[str, str]]:
    order: list[str] = []
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current or current in records:
                    raise ValueError(f"invalid or duplicate FASTA ID at line {line_number}")
                order.append(current)
                records[current] = [line]
            elif current is None:
                if line.strip():
                    raise ValueError(f"sequence before first FASTA header at line {line_number}")
            else:
                records[current].append(line)
    return order, {key: "\n".join(value) + "\n" for key, value in records.items()}


def write_filtered_fasta(input_path: Path, output_path: Path, kept: set[str], expected: set[str]) -> None:
    order, records = read_fasta(input_path)
    ids = set(records)
    if ids != expected:
        raise ValueError(
            f"FASTA/GFF ID mismatch for {input_path}: "
            f"missing={len(expected - ids)}, extra={len(ids - expected)}")
    with output_path.open("w") as handle:
        for model_id in order:
            if model_id in kept:
                handle.write(records[model_id])


def reason(loser: Model, winner: Model) -> str:
    if loser.class_rank != winner.class_rank:
        return "higher_NB-LRR_class_rank"
    if loser.architecture_score != winner.architecture_score:
        return "more_complete_domain_architecture"
    if loser.lrr_hits != winner.lrr_hits:
        return "more_LRR_evidence"
    if loser.genblastg_score != winner.genblastg_score:
        return "higher_GenBlastG_score"
    if loser.genblastg_rank != winner.genblastg_rank:
        return "better_GenBlastG_rank"
    if loser.length_aa != winner.length_aa:
        return "longer_protein"
    return "deterministic_ID_tiebreak"


def _filtering_main() -> int:
    parser = argparse.ArgumentParser(description="Filter same-strand overlapping GenBlastG predictions.")
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--proteins", required=True, type=Path)
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--query-classification", type=Path,
                        help="Classification of the original queries, used to detect gene fusions")
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()

    preamble, models = parse_gff(args.gff)
    add_classification(args.classification, models)

    non_nlr = [model for model in models.values() if model.nlr_class == "non_NLR"]

    short_proteins = [model for model in models.values()
                      if model.nlr_class != "non_NLR"
                      and model.length_aa < MIN_PROTEIN_LENGTH_AA]
    ordinary_eligible = {
        model_id: model for model_id, model in models.items()
        if model.nlr_class != "non_NLR"
        and model.length_aa >= MIN_PROTEIN_LENGTH_AA
        and model.end - model.start + 1 <= MAX_GENE_MODEL_SPAN_BP
    }
    extended_complete = {
        model_id: model for model_id, model in models.items()
        if is_complete_nlr(model)
        and model.length_aa >= MIN_PROTEIN_LENGTH_AA
        and MAX_GENE_MODEL_SPAN_BP < model.end - model.start + 1
        <= MAX_COMPLETE_NLR_SPAN_BP
        and not any(overlaps(model, ordinary)
                    for ordinary in ordinary_eligible.values())
    }
    rescued_extended, _ = select_models(extended_complete)
    eligible_models = dict(ordinary_eligible)
    eligible_models.update({identifier: extended_complete[identifier]
                            for identifier in rescued_extended})
    long_models = [model for model in models.values()
                   if model.nlr_class != "non_NLR"
                   and model.length_aa >= MIN_PROTEIN_LENGTH_AA
                   and model.end - model.start + 1 > MAX_GENE_MODEL_SPAN_BP
                   and model.model_id not in rescued_extended]
    baseline_kept, _ = select_models(eligible_models)
    fusion_queries = fusion_query_ids(args.query_classification)
    fusion_bridges, fusion_decisions = fusion_bridge_models(
        eligible_models, fusion_queries, baseline_kept)
    protected_replacements = {
        model.model_id
        for _, left, right, _ in fusion_decisions
        for model in (left, right)
    }
    allowed_replacement_overlaps = {
        frozenset({left.model_id, right.model_id})
        for _, left, right, _ in fusion_decisions
        if overlaps(left, right)
    }
    overlap_candidates = {identifier: model for identifier, model in eligible_models.items()
                          if identifier not in fusion_bridges}
    kept, removed = select_models(overlap_candidates, protected_replacements,
                                  allowed_replacement_overlaps)

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    gff_output = Path(str(prefix) + ".gff")
    protein_output = Path(str(prefix) + ".pro")
    report_output = Path(str(prefix) + ".overlap_report.tsv")
    non_nlr_output = Path(str(prefix) + ".excluded_non_NB-LRR.tsv")
    thresholds_output = Path(str(prefix) + ".excluded_thresholds.tsv")
    loci_output = Path(str(prefix) + ".retained_loci.tsv")

    with gff_output.open("w") as handle:
        for line in preamble:
            handle.write(line + "\n")
        for model in models.values():
            if model.model_id in kept:
                handle.write("\n".join(model.gff_lines) + "\n")
    write_filtered_fasta(args.proteins, protein_output, kept, set(models))

    with report_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["removed_id", "retained_id", "seqid", "strand", "overlap_bp",
                         "removed_class", "retained_class", "removed_GenBlastG_score",
                         "retained_GenBlastG_score", "removed_GenBlastG_rank",
                         "retained_GenBlastG_rank", "selection_reason",
                         "replacement_overlap_bp", "replacement_overlap_pct_model_1",
                         "replacement_overlap_pct_model_2"])
        for loser, winner, overlap_bp in sorted(removed, key=lambda item: item[0].model_id):
            writer.writerow([loser.model_id, winner.model_id, loser.seqid, loser.strand, overlap_bp,
                             loser.nlr_class, winner.nlr_class, loser.genblastg_score,
                             winner.genblastg_score, loser.genblastg_rank,
                             winner.genblastg_rank, reason(loser, winner), "", "", ""])
        for bridge, left, right, overlap_bp in fusion_decisions:
            replacement_overlap = genomic_overlap_bp(left, right)
            left_pct = 100.0 * replacement_overlap / (left.end - left.start + 1)
            right_pct = 100.0 * replacement_overlap / (right.end - right.start + 1)
            bridge_signature = fusion_signature(bridge.domain_order)
            fusion_reason = {
                "MULTI_NB_BRIDGE":
                    "MULTI_NB_BRIDGE_resolved_by_distinct_NB-LRR_loci",
                "REPEATED_NB-LRR_CORE":
                    "REPEATED_NB-LRR_CORE_resolved_by_distinct_NB-LRR_loci",
                "N/C-TERMINAL_BRIDGE":
                    "N/C-TERMINAL_BRIDGE_resolved_by_compatible_loci",
            }.get(bridge_signature,
                  "fusion_model_resolved_by_distinct_NB-LRR_loci")
            writer.writerow([bridge.model_id, f"{left.model_id},{right.model_id}", bridge.seqid,
                             bridge.strand, overlap_bp, bridge.nlr_class,
                             f"{left.nlr_class},{right.nlr_class}", bridge.genblastg_score,
                             f"{left.genblastg_score},{right.genblastg_score}",
                             bridge.genblastg_rank,
                             f"{left.genblastg_rank},{right.genblastg_rank}",
                             fusion_reason,
                             replacement_overlap, f"{left_pct:.4f}", f"{right_pct:.4f}"])

    with non_nlr_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["excluded_id", "chromosome", "start", "end", "strand",
                         "classification", "exclusion_reason"])
        for model in sorted(non_nlr,
                            key=lambda item: (item.seqid, item.start, item.end, item.strand)):
            writer.writerow([model.model_id, model.seqid, model.start, model.end, model.strand,
                             model.nlr_class, "no_NB-LRR_associated_domain"])

    with thresholds_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["excluded_id", "chromosome", "start", "end", "strand",
                         "nlr_class", "length_aa", "gene_model_span_bp",
                         "exclusion_reason"])
        threshold_exclusions = [
            (model, "protein_shorter_than_50_aa") for model in short_proteins
        ] + [
            (model, ("complete_NB-LRR_20001_to_25000_bp_overlaps_existing_NB-LRR_locus"
                     if is_complete_nlr(model)
                     and model.end - model.start + 1 <= MAX_COMPLETE_NLR_SPAN_BP
                     else "gene_model_longer_than_allowed_span"))
            for model in long_models
        ]
        for model, exclusion_reason in sorted(
                threshold_exclusions,
                key=lambda item: (item[0].seqid, item[0].start, item[0].end, item[0].strand)):
            writer.writerow([model.model_id, model.seqid, model.start, model.end, model.strand,
                             model.nlr_class, model.length_aa,
                             model.end - model.start + 1, exclusion_reason])

    with loci_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "chromosome", "start", "end", "strand",
                         "nlr_class", "architecture", "length_aa", "GenBlastG_score",
                         "GenBlastG_rank"])
        retained_models = (model for model in models.values() if model.model_id in kept)
        for model in sorted(retained_models,
                            key=lambda item: (item.seqid, item.start, item.end, item.strand)):
            writer.writerow([model.model_id, model.seqid, model.start, model.end, model.strand,
                             model.nlr_class, model.architecture, model.length_aa,
                             model.genblastg_score, model.genblastg_rank])

    print(f"Input models: {len(models)}", file=sys.stderr)
    print(f"Excluded non-NB-LRR models: {len(non_nlr)}", file=sys.stderr)
    print(f"Excluded proteins shorter than 50 aa: {len(short_proteins)}", file=sys.stderr)
    print(f"Complete 20-25 kb models rescued in isolated regions: {len(rescued_extended)}",
          file=sys.stderr)
    print(f"Excluded gene models exceeding their permitted span: {len(long_models)}",
          file=sys.stderr)
    print(f"Retained models: {len(kept)}", file=sys.stderr)
    print(f"Removed overlapping models: {len(removed)}", file=sys.stderr)
    print(f"Fusion-spanning models replaced by two loci: {len(fusion_decisions)}", file=sys.stderr)
    print(f"Filtered GFF: {gff_output}", file=sys.stderr)
    print(f"Filtered proteins: {protein_output}", file=sys.stderr)
    print(f"Decision report: {report_output}", file=sys.stderr)
    print(f"Non-NB-LRR exclusions: {non_nlr_output}", file=sys.stderr)
    print(f"Threshold exclusions: {thresholds_output}", file=sys.stderr)
    print(f"Retained loci table: {loci_output}", file=sys.stderr)
    return 0


"""Recover annotated partial NB-LRR genes missed by GenBlastG."""


import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path


FULL_LENGTH_CLASSES = frozenset({"CNL", "TNL", "RNL", "NL"})
PARTIAL_NLR_CLASSES = frozenset({
    "partial_CN", "partial_N", "partial_CL", "partial_C",
    "partial_TN", "partial_TL", "partial_T",
    "partial_RN", "partial_RL", "partial_R", "partial_L",
})
MATCH_KEYS = ("ID", "Name", "protein_id", "transcript_id")
GFF_ID_PREFIXES = frozenset({
    "gene", "mrna", "transcript", "rna", "protein", "polypeptide",
})


@dataclass(frozen=True)
class Feature:
    index: int
    raw: str
    seqid: str
    feature_type: str
    start: int
    end: int
    strand: str
    attributes: dict[str, str]


@dataclass(frozen=True)
class GFFAnomaly:
    line_number: int
    raw: str
    seqid: str
    feature_type: str
    start: str
    end: str
    strand: str
    action: str
    reason: str


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(None, 1)
            value = value.strip().strip('"')
        else:
            continue
        attributes[key.strip()] = value.strip()
    return attributes


def read_gff(path: Path) -> tuple[list[str], list[Feature], list[GFFAnomaly]]:
    """Read valid GFF3 features and record malformed or normalized entries for auditing."""
    headers: list[str] = []
    features: list[Feature] = []
    anomalies: list[GFFAnomaly] = []
    invalid_root_ids: set[str] = set()
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                headers.append(line)
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                anomalies.append(GFFAnomaly(
                    line_number, line, fields[0] if fields else "",
                    fields[2] if len(fields) > 2 else "",
                    fields[3] if len(fields) > 3 else "",
                    fields[4] if len(fields) > 4 else "",
                    fields[6] if len(fields) > 6 else "",
                    "skipped_feature", "expected_9_columns",
                ))
                continue
            attributes = parse_attributes(fields[8])
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                feature_type = fields[2].lower()
                action = ("skipped_root_and_descendants"
                          if feature_type in {"gene", "mrna", "transcript"}
                          else "skipped_feature")
                anomalies.append(GFFAnomaly(
                    line_number, line, fields[0], fields[2], fields[3],
                    fields[4], fields[6], action, "non_integer_coordinates",
                ))
                if action == "skipped_root_and_descendants":
                    invalid_root_ids.update(
                        value for value in attributes.get("ID", "").split(",")
                        if value)
                continue
            if start > end:
                anomalies.append(GFFAnomaly(
                    line_number, line, fields[0], fields[2], fields[3],
                    fields[4], fields[6], "normalized_coordinate_order",
                    "start_exceeds_end",
                ))
                start, end = end, start
                fields[3], fields[4] = str(start), str(end)
                line = "\t".join(fields)
            features.append(Feature(len(features), line, fields[0], fields[2], start, end,
                                    fields[6], attributes))
    blocked_ids = set(invalid_root_ids)
    retained = features
    while blocked_ids:
        newly_blocked: set[str] = set()
        next_retained: list[Feature] = []
        for feature in retained:
            parents = set(feature.attributes.get("Parent", "").split(",")) - {""}
            if parents & blocked_ids:
                anomalies.append(GFFAnomaly(
                    0, feature.raw, feature.seqid, feature.feature_type,
                    str(feature.start), str(feature.end), feature.strand,
                    "skipped_descendant", "parent_feature_is_invalid",
                ))
                newly_blocked.update(
                    value for value in feature.attributes.get("ID", "").split(",")
                    if value)
            else:
                next_retained.append(feature)
        retained = next_retained
        blocked_ids = newly_blocked
    if not features:
        raise ValueError(f"GFF contains no features: {path}")
    if not retained:
        raise ValueError(f"GFF contains no valid features: {path}")
    return headers, retained, anomalies


def write_gff_anomalies(path: Path,
                        inputs: list[tuple[Path, list[GFFAnomaly]]]) -> None:
    """Write a machine-readable audit of every GFF anomaly found at merge."""
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "source_gff", "line_number", "seqid", "feature_type", "start",
            "end", "strand", "action", "reason", "raw_record",
        ])
        for source, anomalies in inputs:
            for item in anomalies:
                writer.writerow([
                    str(source), item.line_number or ".", item.seqid,
                    item.feature_type, item.start, item.end, item.strand,
                    item.action, item.reason, item.raw,
                ])


def preflight_annotation_gff(path: Path, protein_ids: list[str],
                             audit_path: Path,
                             anomaly_path: Path) -> dict[str, int]:
    _headers, features, anomalies = read_gff(path)
    index = build_exact_index(features)
    counts = {
        "proteins": len(protein_ids),
        "features": len(features),
        "matched": 0,
        "unmatched": 0,
        "ambiguous": 0,
        "anomalies": len(anomalies),
    }
    key_counts: dict[str, int] = defaultdict(int)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    with audit_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "protein_id", "status", "matched_attribute", "match_count",
            "matched_feature_type", "matched_gff_id",
        ])
        for identifier in protein_ids:
            match_key, matches = match_feature(identifier, index)
            if not matches:
                counts["unmatched"] += 1
                writer.writerow([identifier, "unmatched", ".", 0, ".", "."])
                continue
            key_counts[match_key] += 1
            if len(matches) != 1:
                counts["ambiguous"] += 1
                writer.writerow([
                    identifier, "ambiguous", match_key, len(matches),
                    ",".join(sorted({item.feature_type for item in matches})),
                    ",".join(sorted({item.attributes.get("ID", ".")
                                     for item in matches})),
                ])
                continue
            counts["matched"] += 1
            item = matches[0]
            writer.writerow([
                identifier, "matched", match_key, 1, item.feature_type,
                item.attributes.get("ID", "."),
            ])
    write_gff_anomalies(anomaly_path, [(path, anomalies)])
    for key in MATCH_KEYS:
        counts[f"matched_by_{key}"] = key_counts.get(key, 0)
    if counts["matched"] == 0:
        raise ValueError(
            "annotation GFF validation found no unambiguous exact match to "
            "the proteome identifiers using ID, Name, protein_id or "
            f"transcript_id; see {audit_path}")
    return counts


def read_fasta(path: Path) -> tuple[list[str], dict[str, str]]:
    order: list[str] = []
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current or current in records:
                    raise ValueError(f"invalid or duplicate FASTA ID at line {line_number}")
                order.append(current)
                records[current] = [line]
            elif current is None:
                if line.strip():
                    raise ValueError(f"sequence before first FASTA header at line {line_number}")
            else:
                records[current].append(line)
    if not records:
        raise ValueError(f"FASTA contains no sequences: {path}")
    return order, {identifier: "\n".join(lines) + "\n"
                   for identifier, lines in records.items()}


def annotation_candidates(path: Path) -> dict[str, str]:
    """Return every complete or partial NB-LRR annotation candidate."""
    candidates: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "nlr_class"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"classification TSV is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            identifier = row["protein_id"]
            if identifier in candidates:
                raise ValueError(f"duplicate classification ID: {identifier}")
            nlr_class = CLASS_ALIASES.get(row["nlr_class"], row["nlr_class"])
            if nlr_class in FULL_LENGTH_CLASSES | PARTIAL_NLR_CLASSES:
                candidates[identifier] = nlr_class
    return candidates


NB_PARTIAL_CLASSES = frozenset({
    "partial_CN", "partial_TN", "partial_RN", "partial_N",
})
NB_REPLACEMENT_CLASSES = FULL_LENGTH_CLASSES | NB_PARTIAL_CLASSES


def row_class(row: dict[str, str]) -> str:
    return CLASS_ALIASES.get(row.get("nlr_class", ""), row.get("nlr_class", ""))


def row_has_nb(row: dict[str, str]) -> bool:
    explicit = row.get("has_NB_ARC", "").strip().lower()
    if explicit in {"1", "true", "yes"}:
        return True
    return "NB" in domain_tokens(row.get("domain_order", ""))


def feature_root(identifier: str, features: list[Feature]) -> Feature | None:
    """Find the transcript/gene carrying a final protein/model identifier."""
    preferred = {"transcript", "mrna", "gene"}
    for feature in features:
        if (feature.feature_type.lower() in preferred
                and identifier in {feature.attributes.get(key, "") for key in MATCH_KEYS}):
            return feature
    for feature in features:
        if identifier in {feature.attributes.get(key, "") for key in MATCH_KEYS}:
            return feature
    return None


def coding_intervals(root: Feature, features: list[Feature]) -> list[tuple[int, int]]:
    indexes = related_feature_indexes(root, features)
    related = [item for item in features if item.index in indexes]
    for accepted in ({"cds", "coding_exon"}, {"exon"}):
        intervals = sorted((item.start, item.end) for item in related
                           if item.feature_type.lower() in accepted)
        if intervals:
            return intervals
    return [(root.start, root.end)]


def interval_overlap_bp(left: list[tuple[int, int]],
                        right: list[tuple[int, int]]) -> int:
    return sum(max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
               for a_start, a_end in left for b_start, b_end in right)


def valid_fusion_replacements(
        annotated: Feature, annotated_intervals: list[tuple[int, int]],
        predicted_roots: dict[str, Feature], predicted: list[Feature],
        predicted_rows: dict[str, dict[str, str]]) -> list[str]:
    """Select the best pair of GenBlastG models that resolves an annotated NB-LRR gene fusion."""
    eligible: list[tuple[str, Feature, int, bool]] = []
    annotated_span = (annotated.seqid, annotated.start, annotated.end,
                      annotated.strand)
    for identifier, root in predicted_roots.items():
        row = predicted_rows.get(identifier, {})
        nlr_class = row_class(row)
        if nlr_class not in NB_REPLACEMENT_CLASSES or not row_has_nb(row):
            continue
        span = (root.seqid, root.start, root.end, root.strand)
        if not span_overlaps(annotated_span, span):
            continue
        overlap = interval_overlap_bp(
            annotated_intervals, coding_intervals(root, predicted))
        if overlap:
            eligible.append((identifier, root, overlap,
                             nlr_class in FULL_LENGTH_CLASSES))
    pairs: list[tuple[tuple[int, int, int, str, str], list[str]]] = []
    for pos, left in enumerate(eligible):
        for right in eligible[pos + 1:]:
            left_id, left_root, left_overlap, left_full = left
            right_id, right_root, right_overlap, right_full = right
            if (left_root.seqid != right_root.seqid
                    or left_root.strand != right_root.strand):
                continue
            if not (left_root.end < right_root.start
                    or right_root.end < left_root.start):
                continue
            if not (left_full or right_full):
                continue
            identifiers = sorted([left_id, right_id])
            score = (int(left_full) + int(right_full),
                     left_overlap + right_overlap,
                     min(left_overlap, right_overlap),
                     identifiers[0], identifiers[1])
            pairs.append((score, identifiers))
    return max(pairs, default=((0, 0, 0, "", ""), []), key=lambda item: item[0])[1]


def valid_nc_terminal_replacements(
        annotated: Feature, annotated_intervals: list[tuple[int, int]],
        terminal_domain: str, predicted_roots: dict[str, Feature],
        predicted: list[Feature],
        predicted_rows: dict[str, dict[str, str]]) -> list[str]:
    """Resolve X--NB--LRR--X as a canonical XNL plus a terminal X locus."""
    required_class = TERMINAL_DOMAIN_CLASS[terminal_domain]
    compatible_classes = TERMINAL_PARTIAL_CLASSES[terminal_domain]
    annotated_span = (annotated.seqid, annotated.start, annotated.end,
                      annotated.strand)
    eligible: list[tuple[str, Feature, int, str]] = []
    for identifier, root in predicted_roots.items():
        row = predicted_rows.get(identifier, {})
        nlr_class = row_class(row)
        if nlr_class not in compatible_classes:
            continue
        if not row_has_terminal_domain(row, terminal_domain):
            continue
        span = (root.seqid, root.start, root.end, root.strand)
        if not span_overlaps(annotated_span, span):
            continue
        overlap = interval_overlap_bp(
            annotated_intervals, coding_intervals(root, predicted))
        if overlap:
            eligible.append((identifier, root, overlap, nlr_class))
    pairs: list[tuple[tuple[int, int, int, str, str], list[str]]] = []
    for pos, left in enumerate(eligible):
        for right in eligible[pos + 1:]:
            left_id, left_root, left_overlap, left_class = left
            right_id, right_root, right_overlap, right_class = right
            if (left_root.seqid != right_root.seqid
                    or left_root.strand != right_root.strand
                    or not (left_root.end < right_root.start
                            or right_root.end < left_root.start)):
                continue
            genomic = sorted((left, right), key=lambda item: (item[1].start,
                                                               item[1].end))
            core, terminal = (genomic if annotated.strand == "+"
                              else list(reversed(genomic)))
            if core[3] != required_class:
                continue
            identifiers = [core[0], terminal[0]]
            score = (int(terminal[3] == required_class),
                     core[2] + terminal[2], min(core[2], terminal[2]),
                     identifiers[0], identifiers[1])
            pairs.append((score, identifiers))
    return max(pairs, default=((0, 0, 0, "", ""), []),
               key=lambda item: item[0])[1]


CLASS_TERMINAL_DOMAIN = {"CNL": "CC", "TNL": "TIR", "RNL": "RPW8"}
TERMINAL_DOMAIN_CLASS = {value: key for key, value in CLASS_TERMINAL_DOMAIN.items()}
TERMINAL_PARTIAL_CLASSES = {
    "CC": frozenset({"partial_C", "partial_CN", "partial_CL", "CNL"}),
    "TIR": frozenset({"partial_T", "partial_TN", "partial_TL", "TNL"}),
    "RPW8": frozenset({"partial_R", "partial_RN", "partial_RL", "RNL"}),
}


def row_has_terminal_domain(row: dict[str, str], domain: str) -> bool:
    field = {"CC": "has_CC", "TIR": "has_TIR", "RPW8": "has_RPW8"}[domain]
    return row.get(field, "").strip().lower() in {"1", "true", "yes"}


def canonical_specific_replacements(
        annotated: Feature, annotated_intervals: list[tuple[int, int]],
        predicted_roots: dict[str, Feature], predicted: list[Feature],
        predicted_rows: dict[str, dict[str, str]]) -> list[str]:
    """Prefer a canonical predicted CNL, TNL or RNL model over a less specific annotated NL."""
    annotated_span = (annotated.seqid, annotated.start, annotated.end,
                      annotated.strand)
    retained: list[tuple[int, int, str]] = []
    for identifier, root in predicted_roots.items():
        row = predicted_rows.get(identifier, {})
        nlr_class = row_class(row)
        terminal = CLASS_TERMINAL_DOMAIN.get(nlr_class)
        if terminal is None or not row_has_nb(row):
            continue
        if not row_has_terminal_domain(row, terminal):
            continue
        try:
            if int(row.get("LRR_hits", "0") or 0) < 1:
                continue
        except ValueError:
            continue
        domains = domain_tokens(row.get("domain_order", ""))
        try:
            terminal_position = domains.index(terminal)
            nb = next(index for index, value in enumerate(domains)
                      if value == "NB" and index > terminal_position)
            next(index for index, value in enumerate(domains)
                 if value == "LRR" and index > nb)
        except (ValueError, StopIteration):
            continue
        if fusion_signature(row.get("domain_order", "")) is not None:
            continue
        span = (root.seqid, root.start, root.end, root.strand)
        if not span_overlaps(annotated_span, span):
            continue
        coding_overlap = interval_overlap_bp(
            annotated_intervals, coding_intervals(root, predicted))
        if coding_overlap:
            retained.append((coding_overlap, root.end - root.start + 1,
                             identifier))
    retained.sort(reverse=True)
    return [retained[0][2]] if retained else []


def build_exact_index(features: list[Feature]) -> dict[str, dict[str, list[Feature]]]:
    """Index exact GFF identifiers, including standard type-prefixed aliases."""
    index = {key: {} for key in MATCH_KEYS}
    for feature in features:
        for key in MATCH_KEYS:
            for value in feature.attributes.get(key, "").split(","):
                if value:
                    aliases = {value}
                    if ":" in value:
                        prefix, unprefixed = value.split(":", 1)
                        if prefix.lower() in GFF_ID_PREFIXES and unprefixed:
                            aliases.add(unprefixed)
                    for alias in aliases:
                        bucket = index[key].setdefault(alias, [])
                        if feature not in bucket:
                            bucket.append(feature)
    return index


def match_feature(identifier: str, index: dict[str, dict[str, list[Feature]]]) -> tuple[str, list[Feature]]:
    """Return the first exact, unique match by key priority."""
    for key in MATCH_KEYS:
        matches = index[key].get(identifier, [])
        if matches:
            return key, matches
    return "", []


def transcript_spans(features: list[Feature]) -> list[tuple[str, int, int, str]]:
    preferred = [feature for feature in features
                 if feature.feature_type.lower() in {"transcript", "mrna", "gene"}]
    source = preferred or features
    return [(feature.seqid, feature.start, feature.end, feature.strand) for feature in source]


def span_overlaps(span: tuple[str, int, int, str],
                  other: tuple[str, int, int, str]) -> bool:
    seqid, start, end, strand = span
    other_seqid, other_start, other_end, other_strand = other
    return (seqid == other_seqid and strand == other_strand
            and start <= other_end and other_start <= end)


def related_feature_indexes(root: Feature, features: list[Feature]) -> set[int]:
    """Include the matched record, its ancestors and its descendants."""
    selected = {root.index}
    by_id = {feature.attributes.get("ID", ""): feature for feature in features
             if feature.attributes.get("ID")}
    current = root
    seen: set[str] = set()
    while current.attributes.get("Parent"):
        parent_ids = [value for value in current.attributes["Parent"].split(",") if value]
        if len(parent_ids) != 1 or parent_ids[0] in seen or parent_ids[0] not in by_id:
            break
        seen.add(parent_ids[0])
        current = by_id[parent_ids[0]]
        selected.add(current.index)
    identifiers = set(root.attributes.get("ID", "").split(",")) - {""}
    changed = True
    while changed and identifiers:
        changed = False
        for feature in features:
            parents = set(feature.attributes.get("Parent", "").split(",")) - {""}
            if feature.index not in selected and parents & identifiers:
                selected.add(feature.index)
                identifiers.update(set(feature.attributes.get("ID", "").split(",")) - {""})
                changed = True
    return selected


def gene_identity(root: Feature, features: list[Feature]) -> tuple[str, Feature]:
    """Resolve the gene ancestor without changing its original GFF identifier."""
    by_id = {feature.attributes.get("ID", ""): feature for feature in features
             if feature.attributes.get("ID")}
    current = root
    seen: set[str] = set()
    while current.feature_type.lower() != "gene" and current.attributes.get("Parent"):
        parents = [value for value in current.attributes["Parent"].split(",") if value]
        if len(parents) != 1 or parents[0] in seen or parents[0] not in by_id:
            break
        seen.add(parents[0])
        current = by_id[parents[0]]
    return current.attributes.get("ID", root.attributes.get("ID", "")), current


def classification_rows(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "protein_id" not in reader.fieldnames:
            raise ValueError(f"classification TSV lacks protein_id: {path}")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            identifier = row["protein_id"]
            if identifier in rows:
                raise ValueError(f"duplicate classification ID: {identifier}")
            rows[identifier] = row
        return list(reader.fieldnames), rows


def _merging_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Recover annotated NB-LRR loci with fusion-aware GenBlastG replacement.")
    parser.add_argument("--annotation-gff", required=True, type=Path)
    parser.add_argument("--proteome", required=True, type=Path)
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--filtered-genblastg-gff", required=True, type=Path)
    parser.add_argument("--filtered-genblastg-proteins", type=Path,
                        help="Filtered GenBlastG proteins; enables merged outputs")
    parser.add_argument("--predicted-classification", type=Path,
                        help="Classification TSV for all predicted GenBlastG proteins")
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args(argv)

    headers, annotation, annotation_anomalies = read_gff(args.annotation_gff)
    _, predicted, predicted_anomalies = read_gff(args.filtered_genblastg_gff)
    order, proteins = read_fasta(args.proteome)
    candidates = annotation_candidates(args.classification)
    unknown = set(candidates) - set(proteins)
    if unknown:
        raise ValueError(f"annotation candidate absent from proteome: {sorted(unknown)[0]}")

    if args.filtered_genblastg_proteins and args.predicted_classification is None:
        raise ValueError("--predicted-classification is required for merged outputs")
    original_fields, original_rows = classification_rows(args.classification)
    predicted_fields: list[str] = []
    predicted_rows: dict[str, dict[str, str]] = {}
    predicted_order: list[str] = []
    predicted_protein_records: dict[str, str] = {}
    if args.predicted_classification:
        predicted_fields, predicted_rows = classification_rows(args.predicted_classification)
    if args.filtered_genblastg_proteins:
        predicted_order, predicted_protein_records = read_fasta(
            args.filtered_genblastg_proteins)
    predicted_roots = {identifier: root for identifier in predicted_rows
                       if (root := feature_root(identifier, predicted)) is not None}
    retained_predicted = set(predicted_order or predicted_roots)

    index = build_exact_index(annotation)
    recovered: dict[str, tuple[str, Feature, str, set[int]]] = {}
    decisions: list[tuple[str, str, str, str, str, str, str, str]] = []

    for identifier in order:
        if candidates.get(identifier) not in FULL_LENGTH_CLASSES:
            continue
        match_key, matches = match_feature(identifier, index)
        if not matches:
            decisions.append((identifier, candidates[identifier], "not_recovered",
                              "no_exact_annotation_match", "", "", "", ""))
            continue
        if len(matches) != 1:
            decisions.append((identifier, candidates[identifier], "not_recovered",
                              "ambiguous_exact_annotation_match", match_key, "", "", ""))
            continue
        feature = matches[0]
        gene_id, locus_feature = gene_identity(feature, annotation)
        span = (locus_feature.seqid, locus_feature.start, locus_feature.end,
                locus_feature.strand)
        signature = fusion_signature(original_rows[identifier].get("domain_order", ""))
        replacements: list[str] = []
        if candidates[identifier] == "NL":
            replacements = canonical_specific_replacements(
                locus_feature, coding_intervals(feature, annotation),
                predicted_roots, predicted, predicted_rows)
            if replacements:
                preferred_class = row_class(predicted_rows.get(replacements[0], {}))
                decisions.append((identifier, candidates[identifier], "not_recovered",
                                  "canonical_specific_GenBlastG_model_preferred",
                                  match_key, gene_id,
                                  f"{preferred_class}_OVER_NL",
                                  ",".join(replacements)))
                continue
        if signature in {"REPEATED_NB-LRR_CORE", "MULTI_NB_BRIDGE"}:
            replacements = valid_fusion_replacements(
                locus_feature, coding_intervals(feature, annotation),
                predicted_roots, predicted, predicted_rows)
        elif signature == "N/C-TERMINAL_BRIDGE":
            terminal_domain = nc_terminal_bridge_domain(
                original_rows[identifier].get("domain_order", ""))
            if terminal_domain is not None:
                replacements = valid_nc_terminal_replacements(
                    locus_feature, coding_intervals(feature, annotation),
                    terminal_domain, predicted_roots, predicted, predicted_rows)
        if replacements:
            decisions.append((identifier, candidates[identifier], "not_recovered",
                              "fusion_resolved_by_GenBlastG", match_key, gene_id,
                              signature or "", ",".join(replacements)))
            continue

        removed = sorted(identifier_ for identifier_, root in predicted_roots.items()
                         if identifier_ in retained_predicted
                         and span_overlaps(span, (root.seqid, root.start, root.end,
                                                 root.strand)))
        retained_predicted.difference_update(removed)
        recovered[identifier] = (match_key, feature, gene_id,
                                 related_feature_indexes(feature, annotation))
        decisions.append((identifier, candidates[identifier], "recovered",
                          ("official_full_length_preferred" if removed
                           else "official_full_length_retained"),
                          match_key, gene_id, signature or "", ",".join(removed)))

    retained_feature_indexes: set[int] = set()
    for identifier in retained_predicted:
        root = predicted_roots.get(identifier)
        if root is not None:
            retained_feature_indexes.update(related_feature_indexes(root, predicted))
    occupied = transcript_spans(
        [feature for feature in predicted if feature.index in retained_feature_indexes])
    for _match_key, feature, _gene_id, _indexes in recovered.values():
        _id, locus = gene_identity(feature, annotation)
        occupied.append((locus.seqid, locus.start, locus.end, locus.strand))

    for identifier in order:
        if candidates.get(identifier) not in PARTIAL_NLR_CLASSES:
            continue
        match_key, matches = match_feature(identifier, index)
        if not matches:
            decisions.append((identifier, candidates[identifier], "not_recovered",
                              "no_exact_annotation_match", "", "", "", ""))
            continue
        if len(matches) != 1:
            decisions.append((identifier, candidates[identifier], "not_recovered",
                              "ambiguous_exact_annotation_match", match_key, "", "", ""))
            continue
        feature = matches[0]
        gene_id, locus_feature = gene_identity(feature, annotation)
        span = (locus_feature.seqid, locus_feature.start, locus_feature.end,
                locus_feature.strand)
        if any(span_overlaps(span, prediction) for prediction in occupied):
            decisions.append((identifier, candidates[identifier], "not_recovered",
                              "overlaps_retained_NB-LRR_locus", match_key, gene_id, "", ""))
            continue
        recovered[identifier] = (match_key, feature, gene_id,
                                 related_feature_indexes(feature, annotation))
        occupied.append(span)
        decisions.append((identifier, candidates[identifier], "recovered",
                          "unoccupied_genomic_region", match_key, gene_id, "", ""))

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    recovered_gff_output = Path(str(prefix) + ".recovered_annotation.gff3")
    recovered_fasta_output = Path(str(prefix) + ".recovered_annotation.fasta")
    report_output = Path(str(prefix) + ".report.tsv")
    merged_gff_output = Path(str(prefix) + ".gff3")
    merged_fasta_output = Path(str(prefix) + ".fasta")
    merged_tsv_output = Path(str(prefix) + ".tsv")
    anomaly_output = Path(str(prefix) + "_gff_anomalies.tsv")
    write_gff_anomalies(anomaly_output, [
        (args.annotation_gff, annotation_anomalies),
        (args.filtered_genblastg_gff, predicted_anomalies),
    ])

    selected_indexes = {index for _, _, _, indexes in recovered.values() for index in indexes}
    with recovered_gff_output.open("w") as handle:
        version_headers = [line for line in headers if line.startswith("##gff-version")]
        for line in version_headers or ["##gff-version 3"]:
            handle.write(line + "\n")
        for feature in annotation:
            if feature.index in selected_indexes:
                handle.write(feature.raw + "\n")
    with recovered_fasta_output.open("w") as handle:
        for identifier in order:
            if identifier in recovered:
                handle.write(proteins[identifier])
    with report_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id", "nlr_class", "status", "reason",
                         "matched_attribute", "gene_id", "fusion_signature",
                         "GenBlastG_models"])
        writer.writerows(decisions)

    if args.filtered_genblastg_proteins:
        duplicate_ids = retained_predicted & set(recovered)
        if duplicate_ids:
            raise ValueError(f"duplicate protein ID while merging: {sorted(duplicate_ids)[0]}")
        with merged_gff_output.open("w") as handle:
            handle.write("##gff-version 3\n")
            for feature in predicted:
                if feature.index in retained_feature_indexes:
                    handle.write(feature.raw + "\n")
            for feature in annotation:
                if feature.index in selected_indexes:
                    handle.write(feature.raw + "\n")
        with merged_fasta_output.open("w") as handle:
            for identifier in predicted_order:
                if identifier in retained_predicted:
                    handle.write(predicted_protein_records[identifier])
            for identifier in order:
                if identifier in recovered:
                    handle.write(proteins[identifier])
        if predicted_fields != original_fields:
            raise ValueError("proteome and predicted-model classification columns differ")
        data_fields = [field for field in predicted_fields if field != "protein_id"]
        with merged_tsv_output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n",
                                    fieldnames=["gene_id", "protein_id", "source"] + data_fields)
            writer.writeheader()
            for identifier in predicted_order:
                if identifier not in retained_predicted:
                    continue
                if identifier not in predicted_rows:
                    raise ValueError(f"filtered model lacks classification: {identifier}")
                row = predicted_rows[identifier]
                writer.writerow({"gene_id": identifier, "protein_id": identifier,
                                 "source": "GenBlastG", **{key: row[key] for key in data_fields}})
            for identifier in order:
                if identifier in recovered:
                    gene_id = recovered[identifier][2]
                    row = original_rows[identifier]
                    writer.writerow({"gene_id": gene_id, "protein_id": identifier,
                                     "source": "original_annotation",
                                     **{key: row[key] for key in data_fields}})

    complete_candidates = sum(value in FULL_LENGTH_CLASSES for value in candidates.values())
    recovered_complete = sum(candidates[item] in FULL_LENGTH_CLASSES for item in recovered)
    resolved_fusions = sum(item[2] == "not_recovered" and
                           item[3] == "fusion_resolved_by_GenBlastG" for item in decisions)
    preserved_specific = sum(
        item[2] == "not_recovered" and
        item[3] == "canonical_specific_GenBlastG_model_preferred"
        for item in decisions)
    print(f"Full-length NB-LRR annotation candidates: {complete_candidates}", file=sys.stderr)
    print(f"Recovered full-length annotation loci: {recovered_complete}", file=sys.stderr)
    print(f"Fusion loci retained as GenBlastG models: {resolved_fusions}", file=sys.stderr)
    print(f"Canonical CNL/TNL/RNL models preferred over official NL: "
          f"{preserved_specific}",
          file=sys.stderr)
    print(f"Recovered annotation loci (all classes): {len(recovered)}", file=sys.stderr)
    print(f"Annotation GFF anomalies recorded: {len(annotation_anomalies)}", file=sys.stderr)
    print(f"GenBlastG GFF anomalies recorded: {len(predicted_anomalies)}", file=sys.stderr)
    print(f"GFF anomaly report: {anomaly_output}", file=sys.stderr)
    print(f"Recovered annotation: {recovered_gff_output}", file=sys.stderr)
    print(f"Recovered proteins: {recovered_fasta_output}", file=sys.stderr)
    print(f"Recovery report: {report_output}", file=sys.stderr)
    if args.filtered_genblastg_proteins:
        print(f"Merged annotation: {merged_gff_output}", file=sys.stderr)
        print(f"Merged proteins: {merged_fasta_output}", file=sys.stderr)
        print(f"Merged classification: {merged_tsv_output}", file=sys.stderr)
    return 0


"""Annotate full-length NLR proteins and use them to run GenBlastG."""


import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


VERSION = "HRPv2"
INTERPRO_APPLICATIONS = (
    "Pfam,SUPERFAMILY,Coils,Gene3D,SMART,PANTHER,CDD,FunFam,"
    "PRINTS,ProSiteProfiles,Phobius"
)
FULL_LENGTH_CLASSES = frozenset({"CNL", "TNL", "RNL", "NL"})
ANSI_ENABLED = True


def styled(text: str, *codes: str) -> str:
    """Apply terminal styling only for an interactive ANSI-capable stream."""
    if not ANSI_ENABLED:
        return text
    return f"\033[{';'.join(codes)}m{text}\033[0m"


def print_section(number: int, title: str, *, top_border: bool = True) -> None:
    """Display one of the six main workflow sections."""
    if top_border:
        print(styled("=" * 88, "36"), file=sys.stderr)
        print(file=sys.stderr)
    print(styled(f"[{number}/6]", "1", "35") + " "
          + styled(title, "1", "37"), file=sys.stderr)
    print(styled("=" * 88, "36"), file=sys.stderr)


def print_reused_outputs(label: str) -> None:
    """Report an explicitly supplied, previously generated analysis output."""
    print("  " + styled("✓", "1", "32") + "  "
          + styled(f"Existing {label} outputs reused", "1", "32"),
          file=sys.stderr)


def is_nb_rescue_integrated_tsv(path: Path) -> bool:
    """Return whether a supplied TSV carries the HRP HMMER/MAST integration."""
    if "with_nb-arc_rescue" in path.name.lower():
        return True
    with path.open(errors="replace") as handle:
        return any("HRP_NB_ARC" in line or "HRP-HMMER" in line for line in handle)


def read_fasta(path: Path) -> tuple[list[str], dict[str, str]]:
    """Read a FASTA while preserving headers, order and sequence wrapping."""
    order: list[str] = []
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open() as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current:
                    raise ValueError(f"empty FASTA identifier at line {number}")
                if current in records:
                    raise ValueError(f"duplicate FASTA identifier: {current}")
                order.append(current)
                records[current] = [line]
            elif current is None:
                if line.strip():
                    raise ValueError(f"sequence before first FASTA header at line {number}")
            else:
                sequence = "".join(line.split())
                if sequence:
                    records[current].append(sequence)
    if not records:
        raise ValueError(f"no FASTA records found in {path}")
    for identifier, record in records.items():
        if len(record) == 1:
            raise ValueError(f"empty FASTA sequence: {identifier}")
    return order, {identifier: "\n".join(record) + "\n"
                   for identifier, record in records.items()}


def remove_terminal_stops(input_fasta: Path, output_fasta: Path) -> int:
    """Remove terminal '*' symbols and reject, rather than erase, internal stops."""
    order, records = read_fasta(input_fasta)
    corrected = 0
    with output_fasta.open("w") as handle:
        for identifier in order:
            lines = records[identifier].rstrip("\n").splitlines()
            header = lines[0]
            sequence = "".join(lines[1:])
            cleaned = sequence.rstrip("*")
            if cleaned != sequence:
                corrected += 1
            if "*" in cleaned:
                raise ValueError(
                    "protein contains an internal stop symbol and cannot be "
                    f"submitted to InterProScan: {identifier}")
            if not cleaned:
                raise ValueError(f"protein contains only terminal stop symbols: {identifier}")
            handle.write(header + "\n")
            for start in range(0, len(cleaned), 60):
                handle.write(cleaned[start:start + 60] + "\n")
    return corrected


def select_full_length_queries(proteome: Path, classification: Path,
                               output_fasta: Path, output_table: Path) -> dict[str, int]:
    """Extract proteins with the required NB-ARC/LRR core architecture."""
    order, records = read_fasta(proteome)
    proteome_ids = set(records)
    selected_rows: dict[str, dict[str, str]] = {}
    counts: dict[str, int] = {name: 0 for name in sorted(FULL_LENGTH_CLASSES)}

    with classification.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "nlr_class", "has_NB_ARC", "LRR_hits"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"classification TSV lacks required columns: {', '.join(missing)}")
        fieldnames = list(reader.fieldnames)
        seen: set[str] = set()
        for row in reader:
            identifier = row["protein_id"]
            if identifier in seen:
                raise ValueError(f"duplicate classification identifier: {identifier}")
            seen.add(identifier)
            if identifier not in proteome_ids:
                raise ValueError(f"classified protein absent from proteome: {identifier}")
            has_nb_arc = row["has_NB_ARC"] == "1"
            try:
                lrr_hits = int(row["LRR_hits"])
            except ValueError as exc:
                raise ValueError(f"invalid LRR_hits for {identifier}") from exc
            if row["nlr_class"] in FULL_LENGTH_CLASSES and has_nb_arc and lrr_hits >= 1:
                selected_rows[identifier] = row
                counts[row["nlr_class"]] += 1

    if not selected_rows:
        raise ValueError("no full-length CNL, TNL, RNL or NL proteins were detected")

    with output_fasta.open("w") as handle:
        for identifier in order:
            if identifier in selected_rows:
                handle.write(records[identifier])
    with output_table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames,
                                lineterminator="\n")
        writer.writeheader()
        for identifier in order:
            if identifier in selected_rows:
                writer.writerow(selected_rows[identifier])
    return counts


def executable(name: str, dry_run: bool) -> str:
    resolved = shutil.which(name)
    if resolved:
        return resolved
    if dry_run:
        return name
    raise RuntimeError(f"required executable not found in PATH: {name}")


LOG_RULE = "=" * 88
LOG_SUBRULE = "-" * 88


def initialize_log(log: Path, proteome: Path, genome: Path | None,
                   annotation: Path | None, threads: int) -> None:
    """Create a fresh, consistently formatted log for one HRPv2 run."""
    with log.open("w") as handle:
        handle.write(f"{LOG_RULE}\nHRPv2 RUN LOG\n")
        handle.write(f"Started: {datetime.now(timezone.utc).isoformat()}\n")
        handle.write(f"Version: {VERSION}\n")
        handle.write(f"Threads: {threads}\n{LOG_RULE}\n")


def log_stage(log: Path, number: int, title: str) -> None:
    with log.open("a") as handle:
        handle.write(f"\n[{number}/6] {title}\n{LOG_RULE}\n")


def finalize_log(log: Path, classification: Path, proteins: Path, gff: Path) -> None:
    with log.open("a") as handle:
        handle.write(f"\nFINAL OUTPUTS\nClassification: {classification}\n")
        handle.write(f"Proteins: {proteins}\nAnnotation: {gff}\n\n")
        handle.write(f"Completed: {datetime.now(timezone.utc).isoformat()}\n")
        handle.write(f"Status: SUCCESS\n{LOG_RULE}\n")


def run_command(label: str, command: list[str], log: Path, cwd: Path | None = None,
                dry_run: bool = False, display: bool = True,
                completed_label: str | None = None) -> None:
    rendered = shlex.join(command)
    if dry_run:
        if not display:
            return
        print("  " + styled("○", "1", "33") + f"  Planned: {label}", file=sys.stderr)
        return
    if display:
        print("  " + styled("▶", "1", "33") + "  "
              + styled(f"Running: {label}...", "1", "33"), file=sys.stderr)
    with log.open("a") as handle:
        handle.write(f"\nCommand ({label}):\n{rendered}\n{LOG_SUBRULE}\n")
        handle.write(f"{label} output:\n")
        completed = subprocess.run(command, cwd=cwd, text=True,
                                   stdout=handle, stderr=subprocess.STDOUT)
        handle.write(f"{LOG_SUBRULE}\n")
    if completed.returncode:
        with log.open("a") as handle:
            handle.write(f"Status: FAILED ({label}, exit code {completed.returncode})\n")
        raise RuntimeError(f"{label} failed; see {log}")
    if display:
        print("  " + styled("✓", "1", "32") + "  "
              + styled(f"Completed: {completed_label or label}", "1", "32"),
              file=sys.stderr)


def run_interproscan_command(label: str, command: list[str], log: Path,
                             *, dry_run: bool = False,
                             completed_label: str | None = None) -> bool:
    """Run InterProScan, retrying without Phobius only when Phobius is unavailable."""
    log_offset = log.stat().st_size if log.exists() else 0
    try:
        run_command(label, command, log, dry_run=dry_run,
                    completed_label=completed_label)
        return True
    except RuntimeError:
        with log.open("rb") as handle:
            handle.seek(log_offset)
            failure_output = handle.read().decode(errors="replace")
        phobius_deactivated = (
            "Analysis Phobius" in failure_output
            and ("does not exist or is deactivated" in failure_output
                 or "is deactivated" in failure_output)
        )
        if not phobius_deactivated or "-appl" not in command:
            raise
        application_index = command.index("-appl") + 1
        applications = command[application_index].split(",")
        if not any(application.lower() == "phobius" for application in applications):
            raise
        fallback = list(command)
        fallback[application_index] = ",".join(
            application for application in applications
            if application.lower() != "phobius"
        )
        warning = (
            "Phobius is deactivated in this InterProScan installation; "
            "continuing without the optional Phobius membrane filter."
        )
        with log.open("a") as handle:
            handle.write(f"Warning: {warning}\n")
        run_command(f"{label} (without Phobius)", fallback, log,
                    dry_run=dry_run, display=False)
        if not dry_run:
            print("  " + styled("✓", "1", "32") + "  "
                  + styled(f"Completed: {completed_label or label}", "1", "32"),
                  file=sys.stderr)
        return False




@dataclass(frozen=True)
class RescueDomain:
    protein_id: str
    protein_length: int
    sequence_start: int
    sequence_end: int
    hmm_score: float
    domain_evalue: float
    hmm_coverage: float


def rescue_command(environment: str, executable_name: str,
                   *arguments: object) -> list[str]:
    command = [executable_name, *(str(value) for value in arguments)]
    if environment:
        return [shutil.which("conda") or "conda", "run", "--no-capture-output",
                "-n", environment, *command]
    return [shutil.which(executable_name) or executable_name,
            *(str(value) for value in arguments)]


def verify_rescue_environment(environment: str, *, dry_run: bool) -> None:
    """Validate the modern rescue toolchain before a long HRP run."""
    if dry_run:
        return
    required = ("hmmsearch", "hmmbuild", "meme", "meme2meme", "mast", "mafft")
    if environment:
        conda = shutil.which("conda")
        if not conda:
            raise RuntimeError("conda is required to access the rescue environment")
        check = " && ".join(f"command -v {name} >/dev/null" for name in required)
        completed = subprocess.run(
            [conda, "run", "-n", environment, "sh", "-c", check],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if completed.returncode:
            raise RuntimeError(
                f"rescue environment '{environment}' is missing or incomplete")
    else:
        missing = [name for name in required if not shutil.which(name)]
        if missing:
            raise RuntimeError("missing rescue executables: " + ", ".join(missing))


def run_captured_command(label: str, command: list[str], output: Path,
                         log: Path, *, dry_run: bool,
                         display: bool = True) -> None:
    """Run a command whose standard output is a required result file."""
    if dry_run:
        if display:
            print("  " + styled("○", "1", "33") + f"  Planned: {label}",
                  file=sys.stderr)
        return
    if display:
        print("  " + styled("▶", "1", "33") + "  "
              + styled(f"Running: {label}...", "1", "33"), file=sys.stderr)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as stdout_handle, log.open("a") as log_handle:
        log_handle.write(f"\nCommand ({label}):\n{shlex.join(command)}\n{LOG_SUBRULE}\n")
        completed = subprocess.run(command, text=True, stdout=stdout_handle,
                                   stderr=log_handle)
        log_handle.write(f"{LOG_SUBRULE}\n")
    if completed.returncode:
        raise RuntimeError(f"{label} failed; see {log}")
    if display:
        print("  " + styled("✓", "1", "32") + "  "
              + styled(f"Completed: {label}", "1", "32"), file=sys.stderr)


def fasta_sequences(path: Path) -> tuple[list[str], dict[str, str]]:
    order, records = read_fasta(path)
    return order, {
        identifier: "".join(records[identifier].splitlines()[1:]).replace("*", "")
        for identifier in order
    }


def extract_pfam_nb_arc_seeds(proteome: Path, interpro_tsv: Path,
                              output: Path, min_length: int,
                              max_length: int) -> tuple[int, set[str]]:
    """Extract non-redundant primary PF00931 regions for model training."""
    order, sequences = fasta_sequences(proteome)
    intervals: dict[str, set[tuple[int, int]]] = defaultdict(set)
    with interpro_tsv.open() as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\n").split("\t")
            if len(columns) < 8:
                raise ValueError(
                    f"InterProScan TSV line {number}: fewer than eight columns")
            if columns[3].lower() != "pfam" or columns[4].upper() != "PF00931":
                continue
            identifier = columns[0]
            if identifier not in sequences:
                raise ValueError(f"PF00931 protein absent from FASTA: {identifier}")
            start, end = int(columns[6]), int(columns[7])
            if start < 1 or end > len(sequences[identifier]) or start > end:
                raise ValueError(
                    f"invalid PF00931 coordinates for {identifier}: {start}-{end}")
            intervals[identifier].add((start, end))
    seen: set[str] = set()
    count = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for identifier in order:
            for index, (start, end) in enumerate(sorted(intervals.get(identifier, set())), 1):
                sequence = sequences[identifier][start - 1:end]
                if not min_length <= len(sequence) <= max_length:
                    continue
                if sequence in seen:
                    continue
                seen.add(sequence)
                count += 1
                handle.write(f">{identifier}|PF00931|{start}-{end}|seed{index}\n")
                for offset in range(0, len(sequence), 60):
                    handle.write(sequence[offset:offset + 60] + "\n")
    return count, set(intervals)


def extract_extended_pfam_nb_arc_training_regions(
        proteome: Path, interpro_tsv: Path, output: Path,
        min_domain_length: int, max_domain_length: int,
        c_terminal_extension: int) -> int:
    """Extract PF00931 regions and their controlled C-terminal extensions for NB-ARC training."""
    order, sequences = fasta_sequences(proteome)
    intervals: dict[str, set[tuple[int, int]]] = defaultdict(set)
    with interpro_tsv.open() as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\n").split("\t")
            if len(columns) < 8:
                raise ValueError(
                    f"InterProScan TSV line {number}: fewer than eight columns")
            if columns[3].lower() != "pfam" or columns[4].upper() != "PF00931":
                continue
            identifier = columns[0]
            if identifier not in sequences:
                raise ValueError(f"PF00931 protein absent from FASTA: {identifier}")
            start, end = int(columns[6]), int(columns[7])
            if start < 1 or end > len(sequences[identifier]) or start > end:
                raise ValueError(
                    f"invalid PF00931 coordinates for {identifier}: {start}-{end}")
            if min_domain_length <= end - start + 1 <= max_domain_length:
                intervals[identifier].add((start, end))

    seen: set[str] = set()
    count = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for identifier in order:
            sequence = sequences[identifier]
            for index, (start, end) in enumerate(
                    sorted(intervals.get(identifier, set())), 1):
                extended_end = min(len(sequence), end + c_terminal_extension)
                region = sequence[start - 1:extended_end]
                if region in seen:
                    continue
                seen.add(region)
                count += 1
                handle.write(
                    f">{identifier}|PF00931_EXTENDED|{start}-{extended_end}|seed{index}\n")
                for offset in range(0, len(region), 60):
                    handle.write(region[offset:offset + 60] + "\n")
    return count


NB_ARC_BLOCKS = (
    ("P_LOOP", 1),
    ("RNBS_A", 1),
    ("KINASE_2", 1),
    ("RNBS_BC", 2),
    ("GLPL", 1),
    ("RNBS_D", 1),
    ("MHD", 1),
)


def _first_anchor(sequence: str, pattern: str, start: int, end: int
                  ) -> tuple[int, int] | None:
    """Return the first conserved anchor inside a constrained sequence span."""
    match = re.search(pattern, sequence[max(start, 0):max(end, 0)])
    if match is None:
        return None
    offset = max(start, 0)
    return offset + match.start(), offset + match.end()


def _anchor_window(sequence: str, span: tuple[int, int], width: int) -> str:
    """Extract a fixed-width sequence window centred on an anchor."""
    centre = (span[0] + span[1]) // 2
    left = max(0, centre - width // 2)
    right = min(len(sequence), left + width)
    left = max(0, right - width)
    return sequence[left:right]


def extract_nb_arc_training_blocks(
        seeds: Path, output_dir: Path, min_width: int, max_width: int
        ) -> tuple[dict[str, Path], dict[str, int], Path]:
    """Partition PF00931 training regions using the P-loop, Kinase-2, GLPL and MHD anchors."""
    order, sequences = fasta_sequences(seeds)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {name: output_dir / f"{name}.fasta" for name, _count in NB_ARC_BLOCKS}
    records: dict[str, list[tuple[str, str]]] = defaultdict(list)
    audit = output_dir / "PF00931_anchor_audit.tsv"
    with audit.open("w", newline="") as audit_handle:
        writer = csv.writer(audit_handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "seed_id", "length", "p_loop", "kinase_2", "glpl", "mhd",
            "training_status",
        ])
        for identifier in order:
            sequence = sequences[identifier].upper()
            length = len(sequence)
            p_loop = _first_anchor(sequence, r"G[A-Z]{2,5}GK[ST]", 0,
                                   int(length * 0.45))
            kinase_2 = (_first_anchor(
                sequence, r"[LIVMFY]{2,5}D[DNE]?[LIVMFYW]",
                p_loop[1] + 10, int(length * 0.72)) if p_loop else None)
            glpl = (_first_anchor(
                sequence, r"G[LIVMF]P[LIVMF]", kinase_2[1] + 10,
                int(length * 0.92)) if kinase_2 else None)
            mhd = (_first_anchor(
                sequence, r"[LIVM]H[DNE]", glpl[1] + 5, length)
                if glpl else None)
            spans = (p_loop, kinase_2, glpl, mhd)
            valid = all(spans)
            writer.writerow([
                identifier, length,
                *(f"{span[0] + 1}-{span[1]}" if span else "." for span in spans),
                "accepted_complete_anchor_set" if valid else "excluded_missing_anchor",
            ])
            if not valid:
                continue
            assert p_loop and kinase_2 and glpl and mhd
            anchor_width = max_width
            candidates = {
                "P_LOOP": _anchor_window(sequence, p_loop, anchor_width),
                "RNBS_A": sequence[p_loop[1]:kinase_2[0]],
                "KINASE_2": _anchor_window(sequence, kinase_2, anchor_width),
                "RNBS_BC": sequence[kinase_2[1]:glpl[0]],
                "GLPL": _anchor_window(sequence, glpl, anchor_width),
                "RNBS_D": sequence[glpl[1]:mhd[0]],
                "MHD": _anchor_window(sequence, mhd, anchor_width),
            }
            for name, candidate in candidates.items():
                candidate = candidate.strip("-")
                if len(candidate) >= min_width:
                    records[name].append((identifier, candidate))

    counts: dict[str, int] = {}
    for name, path in paths.items():
        counts[name] = len(records[name])
        with path.open("w") as handle:
            for identifier, sequence in records[name]:
                handle.write(f">{identifier}|{name}\n")
                for offset in range(0, len(sequence), 60):
                    handle.write(sequence[offset:offset + 60] + "\n")
    return paths, counts, audit


def primary_pf00931_proteins(interpro_tsv: Path) -> set[str]:
    identifiers: set[str] = set()
    with interpro_tsv.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.rstrip("\n").split("\t")
            if len(columns) >= 5 and columns[3].lower() == "pfam" \
                    and columns[4].upper() == "PF00931":
                identifiers.add(columns[0])
    return identifiers


def parse_hmmer_domtblout(path: Path, max_domain_evalue: float,
                          min_hmm_coverage: float) -> dict[str, RescueDomain]:
    best: dict[str, RescueDomain] = {}
    with path.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.split(maxsplit=22)
            if len(columns) < 22:
                continue
            identifier = columns[0]
            protein_length, hmm_length = int(columns[2]), int(columns[5])
            domain_evalue, score = float(columns[12]), float(columns[13])
            hmm_from, hmm_to = int(columns[15]), int(columns[16])
            ali_from, ali_to = int(columns[17]), int(columns[18])
            coverage = (abs(hmm_to - hmm_from) + 1) / hmm_length
            candidate = RescueDomain(
                identifier, protein_length, min(ali_from, ali_to),
                max(ali_from, ali_to), score, domain_evalue, coverage)
            if domain_evalue > max_domain_evalue or coverage < min_hmm_coverage:
                continue
            previous = best.get(identifier)
            if previous is None or (candidate.domain_evalue, -candidate.hmm_score) < (
                    previous.domain_evalue, -previous.hmm_score):
                best[identifier] = candidate
    return best


def parse_mast_hit_list(path: Path, max_site_pvalue: float
                        ) -> dict[str, list[tuple[str, int, int, float]]]:
    """Parse the documented MAST protein hit-list columns."""
    hits: dict[str, list[tuple[str, int, int, float]]] = defaultdict(list)
    with path.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            columns = raw.split()
            if len(columns) < 8:
                continue
            try:
                start, end = int(columns[-4]), int(columns[-3])
                pvalue = float(columns[-1])
            except ValueError:
                continue
            if pvalue <= max_site_pvalue:
                hits[columns[0]].append(
                    (columns[1].lstrip("+-"), min(start, end), max(start, end), pvalue))
    return hits


def motif_order_from_seeds(
        seed_hits: dict[str, list[tuple[str, int, int, float]]],
        seed_count: int, min_prevalence: float) -> list[str]:
    """Learn motif order from prevalent sites and their median seed position."""
    positions: dict[str, list[float]] = defaultdict(list)
    proteins_per_motif: dict[str, set[str]] = defaultdict(set)
    for identifier, hits in seed_hits.items():
        extent = max((end for _motif, _start, end, _pvalue in hits), default=0)
        if not extent:
            continue
        for motif, start, end, _pvalue in hits:
            positions[motif].append(((start + end) / 2) / extent)
            proteins_per_motif[motif].add(identifier)

    def median(values: list[float]) -> float:
        ordered = sorted(values)
        middle = len(ordered) // 2
        return ordered[middle] if len(ordered) % 2 else (
            ordered[middle - 1] + ordered[middle]) / 2

    retained = [motif for motif in positions
                if len(proteins_per_motif[motif]) / max(seed_count, 1) >= min_prevalence]
    return sorted(retained, key=lambda motif: (median(positions[motif]), motif))


def ordered_motif_support(
        hits: list[tuple[str, int, int, float]], domain: RescueDomain,
        expected_order: list[str], motif_labels: dict[str, str] | None = None,
        n_terminal_flank: int = 25, c_terminal_flank: int = 200,
        max_cluster_span: int = 300,
        ) -> tuple[int, str, set[str], int | None, int | None]:
    """Identify an ordered and spatially coherent MAST motif cluster associated with a HMMER hit."""
    ranks = {motif: index for index, motif in enumerate(expected_order)}
    window_start = max(1, domain.sequence_start - n_terminal_flank)
    window_end = min(domain.protein_length,
                     domain.sequence_end + c_terminal_flank)
    contained = sorted(
        (hit for hit in hits if hit[0] in ranks
         and hit[1] >= window_start and hit[2] <= window_end),
        key=lambda hit: (hit[1], hit[2], hit[3]))
    paths: list[list[tuple[str, int, int, float]]] = []
    for hit in contained:
        prefix: list[tuple[str, int, int, float]] = []
        for path in paths:
            prospective_span = hit[2] - path[0][1] + 1
            if (ranks[path[-1][0]] < ranks[hit[0]]
                    and prospective_span <= max_cluster_span
                    and len(path) > len(prefix)):
                prefix = path
        paths.append([*prefix, hit])
    overlapping = [
        path for path in paths
        if path[0][1] <= domain.sequence_end
        and path[-1][2] >= domain.sequence_start
        and path[-1][2] - path[0][1] + 1 <= max_cluster_span
    ]
    best = max(overlapping, key=lambda path: (len(path), -(
        path[-1][2] - path[0][1] + 1)), default=[])
    motif_ids = [hit[0] for hit in best]
    labels = {motif_labels.get(motif, motif) for motif in motif_ids} \
        if motif_labels else set(motif_ids)
    cluster_start = best[0][1] if best else None
    cluster_end = best[-1][2] if best else None
    return (len(set(motif_ids)), ",".join(motif_ids) or ".", labels,
            cluster_start, cluster_end)


def append_rescue_hits(interpro_tsv: Path, output_tsv: Path,
                       proteins: dict[str, int],
                       accepted: dict[str, RescueDomain]) -> None:
    """Append accepted calls as auditable synthetic InterProScan-style rows."""
    shutil.copyfile(interpro_tsv, output_tsv)
    with output_tsv.open("a") as handle:
        for identifier in sorted(accepted):
            hit = accepted[identifier]
            columns = [
                identifier, "-", str(proteins[identifier]), "HRP-HMMER",
                "HRP_NB_ARC", "Species-specific NB-ARC rescue",
                str(hit.sequence_start), str(hit.sequence_end),
                f"{hit.domain_evalue:.3g}", "T",
                datetime.now(timezone.utc).strftime("%d-%m-%Y"),
                "IPR002182", "NB-ARC", "-", "-",
            ]
            handle.write("\t".join(columns) + "\n")


def screen_nb_arc_rescue(
        proteome: Path, interpro_tsv: Path, output_tsv: Path,
        model_dir: Path, scan_dir: Path, args: argparse.Namespace, log: Path,
        *, build_models: bool, dry_run: bool) -> dict[str, object]:
    """Build/reuse species models, scan proteins and integrate accepted calls."""
    running_label = "Specie-specific NB domain detection by HMMER & MEME\\MAST"
    completed_label = "Specie-specific NB domain detection by HMMER & MEME\\MAST"
    if dry_run:
        print("  " + styled("○", "1", "33") + f"  Planned: {running_label}",
              file=sys.stderr)
    else:
        print("  " + styled("▶", "1", "33") + "  "
              + styled(f"Running: {running_label}...", "1", "33"),
              file=sys.stderr)
    if not dry_run:
        model_dir.mkdir(parents=True, exist_ok=True)
        scan_dir.mkdir(parents=True, exist_ok=True)
    seeds = model_dir / "PF00931_validated_domains.fasta"
    extended_seeds = model_dir / "PF00931_extended_training_regions.fasta"
    alignment = model_dir / "PF00931_species_alignment.fasta"
    hmm_model = model_dir / "PF00931_species-specific.hmm"
    meme_dir = model_dir / "MEME_NB-ARC_motifs"
    motif_file = meme_dir / "meme.txt"
    seed_mast = model_dir / "PF00931_seed_MAST_hits.tsv"
    motif_order_file = model_dir / "NB-ARC_motif_order.txt"

    if build_models:
        if dry_run:
            seed_count = 0
        else:
            seed_count, _primary_ids = extract_pfam_nb_arc_seeds(
                proteome, interpro_tsv, seeds, args.seed_min_length,
                args.seed_max_length)
            if seed_count < args.rescue_min_seeds:
                raise RuntimeError(
                    f"only {seed_count} non-redundant PF00931 seeds; "
                    f"at least {args.rescue_min_seeds} are required")
            extended_seed_count = extract_extended_pfam_nb_arc_training_regions(
                proteome, interpro_tsv, extended_seeds,
                args.seed_min_length, args.seed_max_length,
                args.nb_arc_c_terminal_extension)
            if extended_seed_count < args.rescue_min_seeds:
                raise RuntimeError(
                    f"only {extended_seed_count} extended PF00931 training "
                    "regions were retained")
        run_captured_command(
            "MAFFT alignment of primary PF00931 domains",
            rescue_command(args.rescue_env, "mafft", "--auto", "--thread",
                           args.threads, seeds), alignment, log,
            dry_run=dry_run, display=False)
        run_command(
            "Species-specific NB-ARC HMM construction",
            rescue_command(args.rescue_env, "hmmbuild", "--amino", "--cpu",
                           args.threads, hmm_model, alignment),
            log, dry_run=dry_run, display=False)
        block_fastas: dict[str, Path] = {}
        block_counts: dict[str, int] = {}
        anchor_audit = model_dir / "PF00931_training_blocks" / "PF00931_anchor_audit.tsv"
        if not dry_run:
            meme_dir.mkdir(parents=True, exist_ok=True)
            block_fastas, block_counts, anchor_audit = extract_nb_arc_training_blocks(
                extended_seeds, model_dir / "PF00931_training_blocks",
                args.meme_min_width, args.meme_max_width)
            insufficient = [name for name, count in block_counts.items()
                            if count < args.rescue_min_seeds]
            if insufficient:
                details = ", ".join(
                    f"{name}={block_counts[name]}" for name in insufficient)
                raise RuntimeError(
                    "insufficient complete PF00931 training sequences for "
                    f"NB-ARC blocks: {details}; see {anchor_audit}")
        motif_sources: list[Path] = []
        for block_name, motif_number in NB_ARC_BLOCKS:
            block_dir = meme_dir / block_name
            block_fasta = (block_fastas.get(block_name)
                           or model_dir / "PF00931_training_blocks" / f"{block_name}.fasta")
            run_command(
                f"MEME training of {block_name}",
                rescue_command(
                    args.rescue_env, "meme", "-protein",
                    "-mod", "oops" if block_name in {
                        "P_LOOP", "KINASE_2", "GLPL", "MHD"} else "zoops",
                    "-nmotifs", motif_number,
                    "-minw", args.meme_min_width,
                    "-maxw", args.meme_max_width,
                    "-oc", block_dir, block_fasta),
                log, dry_run=dry_run, display=False)
            motif_sources.append(block_dir / "meme.txt")
        run_captured_command(
            "Assembly of eight species-specific NB-ARC motif profiles",
            rescue_command(args.rescue_env, "meme2meme", *motif_sources),
            motif_file, log, dry_run=dry_run, display=False)
        if not dry_run:
            motif_total = sum(
                1 for line in motif_file.read_text(errors="replace").splitlines()
                if line.startswith("MOTIF "))
            if motif_total != 8:
                raise RuntimeError(
                    f"expected eight species-specific NB-ARC profiles, found "
                    f"{motif_total} in {motif_file}")
            summary_file = model_dir / "NB-ARC_training_block_summary.tsv"
            with summary_file.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["training_block", "requested_motifs", "seed_sequences"])
                for block_name, motif_number in NB_ARC_BLOCKS:
                    writer.writerow([block_name, motif_number,
                                     block_counts.get(block_name, 0)])
        run_captured_command(
            "MAST calibration on extended primary PF00931 regions",
            rescue_command(args.rescue_env, "mast", "-hit_list",
                           motif_file, extended_seeds),
            seed_mast, log, dry_run=dry_run, display=False)
        if not dry_run:
            seed_hits = parse_mast_hit_list(seed_mast, float("inf"))
            motif_order = motif_order_from_seeds(
                seed_hits, min(block_counts.values()), 0.0)
            if len(motif_order) != 8:
                raise RuntimeError(
                    "the eight ordered NB-ARC profiles could not all be calibrated; "
                    f"found {len(motif_order)} prevalent profiles")
            biological_order = (
                "P_LOOP", "RNBS_A", "KINASE_2", "RNBS_B",
                "RNBS_C", "GLPL", "RNBS_D", "MHD")
            motif_order_file.write_text("\n".join(
                f"{motif}\t{label}"
                for motif, label in zip(motif_order, biological_order)) + "\n")
            calibration = model_dir / "NB-ARC_motif_calibration.tsv"
            with calibration.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    "mast_motif_id", "biological_motif", "seed_hits",
                    "extended_seed_prevalence", "used_in_model",
                ])
                for motif, label in zip(motif_order, biological_order):
                    carrying = sum(
                        motif in {hit[0] for hit in hits}
                        for hits in seed_hits.values())
                    writer.writerow([
                        motif, label, carrying,
                        f"{carrying / max(extended_seed_count, 1):.6f}", 1,
                    ])
    else:
        seed_count = None
        if not dry_run:
            for required in (hmm_model, motif_file, motif_order_file):
                if not required.is_file():
                    raise RuntimeError(f"NB-ARC rescue model is missing: {required}")

    domtblout = scan_dir / "NB-ARC_hmmsearch.domtblout"
    mast_hits_file = scan_dir / "NB-ARC_MAST_hits.tsv"
    run_command(
        "Species-specific HMMER NB-ARC screening",
        rescue_command(args.rescue_env, "hmmsearch", "--cpu", args.threads,
                       "--domtblout", domtblout, hmm_model, proteome),
        log, dry_run=dry_run, display=False)
    run_captured_command(
        "MAST NB-ARC motif screening",
        rescue_command(args.rescue_env, "mast", "-hit_list",
                       motif_file, proteome),
        mast_hits_file, log, dry_run=dry_run, display=False)

    audit = scan_dir / "NB-ARC_rescue_evidence.tsv"
    accepted: dict[str, RescueDomain] = {}
    hmm_candidates = 0
    if not dry_run:
        proteins = fasta_lengths(proteome)
        domains = parse_hmmer_domtblout(
            domtblout, float("inf"), 0.0)
        mast_hits = parse_mast_hit_list(mast_hits_file, float("inf"))
        motif_rows = [line.split("\t", 1)
                      for line in motif_order_file.read_text().splitlines()
                      if line.strip()]
        motif_order = [row[0] for row in motif_rows]
        motif_labels = {row[0]: row[1] if len(row) > 1 else row[0]
                        for row in motif_rows}
        primary_ids = primary_pf00931_proteins(interpro_tsv)
        with audit.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "protein_id", "hmm_start", "hmm_end", "hmm_score",
                "domain_evalue", "hmm_coverage", "ordered_motif_count",
                "ordered_motifs", "functional_motifs", "motif_cluster_start",
                "motif_cluster_end", "motif_cluster_span", "hmm_motif_overlap",
                "primary_PF00931", "rescue_status",
            ])
            for identifier in sorted(domains):
                domain = domains[identifier]
                is_primary = identifier in primary_ids
                motif_count, ordered, functional, cluster_start, cluster_end = \
                    ordered_motif_support(
                    mast_hits.get(identifier, []), domain, motif_order,
                    motif_labels, args.mast_n_terminal_flank,
                    args.mast_c_terminal_flank, args.mast_max_cluster_span)
                motif_valid = motif_count >= args.mast_min_motifs
                if is_primary:
                    status = "primary_PF00931_not_rescued"
                elif motif_valid:
                    status = "HMMER_MAST_accepted"
                    accepted[identifier] = domain
                else:
                    status = "HMMER_only_rejected"
                if not is_primary:
                    hmm_candidates += 1
                writer.writerow([
                    identifier, domain.sequence_start, domain.sequence_end,
                    f"{domain.hmm_score:.3f}", f"{domain.domain_evalue:.3g}",
                    f"{domain.hmm_coverage:.3f}", motif_count, ordered,
                    ",".join(sorted(functional)) or ".",
                    cluster_start if cluster_start is not None else ".",
                    cluster_end if cluster_end is not None else ".",
                    (cluster_end - cluster_start + 1
                     if cluster_start is not None and cluster_end is not None else "."),
                    int(cluster_start is not None), int(is_primary), status,
                ])
        append_rescue_hits(interpro_tsv, output_tsv, proteins, accepted)
    result = {
        "seed_count": seed_count,
        "hmm_candidates_without_primary_PF00931": hmm_candidates,
        "accepted_novel_NB_ARC_hits": len(accepted),
        "hmm_model": str(hmm_model),
        "meme_motifs": str(motif_file),
        "audit": str(audit),
        "integrated_interpro_tsv": str(output_tsv),
    }
    if not dry_run:
        print("  " + styled("✓", "1", "32") + "  "
              + styled(f"Completed: {completed_label}", "1", "32"),
              file=sys.stderr)
    return result


def build_genblastg_command(executable_path: str, query_fasta: Path, genome: Path,
                            output_name: str, args: argparse.Namespace) -> list[str]:
    """Build the GenBlastG command while preserving its default prediction parameters."""
    return [
        executable_path, "-q", str(query_fasta), "-t", str(genome),
        "-gff", "-pro", "-o", output_name,
    ]


def stage_genblastg_input(source: Path, workdir: Path, local_name: str) -> Path:
    source = source.resolve()
    destination = workdir / local_name
    if destination.is_symlink():
        destination.unlink()
    elif destination.exists():
        raise RuntimeError(
            f"cannot stage GenBlastG input because {destination} already exists")
    destination.symlink_to(source)
    return Path(local_name)


def find_single_genblastg_output(directory: Path, output_name: str, suffix: str) -> Path:
    """Locate the single GenBlastG result carrying the requested suffix."""
    matches = sorted(path for path in directory.glob(f"{output_name}*{suffix}")
                     if path.is_file())
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one GenBlastG {suffix} file matching {output_name}*, found {len(matches)}")
    return matches[0]


def write_predicted_only_catalog(gff: Path, proteins: Path, classification: Path,
                                 prefix: Path) -> None:
    """Create the three final outputs when no annotation GFF3 is supplied."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(gff, Path(str(prefix) + ".gff3"))
    shutil.copyfile(proteins, Path(str(prefix) + ".fasta"))
    order, _ = read_fasta(proteins)
    with classification.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "protein_id" not in reader.fieldnames:
            raise ValueError("predicted-model classification lacks protein_id")
        fields = list(reader.fieldnames)
        rows = {row["protein_id"]: row for row in reader}
    data_fields = [field for field in fields if field != "protein_id"]
    with Path(str(prefix) + ".tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n",
                                fieldnames=["gene_id", "protein_id", "source"] + data_fields)
        writer.writeheader()
        for identifier in order:
            if identifier not in rows:
                raise ValueError(f"filtered model lacks classification: {identifier}")
            row = rows[identifier]
            writer.writerow({"gene_id": identifier, "protein_id": identifier,
                             "source": "GenBlastG",
                             **{key: row[key] for key in data_fields}})


def summarize_final_catalog(path: Path) -> tuple[int, int, int]:
    """Count total, full-length and partial NLRs in the final catalogue."""
    full_length = 0
    partial = 0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "nlr_class" not in reader.fieldnames:
            raise ValueError(f"final classification TSV lacks nlr_class: {path}")
        for row in reader:
            nlr_class = row["nlr_class"]
            if nlr_class == "non_NLR":
                continue
            if nlr_class in FULL_LENGTH_CLASSES:
                full_length += 1
            else:
                partial += 1
    return full_length + partial, full_length, partial


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=("Welcome to The full-length Homology-based R-gene Prediction "
                     "(HRPv2)"))
    ap.add_argument("--version", action="version", version="HRPv2")
    ap.add_argument("--proteome", type=Path,
                    help="Protein FASTA for the species under analysis")
    ap.add_argument("--genome", type=Path,
                    help="Genome FASTA from the same species")
    ap.add_argument("--annotation-gff", type=Path,
                    help="Optional annotation GFF3 corresponding to the proteome")
    ap.add_argument(
        "--workdir", type=Path,
        help=("Optional output root for the six numbered HRPv2 directories; "
              "when omitted, they are created in the current directory"),
    )
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--interproscan-bin", default="interproscan.sh")
    ap.add_argument("--classifier", type=Path,
                    default=Path(__file__))
    ap.add_argument("--filter-script", type=Path,
                    default=Path(__file__),
                    help="Module used to filter overlapping GenBlastG gene models")
    ap.add_argument("--recovery-script", type=Path,
                    default=Path(__file__),
                    help="Module used to recover annotated partial NB-LRR loci")
    ap.add_argument("--interpro-tsv", type=Path,
                    help="Reuse an existing InterProScan TSV instead of running InterProScan")
    ap.add_argument("--genblastg-gff", type=Path,
                    help="Reuse an existing GenBlastG GFF and skip step-3 execution")
    ap.add_argument("--genblastg-proteins", type=Path,
                    help="Protein FASTA/.pro paired with --genblastg-gff")
    ap.add_argument("--predicted-interpro-tsv", type=Path,
                    help="Reuse the InterProScan TSV for the supplied GenBlastG proteins")
    ap.add_argument("--genblastg-bin", default="run_genblastG")
    ap.add_argument(
        "--rescue-env", default="hrp_rescue_env",
        help=("Conda environment containing HMMER, MEME Suite and MAFFT; "
              "use an empty value only when all tools are already in PATH"),
    )
    ap.add_argument("--skip-nb-rescue", action="store_true",
                    help="Disable the species-specific HMMER plus MEME/MAST rescue")
    ap.add_argument("--rescue-min-seeds", type=int, default=10,
                    help="Minimum number of non-redundant PF00931 domains")
    ap.add_argument("--seed-min-length", type=int, default=130,
                    help="Minimum length of a complete PF00931 training region")
    ap.add_argument("--seed-max-length", type=int, default=350,
                    help="Maximum length of a complete PF00931 training region")
    ap.add_argument(
        "--nb-arc-c-terminal-extension", type=int, default=200,
        help=("C-terminal residues added only for GLPL/MHD motif training "
              "(default: 200 aa)"),
    )
    ap.add_argument("--meme-motifs", type=int, default=8,
                    help=("Number of block-constrained NB-ARC profiles; the "
                          "Meyers architecture requires exactly eight"))
    ap.add_argument("--meme-min-width", type=int, default=6)
    ap.add_argument("--meme-max-width", type=int, default=12,
                    help="Maximum MEME motif width (default: 12 aa)")
    ap.add_argument("--mast-min-motifs", type=int, default=3,
                    help=("Minimum distinct, ordered motifs in a local cluster "
                          "spatially compatible with the HMMER hit"))
    ap.add_argument(
        "--mast-n-terminal-flank", type=int, default=25,
        help="Allowed MAST-cluster extension upstream of the HMMER hit",
    )
    ap.add_argument(
        "--mast-c-terminal-flank", type=int, default=200,
        help="Allowed MAST-cluster extension downstream of the HMMER hit",
    )
    ap.add_argument(
        "--mast-max-cluster-span", type=int, default=300,
        help="Maximum amino-acid span of an ordered MAST motif cluster",
    )
    ap.add_argument("--skip-genblastg", action="store_true",
                    help="Stop after creating Annotated_full-length_NB-LRRs.fasta")
    ap.add_argument("--dry-run", action="store_true",
                    help="Validate inputs and display external commands without running them")
    return ap


def _pipeline_main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    interactive_inputs = args.proteome is None or (not args.skip_genblastg and args.genome is None)
    if args.proteome is None:
        value = input("Proteome FASTA file: ").strip().strip("'\"")
        if not value:
            raise ValueError("the proteome path cannot be empty")
        args.proteome = Path(value).expanduser()
    if not args.skip_genblastg and args.genome is None:
        value = input("Genome FASTA file: ").strip().strip("'\"")
        if not value:
            raise ValueError("the genome path cannot be empty")
        args.genome = Path(value).expanduser()
    if interactive_inputs and not args.skip_genblastg and args.annotation_gff is None:
        value = input("Genome annotation GFF3 file (optional; press Enter to skip): ").strip().strip("'\"")
        if value:
            args.annotation_gff = Path(value).expanduser()

    proteome = args.proteome.resolve()
    genome = args.genome.resolve() if args.genome else None
    classifier = args.classifier.resolve()
    filter_script = args.filter_script.resolve()
    recovery_script = args.recovery_script.resolve()
    annotation_gff = args.annotation_gff.resolve() if args.annotation_gff else None
    if not proteome.is_file():
        raise ValueError(f"proteome not found: {proteome}")
    proteome_order, _proteome_records = read_fasta(proteome)
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    if args.rescue_min_seeds < 2:
        raise ValueError("--rescue-min-seeds must be at least 2")
    if not 1 <= args.seed_min_length <= args.seed_max_length:
        raise ValueError("PF00931 seed-length limits are invalid")
    if args.nb_arc_c_terminal_extension < 0:
        raise ValueError("--nb-arc-c-terminal-extension cannot be negative")
    if args.meme_motifs != 8:
        raise ValueError("--meme-motifs must be 8 for the block-constrained model")
    if args.mast_min_motifs < 3:
        raise ValueError("--mast-min-motifs must be at least 3")
    if args.mast_n_terminal_flank < 0 or args.mast_c_terminal_flank < 0:
        raise ValueError("MAST integration flanks cannot be negative")
    if args.mast_max_cluster_span < args.meme_max_width:
        raise ValueError("--mast-max-cluster-span is too short")
    if not 1 <= args.meme_min_width <= args.meme_max_width:
        raise ValueError("MEME motif widths are invalid")
    if not args.skip_nb_rescue:
        verify_rescue_environment(args.rescue_env, dry_run=args.dry_run)
    reuse_genblastg = args.genblastg_gff is not None or args.genblastg_proteins is not None
    if reuse_genblastg and not (args.genblastg_gff and args.genblastg_proteins):
        raise ValueError(
            "--genblastg-gff and --genblastg-proteins must be supplied together")
    if args.predicted_interpro_tsv and not reuse_genblastg:
        raise ValueError(
            "--predicted-interpro-tsv requires reusable GenBlastG GFF and proteins")
    reusable_gff = args.genblastg_gff.resolve() if args.genblastg_gff else None
    reusable_proteins = (
        args.genblastg_proteins.resolve() if args.genblastg_proteins else None)
    reusable_predicted_ips = (
        args.predicted_interpro_tsv.resolve()
        if args.predicted_interpro_tsv else None)
    for label, path in (
            ("GenBlastG GFF", reusable_gff),
            ("GenBlastG proteins", reusable_proteins),
            ("predicted-model InterProScan TSV", reusable_predicted_ips)):
        if path is not None and not path.is_file():
            raise ValueError(f"{label} not found: {path}")
    genome_record_count: int | None = None
    if not args.skip_genblastg:
        if genome is None or not genome.is_file():
            raise ValueError(f"genome not found: {genome}")
        genome_order, _genome_records = read_fasta(genome)
        genome_record_count = len(genome_order)
    if not classifier.is_file():
        raise ValueError(f"classifier not found: {classifier}")
    if not filter_script.is_file():
        raise ValueError(f"filtering module not found: {filter_script}")
    if annotation_gff is not None:
        if not annotation_gff.is_file():
            raise ValueError(f"annotation GFF3 not found: {annotation_gff}")
        if not recovery_script.is_file():
            raise ValueError(f"partial-locus recovery module not found: {recovery_script}")

    workdir = args.workdir.resolve() if args.workdir else Path.cwd().resolve()
    input_dir = workdir / "01_preflight_validation_inputs"
    interpro_dir = workdir / "02_full-length_NB-LRR_annotation"
    query_dir = interpro_dir
    rescue_model_dir = interpro_dir / "NB-ARC_species_rescue"
    genblastg_dir = workdir / "03_gene_model_prediction"
    predicted_interpro_dir = workdir / "04_NB-LRR_model_annotation"
    filtered_dir = workdir / "05_filtered_NB-LRR_models"
    merged_dir = workdir / "06_merged_NB-LRR_genes"
    log = workdir / "HRPv2.log"

    print(styled("#" * 88, "1", "36"), file=sys.stderr)
    print(file=sys.stderr)
    print(styled("The full-length Homology-based R-gene Prediction (HRPv2)",
                 "1", "36"), file=sys.stderr)
    print(styled("Prediction of NB-LRR resistance genes based on full-length sequence homology",
                 "1", "37"), file=sys.stderr)
    print(file=sys.stderr)
    print(styled("#" * 88, "1", "36"), file=sys.stderr)
    print(file=sys.stderr)
    # Step 1/6 - Validate inputs and write preparation reports to 01_preflight_validation_inputs.
    print_section(1, "Preflight validation and preparation of input files",
                  top_border=False)
    print("  " + styled("▶", "1", "33") + "  "
          + styled("Running: Input file validation...", "1", "33"), file=sys.stderr)
    print(f"  Proteome:        {proteome.name}", file=sys.stderr)
    if genome is not None:
        print(f"  Genome:          {genome.name}", file=sys.stderr)
    if annotation_gff is not None:
        print(f"  Annotation:      {annotation_gff.name}", file=sys.stderr)
    if not args.dry_run:
        for directory in (input_dir, interpro_dir, query_dir, genblastg_dir,
                          predicted_interpro_dir, filtered_dir, merged_dir):
            directory.mkdir(parents=True, exist_ok=True)
        initialize_log(log, proteome, genome, annotation_gff, args.threads)
        log_stage(log, 1, "Preflight validation and preparation of input files")
        with log.open("a") as handle:
            handle.write("Input files:\n")
            handle.write(f"  Proteome: {proteome}\n")
            handle.write(f"  Genome: {genome if genome else 'not supplied'}\n")
            handle.write(
                f"  Annotation: {annotation_gff if annotation_gff else 'not supplied'}\n")
            handle.write("Input validation details:\n")
            handle.write(f"  Proteome records: {len(proteome_order)}\n")
            if genome_record_count is not None:
                handle.write(f"  Genome records: {genome_record_count}\n")

    if annotation_gff is not None:
        if args.dry_run:
            print("  " + styled("○", "1", "33")
                  + "  Planned: annotation GFF validation", file=sys.stderr)
        else:
            gff_audit = input_dir / "annotation_gff_identifier_audit.tsv"
            gff_anomalies = input_dir / "annotation_gff_anomalies.tsv"
            gff_counts = preflight_annotation_gff(
                annotation_gff, proteome_order, gff_audit, gff_anomalies)
            with log.open("a") as handle:
                handle.write(
                    "Annotation validation details:\n"
                    f"  Valid GFF features: {gff_counts['features']}\n"
                    f"  Proteome identifiers matched unambiguously: "
                    f"{gff_counts['matched']}\n"
                    f"  Proteome identifiers unmatched: {gff_counts['unmatched']}\n"
                    f"  Proteome identifiers ambiguous: {gff_counts['ambiguous']}\n"
                    f"  GFF anomalies recorded: {gff_counts['anomalies']}\n"
                    f"  Identifier audit: {gff_audit}\n"
                    f"  Anomaly report: {gff_anomalies}\n")
    cleaned_proteome = input_dir / "proteome_without_terminal_stops.fasta"
    if args.dry_run:
        terminal_stops_removed = 0
        print("  " + styled("○", "1", "33") + "  Planned: protein input cleaning",
              file=sys.stderr)
    else:
        terminal_stops_removed = remove_terminal_stops(proteome, cleaned_proteome)
        with log.open("a") as handle:
            handle.write("Input preparation results:\n")
            handle.write(f"  Terminal stop symbols removed: {terminal_stops_removed}\n")
            handle.write(f"  Cleaned proteome: {cleaned_proteome}\n")
            handle.write("  Completed: Input file validation\n")
            handle.write(f"{LOG_RULE}\n")
    if not args.dry_run:
        print("  " + styled("✓", "1", "32") + "  "
              + styled("Completed: Input file validation", "1", "32"),
              file=sys.stderr)
        print(styled("=" * 88, "36"), file=sys.stderr)
        print(file=sys.stderr)

    interpro_prefix = interpro_dir / "proteome_NB-LRR"
    classification = query_dir / "proteome_NB-LRR.classification.tsv"
    query_fasta = query_dir / "annotated_full-length_NB-LRRs.fasta"
    query_table = query_dir / "annotated_full-length_NB-LRRs.tsv"

    # Step 2/6 - Annotate proteins and write the full-length NB-LRR query set to 02_full-length_NB-LRR_annotation.
    print_section(2, "Full-length NB-LRR annotation using a protein motif/domain-based search (PDS)",
                  top_border=False)
    if not args.dry_run:
        log_stage(log, 2, "Full-length NB-LRR annotation using PDS")
    phobius_used: bool | None = None
    if args.interpro_tsv:
        interpro_tsv = args.interpro_tsv.resolve()
        if not interpro_tsv.is_file():
            raise ValueError(f"InterProScan TSV not found: {interpro_tsv}")
        if args.dry_run:
            print("  " + styled("○", "1", "33")
                  + "  Planned: InterProScan protein evidence loading",
                  file=sys.stderr)
        if not args.dry_run:
            with log.open("a") as handle:
                handle.write(f"Existing InterProScan TSV reused: {interpro_tsv}\n")
            print_reused_outputs("InterProScan")
            if (args.skip_nb_rescue
                    and is_nb_rescue_integrated_tsv(interpro_tsv)):
                print_reused_outputs("HMMER\\MAST")
    else:
        interproscan = executable(args.interproscan_bin, args.dry_run)
        interpro_tsv = Path(str(interpro_prefix) + ".tsv")
        command = [interproscan, "-i", str(cleaned_proteome), "-appl", INTERPRO_APPLICATIONS,
                   "-f", "TSV,GFF3", "-b", str(interpro_prefix),
                   "-cpu", str(args.threads), "-dp"]
        phobius_used = run_interproscan_command(
            "Protein domain annotation by InterProScan", command, log,
            dry_run=args.dry_run,
            completed_label="Protein domain annotation by InterProScan")
        if not args.dry_run and not interpro_tsv.is_file():
            raise RuntimeError(f"InterProScan did not produce {interpro_tsv}")

    proteome_rescue_summary: dict[str, object] | None = None
    classification_tsv = interpro_tsv
    if not args.skip_nb_rescue:
        integrated_tsv = interpro_dir / "proteome_NB-LRR.with_NB-ARC_rescue.tsv"
        proteome_rescue_summary = screen_nb_arc_rescue(
            cleaned_proteome, interpro_tsv, integrated_tsv,
            rescue_model_dir, rescue_model_dir / "proteome_scan", args, log,
            build_models=True, dry_run=args.dry_run)
        classification_tsv = integrated_tsv

    classifier_command = [sys.executable, str(classifier), "__classify__",
                          "--proteins", str(cleaned_proteome),
                          "--interpro-tsv", str(classification_tsv), "--output", str(classification),
                          "--quiet", "--summary-only", "--no-authorship"]
    run_command("NB-LRR protein classification", classifier_command, log,
                dry_run=args.dry_run)

    if args.dry_run:
        print("  " + styled("○", "1", "33")
              + "  Planned: full-length NB-LRR query selection", file=sys.stderr)
        counts: dict[str, int] = {}
    else:
        counts = select_full_length_queries(cleaned_proteome, classification,
                                            query_fasta, query_table)
        total = sum(counts.values())
        print(file=sys.stderr)
        print(f"Full-length NB-LRR queries: {total}", file=sys.stderr)
        for name in ("CNL", "TNL", "RNL", "NL"):
            print(f"  {name}: {counts[name]}", file=sys.stderr)

    # Step 3/6 - Predict genomic NB-LRR gene models and save GenBlastG outputs in 03_gene_model_prediction.
    print_section(3, "Prediction of NB-LRR gene models")
    if not args.dry_run:
        log_stage(log, 3, "Prediction of NB-LRR gene models")
    genblastg_command: list[str] | None = None
    predicted_gff: Path | None = None
    predicted_proteins: Path | None = None
    predicted_classification: Path | None = None
    filtered_prefix: Path | None = None
    if args.skip_genblastg:
        print("GenBlastG skipped by --skip-genblastg", file=sys.stderr)
    else:
        output_name = "genblastg_Annotated_full-length_NB-LRRs"
        if reuse_genblastg:
            predicted_gff = reusable_gff
            raw_predicted_proteins = reusable_proteins
            if not args.dry_run:
                with log.open("a") as handle:
                    handle.write(f"Existing GenBlastG GFF reused: {predicted_gff}\n")
                    handle.write(
                        f"Existing GenBlastG proteins reused: {raw_predicted_proteins}\n")
                print_reused_outputs("GenBlastG")
        else:
            genblastg = executable(args.genblastg_bin, args.dry_run)
            assert genome is not None
            local_query = Path("full_length_queries.fasta")
            local_genome = Path("target_genome.fasta")
            staged_links: list[Path] = []
            if not args.dry_run:
                local_query = stage_genblastg_input(
                    query_fasta, genblastg_dir, local_query.name)
                staged_links.append(genblastg_dir / local_query)
                local_genome = stage_genblastg_input(
                    genome, genblastg_dir, local_genome.name)
                staged_links.append(genblastg_dir / local_genome)
            genblastg_command = build_genblastg_command(
                genblastg, local_query, local_genome, output_name, args)
            try:
                run_command("Gene model prediction by GenBlastG", genblastg_command,
                            log, cwd=genblastg_dir, dry_run=args.dry_run)
            finally:
                for staged_link in staged_links:
                    if staged_link.is_symlink():
                        staged_link.unlink()

        predicted_interpro_prefix = predicted_interpro_dir / "predicted_gene_models"
        predicted_interpro_tsv = Path(str(predicted_interpro_prefix) + ".tsv")
        predicted_classification = predicted_interpro_dir / "predicted_gene_models.classification.tsv"
        filtered_prefix = filtered_dir / "filtered_NB-LRRs"
        if reuse_genblastg:
            assert predicted_gff is not None and raw_predicted_proteins is not None
        elif args.dry_run:
            predicted_gff = genblastg_dir / f"{output_name}*.gff"
            raw_predicted_proteins = genblastg_dir / f"{output_name}*.pro"
        else:
            predicted_gff = find_single_genblastg_output(genblastg_dir, output_name, ".gff")
            raw_predicted_proteins = find_single_genblastg_output(
                genblastg_dir, output_name, ".pro")

        predicted_proteins = (
            predicted_interpro_dir /
            "predicted_gene_models_without_terminal_stops.fasta"
        )
        if not args.dry_run:
            predicted_stops_removed = remove_terminal_stops(
                raw_predicted_proteins, predicted_proteins)
            with log.open("a") as handle:
                handle.write(
                    "GenBlastG terminal stop symbols removed: "
                    f"{predicted_stops_removed}\n")

        # Step 4/6 - Annotate and classify predicted proteins in 04_NB-LRR_model_annotation.
        print_section(4, "Annotation and classification of gene models")
        if not args.dry_run:
            log_stage(log, 4, "Annotation and classification of NB-LRR gene models")
        if reusable_predicted_ips is not None:
            predicted_interpro_tsv = reusable_predicted_ips
            predicted_phobius_used = None
            if args.dry_run:
                print("  " + styled("○", "1", "33")
                      + "  Planned: InterProScan gene-model evidence loading",
                      file=sys.stderr)
            if not args.dry_run:
                with log.open("a") as handle:
                    handle.write(
                        "Existing predicted-model InterProScan TSV reused: "
                        f"{predicted_interpro_tsv}\n")
                print_reused_outputs("InterProScan")
                if (args.skip_nb_rescue
                        and is_nb_rescue_integrated_tsv(predicted_interpro_tsv)):
                    print_reused_outputs("HMMER\\MAST")
        else:
            predicted_ips_command = [
                executable(args.interproscan_bin, args.dry_run),
                "-i", str(predicted_proteins), "-appl", INTERPRO_APPLICATIONS,
                "-f", "TSV,GFF3", "-b", str(predicted_interpro_prefix),
                "-cpu", str(args.threads), "-dp",
            ]
            predicted_phobius_used = run_interproscan_command(
                "Protein domain annotation by InterProScan", predicted_ips_command,
                log, dry_run=args.dry_run,
                completed_label="Protein domain annotation by InterProScan")
        predicted_rescue_summary: dict[str, object] | None = None
        predicted_classification_tsv = predicted_interpro_tsv
        if not args.skip_nb_rescue:
            predicted_integrated_tsv = (
                predicted_interpro_dir /
                "predicted_gene_models.with_NB-ARC_rescue.tsv"
            )
            predicted_rescue_summary = screen_nb_arc_rescue(
                predicted_proteins, predicted_interpro_tsv,
                predicted_integrated_tsv,
                rescue_model_dir,
                predicted_interpro_dir / "NB-ARC_species_rescue", args, log,
                build_models=False, dry_run=args.dry_run)
            predicted_classification_tsv = predicted_integrated_tsv

        predicted_classifier_command = [
            sys.executable, str(classifier), "__classify__",
            "--proteins", str(predicted_proteins),
            "--interpro-tsv", str(predicted_classification_tsv),
            "--output", str(predicted_classification), "--quiet", "--summary-only",
            "--no-authorship",
        ]
        run_command("NB-LRR protein classification", predicted_classifier_command,
                    log, dry_run=args.dry_run)
        # Step 5/6 - Remove implausible or redundant models and write retained predictions to 05_filtered_NB-LRR_models.
        print_section(5, "Filtering of redundant gene models")
        if not args.dry_run:
            log_stage(log, 5, "Filtering of redundant NB-LRR gene models")
        filtering_command = [
            sys.executable, str(filter_script), "__filter__", "--gff", str(predicted_gff),
            "--proteins", str(predicted_proteins),
            "--classification", str(predicted_classification),
            "--query-classification", str(classification),
            "--output-prefix", str(filtered_prefix),
        ]
        run_command("Filtering of redundant NB-LRR gene models", filtering_command,
                    log, dry_run=args.dry_run)

        # Step 6/6 - Merge eligible annotated and predicted loci into the final catalogue in 06_merged_NB-LRR_genes.
        print_section(6, "Full NB-LRR resistance gene repertoire")
        if not args.dry_run:
            log_stage(log, 6, "Full NB-LRR resistance gene repertoire")
        recovery_prefix = merged_dir / "final_NB-LRR-prediction"
        filtered_gff = Path(str(filtered_prefix) + ".gff")
        filtered_proteins = Path(str(filtered_prefix) + ".pro")
        if annotation_gff is None:
            print("Skipped: no genome annotation GFF3 was provided.", file=sys.stderr)
            if args.dry_run:
                print("  " + styled("○", "1", "33")
                      + "  Planned: final catalogue generation", file=sys.stderr)
            else:
                write_predicted_only_catalog(filtered_gff, filtered_proteins,
                                             predicted_classification, recovery_prefix)
                with log.open("a") as handle:
                    handle.write("Merging skipped: no annotation GFF3 was supplied.\n")
        else:
            if not args.dry_run:
                merged_dir.mkdir(parents=True, exist_ok=True)
            recovery_command = [
                sys.executable, str(recovery_script), "__merge__",
                "--annotation-gff", str(annotation_gff),
                "--proteome", str(cleaned_proteome),
                "--classification", str(classification),
                "--filtered-genblastg-gff", str(filtered_gff),
                "--filtered-genblastg-proteins", str(filtered_proteins),
                "--predicted-classification", str(predicted_classification),
                "--output-prefix", str(recovery_prefix),
            ]
            run_command("Merging of annotated and predicted gene models", recovery_command, log,
                        dry_run=args.dry_run)

    if not args.dry_run:
        manifest = {
            "program": "HRPv2",
            "version": VERSION,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "proteome": str(proteome),
            "cleaned_proteome": str(cleaned_proteome),
            "proteins_with_terminal_stop_removed": terminal_stops_removed,
            "genome": str(genome) if genome else None,
            "annotation_gff": str(annotation_gff) if annotation_gff else None,
            "interproscan_tsv": str(interpro_tsv),
            "classification_evidence_tsv": str(classification_tsv),
            "species_specific_nb_arc_rescue": proteome_rescue_summary,
            "phobius_used_for_proteome": phobius_used,
            "phobius_used_for_predicted_models": (
                predicted_phobius_used if not args.skip_genblastg else None),
            "classification": str(classification),
            "full_length_query_fasta": str(query_fasta),
            "full_length_query_table": str(query_table),
            "full_length_class_counts": counts,
            "genblastg_command": genblastg_command,
            "predicted_genblastg_gff": str(predicted_gff) if predicted_gff else None,
            "predicted_genblastg_proteins": str(predicted_proteins) if predicted_proteins else None,
            "predicted_models_classification": (
                str(predicted_classification) if predicted_classification else None),
            "predicted_models_nb_arc_rescue": (
                predicted_rescue_summary if not args.skip_genblastg else None),
            "filtered_output_prefix": str(filtered_prefix) if filtered_prefix else None,
            "genblastg_outputs": sorted(
                path.name for path in genblastg_dir.glob("genblastg_Annotated_full-length_NB-LRRs*")
                if path.is_file()),
        }
        (workdir / "HRPv2_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        if not args.skip_genblastg:
            final_tsv = Path(str(recovery_prefix) + ".tsv")
            total_nlrs, full_length_nlrs, partial_nlrs = summarize_final_catalog(final_tsv)
            print(file=sys.stderr)
            print(f"Identified NB-LRR genes: {total_nlrs}", file=sys.stderr)
            print(f"  Full-length NB-LRRs: {full_length_nlrs}", file=sys.stderr)
            print(f"  Partial NB-LRRs: {partial_nlrs}", file=sys.stderr)
            finalize_log(
                log,
                final_tsv,
                Path(str(recovery_prefix) + ".fasta"),
                Path(str(recovery_prefix) + ".gff3"),
            )
        print(styled("=" * 88, "36"), file=sys.stderr)
        print(file=sys.stderr)
        print(styled(AUTHORSHIP, "1", "36"), file=sys.stderr)
        print(styled("=" * 88, "36"), file=sys.stderr)
    return 0



def _dispatch() -> int:
    """Route private component calls or start the public HRPv2 workflow."""
    if len(sys.argv) > 1 and sys.argv[1] in {"__classify__", "__filter__", "__merge__"}:
        action = sys.argv.pop(1)
        if action == "__classify__":
            return _classification_main()
        if action == "__filter__":
            return _filtering_main()
        return _merging_main()
    return _pipeline_main()


AUTHORSHIP = (
    'Giuseppe Andolfo, Ph.D.\n'
    'University of Naples "Federico II", Italy\n'
)


if __name__ == "__main__":
    try:
        raise SystemExit(_dispatch())
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
