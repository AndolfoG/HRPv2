#!/usr/bin/env python3
"""Generate the alignscore.txt file required by GenBlastG.

Data source: NCBI BLOSUM62 matrix by Aaron Ucko and Mike Gertz:
https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/util/tables/sm_blosum62.c
"""

import argparse
import difflib
from pathlib import Path
import sys


SYMBOLS = "ARNDCQEGHILKMFPSTWYVBJZX*"
STANDARD = "ARNDCQEGHILKMFPSTWYV"
AMBIGUOUS = "BJZ"

BLOSUM62 = (
    (  4,  -1,  -2,  -2,   0,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   0,  -3,  -2,   0,  -2,  -1,  -1,  -1,  -4),
    ( -1,   5,   0,  -2,  -3,   1,   0,  -2,   0,  -3,  -2,   2,  -1,  -3,  -2,  -1,  -1,  -3,  -2,  -3,  -1,  -2,   0,  -1,  -4),
    ( -2,   0,   6,   1,  -3,   0,   0,   0,   1,  -3,  -3,   0,  -2,  -3,  -2,   1,   0,  -4,  -2,  -3,   4,  -3,   0,  -1,  -4),
    ( -2,  -2,   1,   6,  -3,   0,   2,  -1,  -1,  -3,  -4,  -1,  -3,  -3,  -1,   0,  -1,  -4,  -3,  -3,   4,  -3,   1,  -1,  -4),
    (  0,  -3,  -3,  -3,   9,  -3,  -4,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -3,  -1,  -1,  -2,  -2,  -1,  -3,  -1,  -3,  -1,  -4),
    ( -1,   1,   0,   0,  -3,   5,   2,  -2,   0,  -3,  -2,   1,   0,  -3,  -1,   0,  -1,  -2,  -1,  -2,   0,  -2,   4,  -1,  -4),
    ( -1,   0,   0,   2,  -4,   2,   5,  -2,   0,  -3,  -3,   1,  -2,  -3,  -1,   0,  -1,  -3,  -2,  -2,   1,  -3,   4,  -1,  -4),
    (  0,  -2,   0,  -1,  -3,  -2,  -2,   6,  -2,  -4,  -4,  -2,  -3,  -3,  -2,   0,  -2,  -2,  -3,  -3,  -1,  -4,  -2,  -1,  -4),
    ( -2,   0,   1,  -1,  -3,   0,   0,  -2,   8,  -3,  -3,  -1,  -2,  -1,  -2,  -1,  -2,  -2,   2,  -3,   0,  -3,   0,  -1,  -4),
    ( -1,  -3,  -3,  -3,  -1,  -3,  -3,  -4,  -3,   4,   2,  -3,   1,   0,  -3,  -2,  -1,  -3,  -1,   3,  -3,   3,  -3,  -1,  -4),
    ( -1,  -2,  -3,  -4,  -1,  -2,  -3,  -4,  -3,   2,   4,  -2,   2,   0,  -3,  -2,  -1,  -2,  -1,   1,  -4,   3,  -3,  -1,  -4),
    ( -1,   2,   0,  -1,  -3,   1,   1,  -2,  -1,  -3,  -2,   5,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2,   0,  -3,   1,  -1,  -4),
    ( -1,  -1,  -2,  -3,  -1,   0,  -2,  -3,  -2,   1,   2,  -1,   5,   0,  -2,  -1,  -1,  -1,  -1,   1,  -3,   2,  -1,  -1,  -4),
    ( -2,  -3,  -3,  -3,  -2,  -3,  -3,  -3,  -1,   0,   0,  -3,   0,   6,  -4,  -2,  -2,   1,   3,  -1,  -3,   0,  -3,  -1,  -4),
    ( -1,  -2,  -2,  -1,  -3,  -1,  -1,  -2,  -2,  -3,  -3,  -1,  -2,  -4,   7,  -1,  -1,  -4,  -3,  -2,  -2,  -3,  -1,  -1,  -4),
    (  1,  -1,   1,   0,  -1,   0,   0,   0,  -1,  -2,  -2,   0,  -1,  -2,  -1,   4,   1,  -3,  -2,  -2,   0,  -2,   0,  -1,  -4),
    (  0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   5,  -2,  -2,   0,  -1,  -1,  -1,  -1,  -4),
    ( -3,  -3,  -4,  -4,  -2,  -2,  -3,  -2,  -2,  -3,  -2,  -3,  -1,   1,  -4,  -3,  -2,  11,   2,  -3,  -4,  -2,  -2,  -1,  -4),
    ( -2,  -2,  -2,  -3,  -2,  -1,  -2,  -3,   2,  -1,  -1,  -2,  -1,   3,  -3,  -2,  -2,   2,   7,  -1,  -3,  -1,  -2,  -1,  -4),
    (  0,  -3,  -3,  -3,  -1,  -2,  -2,  -3,  -3,   3,   1,  -2,   1,  -1,  -2,  -2,   0,  -3,  -1,   4,  -3,   2,  -2,  -1,  -4),
    ( -2,  -1,   4,   4,  -3,   0,   1,  -1,   0,  -3,  -4,   0,  -3,  -3,  -2,   0,  -1,  -4,  -3,  -3,   4,  -3,   0,  -1,  -4),
    ( -1,  -2,  -3,  -3,  -1,  -2,  -3,  -4,  -3,   3,   3,  -3,   2,   0,  -3,  -2,  -1,  -2,  -1,   2,  -3,   3,  -3,  -1,  -4),
    ( -1,   0,   0,   1,  -3,   4,   4,  -2,   0,  -3,  -3,   1,  -1,  -3,  -1,   0,  -1,  -2,  -2,  -2,   0,  -3,   4,  -1,  -4),
    ( -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -4),
    ( -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,   1),
)


def validate_matrix() -> None:
    if len(BLOSUM62) != len(SYMBOLS) or any(
        len(row) != len(SYMBOLS) for row in BLOSUM62
    ):
        raise ValueError("BLOSUM62 must be a 25 x 25 matrix")
    if any(
        BLOSUM62[i][j] != BLOSUM62[j][i]
        for i in range(25)
        for j in range(25)
    ):
        raise ValueError("BLOSUM62 must be symmetric")
    x = SYMBOLS.index("X")
    stop = SYMBOLS.index("*")
    if any(BLOSUM62[x][j] != -1 for j in range(24)):
        raise ValueError("GenBlastG's Xx wildcard requires X scores of -1")
    if any(BLOSUM62[stop][j] != -4 for j in range(24)):
        raise ValueError("GenBlastG's *x wildcard requires stop scores of -4")
    if BLOSUM62[stop][stop] != 1:
        raise ValueError("The stop/stop score must be 1")


def generate() -> str:
    validate_matrix()
    index = {symbol: i for i, symbol in enumerate(SYMBOLS)}
    lines = ["-o:-11", "-e:-1"]

    for i, first in enumerate(STANDARD):
        for second in STANDARD[i:]:
            lines.append(f"{first}{second}:{BLOSUM62[index[first]][index[second]]}")
    for first in STANDARD:
        for second in AMBIGUOUS:
            lines.append(f"{first}{second}:{BLOSUM62[index[first]][index[second]]}")
    for i, first in enumerate(AMBIGUOUS):
        for second in AMBIGUOUS[i:]:
            lines.append(f"{first}{second}:{BLOSUM62[index[first]][index[second]]}")

    lines.extend(("Xx:-1", "*x:-4", "**:1"))
    if len(lines) != 281:
        raise AssertionError(f"Expected 281 lines, got {len(lines)}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("alignscore.txt"),
        help="output path (default: ./alignscore.txt)",
    )
    parser.add_argument(
        "--compare", type=Path, metavar="REFERENCE",
        help="also compare the generated file byte for byte with REFERENCE",
    )
    args = parser.parse_args()
    if args.compare is not None and args.output.resolve() == args.compare.resolve():
        parser.error("--compare must differ from --output")
    output = generate().encode("ascii")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(output)
    print(f"Wrote {args.output} ({len(output)} bytes)", file=sys.stderr)

    if args.compare is not None:
        reference = args.compare.read_bytes()
        if output != reference:
            diff = difflib.unified_diff(
                reference.decode("ascii").splitlines(keepends=True),
                output.decode("ascii").splitlines(keepends=True),
                fromfile=str(args.compare), tofile=str(args.output),
            )
            sys.stderr.writelines(diff)
            print("Files differ (including possible line endings).", file=sys.stderr)
            return 1
        print("Byte-for-byte identical to reference.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# Author: Giuseppe Andolfo, Ph.D., University of Naples “Federico II”, Italy
