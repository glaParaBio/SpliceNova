#!/usr/bin/env python3

import argparse
import re
import sys

from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(
    description="Add transcript id to selected feature types"
)
parser.add_argument(
    "gff", type=str, help="Input GFF file [%(default)s]", default="-", nargs="?"
)
parser.add_argument(
    "--type",
    "-t",
    type=str,
    help="List of feature types (column 3) identifying a transcript child %(default)s",
    default=["exon", "CDS"],
    nargs="+",
)
parser.add_argument(
    "--transcript-key",
    "-k",
    help="Attribute for transcript [%(default)s]",
    default="transcript_id",
)
parser.add_argument(
    "--overwrite-id",
    "-o",
    help="Overwrite attribute if it exists. Default is to exit with error",
    action="store_true",
)
parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")

args = parser.parse_args()


def key_exists(attributes, key):
    found = [x for x in attributes if x.startswith(f"{key}=")]
    return len(found) > 0


if args.gff == "-":
    fin = sys.stdin
else:
    fin = open(args.gff)

for line in fin:
    if line.startswith("#"):
        sys.stdout.write(line)
        continue
    line = line.strip().split("\t")
    if line[2] in args.type:
        attrs = line[8].split(";")
        if not args.overwrite_id and key_exists(attrs, args.transcript_key):
            sys.stderr.write(f'Key "{args.transcript_key}" already present\n')
            sys.exit(1)
        attrs = [x for x in attrs if not x.startswith(f"{args.transcript_key}=")]
        parent = [x for x in attrs if x.startswith("Parent=")]
        if len(parent) != 1:
            sys.stderr.write("Invalid Parent\n")
            sys.exit(1)
        parent = re.sub("^Parent=", "", parent[0])
        attrs.append(f"{args.transcript_key}={parent}")
        line[8] = ";".join(attrs)
    print("\t".join(line))
