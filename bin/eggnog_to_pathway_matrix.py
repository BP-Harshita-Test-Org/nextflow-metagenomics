#!/usr/bin/env python3
"""
Convert eggNOG-mapper annotations to a KEGG pathway abundance table.

eggNOG output columns (tab-separated, header starts with #):
  query  seed_ortholog  evalue  score  eggNOG_OGs  max_annot_lvl  COG_category
  Description  Preferred_name  GOs  EC  KEGG_ko  KEGG_Pathway  KEGG_Module
  KEGG_Reaction  KEGG_rclass  BRITE  KEGG_TC  CAZy  BiGG_Reaction  PFAMs

Output columns (same structure as Bracken, so format_lefse_input.py works unchanged):
  name                  KEGG pathway ID (e.g. map00010)
  fraction_total_reads  Relative abundance (genes in pathway / total annotated genes)
  count                 Number of genes assigned to this pathway
"""

import argparse
import csv
import sys
from collections import defaultdict


KEGG_PATHWAY_COL = "KEGG_Pathway"
QUERY_COL = "query"


def parse_eggnog_annotations(filepath: str) -> dict:
    """Read eggNOG annotations and count genes per KEGG pathway."""
    pathway_counts: dict[str, int] = defaultdict(int)
    total_annotated = 0

    with open(filepath) as f:
        # Skip comment lines until we hit the header
        header_line = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                # This is the column header line — strip the leading #
                header_line = line.lstrip("#").strip()
                break

        if header_line is None:
            print("ERROR: Could not find header line in eggNOG annotations file.", file=sys.stderr)
            sys.exit(1)

        columns = [c.strip() for c in header_line.split("\t")]

        if KEGG_PATHWAY_COL not in columns:
            print(f"ERROR: Column '{KEGG_PATHWAY_COL}' not found in eggNOG output.", file=sys.stderr)
            print(f"       Found columns: {columns}", file=sys.stderr)
            sys.exit(1)

        pathway_idx = columns.index(KEGG_PATHWAY_COL)

        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) <= pathway_idx:
                continue

            total_annotated += 1
            pathway_field = row[pathway_idx].strip()

            # KEGG_Pathway field contains comma-separated pathway IDs
            # e.g. "map00010,map00020,map01100"
            # Entries with no pathway are "-" or empty
            if pathway_field and pathway_field != "-":
                for pathway in pathway_field.split(","):
                    pathway = pathway.strip()
                    if pathway:
                        pathway_counts[pathway] += 1

    return pathway_counts, total_annotated


def write_pathway_table(pathway_counts: dict, total_annotated: int, output_path: str) -> None:
    """
    Write pathway abundance table in Bracken-compatible format so that
    format_lefse_input.py can consume it unchanged.
    """
    if total_annotated == 0:
        print("WARNING: No annotated genes found. Output will be empty.", file=sys.stderr)

    rows = []
    for pathway, count in pathway_counts.items():
        rel_abund = count / total_annotated if total_annotated > 0 else 0.0
        rows.append({
            "name": pathway,
            "count": count,
            "fraction_total_reads": round(rel_abund, 8),
        })

    # Sort by relative abundance descending
    rows.sort(key=lambda r: r["fraction_total_reads"], reverse=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "count", "fraction_total_reads"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written {len(rows)} KEGG pathways from {total_annotated} annotated genes → {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert eggNOG-mapper annotations to KEGG pathway abundance table"
    )
    parser.add_argument("--annotations", required=True,
                        help="eggNOG-mapper .emapper.annotations file")
    parser.add_argument("--output", required=True,
                        help="Output pathway abundance TSV")
    parser.add_argument("--sample_id", default="sample",
                        help="Sample identifier (for logging)")
    args = parser.parse_args()

    print(f"Processing eggNOG annotations for sample: {args.sample_id}")
    print(f"Input: {args.annotations}")

    pathway_counts, total_annotated = parse_eggnog_annotations(args.annotations)

    print(f"Total annotated genes : {total_annotated}")
    print(f"Unique KEGG pathways  : {len(pathway_counts)}")

    write_pathway_table(pathway_counts, total_annotated, args.output)


if __name__ == "__main__":
    main()
