#!/usr/bin/env python3
"""Merge individual Bracken/pathway abundance TSVs into a samples x features matrix."""

import argparse
import csv
import os


def main():
    parser = argparse.ArgumentParser(description="Merge per-sample abundance tables")
    parser.add_argument("--input_files", nargs="+", required=True,
                        help="Per-sample abundance TSV files")
    parser.add_argument("--output", required=True, help="Output merged matrix TSV")
    parser.add_argument("--type", choices=["taxa", "pathways"], default="taxa",
                        help="Type of abundance data")
    args = parser.parse_args()

    all_features = {}   # feature_name -> {sample_id: abundance}
    sample_ids = []

    for filepath in args.input_files:
        basename = os.path.basename(filepath)
        # Extract sample_id from filename conventions:
        #   SAMPLE1_bracken_species.tsv  -> SAMPLE1
        #   SAMPLE1_pathway_abundance.tsv -> SAMPLE1
        if "_bracken_" in basename:
            sample_id = basename.split("_bracken_")[0]
        elif "_pathway_" in basename:
            sample_id = basename.split("_pathway_")[0]
        else:
            sample_id = os.path.splitext(basename)[0]

        sample_ids.append(sample_id)

        with open(filepath) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                name = row.get("name", "")
                abundance = float(row.get("fraction_total_reads", 0))
                if name:
                    if name not in all_features:
                        all_features[name] = {}
                    all_features[name][sample_id] = abundance

    # Write merged matrix: rows = features, columns = samples
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["feature"] + sample_ids)
        for feature in sorted(all_features.keys()):
            row = [feature]
            for sid in sample_ids:
                row.append(all_features[feature].get(sid, 0.0))
            writer.writerow(row)

    print(f"Merged {len(sample_ids)} samples, {len(all_features)} features -> {args.output}")


if __name__ == "__main__":
    main()
