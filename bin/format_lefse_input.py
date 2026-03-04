#!/usr/bin/env python3
"""Format Bracken/eggNOG abundance tables into LEfSe-compatible input."""

import argparse
import csv
import sys


def main():
    parser = argparse.ArgumentParser(description="Format abundance data for LEfSe")
    parser.add_argument("--abundance", required=True, help="Abundance table (TSV)")
    parser.add_argument("--metadata", required=True, help="Sample metadata (CSV)")
    parser.add_argument("--output", required=True, help="Output LEfSe-formatted file")
    args = parser.parse_args()

    # Read metadata to get sample-group mapping
    sample_groups = {}
    with open(args.metadata) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = row.get("sample_id", "")
            group = row.get("group", row.get("diagnosis", "unknown"))
            sample_groups[sid] = group

    # Read abundance table
    features = {}
    with open(args.abundance, delimiter="\t") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row.get("name", "")
            abundance = row.get("fraction_total_reads", "0")
            if name:
                features[name] = abundance

    # Write LEfSe format:
    # Row 1: class labels (ASD / neurotypical)
    # Row 2+: feature\tabundance_per_sample
    with open(args.output, "w") as f:
        # For single-sample, write a simple two-column format
        # For multi-sample, this would be expanded to columns per sample
        if sample_groups:
            samples = list(sample_groups.keys())
            groups = [sample_groups[s] for s in samples]
            f.write("class\t" + "\t".join(groups) + "\n")
        else:
            f.write("class\tsample\n")

        for feature, abund in features.items():
            f.write(f"{feature}\t{abund}\n")

    print(f"LEfSe input written to {args.output} ({len(features)} features)")


if __name__ == "__main__":
    main()
