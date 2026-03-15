#!/usr/bin/env python3
"""Format merged abundance matrix + metadata into LEfSe-compatible input."""

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description="Format abundance data for LEfSe")
    parser.add_argument("--abundance", required=True, help="Merged abundance matrix TSV")
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

    # Read merged abundance matrix (feature \t sample1 \t sample2 ...)
    with open(args.abundance) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)  # ["feature", "sample1", "sample2", ...]
        sample_ids = header[1:]

        features = {}
        for row in reader:
            feature_name = row[0]
            abundances = row[1:]
            features[feature_name] = abundances

    # Write LEfSe format
    with open(args.output, "w") as f:
        # Row 1: class labels (ASD / neurotypical)
        class_labels = [sample_groups.get(s, "unknown") for s in sample_ids]
        f.write("class\t" + "\t".join(class_labels) + "\n")

        # Row 2: subject IDs
        f.write("subject\t" + "\t".join(sample_ids) + "\n")

        # Feature rows
        for feature, abundances in features.items():
            f.write(f"{feature}\t" + "\t".join(abundances) + "\n")

    print(f"LEfSe input written to {args.output} ({len(features)} features, {len(sample_ids)} samples)")


if __name__ == "__main__":
    main()
