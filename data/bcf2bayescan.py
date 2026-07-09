#!/usr/bin/env python3

import argparse
from collections import defaultdict
import pysam
import sys

def load_popmap(filename):
    popmap = {}
    with open(filename) as fh:
        for line_number, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(
                    f"{filename}:{line_number}: expected two columns "
                    "(sample population)"
                )
            sample, population = fields
            if sample in popmap:
                raise ValueError(
                    f"Duplicate sample '{sample}' in popmap."
                )
            popmap[sample] = population
    return popmap


def main():
    parser = argparse.ArgumentParser(
        description="Convert a BCF/VCF into BayeScan format."
    )
    parser.add_argument("bcf")
    parser.add_argument("popmap")
    #parser.add_argument("output")

    args = parser.parse_args()
    popmap = load_popmap(args.popmap)
    vcf = pysam.VariantFile(args.bcf)
    samples = list(vcf.header.samples)
    missing = [s for s in samples if s not in popmap]
    if missing:
        raise RuntimeError(
            "The following samples are missing from the population map:\n"
            + "\n".join(missing)
        )
    extra = set(popmap) - set(samples)
    if extra:
        raise RuntimeError(
            "The following samples are present in the population map "
            "but not in the BCF:\n"
            + "\n".join(sorted(extra))
        )

    populations = sorted(set(popmap.values()))
    pop_number = {p: i + 1 for i, p in enumerate(populations)}
    sample_indices = defaultdict(list)
    for idx, sample in enumerate(samples):
        sample_indices[popmap[sample]].append(idx)

    # counts[population] = list of (ref_count, alt_count)
    counts = {p: [] for p in populations}
    nloci = 0

    for record in vcf:
        # Skip anything except biallelic SNPs.
        if len(record.alleles) != 2:
            continue
        if len(record.ref) != 1:
            continue
        if len(record.alts[0]) != 1:
            continue
        nloci += 1
        for pop in populations:
            ref = 0
            alt = 0
            for idx in sample_indices[pop]:
                gt = record.samples[idx]["GT"]
                if gt is None:
                    continue
                for allele in gt:
                    if allele is None:
                        continue
                    if allele == 0:
                        ref += 1
                    elif allele == 1:
                        alt += 1
            counts[pop].append((ref, alt))

    #with open(args.output, "w") as out:
    sys.stdout.write(f"[loci]={nloci}\n\n")
    sys.stdout.write(f"[populations]={len(populations)}\n\n")
    for pop in populations:
        sys.stdout.write(f"[pop]={pop_number[pop]}\n")
        for locus, (ref, alt) in enumerate(counts[pop], start=1):
            nchrom = ref + alt
            sys.stdout.write(f" {locus}\t{nchrom}\t2\t{alt} {ref}\n")
            #sys.stdout.write(f"{locus} 2 {ref} {alt}\n")
        sys.stdout.write("\n")

if __name__ == "__main__":
    main()
