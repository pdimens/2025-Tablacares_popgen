#! /usr/bin/env bash

# Extract outlier markers
bcftools view -Ob -l9 -T ../outlier_positions.txt ../YFT.kinless.bcf > ../YFT.outlier.kinless.bcf

# Extract neutral markers (by excluding outliers)
bcftools view -Ob -l9 -T ../neutral_positions.txt ../YFT.kinless.bcf > ../YFT.neutral.kinless.bcf
