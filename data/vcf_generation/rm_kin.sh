#! /usr/bin/env bash

bcftools view -Ob -l9 -S ^../../kinship/kin_to_rm.txt ../out.14.bcf > ../YFT.kinless.bcf