freebayes -b cat-RRG.bam -v YFT_nonlarvae.vcf -l YFT_TXL.4.recode.vcf -f reference.fasta -m 5 -q 5 -E 3 --min-repeat-entropy 1 -V --populations popmap -n 10 -F 0.1

freebayes-parallel <(fasta_generate_regions.py reference.fasta.fai 100000) 20 -b cat-RRG.bam -l YFT_TXL.4.recode.vcf -f reference.fasta -m 5 -q 5 -E 3 --min-repeat-entropy 1 -V --populations popmap -n 10 -F 0.1 > YFT_nonlarvae.vcf