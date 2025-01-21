from snakemake.utils import R

rule all:
    input: "out.14.recode.vcf", "out.12.clones"
    output: "filter.summary"
    shell:
        """
        echo "filename ind_after ind_before sites_after sites_before" > {output}.tmp

        for i in *.log; do
            SITES=$(grep "Sites$" $i | cut -d" " -f9,4)
            INDV=$(grep "Individuals$" $i | cut -d" " -f7,4)
            echo "$i $INDV $SITES" >> {output}.tmp
        done
        cat {output}.tmp | column -t > {output} && rm {output}.tmp
        """

rule compress_VCF:
    input: "TotalRawSNPs.vcf"
    output: "TotalRawSNPs.vcf.gz"
    message: "bgzip compressing {input}"
    shell: "bgzip {input} && tabix {output}"

rule missing_50:
    input: "TotalRawSNPs.vcf.gz"
    output: "out.1.recode.vcf.gz"
    log: report("out.1.log")
    message: "Filter sites with >50% missing data"
    params:
        param1 = "--max-missing 0.5",
        out = "--out out.1"
    shell: 
        """
        vcftools --gzvcf {input} {params} --recode --recode-INFO-all --stdout 2> {log}.tmp | bgzip > {output} 
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        tabix {output}
        """

rule depth10_quality20:
    input: "out.1.recode.vcf.gz"
    output: "out.2.recode.vcf.gz"
    log: report("out.2.log")
    message: "Filter genotypes with depth < 10 and genotype quality < 20"
    params:
        param1 = "--minDP 10",
        param2 = "--minGQ 20",
        out = "--out out.2"
    shell: 
        """
        vcftools --gzvcf {input} {params} --recode --recode-INFO-all --stdout 2> {log}.tmp | bgzip > {output} 
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        tabix {output}
        """

rule purge_monomorphic:
    input: "out.2.recode.vcf.gz"
    output: "out.3.recode.vcf.gz"
    message: "Filter out sites that were made monomorphic by the previous filter"
    params:
        param1 = "-q 0.001:minor"
    shell: 
        """
        bcftools view -Oz -q 0.001:minor {input} > {output}
    	tabix {output}
        """

rule calculate_missing_indiv:
    input: "out.3.recode.vcf.gz"
    output: report("out.3.imiss")
    message: "Produce a file with missingness per individual"
    params:
        out = "--out out.3",
        param1 = "--missing-indv"
    shell: "vcftools --gzvcf {input} {params} 2> /dev/null"


rule find_missing70_indiv:
    input: "out.3.imiss"
    output: report("remove.3.indivs")
    log: report("out.3.imiss.png")
    message: "finding individuals with more than 70% missing data"
    run:
        R("""
        suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
        library(stringr)
        
        out_3_imiss <- read_tsv("{input}") %>% suppressMessages()
        png(filename = "out.3.imiss.png", width = 480, height = 480, units = "px")
        
        # Plot a quick histogram of the data
        qplot(out_3_imiss$F_MISS)
        dev.off()
        
        # Select individuals with more than 70% missing data
        miss_70 <- filter(out_3_imiss, F_MISS > 0.7) %>%
        select(INDV)

        write_delim(miss_70, "{output}", col_names = FALSE)
        """)

rule missing70:
    input: 
        vcf = "out.3.recode.vcf.gz",
        inds = "remove.3.indivs"
    output: 
        vcf = "out.4.recode.vcf",
    log: report("out.4.log")
    message: "Remove individuals with >70% missing data"
    params:
        out = "--out out.4"
    shell: 
        """
        vcftools --gzvcf {input.vcf} {params} --remove {input.inds} --recode --recode-INFO-all 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule overall_qual20:
    input: "out.4.recode.vcf"
    output: "out.5.vcf"
    log: report("out.5.log")
    message: "Remove loci with overall quality < 20"
    params: "QUAL > 20"
    shell: 
        """
        vcffilter -s -f '{params}' {input} > {output}
        vcftools --vcf {output} 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule dDocent_filters:
    input: "out.5.vcf"
    output:
        balance = temp("out.5.balance.vcf"),
        pairing = temp("out.5.pairing.vcf"),
        overlap = "out.5.ddfilters.vcf"
    message: "Applying dDocent_filters for allele balance, proper pairing, and overlapping"
    log: 
        balance = report("out.5.balance.log"),
        pairing = report("out.5.pairing.log"),
        overlap = report("out.5.ddfilters.log")
    shell:
        """
        echo "Removing loci with extreme allele balance"
        vcffilter -s -f 'AB > 0.2 & AB < 0.8' {input} > {output.balance}
        vcftools --vcf {output.balance} 2> {log.balance}.tmp
        grep -v "^Warning" {log.balance}.tmp > {log.balance} && rm {log.balance}.tmp

        echo "Filtering sites called from improperly paired reads using the defaults from dDocent_filters"
        vcffilter -s -f 'PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05' {output.balance} > {output.pairing}
        vcftools --vcf {output.pairing} 2> {log.pairing}.tmp
        grep -v "^Warning" {log.pairing}.tmp > {log.pairing} && rm {log.pairing}.tmp

        echo "Filtering sites that are called from overlapping forward and reverse reads"
        vcffilter -s -f 'SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 50' {output.pairing} > {output.overlap}
        vcftools --vcf {output.overlap} 2> {log.overlap}.tmp
        grep -v "^Warning" {log.overlap}.tmp > {log.overlap} && rm {log.overlap}.tmp
        """


rule calculate_site_depth:
    input: "out.5.ddfilters.vcf"
    output: 
        depth = report("out.5.ldepth"),
        inds = "out.5.indvcount",
        log = temp(".out.5.count.log")
    message: "Produce a file with depth per site and the number of individuals"
    params:
        out = "--out out.5",
        param1 = "--site-depth"
    shell: 
        """
        vcftools --vcf {input} {params} 2> {output.log} 
        grep "Individuals$"  {output.log} | rev | cut -d" " -f5 | rev > {output.inds}
        """


rule find_depth_threshold:
    input: 
        depth = "out.5.ldepth",
        counts = "out.5.indvcount"
    output: report("remove.5.sites")
    log: 
        report("out.5.ldepth.png"),
        report("out.5.ldepth.post.png")
    message: "finding sites with too much depth"
    run:
        R("""
        suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
        library(stringr)
        
        inds <- as.vector(unlist(read.table("{input.counts}")))
        site_depth <- suppressMessages(read_tsv("{input.depth}")) %>%
            mutate(MEAN_DEPTH = SUM_DEPTH / inds)
        
        # Plot a histogram of the mean site depth per individual
        qplot(site_depth$MEAN_DEPTH)
        ggsave("out.5.ldepth.png", width = 4, height = 4, units = "in")
        
        # Filter out loci with a mean site depth > 95% quantile
        mean_site_depth <- mean(site_depth$MEAN_DEPTH)
        thresh <- quantile(site_depth$MEAN_DEPTH, .95)
        to_keep <- filter(site_depth, MEAN_DEPTH < thresh)
        mean_site_depth_filt <- mean(to_keep$MEAN_DEPTH)
        
        # Plot the distribution again
        qplot(to_keep$MEAN_DEPTH)
        ggsave("out.5.ldepth.post.png", width = 4, height = 4, units = "in")

        # Make a list of the sites to filter
        to_filter <- filter(site_depth, MEAN_DEPTH >= thresh) %>%
            select(CHROM, POS)
        
        write_delim(to_filter, "remove.5.sites", col_names = FALSE)
        """)


rule site_depth:
    input: 
        vcf = "out.5.ddfilters.vcf",
        sites = "remove.5.sites"
    output: 
        vcf = "out.6.recode.vcf",
    log: report("out.6.log")
    message: "Remove really high coverage sites"
    params:
        out = "--out out.6"
    shell: 
        """
        vcftools --vcf {input.vcf} {params} --exclude-positions {input.sites} --recode --recode-INFO-all 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule site_missing75:
    input: "out.6.recode.vcf"
    output: "out.7.recode.vcf"
    log: report("out.7.log")
    message: "Remove sites with >75% missing data"
    params:
        out = "--out out.7",
        param1 = "--max-missing 0.75"
    shell: 
        """
        vcftools --vcf {input} {params} --recode --recode-INFO-all 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule calculate_ind_missing:
    input: "out.7.recode.vcf"
    output: report("out.7.imiss")
    message: "Calculating individual missingness"
    params:
        out = "--out out.7",
        param1 = "--missing-indv"
    shell: "vcftools --vcf {input} {params} 2> /dev/null"

rule find_missing_indivs:
    input: "out.7.imiss"
    output: report("remove.7.indivs")
    log: report("out.7.imiss.png")
    message: "finding sites with too much depth"
    run:
        R("""
        suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
        library(stringr)
        
        out_9_imiss <- suppressMessages(read_tsv("{input}")) %>%
            extract(INDV, "pop", "(\\\w+)_", remove = FALSE)

        out_9_imiss %>% count(pop)

        # Plot a quick histogram of the data
        ggplot(out_9_imiss, aes(x = F_MISS, fill = pop)) +
            geom_histogram()
        ggsave("{log}")

        # Select individuals with more than 40% missing data
        miss_40 <- filter(out_9_imiss, F_MISS > 0.40) %>%
            select(INDV)

        # Write the individuals to remove to a file
        write_delim(miss_40, "{output}", col_names = FALSE)
        """)

rule missing_indiv40:
    input: 
        vcf = "out.7.recode.vcf",
        inds = "remove.7.indivs"
    output: "out.8.recode.vcf"
    message: "Removing individuals with >40% missing genotypes"
    log: report("out.8.log")
    params:
        out = "--out out.8"
    shell: 
        """
        vcftools --vcf {input.vcf} --remove {input.inds} {params} --recode --recode-INFO-all 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule missing_by_population:
    input: 
        vcf = "out.8.recode.vcf",
        popmap = "popmap"
    output: 
        vcf = "out.9.recode.vcf",
        loci = report("out.9.badloci"),
        all_loci = temp("out.9.badloci.all")
    log: report("out.9.log")
    message: "Remove sites with more than 15% missing data in >2 populations"
    params:
        proportion = 0.15,
        cutoff = 2,
        out = "out.9"
    shell: 
        """
        for i in $(cut -f2 {input.popmap} | sort | uniq); do
            grep -w "$i" {input.popmap} | cut -f1 > keep.$i
            echo "Working on population $i"
            vcftools --vcf {input.vcf} --keep keep.$i --missing-site --out $i 2> /dev/null
            mawk '!/CHROM/' $i.lmiss | mawk -v x={params.proportion} '$6 > x' | cut -f1,2 >> {output.all_loci}
            rm keep.$i $i.lmiss
        done

        mawk '!/CH/' {output.all_loci} | perl -e 'while (<>) {{chomp; $z{{$_}}++;}} while(($k,$v) = each(%z)) {{print "$v\t$k\n";}}' | mawk -v x={params.cutoff} '$1 >= x' | cut -f2,3  > {output.loci}
        echo "excluding $(cat {output.loci} | wc -l) positions"
        vcftools --vcf {input.vcf} --exclude-positions {output.loci} --recode --recode-INFO-all --out {params.out} 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule decompose_mnp:
    input: "out.9.recode.vcf"
    output: "out.10.snp.vcf"
    log: report("out.10.log")
    message: "Decomposing MNPs into SNPs"
    shell: 
        """
        vcfallelicprimitives -k -g {input} > {output}
        vcftools --vcf {output} 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule biallelic_noindel:
    input: "out.10.snp.vcf"
    output: "out.11.recode.vcf"
    log: report("out.11.log")
    message: "Restricting data to biallelic SNPs and removing indels"
    params:
        out = "--out out.11",
        param1 = "--remove-indels",
        param2 = "--max-alleles 2",
        param3 = "--min-alleles 2"
    shell: 
        """
        vcftools --vcf {input} {params} --recode --recode-INFO-all 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule minor_allele_05:
    input: "out.11.recode.vcf"
    output: "out.12.recode.vcf"
    message: "Filter for minor allele frequency <0.05"
    log: report("out.12.log")
    params:
        out = "--out out.12",
        param1 = "--maf 0.05"
    shell: 
        """
        vcftools --vcf {input} {params} --recode --recode-INFO-all 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule calculate_relatedness:
    input: "out.12.recode.vcf"
    output: "out.12.relatedness"
    message: "Calculating relatedness"
    params:
        out = "--out out.12"
    shell: "vcftools --vcf {input} {params} --relatedness 2> /dev/null"

rule find_duplicates:
    input: "out.12.relatedness"
    output: "out.12.clones"
    message: "Finding duplicates / clones"
    params:
        clone_thresh = 0.8
    run:
        R("""
        suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
        library(stringr)
        
        out_relate <- read_tsv("{input}") %>% supppressMessages()

        # Check for abnormally high relatedness
        relatives <- filter(out_relate, INDV1 != INDV2) %>% 
            filter(RELATEDNESS_AJK > {params.clone_thresh})

        # Make a list of individuals that may be contaminated or clones
        inds_to_remove <- rbind(relatives$INDV1, relatives$INDV2) %>% as.data.frame()

        # Write the individuals to remove to a file
        write_delim(inds_to_remove, "{output}", col_names = FALSE)
        """)

rule unphase:
    input: "out.12.recode.vcf"
    output: 
        unphased = "out.12.unphased.vcf",
        missrecode = "out.12.missrecode.vcf"
    message: "Unphase genotypes and recode missing values"
    shell:
        """
        sed '/^##/! s/|/\//g' {input} > {output.unphased}
        cat {output.unphased} | perl -pe 's/\s\.:/\t.\/.:/g' > {output.missrecode}
        """

rule hardyweinberg:
    input: 
        vcf = "out.12.missrecode.vcf",
        popmap = "popmap"
    output: "out.13.recode.vcf"
    message: "Filter based on HWE in each population"
    params:
        out = "--out out.13",
        h = "--hwe 0.001",
        c = "--cutoff 0.33"
    shell: 
        """
        perl $(which filter_hwe_by_pop.pl) -v {input.vcf} -p {input.popmap} {params}
        rm *.inds *.hwe
        """

rule missing90:
    input: "out.13.recode.vcf"
    output: "out.14.recode.vcf"
    log: report("out.14.log")
    message: "Filter sites with >10% missing data"
    params:
        out = "--out out.14",
        param1 = "--max-missing 0.9"
    shell: 
        """
        vcftools --vcf {input} {params} --recode --recode-INFO-all 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """

rule thin:
    input: "out.14.recode.vcf"
    output: 
        vcf = "out.15.vcf",
        sites = "out.15.sites"
    message: "Thinning sites by 1000bp"
    params:
        thinning_val = "1000"
    shell:
        """
        bcftools +prune -w {params.thinning_val}bp -n 1 {input} -o {output.vcf}

        # Get a list of the scaffolds and positions
        grep -v '#' {output.vcf} | cut -f1,2 > {output.sites}
        """

rule thin_to_keep:
    input: "out.15.sites"
    output: "out.15.keep"
    message: "Finding sites to keep"
    params:
        prefix = "HiC_scaffold_"
    run:
        R("""
        suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
        library(stringr)
        read_tsv("{input}", col_names = FALSE) %>%
            extract(X1, "chrom", "{params.prefix}(\\d+)_", remove = FALSE) %>%
            mutate(chrom = as.integer(chrom)) %>%
            filter(chrom <= 19) %>%
            select(X1, X2) %>%
        
        write_tsv("{output}"), col_names = FALSE)
        """)


rule thin_keep:
    input:
        vcf = "out.15.vcf",
        sites = "out.15.keep"
    output: "out.16.recode.vcf"
    log: "out.16.log"
    message: "Retaining only the sites listed in {input.sites}"
    params:
        out = "--out out.16"
    shell:
        """
        vcftools --vcf {input.vcf} {params} --positions {input.sites} --recode --recode-INFO-all 2> {log}.tmp
        grep -v "^Warning" {log}.tmp > {log} && rm {log}.tmp
        """