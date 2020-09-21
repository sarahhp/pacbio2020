"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2020-09-02"

FASTQ_DIR = config["fastq_dir"]

import pandas
import get_fastq_sum as tools


rule cupcake_all:
    input:
        expand("cupcake/{sample}/filt-f.gff",
                sample = config["samples"]),
        expand("cupcake/{sample}/filt.min_fl_50-f.gff",
                sample =  config["samples"]),
        expand("cupcake/{sample}/filt.min_fl_100-f.gff",
                sample = config["samples"])

rule raw_abundance:
    input:
        expand("cupcake/{sample}/A.collapsed.abundance.txt",
                sample = config["samples"])

rule cupcake_collapse:
    input: 
        fq = "filtered_fastp/{sample}_hq.f.fastq",
        sam = "hotair_selection/{sample}_hq.sam"
    params:
        out_prefix = lambda wildcards, output: output.gff[:-14],
        extra = "--dun-merge-5-shorter"
    output:
        groups = "cupcake/{sample}/A.collapsed.group.txt",
        gff = "cupcake/{sample}/A.collapsed.gff",
        fq = "cupcake/{sample}/A.collapsed.rep.fq"
    shell:
        "collapse_isoforms_by_sam.py --input {input.fq} --fq "
            "-s {input.sam} {params.extra} -o {params.out_prefix} "

rule custom_filter:
    input:
        gff = "cupcake/{sample}/A.collapsed.gff",
        groups = "cupcake/{sample}/A.collapsed.group.txt",
        fq = "cupcake/{sample}/A.collapsed.rep.fq"
    params:
        position = "53965961"
    output:
        included = "cupcake/{sample}/included.txt",
        filtered = "cupcake/{sample}/filt.gff",
        groups = "cupcake/{sample}/filt.group.txt",
        fq = "cupcake/{sample}/filt.rep.fq"
    shell:
        "grep 'transcript[[:space:]]' {input.gff} | " #select transcript lines
            "grep '[[:space:]]-[[:space:]]' | "       #select - strand
            "awk '$5 > {params.position}' | "         #select by start site
            "cut -f 9 | cut -f 4 -d '\"' > {output.included}; "
        "grep -wf {output.included} {input.gff} > {output.filtered}; " #filter 
        "grep -wf {output.included} {input.groups} > {output.groups}; "
        "grep -A 3 --no-group-separator -wf {output.included} {input.fq} > {output.fq} "

rule filter_all_sams:
    input:
        expand("hotair_selection/ccs_tracks/{prfx}{sample}{sfx}.cfilt.bam",
                prfx = config["datatype_prefix"],
                sample = config["samples"],
                sfx = config["datatype_suffix"])

rule sam_to_bed:
    input:
        "{sample}.sam"
    output:
        "{sample}.bed"
    shell:
        "bedtools bamtobed -bed12 -i {input} > {output}"

rule custom_bed_filter:
    input:
        bed = "{dir}/{sample}.bed",
        sam = "{dir}/{sample}.sam"
    output:
        included = "{dir}/ccs_tracks/{sample}.included.txt",
        bam = "{dir}/ccs_tracks/{sample}.cfilt.bam"
    shell:
        "awk '$6 == \"-\"' {input.bed} | "   #select - strand
            "awk '$10 > 1' | "               #must have > 1 exon
            "cut -f 4 > {output.included}; " 
        "picard FilterSamReads I={input.sam} O={output.bam} "
            "READ_LIST_FILE={output.included} FILTER=includeReadList "
            "CREATE_INDEX=true "

#rule bed_to_bam:
#    input:
#        genome = expand("{dir}/{name}.chrom.sizes",
#                        dir = config["genome_dir"],
#                        name = config["genome_name"]),
#        bed = "{sample}.bed"
#    output:
#        "{sample}.bam"
#    shell:
#        "bedtools bedtobam -bed12 -i {input.bed} -g {input.genome} | "
#            "samtools sort - -o {output}"

rule generate_counts:
    input:
       report = "new_transcripts/reports/{sample}_polished.cluster_report.csv",
       groups = "cupcake/{sample}/{id}.group.txt"
   # wildcard_constraints:
    #    id = "(?!\d)" #may not contain numbers
    params:
        inprefix = lambda wildcards, input: input.groups[:-10]
    output:
        "cupcake/{sample}/{id}.abundance.txt"
    shell:
        "get_abundance_post_collapse.py {params.inprefix} {input.report}"

rule filter_by_cov:
    input:
        gff = "cupcake/{sample}/{id}.gff",
        counts = "cupcake/{sample}/{id}.abundance.txt",
        fq = "cupcake/{sample}/{id}.rep.fq"
    params: 
        inprefix = lambda wildcards, input: input.gff[:-4]
    wildcard_constraints:
        threshold = "\d+"
    output:
        filtered = "cupcake/{sample}/{id}.min_fl_{threshold}.gff",
        fq = "cupcake/{sample}/{id}.min_fl_{threshold}.rep.fq",
        counts = "cupcake/{sample}/{id}.min_fl_{threshold}.abundance.txt"
    shell:
        "filter_by_count.py --min_count {wildcards.threshold} --dun_use_group_count "
            "{params.inprefix} "

rule add_cov_to_gff:
    input:
        gff = "cupcake/{sample}/{id}.gff",
        counts = "cupcake/{sample}/{id}.abundance.txt"
    output:
        gff = "cupcake/{sample}/{id}-f.gff"
    run:
        tools.addCoverageToGFF(input.counts, input.gff, output.gff)
