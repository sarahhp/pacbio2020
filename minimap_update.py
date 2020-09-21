"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2020-09-07"

FASTQ_DIR = config["fastq_dir"]

import os
import pandas
import get_fastq_sum as tools

rule all:
    input:
        expand("hotair_selection/{sample}{sfx}.mms.sam",
                sample=config["samples"],
                sfx = config["datatype_suffix"]),


rule format_fastq_for_sam:
    input: 
         FASTQ_DIR + "/{sample}.fastq"
    output:
         FASTQ_DIR + "/{sample}.f.fastq"
    shell:
        "sed 's/full_length_coverage=/fc:i:/g' {input} "
            "| sed 's/;length=.*;num_subreads=/ SR:i:/g' "
            "> {output} "


rule minimap:
    input:
        ccs = FASTQ_DIR + "/{sample}.fastq",
        index = config["genome_dir"] + "/minimap_indexes/" + \
                config["genome_file"][:-2] + "splice:hq.mmi"
    params:
        extra = config["params"]["minimap2"]
    threads: 12 
    output:
        sam = "minimap_output/{sample}.sam",
        tmpdir = temp(directory(os.environ["TMPDIR"] + "{sample}"))
    shell:
        "minimap2 -ax splice:hq {params.extra} -t {threads} "
            "--split-prefix {output.tmp_dir} "
            "{input.index} {input.ccs} -o {output.sam} "

rule sam_to_bam:
    input:
        "{sample}.sam"
    threads: 2
    params: 
        extra = config["params"]["bam_filters"]
    output:
        "{sample}.filt.bam"
    shell:
        "samtools view -b -@ {threads} {params.extra} {input} "
            "| samtools sort - -o {output} "

rule index_bam:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule minimap_all:
    input:
        expand("minimap_output/{prfx}{sample}{sfx}.filt.bam",
                prfx = config["datatype_prefix"], #if any
                sample = config["samples"],
                sfx= config["datatype_suffix"]), #if any
        expand("minimap_output/{prfx}{sample}{sfx}.filt.bam.bai",
                prfx = config["datatype_prefix"],#if any
                sample = config["samples"],
                sfx = config["datatype_suffix"]) #if any

rule select_locus:
    input:
        "minimap_output/{sample}.filt.bam"
    params:
        locus = "12:53962308-53974956"
    output:
        "hotair_selection/{sample}.sam"
    shell:
        "samtools view -h {input} {params.locus}  | "
            "samtools sort - -o {output} "

rule list_mulitmappers:
    input:
        "minimap_output/{sample}.sam"
    output:
        "minimap_output/multimapping/{sample}.txt"
    shell:
        "samtools view -f 0x100 {input} | cut -f 1 > {output}; "
        "samtools view -f 0x800 {input} | cut -f 1 >> {output}"

rule select_multimappers:
    input:
        sam = "hotair_selection/{sample}.sam",
        mms = "minimap_output/multimapping/{sample}.txt"
    output:
        "hotair_selection/{sample}.mms.sam"
    shell:
        "picard FilterSamReads I={input.sam} O={output} "
            "READ_LIST_FILE={input.mms} FILTER=includeReadList "

rule cupcake_collapse:
    input: 
        fq = "filtered_fastp/{dataset}_{sample}_hq.f.fastq",
        sam = "hotair_selection/{dataset}_{sample}_hq.sam"
    params:
        out_prefix = lambda wildcards, output: output.gff[:-14],
        extra = "--dun-merge-5-shorter"
    output:
        groups = "cupcake/{dataset}_{sample}/{sample}.collapsed.group.txt",
        gff = "cupcake/{dataset}_{sample}/{sample}.collapsed.gff",
        fq = "cupcake/{dataset}_{sample}/{sample}.collapsed.rep.fq"
    shell:
        "collapse_isoforms_by_sam.py --input {input.fq} --fq "
            "-s {input.sam} {params.extra} -o {params.out_prefix} "

#rule generate_counts:
#    input:
#       report = "new_transcripts/reports/{dataset}_{sample}_polished.cluster_report.csv",
#       groups = "cupcake/{dataset}_{sample}/{sample}.collapsed.group.txt"
#    params:
#        inprefix = lambda wildcards, output: output[0][:-14]
#    output:
#        "cupcake/{dataset}_{sample}/{sample}.collapsed.abundance.txt"
#    shell:
#        "get_abundance_post_collapse.py {params.inprefix} {input.report}"

rule custom_filter:
    input:
        gff = "cupcake/{dataset}_{sample}/{sample}.collapsed.gff",
        groups = "cupcake/{dataset}_{sample}/{sample}.collapsed.group.txt",
        fq = "cupcake/{dataset}_{sample}/{sample}.collapsed.rep.fq"
    params:
        position = "53965961"
    output:
        included = "cupcake/{dataset}_{sample}/{sample}.included.txt",
        filtered = "cupcake/{dataset}_{sample}/{sample}.filt.gff",
        groups = "cupcake/{dataset}_{sample}/{sample}.filt.group.txt",
        fq = "cupcake/{dataset}_{sample}/{sample}.filt.rep.fq"
    shell:
        "grep 'transcript[[:space:]]' {input.gff} | " #select transcript lines
            "grep '[[:space:]]-[[:space:]]' | "       #select - strand
            "awk '$5 > {params.position}' | "         #select by start site
            "cut -f 9 | cut -f 4 -d '\"' > {output.included}; "
        "grep -wf {output.included} {input.gff} > {output.filtered}; " #filter 
        "grep -wf {output.included} {input.groups} > {output.groups}; "
        "grep -A 3 --no-group-separator -wf {output.included} {input.fq} > {output.fq} "

rule add_cov_to_gff:
    input:
        gff = "cupcake/{dataset}_{sample}/{sample}.filt.gff",
        counts = "cupcake/{dataset}_{sample}/{sample}.collapsed.abundance.txt"
    output:
        gff = "cupcake/{dataset}_{sample}/{sample}.filt.f.gff"
    run:
        tools.addCoverageToGFF(input.counts, input.gff, output.gff)

#rule generate_counts_f:
#    input:
#       report = "new_transcripts/reports/{dataset}_{sample}_polished.cluster_report.csv",
#       groups = "cupcake/{dataset}_{sample}/{sample}.filt.group.txt"
#    params:
#        inprefix = lambda wildcards, input: input.groups[:-10]
#    output:
#        "cupcake/{dataset}_{sample}/{sample}.filt.abundance.txt"
#    shell:
#        "get_abundance_post_collapse.py {params.inprefix} {input.report}"

rule generate_counts_g:
    input:
       report = "new_transcripts/reports/{dataset}_{sample}_polished.cluster_report.csv",
       groups = "cupcake/{dataset}_{sample}/{sample}.{inprefix}.group.txt"
    params:
        inprefix = lambda wildcards, input: input.gff[:-10]
    output:
        "cupcake/{dataset}_{sample}/{sample}.{inprefix}.abundance.txt"
    shell:
        "get_abundance_post_collapse.py {params.inprefix} {input.report}"

rule filter_by_cov:
    input:
        gff = "cupcake/{dataset}_{sample}/{sample}.{id}.gff",
        groups = "cupcake/{dataset}_{sample}/{sample}.{id}.group.txt",
        fq = "cupcake/{dataset}_{sample}/{sample}.{id}.rep.fq"
    params: 
        inprefix = lambda wildcards, input: input.gff[:-4]
    output:
        filtered = "cupcake/{dataset}_{sample}/{sample}.{id}.min_fl_{threshold}.gff",
        groups = "cupcake/{dataset}_{sample}/{sample}.{id}.min_fl_{threshold}.group.txt",
        fq = "cupcake/{dataset}_{sample}/{sample}.{id}.min_fl_{threshold}.rep.fq"
    shell:
        "filter_by_count.py --min_count {wildcards.threshold} --dun_use_group_count "
            "{params.inprefix} "
