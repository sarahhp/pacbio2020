"""Align PacBio isoseq or ccs reads to genome with minimap.
And if necessary filter the output sam file by quality or by clipped bases,
select a particular locus to inspect,
or examine multimapping reads.  
   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2020-09-07"

FASTQ_DIR = config["fastq_dir"]

import pandas
import get_fastq_sum as tools

rule all:
    input:
        expand("hotair_selection/{prfx}{sample}{sfx}.mms.sam",
                prfx = config["datatype_prefix"],
                sample=config["samples"],
                sfx = config["datatype_suffix"]),

rule filter_and_select_hotair:
    input:
        expand("hotair_selection/{prfx}{sample}{sfx}.sam",
                prfx = config["datatype_prefix"],
                sample=config["samples"],
                sfx = config["datatype_suffix"]),
        expand("minimap_output/{prfx}{sample}{sfx}.filt.bam.bai",
                prfx = config["datatype_prefix"],
                sample=config["samples"],
                sfx = config["datatype_suffix"])

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
        extra = config["params"]["minimap2"],
        tmp = lambda wildcards: os.environ["TMPDIR"] + wildcards.sample
    threads: 12
    output:
        sam = "minimap_output/{sample}.sam",
    shell:
        "minimap2 -ax splice:hq {params.extra} -t {threads} "
            "--split-prefix {params.tmp} "
            "{input.index} {input.ccs} -o {output.sam} "

rule sam_to_filtered_bam:
    input:
        sam = "{sample}.sam",
        genome = config["genome_dir"] + "/" + config["genome_file"]
    threads: 2
    params: 
        bam_filters = config["params"]["bam_filters"],
        max_clip = config["params"]["max_clip"]
    output:
        "{sample}.filt.bam"
    shell:
        "samclip --ref {input.genome} --max {params.max_clip} < {input.sam} "
            "| samtools view - -b -@ {threads} {params.bam_filters}  "
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
        bam = "minimap_output/{sample}.filt.bam",
        index = "minimap_output/{sample}.filt.bam.bai"
    params:
        locus = "12:53962308-53974956"
    output:
        "hotair_selection/{sample}.sam"
    shell:
        "samtools view -h {input.bam} {params.locus}  | "
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
